
/* STATEMENT:

 Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 Geotop 2.0.0 - 31 Oct 2013

 Copyright (c), 2013 - Stefano Endrizzi

 This file is part of Geotop 2.0.0

 Geotop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE

 Geotop 2.0.0  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the Geotop model. Any feedback will be highly appreciated.

 If you have satisfactorily used the code, please acknowledge the authors.

 */

#include "constants.h"
#include "struct.geotop.h"
#include "meteo.h"
#include "meteodata.h"
#include "meteodistr.h"
#include "tabs.h"
#include "times.h"
#include "snow.h"
#include "math.optim.h"

extern long number_novalue, number_absent;
extern T_INIT *UV;
extern const char *WORKING_DIRECTORY;
extern char *logfile;
extern long Nl, Nr, Nc;
extern long i_sim;

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void meteo_distr(long *line, long lineLR, METEO *met, WATER *wat, TOPO *top,
                 PAR *par, double JD0, double JDbeg, double JDend)
{

  long i,r,c,year,j;
  double JD;

  //INTERPOLATION OF METEO VARIABLES
  for (i=1; i<=met->st->Z->nh; i++)
    {
      if ( (*par->linear_interpolation_meteo)(i) == 1)
        {
          time_interp_linear(JD0, JDbeg, JDend, met->var[i-1], met->data[i-1],
                             met->numlines[i-1], nmet, iJDfrom0, 1, &(line[i-1]));
        }
      else
        {
          time_interp_constant(JD0, JDbeg, JDend, met->var[i-1], met->data[i-1],
                               met->numlines[i-1], nmet, iJDfrom0, 1, &(line[i-1]));
        }
    }

  //LAPSE RATES
  if (par->LRflag==1)
    {
      time_interp_linear(JD0, JDbeg, JDend, met->LRv, met->LRs, met->LRsnr, nlstot,
                         0, 0, &lineLR);
    }
  else
    {
      convert_JDfrom0_JDandYear ( (JDbeg+JDend)/2., &JD, &year);
      for (i=0; i<nlstot; i++)
        {
          j = floor( JD * met->LRcnc[i] / (365. + is_leap(year)) );
          met->LRv[i] = met->LRc[i][j];
        }
    }

  for (i=0; i<nlstot; i++)
    {
      if ( (long)met->LRv[i] == number_novalue) met->LRv[i] = met->LRd[i];
    }

  //DISTRIBUTION OF METEROLOGICAL VARIABLES FROM MEASUREMENTS IN SOME STATIONS
  Meteodistr((*UV->U)(2), (*UV->U)(1), top->East.get(), top->North.get(), top->Z0.get(), top->curvature1.get(), top->curvature2.get(),
             top->curvature3.get(), top->curvature4.get(), top->slope.get(), top->aspect.get(), met, par->slopewtD, par->curvewtD,
             par->slopewtI, par->curvewtI,
             par->Vmin, par->RHmin, par->dn, par->iobsint, iT, iTdew, iWsx, iWsy, iWs, iPrecInt, met->Tgrid.get(),
             met->RHgrid.get(), met->Vgrid.get(), met->Vdir.get(), met->Pgrid.get(), wat->PrecTot.get(), met->LRv[ilsTa],
             met->LRv[ilsTdew], met->LRv[ilsPrec], par->MaxIncrFactWithElev, par->MinIncrFactWithElev, par->dew,
             par->T_rain, par->T_snow,
             par->snowcorrfact, par->raincorrfact);

  if (par->en_balance==0)
    {
      for (r=1; r<=Nr; r++)
        {
          for (c=1; c<=Nc; c++)
            {
              (*wat->Pnet)(r,c) = par->raincorrfact*(*wat->PrecTot)(r,c)*((JDend-JDbeg)*24.)*cos((*top->slope)(r,c)*GTConst::Pi/180.); //from [mm/h] to [mm]
              if (par->point_sim==1) (*wat->Pnet)(r,c) *= cos((*top->slope)(r,c)*GTConst::Pi/180.);
            }
        }
    }
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double pressure(double Z)
{
  double scale_ht = 8500.0;
  return GTConst::Pa0 * exp(-Z/scale_ht);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double temperature(double Z, double Z0, double T0, double gamma)
{
  return T0 - (Z-Z0)*gamma*1.E-3;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void part_snow(double prec_total, double *prec_rain, double *prec_snow,
               double temperature, double t_rain, double t_snow)

/*  Inputs: prec_total:   matrix of total precipitation
      temperature:  matrix of air temperature
      t_rain:     temperature above wich all precipitation is rain
      t_snow:     temperature below wich all precipitation is snow
  Outputs:prec_snow:    matrix of solid precipitation
      prec_rain:    matrix of liquid precipitation
*/

{

  if (temperature<=t_snow)
    {
      *prec_snow=prec_total;
      *prec_rain=0.0;
    }
  else if (temperature>=t_rain)
    {
      *prec_snow=0.0;
      *prec_rain=prec_total;
    }
  else
    {
      *prec_snow=prec_total*(t_rain-temperature)/(t_rain-t_snow);
      *prec_rain=prec_total*(temperature-t_snow)/(t_rain-t_snow);
    }

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

/*
In the following subroutines:
e: vapour pressure [mbar]
P: atmospheric pressure [mbar]
T: air temperature in [C] 
*/

//Saturated vapour pressure
double SatVapPressure(double T, double P)
{
  double A, b, c,e;
  A=6.1121*(1.0007+3.46E-6*P);
  b=17.502;
  c=240.97; 
  // EC 20200902 why may specific humidity became negative?  // CHECK TEMPERATURE RANGE VALIDITY! https://en.wikipedia.org/wiki/Vapour_pressure_of_water
  // https://en.wikipedia.org/wiki/Clausius%E2%80%93Clapeyron_relation#cite_note-11
  // https://watermark.silverchair.com/bams-86-2-225.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAApEwggKNBgkqhkiG9w0BBwagggJ-MIICegIBADCCAnMGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMfZbADrzR64Sek9raAgEQgIICRGWLRcBBnY2qFIWY3bSyCIUJshRez25alWui0C8DIn_d8Z8-w-F_4_lngTuI5Ak6-438_kHY7RmMBqVE8WPti7tloGd12kyAGa6NjUPAtuIyB_qHrljYMPLciUN9Y9QpaETsxEm6mneoY8uZZhaBjJlM2BISbaw0AkTe7Kgz1hyIEpdz6DL0IideUGgjLdG_baefNK_PYp2SHXpE4IkqrgmEqP1Zc5W_1lw4Gu0cXfboS7WG7Tv1pPnLQZ72Ni6unT4Z8p5OTtoiBHDdPVxcdCDR5B0FnqkesgRCyuqARlLYn8975W2DV3Yt0gZAICX0QUGy9a8iT_o7ONXbettJytGoQD5f0UiAeZK-p0Lc829rbjXuRNRUuXT-pJ0i_u4D1MGDpjYR1w28irA2qfLuj6WVxH2H7nRtIPM0OOwkGIk2F3xdqhCo2mGo6JQym_piFUk-0fHVccS0h4_TVafAXIupX6Sf8FFPXPoOuo6ufRpTdQvRlu1g4-nnnbQcd5_iMo7t000Xt725YcbRB-_dfeR-hQlVueK-uPjYiz8GNFJsY2ufB2yUMckpOtmUXoO2gvJGIHoeMk7M08gdxa6036Wcfzil2RpJStIk4Z9yuV_nKP9QusMIiQoiiPdoMCwlVeHDhte-p60XvFHIZDnENJg_L90Jm5mR8IInEljr4WsqPp7SCoNQqrW-V5JT-JbivM3--QJZ1hdgh8P13ZEHXyHqFRSmmBiyBA57hrRjtvG6YJf7NO-gWGrzdeJrDKkObAbIBXo
  // https://meteo-wagenborgen.nl/wp/2019/09/26/corrections-to-the-august-roche-magnus-equation-in-pwsfwi/
  if (T<(-c+1)) {T=-c+1;} // EC 20200902

  //printf("SatVapPressure T=%f\n",T); // EC 20200902
  //printf("SatVapPressure P=%f\n",P); // EC 20200902
  e=A*exp(b*T/(c+T));
  //printf("SatVapPressure e=%f\n",e); // EC 20200902

  return e;
}

//Finds the saturated vapour pressure and its derivative with respect to temperature
void SatVapPressure_2(double *e, double *de_dT, double T, double P)
{
  double A, b, c;
  A=6.1121*(1.0007+3.46E-6*P);
  b=17.502;
  c=240.97;
  if (T<(-c+1)) {T=-c+1;} // EC 20200902
  *e=A*exp(b*T/(c+T));
  *de_dT=(*e)*(b/(c+T)-b*T/pow_2(c+T));
}

//Temperature from saturated water vapour
double TfromSatVapPressure(double e, double P)
{
  double A, b, c,T;
  printf("TfromSatVapPressure e=%f\n",e); // EC 20200911
  printf("TfromSatVapPressure P=%f\n",P); // EC 20200911
  A=6.1121*(1.0007+3.46E-6*P);
  b=17.502;
  c=240.97;

  T=c*log(e/A)/(b-log(e/A)); // EC 20200902
  printf("TfromSatVapPressure T=%f\n",T); // EC 20200911
  if (T<(-c+1)) {T=-c+1;} // EC 20200902
  printf("TfromSatVapPressure after if control T=%f\n",T); // EC 20200911
  return T;
}

double SpecHumidity(double e, double P)
{
  double q;
  printf("SpecHumidity e=%f\n",e); // EC 20200902
  printf("SpecHumidity P=%f\n",P); // EC 20200902
  q=0.622*e/(P-0.378*e);
  printf("SpecHumidity q=%f\n",q); // EC 20200902
  return q;
}

void SpecHumidity_2(double *Q, double *dQ_dT, double RH, double T, double P)
{
  double e, de_dT, dQ_de;
  SatVapPressure_2(&e, &de_dT, T, P);
  *Q=SpecHumidity(RH*e, P);
  dQ_de=0.622/(P-0.378*e)+0.235116*e/pow_2(P-0.378*e);
  *dQ_dT=dQ_de*de_dT;
}

double VapPressurefromSpecHumidity(double Q, double P)
{
  return Q*P/(0.378*Q+0.622);
}

double Tdew(double T, double RH, double Z)
{
  double e, P;
  P = pressure(Z);
  e = SatVapPressure(T, P);
  return TfromSatVapPressure(e*std::min<double>(RH,1.0), P);
}

double RHfromTdew(double T, double Tdew, double Z)
{
  double e, es, P, RH;
  P = pressure(Z);
  e = SatVapPressure(Tdew, P);
  es = SatVapPressure(T, P);
  RH = e/es;
  if (RH>1) RH=1.0;
  if (RH<0) RH=0.0;
  return RH;
}

double air_density(double T, double Q, double P)    //[kg/m3]
{
  double rho;
  //printf("air_density_01: T=%f\n",T); // EC 20200907
  if (T<(-273.15+1)) {T=-273.15+1;} // EC 20200907
  //printf("air_density_02: T=%f\n",T); // EC 20200907
  rho=P*100/(287.04*(T+273.15))*(1- (Q * P/(0.622+0.368*Q) ) / P*(1-0.622) );
  //printf("air_density_03: rho=%f\n",rho); // EC 20200907

  return (rho);
}

double air_cp(double
              T)   //air specific heat at constant pressure [J/(kg K)] (Garrat,1992)
{
  double cp;

  //printf("air_cp_01: T=%f\n",T); // EC 20200907
  if (T<(-273.15+1)) {T=-273.15+1;} // EC 20200907
  //printf("air_cp_02: T=%f\n",T); // EC 20200907
  cp=1005.00+(T+23.15)*(T+23.15)/3364.0;
  return (cp);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
