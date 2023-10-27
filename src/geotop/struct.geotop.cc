//
// Created by alberto on 5/2/18.
//
#include "struct.geotop.h"
#include "t_datamanipulation.h"

STATEVAR_3D::STATEVAR_3D(double nan, long nl, long nr, long nc) :

        type{new Matrix<short>{nr, nc}},
        lnum{new Matrix<long>{nr, nc}},
        Dzl{new Tensor<double>{nl, nr, nc}},
        w_liq{new Tensor<double>{nl, nr, nc}},
        w_ice{new Tensor<double>{nl, nr, nc}},
        T{new Tensor<double>{nl, nr, nc}}
{
  *type = 2;
  *lnum = 0;
  *Dzl = 0.;
  *T = nan;
  *w_ice = 0.;
  *w_liq = 0.;
}

STATEVAR_3D::~STATEVAR_3D() = default;

// added by EC 20210323

STATEVAR_1D::STATEVAR_1D(double nan, long nl) :


     //   type{new short},
     //  lnum{new long},
        Dzl{new Vector<double>{nl}},
        w_liq{new Vector<double>{nl}},
        w_ice{new Vector<double>{nl}},
        T{new Vector<double>{nl}}
{
  //*type = 2;
  //*lnum = 0;
  *Dzl = 0.;
  *T = nan;
  *w_ice = 0.;
  *w_liq = 0.;
}














SOIL_STATE::SOIL_STATE(const long n, const long nl) :

        P{new Matrix<double>{nl,0,n,1}},
        thi{new Matrix<double>{nl,n}},
        T{new Matrix<double>{nl,n}} {}

