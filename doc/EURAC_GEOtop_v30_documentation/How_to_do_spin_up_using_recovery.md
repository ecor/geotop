# To run a simulation with recovery files:

## Recovery keywords
- ```RecoverSim```
  - simulation recovered
  - [0 - max saving points]
  - default = 0
  - the number of recovered run from which we want to restart the simulation

- ```ContinuousRecovery```
  - number of days after which recovery files are created

- ```SubfolderRecoveryFiles```
  - folder name

- ```SavingPoints```
  - recover every tot days

create the folder rec
Add to the inpts file the following lines:

###################################################################################################
## 1) Using the continuous recovery option
###################################################################################################

You need the settings:
```
SubfolderRecoveryFiles="rec"
ContinuousRecovery = 10 ! number of days after witch recovery files are created
RecoverSim = 0
```
Let it run a simulation with this option and the  folder **rec** will be filled with recovery files.
The time when the recovery is written is saved in the file RecoveryTime.txt.

The number of **ContinuousRecovery** is the time interval in days when recovery files are weritten.

This file is made:
```
Time[s],Time[d],n,i_run,i_sim,cum_time[s],elapsed_time[s],last_recover
2592000.000000,30.000000,3,1,1,2592000.000000,20.000000,1
```
Where
* Time[d], time in days after the simulation start
* i_run, counters of spin up (see later)
* i_sim, counters of spin up (see later)
* cum_time[s], time in s after the first simulation start (you do many recoveries)
* elapsed_time[s], the simulation run time

Only the last recovery is written. The previous is saved with an  extension .old.
In this way, if your simulations is interrupted, you can start from the last recovery.

### To start from a recovery status:

If your simulation crashes for any reason, just restart your simulation
the simulation will start from the last recovery
new timeseries files with the extension **.crec_001** will  be created.
Then you need to manually join the files with the timeseries produced before (.txt) and after recovery (.crec).

### To spin up a simulation manually:

If you need now spin up (i.e. start the simulation again, having as initial state the final state of a previous simulation)
 you need to:
- make sure it writes a recovery at the end of the simulation, setting up a proper # of days for ContinuousRecovery   
- change the file RecoveryTime.txt in the rec folder  in order to make the simulation start at the beginning time and not after the last recovery, or, in other words, change the time stamp of the last recovery to the begining time of the simulation.

**Example:** the file RecoveryTime.txt is as follows (last recovery at day 30 after beginning simulation)
you change from:
```
Time[s],Time[d],n,i_run,i_sim,cum_time[s],elapsed_time[s],last_recover
2592000.000000,30.000000,3,1,1,2592000.000000,20.000000,1
```
to:
```
Time[s],Time[d],n,i_run,i_sim,cum_time[s],elapsed_time[s],last_recover
0.000000,0,1,1,1,0,0,1
```

So the model will use the last recovery as initial state af the time 0 from the beginning simulation.


###################################################################################################
## 2) Using the saving points option
###################################################################################################

In this case, you write recovery files only at specfic times,
defined by and array "savings points" where you specify the days after the start of the simulation
when you wanto to have recovery files

```
ContinuousRecovery = 0
SubfolderRecoveryFiles = "rec"
RecoverSim	=	0
SavingPoints = 20, 50, 1000
```

In this case, in the folder "rec" you will have three sets of recovery files for the days 20, 50, 1000 after the beginning of the simulation.

If you want then to restart a simulation, for example from the third recovery file set (at day 1000 in this case), you write:
```
RecoverSim	=	3
```

###################################################################################################
## 3) Performing spin-up using the NumSimulationTimes option
###################################################################################################

In GEOtop 3.0, Spin-up can be done in the following manner:

```
InitDateDDMMYYYYhhmm        = 01/09/2015 00:00, 31/08/2016 23:00
EndDateDDMMYYYYhhmm        = 31/08/2016 23:00, 31/08/2017 23:00
NumSimulationTimes               = 100, 1
DtPlotPoint 				 = 0,  1
```                    

That's by defining two time periods and generating the output for the final one simulation.
In this case you spin up 100 times the simulation from 01/09/2015 to 31/08/2016 (without writing output for more efficiency)
then you run the simulation one time from 31/08/2016 23:00  to 31/08/2017 23:00 writing the outputs.

This option is particularly useful for long 1d sumulations.

###################################################################################################
## 4) Restart a near real - time simulation using the recovery features
###################################################################################################

We explain  here how to run a simulation in near real time. It means that you run the simulation for a certain period, and then, as soon you have new input meteo data, you continue the simulation starting from the final state of the previous one.  


## Using the continous recovery option

For example we have run a simulation from 30/09/2017 to 30/09/2018 and you write every day a recovery.

```
InitDateDDMMYYYYhhmm = 30/09/2017 23:00
EndDateDDMMYYYYhhmm =  30/09/2018 23:00


SubfolderRecoveryFiles = "rec"
ContinuousRecovery = 1
RecoverSim = 0
```

Once the simulation is finished you want to run the simulation for further 10 days. You need to:

1. Update the meteo input files adding the new week of data
2. change only the final date of the simulation to
```
EndDateDDMMYYYYhhmm =  10/10/2018 23:00
```
3. Restart the simulation

The simulation will run in the same folder, and new ouptut time series maps will appear with the extension **.crec_001**. If you repeat the operation the extension will be **.crec_002** and so on.
The new output maps will be named sequentially starting from the last number of the previous simulation.
