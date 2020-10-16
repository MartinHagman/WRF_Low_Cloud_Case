import numpy as N
import netCDF4 as N4
import matplotlib.pyplot as plt
import os
import pandas as pd
import datetime
import sys
import matplotlib.colors as mplc
from matplotlib.colors import LogNorm
from pylab import *


'''Det viktiga i detta skript är att körningarnas längsd måste överensstämma, dvs har jag tiominutersintervall på utfilerns i WRF och 1h i WRF2, så måste jag plocka ut 72
filer ur WRF [0:73] multiplar av 18+1 om vi tittar på 18h-körningar. Väljer jag ut bara 18 även av WRF kommer denna kortare period att smetas ut över hela 18-timmars-
intervallet.

I SCM väljer man multiplar av 1080 per 6h-period och lägger till 1. I 18h-fallet 3241. Denna period är också 18h som i 3D ovan.

xlim behöver man dock ej ändra på här. Det behövs däremot i FINAL_CROSS_SECTIONS_WRF_EC_OBS'''

###############
# OBS NC-FILE #
#################################################################################################################
print('Nc-obs')


Filobjekt_OBS = N4.Dataset('/home/sm_marha/sodankyla_yopp.nc', mode='r')

LONGWAVE_DOWN_SURFACE_OBS = Filobjekt_OBS.variables['rlds'][:]
LONGWAVE_UP_SURFACE_OBS = Filobjekt_OBS.variables['rlus'][:]

TSK_OBS_K = (LONGWAVE_UP_SURFACE_OBS/5.67e-8)**0.25

TSK_OBS_C = TSK_OBS_K - 273.15

TIMESTEPS_NC_OBS = Filobjekt_OBS.variables['time'][:]
time = Filobjekt_OBS.variables['time']     #OBS, ingen array
dates = N4.num2date(TIMESTEPS_NC_OBS, time.units, time.calendar)

ALL_DATES = []

'''Gör om tiden till ett visst datumformat.'''

for date in dates:
    datum = date.strftime('%y%m%d %Hz')
    ALL_DATES.append(datum)
    
for place, date in enumerate(ALL_DATES):
    if date == '180218 00z':                                #ÄNDRA
        print(place)
        datum_START = place
    if date == '180218 18z':                                #ÄNDRA
        print(place)
        datum_END = place






#####################
# AUTOMATIC STATION #
#################################################################################################################
print('Automatic station')


df = pd.read_csv(r'/home/sm_marha/TEXTFILER/Automatic_Station_jan_feb_all_parameters.txt', encoding='latin-1', delimiter=',')

df = df.loc[:,['DATE TIME','T']]

df_numpy_array = df.values

TIDSVEKTOR_OBS = df_numpy_array[1:, 0]

TIDSVEKTOR_MINUTER_OBS = []

'''Letar upp platserna i arrayen för ändpunkterna på de datum jag är intresserad av.'''


for place, date in enumerate(TIDSVEKTOR_OBS):
    if date == '2018-02-18 00:00:00+00':                #ÄNDRA
        print(place)
        date_START = place
    if date == '2018-02-18 18:00:00+00':                #ÄNDRA
        print(place)
        date_END = place


'''Gör om formatet på datumen.'''

for date in TIDSVEKTOR_OBS[date_START:date_END+1]: 
    date = date[0:16]
    date = datetime.datetime.strptime(date,'%Y-%m-%d %H:%M').strftime('%d/%m %H')
    TIDSVEKTOR_MINUTER_OBS.append(date)

T2_TIDSSERIE_OBS = df_numpy_array[date_START:date_END+1, 1]

TIDSVEKTOR_TIMMAR_OBS = TIDSVEKTOR_MINUTER_OBS[0::6]

T2_GLES_OBS = (T2_TIDSSERIE_OBS[0::6])


'''Strängarna görs om till flytal.'''

for place, temp in enumerate(T2_TIDSSERIE_OBS):
    T2_TIDSSERIE_OBS[place] = float(temp)


for place, temp in enumerate(T2_GLES_OBS):
    T2_GLES_OBS[place] = float(temp)
    





#######
# WRF #
################################################################################################################
print('WRF')



'''Ordnar lista med med de olika WRF-körningarna som senare loopas över. '''


wrfout_list = []

file_list = os.listdir('/nobackup/smhid12/sm_marha/2018-02-18_00z_91_LEVELS_QC_QI_42h_T_to_Td_1e-4/run/') 
#file_list = os.listdir('/nobackup/smhid12/sm_marha/2018-02-27_00z_91_LEVELS_Clouds_Where_Clouds_In_Sounding_YSU_TOPD1_CORRECT_CLW_100_dryer_below250/run/')        #ÄNDRA
#file_list = os.listdir('/nobackup/smhid12/sm_marha/2018-02-26_00z_91_LEVELS_REFERENCE/run/')


for file in file_list:
    if file.startswith('wrfout_d01'):
        wrfout_list.append(file)

wrfout_list.sort()

wrfout_list = wrfout_list[0:19]    #ÄNDRA! SÄTT LÄNGDEN!! 



'''Initialisering av tomma arrays'''

QCLOUD_MATRIX_WRF = N.zeros((90, len(wrfout_list)))

THETA_MATRIX_WRF = N.zeros((90, len(wrfout_list)))

THETA_L_MATRIX_WRF = N.zeros((90, len(wrfout_list)))

TEMP_K_MATRIX_WRF = N.zeros((90, len(wrfout_list)))





GLW_TIME_SERIES_WRF = N.empty_like(wrfout_list)

T2_C_TIME_SERIES_WRF = N.empty_like(wrfout_list)

TIME_SERIES_WRF =  N.arange(len(wrfout_list))

TSK_TIME_SERIES_WRF = N.arange(len(wrfout_list))

'''Loopar över respektive wrfout-fil och lagrar de aktuella variabelvärdena i de tomma arrayerna.'''


for place, file in enumerate(wrfout_list):

        print(file)
        
        Filobjekt_WRF = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-18_00z_91_LEVELS_QC_QI_42h_T_to_Td_1e-4/run/'  + file, mode='r')       #ÄNDRA

        LONGWAVE_DOWN_SURFACE_WRF = Filobjekt_WRF.variables['GLW'][0, 563, 352]

        T2_K_WRF = Filobjekt_WRF.variables['T2'][0, 563, 352]

        T2_C_WRF = T2_K_WRF-273.15

        TSK_WRF = Filobjekt_WRF.variables['TSK'][0, 563, 352]-273.15

        

        TSK_TIME_SERIES_WRF[place] = TSK_WRF

        GLW_TIME_SERIES_WRF[place] =  LONGWAVE_DOWN_SURFACE_WRF

        T2_C_TIME_SERIES_WRF[place] = T2_C_WRF
        

    
        
        THETA_PERT_WRF = Filobjekt_WRF.variables['T'][0, :, 563, 352]
        P_PERT_WRF = Filobjekt_WRF.variables['P'][0, :, 563, 352]
        P_BASE_WRF = Filobjekt_WRF.variables['PB'][0, :, 563, 352]   
        
        P_TOT_WRF= P_BASE_WRF + P_PERT_WRF
        THETA_WRF = THETA_PERT_WRF+300
        
        TEMP_K_WRF = THETA_WRF/((1000/(P_TOT_WRF/100))**0.286)
        TEMP_K_MATRIX_WRF[:, place] = TEMP_K_WRF

        
        TEMP_C_WRF = TEMP_K_WRF-273.15


        PH_WRF = Filobjekt_WRF.variables['PH'][:]
        PHB_WRF = Filobjekt_WRF.variables['PHB'][:]
        MODELLNIVAHOJD_WRF = (PH_WRF+PHB_WRF)/9.81
        MASSLEVELS_WRF = 0.5*(MODELLNIVAHOJD_WRF[0,:-1, 563, 352] + MODELLNIVAHOJD_WRF[0, 1:, 563, 352])
        TERRAIN_HEIGHT_WRF = Filobjekt_WRF.variables['HGT'][0, 563, 352]
        MASSLEVELS_1D_MINUS_TER_WRF = MASSLEVELS_WRF - TERRAIN_HEIGHT_WRF

       
        CLOUD_WATER_WRF = Filobjekt_WRF.variables['QCLOUD'][0, :, 563, 352]
        QCLOUD_MATRIX_WRF[:, place] = CLOUD_WATER_WRF


                
        THETA_MATRIX_WRF[:, place] = THETA_WRF        

THETA_L_MATRIX_WRF = THETA_MATRIX_WRF - ((2.26e6*THETA_MATRIX_WRF)/(1004*TEMP_K_MATRIX_WRF)* QCLOUD_MATRIX_WRF)


QCLOUD_MATRIX_WRF = QCLOUD_MATRIX_WRF*1000
        
        
'''Vid tidssteg 0 har inte strålningsfysiken anropats än och värdet är 0. Sätter det därför till samma värde som nästa tidssteg. Plottningstekniskt.'''

GLW_TIME_SERIES_WRF[0] = GLW_TIME_SERIES_WRF[1]

TIME_WRF, LEVELS = N.meshgrid(TIMESTEPS_NC_OBS[datum_START:datum_END+1], MASSLEVELS_1D_MINUS_TER_WRF)  #Lika många tidssteg som timobservationerna då tidssteget i 3D är 1h


TIMESTEPS_WRF = N.arange(datum_START, datum_END,((datum_END - datum_START)/(len(wrfout_list)-1)))

TIMESTEPS_WRF = N.append(TIMESTEPS_WRF,datum_END)

TIME_WRF, LEVELS_WRF = N.meshgrid(TIMESTEPS_WRF, MASSLEVELS_1D_MINUS_TER_WRF)  #Lika många tidssteg som timobservationerna då tidssteget i 3D är 1h



#QCLOUD_MATRIX_WRF = N.transpose(QCLOUD_MATRIX_WRF)  #Här behöver man inte transponera QCLOUD då vi spar undan i i aaray i stället för  i lista
                                                     #som i Plotting_Crossection_6st_WRF_3D   


#########
# WRF 2 #
#################################################################################################################

print('WRF2')


wrfout_list = []

#file_list = os.listdir('/nobackup/smhid12/sm_marha/2018-02-27_00z_91_LEVELS_Clouds_Where_Clouds_In_Sounding_YSU_TOPD0_100_CORRECT_CLW_AMOUNTS/run/')         #ÄNDRA
file_list = os.listdir('/nobackup/smhid12/sm_marha/2018-02-18_00z_91_LEVELS_42h/run/')


for file in file_list:
    if file.startswith('wrfout_d01'):
        wrfout_list.append(file)

wrfout_list.sort()

wrfout_list = wrfout_list[0:109:6] #[0:19]    #ÄNDRA! SÄTT LÄNGDEN!! Kolla upp om det är filer var tionde minut, eller varje timme.



'''Initialisering av tomma arrays'''

'''Initialisering av tomma arrays'''

QCLOUD_MATRIX_WRF2 = N.zeros((90, len(wrfout_list)))

THETA_MATRIX_WRF2 = N.zeros((90, len(wrfout_list)))

THETA_L_MATRIX_WRF2 = N.zeros((90, len(wrfout_list)))

TEMP_K_MATRIX_WRF2 = N.zeros((90, len(wrfout_list)))





GLW_TIME_SERIES_WRF2 = N.empty_like(wrfout_list)

T2_C_TIME_SERIES_WRF2 = N.empty_like(wrfout_list)

TIME_SERIES_WRF2 =  N.arange(len(wrfout_list))

TSK_TIME_SERIES_WRF2 = N.arange(len(wrfout_list))


'''Loopar över respektive wrfout-fil och lagrar de aktuella variabelvärdena i de tomma arrayerna.'''


for place, file in enumerate(wrfout_list):

        print(file)
        
        Filobjekt_WRF2 = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-18_00z_91_LEVELS_42h/run/'  + file, mode='r')       #ÄNDRA

        LONGWAVE_DOWN_SURFACE_WRF2 = Filobjekt_WRF2.variables['GLW'][0, 563, 352]

        T2_K_WRF2 = Filobjekt_WRF2.variables['T2'][0, 563, 352]

        T2_C_WRF2 = T2_K_WRF2-273.15

        TSK_WRF2 = Filobjekt_WRF2.variables['TSK'][0, 563, 352]-273.15

        GLW_TIME_SERIES_WRF2[place] =  LONGWAVE_DOWN_SURFACE_WRF2

        T2_C_TIME_SERIES_WRF2[place] = T2_C_WRF2

        TSK_TIME_SERIES_WRF2[place] = TSK_WRF2

        

    
        
        THETA_PERT_WRF2 = Filobjekt_WRF2.variables['T'][0, :, 563, 352]
        P_PERT_WRF2 = Filobjekt_WRF2.variables['P'][0, :, 563, 352]
        P_BASE_WRF2 = Filobjekt_WRF2.variables['PB'][0, :, 563, 352]   
        
        P_TOT_WRF2= P_BASE_WRF2 + P_PERT_WRF2
        THETA_WRF2 = THETA_PERT_WRF2+300
        
        TEMP_K_WRF2 = THETA_WRF/((1000/(P_TOT_WRF2/100))**0.286)
        TEMP_K_MATRIX_WRF2[:, place] = TEMP_K_WRF2

        
        TEMP_C_WRF2 = TEMP_K_WRF2-273.15

        
        

        PH_WRF2 = Filobjekt_WRF2.variables['PH'][:]
        PHB_WRF2 = Filobjekt_WRF2.variables['PHB'][:]
        MODELLNIVAHOJD_WRF2 = (PH_WRF2+PHB_WRF2)/9.81
        MASSLEVELS_WRF2 = 0.5*(MODELLNIVAHOJD_WRF2[0,:-1, 563, 352] + MODELLNIVAHOJD_WRF2[0, 1:, 563, 352])
        TERRAIN_HEIGHT_WRF2 = Filobjekt_WRF2.variables['HGT'][0, 563, 352]
        MASSLEVELS_1D_MINUS_TER_WRF2 = MASSLEVELS_WRF2 - TERRAIN_HEIGHT_WRF2

       
        CLOUD_WATER_WRF2 = Filobjekt_WRF2.variables['QCLOUD'][0, :, 563, 352]
        QCLOUD_MATRIX_WRF2[:, place] = CLOUD_WATER_WRF2


                
        THETA_MATRIX_WRF2[:, place] = THETA_WRF2

        
THETA_L_MATRIX_WRF2 = THETA_MATRIX_WRF2 - ((2.26e6*THETA_MATRIX_WRF2)/(1004*TEMP_K_MATRIX_WRF2)* QCLOUD_MATRIX_WRF2)


QCLOUD_MATRIX_WRF2 = QCLOUD_MATRIX_WRF2*1000



        
'''Vid tidssteg 0 har inte strålningsfysiken anropats än och värdet är 0. Sätter det därför till samma värde som nästa tidssteg. Plottningstekniskt.'''

GLW_TIME_SERIES_WRF2[0] = GLW_TIME_SERIES_WRF2[1]

TIMESTEPS_WRF2 = N.arange(datum_START, datum_END,((datum_END - datum_START)/(len(wrfout_list)-1)))

TIMESTEPS_WRF2 = N.append(TIMESTEPS_WRF2,datum_END) 

TIME_WRF2, LEVELS_WRF2 = N.meshgrid(TIMESTEPS_WRF2, MASSLEVELS_1D_MINUS_TER_WRF2)  #Lika många tidssteg som timobservationerna då tidssteget i 3D är 1h





###########
# WRF SCM #
################################################################################################################
print('WRF_SCM')


t=3241
        
Filobjekt_SCM_WRF = N4.Dataset('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfout_d01_2018-02-18_00:00:00_QC_QI_T_to_Td_1e-5', mode='r')       #ÄNDRA

TIME_SCM_WRF = Filobjekt_SCM_WRF.variables['XTIME'][0:t]
time_SCM = Filobjekt_SCM_WRF.variables['XTIME']     #OBS, ingen array
dates_SCM = N4.num2date(TIME_SCM_WRF, time.units)



LONGWAVE_DOWN_SURFACE_SCM_WRF = Filobjekt_SCM_WRF.variables['GLW'][0:t, 1, 1]

T2_K_SCM_WRF = Filobjekt_SCM_WRF.variables['T2'][0:t, 1, 1]

TSK_SCM_WRF = Filobjekt_SCM_WRF.variables['TSK'][0:t, 1, 1]-273.15

T2_C_SCM_WRF = T2_K_SCM_WRF-273.15


#Tvåmeterstemperaturen finns ej vid T0, då den diagnosticeras

T2_C_SCM_WRF[0] = T2_C_SCM_WRF[1]


THETA_PERT_SCM_WRF = Filobjekt_SCM_WRF.variables['T'][0:t, :, 1, 1]
P_PERT_SCM_WRF = Filobjekt_SCM_WRF.variables['P'][0:t, :, 1, 1]
P_BASE_SCM_WRF = Filobjekt_SCM_WRF.variables['PB'][0:t, :, 1, 1]   

P_TOT_SCM_WRF= P_BASE_SCM_WRF + P_PERT_SCM_WRF
THETA_SCM_WRF = THETA_PERT_SCM_WRF+300
TEMP_K_SCM_WRF = THETA_SCM_WRF/((1000/(P_TOT_SCM_WRF/100))**0.286)
TEMP_C_SCM_WRF = TEMP_K_SCM_WRF-273.15     


PH_SCM_WRF = Filobjekt_SCM_WRF.variables['PH'][:]
PHB_SCM_WRF = Filobjekt_SCM_WRF.variables['PHB'][:]
MODELLNIVAHOJD_SCM_WRF = (PH_SCM_WRF+PHB_SCM_WRF)/9.81
MASSLEVELS_SCM_WRF = 0.5*(MODELLNIVAHOJD_SCM_WRF[0:t,:-1, 1, 1] + MODELLNIVAHOJD_SCM_WRF[0:t, 1:, 1, 1])
TERRAIN_HEIGHT_SCM_WRF = Filobjekt_SCM_WRF.variables['HGT'][0, 1, 1]
MASSLEVELS_1D_MINUS_TER_SCM_WRF = MASSLEVELS_SCM_WRF - TERRAIN_HEIGHT_SCM_WRF


CLOUD_WATER_SCM_WRF = Filobjekt_SCM_WRF.variables['QCLOUD'][0:t, :, 1, 1]

        
#Initialisering av tidsmatris med samma shape som MASSLEVELS-matrisen

TIMESTEP_MATRIX_SCM = N.zeros_like(MASSLEVELS_1D_MINUS_TER_SCM_WRF)


#Fyll varje kolumn med tidsvektorn

for column in range(N.shape(TIMESTEP_MATRIX_SCM)[1]):
    TIMESTEP_MATRIX_SCM[:, column] = TIME_SCM_WRF        
        

#Transponera så att antal rade blir antal nivåer i ställer för antal tidssteg

TIMESTEP_MATRIX_SCM = N.transpose(TIMESTEP_MATRIX_SCM)
MASSLEVELS_1D_MINUS_TER_SCM_WRF = N.transpose(MASSLEVELS_1D_MINUS_TER_SCM_WRF)


THETA_L_SCM_WRF = THETA_SCM_WRF - ((2.26e6*THETA_SCM_WRF)/(1004*TEMP_K_SCM_WRF)* CLOUD_WATER_SCM_WRF)

THETA_SCM_WRF = N.transpose(THETA_SCM_WRF)
THETA_L_SCM_WRF = N.transpose(THETA_L_SCM_WRF)
CLOUD_WATER_SCM_WRF = N.transpose(CLOUD_WATER_SCM_WRF)

#Lägger ut tidsstegen mellan datum_START och datum_END, för att det ska gå att plotta på samma axel för temperatur och longwave-graferna.
#OBS!!!! Fungerar bara om SCM-körningen är lika lång som det timintervall man valt i START och END

TIMESTEPS_SCM_WRF = N.arange(datum_START, datum_END,(datum_END-datum_START)/len(TIME_SCM_WRF))



LONGWAVE_DOWN_SURFACE_SCM_WRF[0] = LONGWAVE_DOWN_SURFACE_SCM_WRF[1]



if len(LONGWAVE_DOWN_SURFACE_SCM_WRF) == 3241:
    LONGWAVE_DOWN_SURFACE_SCM_WRF = N.append(LONGWAVE_DOWN_SURFACE_SCM_WRF, LONGWAVE_DOWN_SURFACE_SCM_WRF[-1])


if len(T2_C_SCM_WRF) == 3241:
    T2_C_SCM_WRF = N.append(T2_C_SCM_WRF, T2_C_SCM_WRF[-1])


if len(TSK_SCM_WRF) == 3241:
    TSK_SCM_WRF = N.append(TSK_SCM_WRF, TSK_SCM_WRF[-1])


#############
# WRF SCM 2 #
#################################################################################################################
print('SCM2')

t2 = 3241

    
Filobjekt_SCM_WRF2 = N4.Dataset('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfout_d01_2018-02-18_00:00:00_REFERENCE', mode='r')
                               

TIME_SCM_WRF2 = Filobjekt_SCM_WRF2.variables['XTIME'][0:t2]
time_SCM2 = Filobjekt_SCM_WRF2.variables['XTIME']     #OBS, ingen array
dates_SCM2 = N4.num2date(TIME_SCM_WRF2, time.units)
                               

LONGWAVE_DOWN_SURFACE_SCM_WRF2 = Filobjekt_SCM_WRF2.variables['GLW'][0:t2, 1, 1]

T2_K_SCM_WRF2 = Filobjekt_SCM_WRF2.variables['T2'][0:t2, 1, 1]

T2_C_SCM_WRF2 = T2_K_SCM_WRF2-273.15

TSK_SCM_WRF2 = Filobjekt_SCM_WRF2.variables['TSK'][0:t2, 1, 1]-273.15


T2_C_SCM_WRF2[0] = T2_C_SCM_WRF2[1]


THETA_PERT_SCM_WRF2 = Filobjekt_SCM_WRF2.variables['T'][0:t2, :, 1, 1]
P_PERT_SCM_WRF2 = Filobjekt_SCM_WRF2.variables['P'][0:t2, :, 1, 1]
P_BASE_SCM_WRF2 = Filobjekt_SCM_WRF2.variables['PB'][0:t2, :, 1, 1]   

P_TOT_SCM_WRF2= P_BASE_SCM_WRF2 + P_PERT_SCM_WRF2
THETA_SCM_WRF2 = THETA_PERT_SCM_WRF2+300
TEMP_K_SCM_WRF2 = THETA_SCM_WRF2/((1000/(P_TOT_SCM_WRF2/100))**0.286)
TEMP_C_SCM_WRF2 = TEMP_K_SCM_WRF2-273.15



PH_SCM_WRF2 = Filobjekt_SCM_WRF2.variables['PH'][:]
PHB_SCM_WRF2 = Filobjekt_SCM_WRF2.variables['PHB'][:]
MODELLNIVAHOJD_SCM_WRF2 = (PH_SCM_WRF2+PHB_SCM_WRF2)/9.81
MASSLEVELS_SCM_WRF2 = 0.5*(MODELLNIVAHOJD_SCM_WRF2[0:t2,:-1, 1, 1] + MODELLNIVAHOJD_SCM_WRF2[0:t2, 1:, 1, 1])
TERRAIN_HEIGHT_SCM_WRF2 = Filobjekt_SCM_WRF2.variables['HGT'][0, 1, 1]
MASSLEVELS_1D_MINUS_TER_SCM_WRF2 = MASSLEVELS_SCM_WRF2 - TERRAIN_HEIGHT_SCM_WRF2


CLOUD_WATER_SCM_WRF2 = Filobjekt_SCM_WRF2.variables['QCLOUD'][0:t2, :, 1, 1]
CLOUD_WATER_SCM_WRF2 = N.transpose(CLOUD_WATER_SCM_WRF2)
        
#Initialisering av tidsmatris med samma shape som MASSLEVELS-matrisen

TIMESTEP_MATRIX_SCM2 = N.zeros_like(MASSLEVELS_1D_MINUS_TER_SCM_WRF2)


#Fyll varje kolumn med tidsvektorn

for column in range(N.shape(TIMESTEP_MATRIX_SCM2)[1]):
    TIMESTEP_MATRIX_SCM2[:, column] = TIME_SCM_WRF2        
        

#Transponera så att antal rade blir antal nivåer i ställer för antal tidssteg

TIMESTEP_MATRIX_SCM2 = N.transpose(TIMESTEP_MATRIX_SCM2)
MASSLEVELS_1D_MINUS_TER_SCM_WRF2 = N.transpose(MASSLEVELS_1D_MINUS_TER_SCM_WRF2)

THETA_SCM_WRF2 = N.transpose(THETA_SCM_WRF2)



#Lägger ut tidsstegen mellan datum_START och datum_END, för att det ska gå att plotta på samma axel för temperatur och longwave-graferna.
#OBS!!!! Fungerar bara om SCM-körningen är lika lång som det timintervall man valt i START och END

TIMESTEPS_SCM_WRF2 = N.arange(datum_START, datum_END,(datum_END-datum_START)/len(TIME_SCM_WRF2))  #TIMESTEPS_SCM_WRF_T2 kan lika gärna vara TIMESTEPS_SCM_WRF_LONG_WAVE, då det bara var ett sätt att namnge tidssteget för 2-dimensionella tidsserierna.



LONGWAVE_DOWN_SURFACE_SCM_WRF2[0] = LONGWAVE_DOWN_SURFACE_SCM_WRF2[1]


if len(LONGWAVE_DOWN_SURFACE_SCM_WRF2) == 3241:
    LONGWAVE_DOWN_SURFACE_SCM_WRF2 = N.append(LONGWAVE_DOWN_SURFACE_SCM_WRF2, LONGWAVE_DOWN_SURFACE_SCM_WRF2[-1])


if len(T2_C_SCM_WRF2) == 3241:
    T2_C_SCM_WRF2 = N.append(T2_C_SCM_WRF2, T2_C_SCM_WRF2[-1])

if len(TSK_SCM_WRF2) == 3241:
    TSK_SCM_WRF2 = N.append(TSK_SCM_WRF2, TSK_SCM_WRF2[-1])




#########
# ECMWF #
#################################################################################################################

print('ECMWF')



Filobjekt_GRIB_sfc = N4.Dataset('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1sfc2018021800_M_D.nc' , mode='r')       #ÄNDRA


TIMESTEPS_EC = Filobjekt_GRIB_sfc.variables['time'][:]
time = Filobjekt_GRIB_sfc.variables['time']     #OBS, ingen array
dates = N4.num2date(TIMESTEPS_EC, time.units, time.calendar)
ALLA_DATUM=[]




'''Gör om tiden till ett visst datumformat.'''

for date in dates:
    datum = date.strftime('%d/%m %Hz')
    ALLA_DATUM.append(datum)


ECMWF_TIMESTEPS = N.arange(datum_START, datum_END+1, 3)    #Fixar nya TIMESTEPS som börjar från datum_START för att x-axlarna ska sammanfalla


T2_TEMP_K_EC_TIDSSERIE = Filobjekt_GRIB_sfc.variables['t2m'][0:len(ECMWF_TIMESTEPS) , 27, 192]

T2_TEMP_C_EC_TIDSSERIE = T2_TEMP_K_EC_TIDSSERIE - 273.15




#sys.exit()



############
# PLOTTING #
#################################################################################################################




print('PLOTTING CLOUD WATER')


        
fig1 = plt.figure(1, figsize=(14,12))

fig1.subplots_adjust(right=0.80)

fig1.subplots_adjust(hspace=0.2)

cbar_ax = fig1.add_axes([0.82, 0.53, 0.01, 0.35])

rc('axes', linewidth=1.5)

'''Här kan man välja om man vill ha WRF_SCM eller WRF_3D i första subploten. Vill man byta till WRF2 eller WRF_SCM2 byter man inne i plotkommandona'''

ax1 = plt.subplot(211)
ax1.tick_params(which='major', direction='in')


####Plottar molnvatten##########


# SCM  (Välj SCM eller 3D)

'''
norm = mplc.Normalize(0, 5e-1)
Antal = 5e-1/20
myplot = ax1.contourf(TIMESTEP_MATRIX_SCM, MASSLEVELS_1D_MINUS_TER_SCM_WRF,CLOUD_WATER_SCM_WRF, N.arange(0, 5e-1, Antal), norm = norm, extend = 'both', cmap ='binary')  
myplot = ax1.contourf(TIMESTEP_MATRIX_SCM, MASSLEVELS_1D_MINUS_TER_SCM_WRF, CLOUD_WATER_SCM_WRF, extend = 'both', cmap ='binary')   #decode("utf-8") Denna behövs med scipy, men ej med Netcdf4
cbar = fig1.colorbar(myplot, cax=cbar_ax, format = "%8.1e")
contour_levels = N.arange(0, 330, 1)
myplot_2 = ax1.contour(TIMESTEP_MATRIX_SCM, MASSLEVELS_1D_MINUS_TER_SCM_WRF, THETA_SCM_WRF, contour_levels, colors='red', linestyles='dotted', linewidth='0.2')
myplot_2 = ax1.contour(TIMESTEP_MATRIX_SCM, MASSLEVELS_1D_MINUS_TER_SCM_WRF, THETA_L_SCM_WRF, contour_levels, colors='blue', linestyles='dotted', linewidth='0.2')
#plt.clabel(myplot_2, inline=True, fmt='%1.1f', fontsize=10)
ax1.set_xticks(N.arange(0,(len(TIMESTEPS_SCM_WRF)/3), (len(TIMESTEPS_SCM_WRF)-1)/(len(wrfout_list)-1))) # DELA MED LÄNGDEN AV PROGNOSEN
ax1.set_xticklabels([])
ax1.set_ylim(0, 2000)
ax1.set_ylabel('Height (m)', fontsize=12)
#ax1.set_title('Cloudwater and potential temperature WRF', fontsize=12)

cbar.ax.set_title('                 (g/kg)')

'''


# 3D    (Välj SCM eller 3D)

norm = mplc.Normalize(0, 5e-1)
Antal = 5e-1/20
myplot = ax1.contourf(TIME_WRF, LEVELS_WRF, QCLOUD_MATRIX_WRF, N.arange(0, 5e-1, Antal), norm = norm, extend = 'both', cmap ='binary')   #decode("utf-8") Denna behövs med scipy, men ej med Netcdf4
#myplot = ax1.contourf(TIME_WRF, LEVELS_WRF, QCLOUD_MATRIX_WRF, extend = 'both', cmap ='binary')
#myplot = ax1.scatter(TIME_WRF, LEVELS_WRF, QCLOUD_MATRIX_WRF, cmap ='binary')
cbar = fig1.colorbar(myplot,cax=cbar_ax, format = "%8.1e")
contour_levels = N.arange(0, 330, 1.5)
myplot_2 = ax1.contour(TIME_WRF, LEVELS_WRF, THETA_MATRIX_WRF, contour_levels, colors='red', linestyles='solid', linewidths=0.7)
plt.clabel(myplot_2, contour_levels[::2], inline=True, fmt='%1.0f', fontsize=14)
myplot_2 = ax1.contour(TIME_WRF, LEVELS_WRF, THETA_L_MATRIX_WRF, contour_levels, colors='blue', linestyles='solid', linewidths=1)
plt.clabel(myplot_2, contour_levels[::2], inline=True, fmt='%1.0f', fontsize=14)
ax1.set_xticks(TIMESTEPS_NC_OBS[datum_START:datum_END+1:1])
ax1.set_xticklabels([])
ax1.set_ylim(0, 2000)
ax1.set_xlim(datum_START-0.05,datum_END)
ax1.set_ylabel('Height ($m$)', fontsize=20)
#ax1.set_title('Cloudwater and potential temperature WRF', fontsize=12)
ax1.text(7049.5, 1900, 'a', fontsize = 18, fontweight='bold')
cbar.ax.set_title('                 $(gkg^{-1})$', fontsize = 15)
cbar.ax.tick_params(labelsize = 16)
ax1.tick_params(labelsize=15)

'''Plottar masslevels längs y-axeln''' 
X = N.ones_like(MASSLEVELS_1D_MINUS_TER_WRF)*datum_START
ax1.scatter(X, MASSLEVELS_1D_MINUS_TER_WRF, marker= '>', color='black', s=7)




####Plottar långvåg##########

'''Kombinera fritt vilka fält som ska plottas'''


print('LW')

# Plottar lw-down radiation vid ytan


ax2 = plt.subplot(413)
ax2.tick_params(which='major', direction='in')

ax2.plot(TIMESTEPS_NC_OBS[datum_START:datum_END+1], LONGWAVE_DOWN_SURFACE_OBS[datum_START:datum_END+1], color = 'k', linewidth = 2, label = 'OBS')
ax2.plot(TIMESTEPS_WRF, GLW_TIME_SERIES_WRF, color = 'C1', linewidth = 1.5, label = '3D QC')   #Använder TIMESTEPS från obs nc-file, då dessa har värder 7032-7074
ax2.plot(TIMESTEPS_WRF2, GLW_TIME_SERIES_WRF2, color = 'C2', linewidth = 1.5, label = '3D REF')   #De andra(TIDSVEKTOR_TIMMAR_OBS och TIME_SERIES) starta på 0 och går till 43)
ax2.plot(TIMESTEPS_SCM_WRF, LONGWAVE_DOWN_SURFACE_SCM_WRF, color = 'C3', linewidth = 1.5, label = 'SCM QC')
ax2.plot(TIMESTEPS_SCM_WRF2, LONGWAVE_DOWN_SURFACE_SCM_WRF2, color = 'C4', linewidth = 1.5, label = 'SCM REF')  
#ax2.legend(fontsize = 'large', loc='upper left', bbox_to_anchor=(1, 0.85))
ax2.set_xticks(TIMESTEPS_NC_OBS[datum_START:datum_END+1:1])
#ax2.set_xticklabels(ALL_DATES[datum_START:datum_END+1:3], rotation=45, ha='right', fontsize=12)
ax2.set_xticklabels([])
ax2.set_xlim(datum_START-0.05, datum_END)
ax2.set_ylabel('Radiation  ($Wm^{-2}$)', fontsize=20) 
#ax2.set_xlabel('Time (Hours)')
#ax2.set_title('Longwave radiation down', fontsize=12)
ax2.grid(linestyle="dotted")
ax2.text(7049.5, 258, 'b', fontsize = 18, fontweight='bold')
ax2.tick_params(labelsize=14)






####Plottar Temp2m##########


'''Kombinera fritt vilka fält som ska plottas'''




print('2m temp')

# Plottar 2-meters temperatur

ax3 = plt.subplot(414)
ax3.tick_params(which='major', direction='in')

ax3.plot(TIMESTEPS_NC_OBS[datum_START:datum_END+1], T2_GLES_OBS, color = 'k', linewidth = 2, label = 'OBS')
ax3.plot(TIMESTEPS_NC_OBS[datum_START:datum_END+1], TSK_OBS_C[datum_START:datum_END+1],linestyle = 'dashed', color = 'k', linewidth = 1.3, label = 'OBS SKIN')
ax3.plot(TIMESTEPS_WRF, T2_C_TIME_SERIES_WRF, color = 'C1', linewidth = 1.5, label = '3D QC')
ax3.plot(TIMESTEPS_WRF2, T2_C_TIME_SERIES_WRF2, color = 'C2', linewidth = 1.5, label = '3D REF')
ax3.plot(TIMESTEPS_SCM_WRF, T2_C_SCM_WRF, color = 'C3', linewidth = 1.5, label = 'SCM QC')     
ax3.plot(TIMESTEPS_SCM_WRF2, T2_C_SCM_WRF2, color = 'C4', linewidth = 1.5, label = 'SCM REF')
ax3.plot(TIMESTEPS_WRF, TSK_TIME_SERIES_WRF, color = 'C1',linewidth = 1.5, linestyle = 'dashed', label = '3D QC SKIN')
ax3.plot(TIMESTEPS_WRF2, TSK_TIME_SERIES_WRF2,color = 'C2', linewidth = 1.5, linestyle = 'dashed', label = '3D REF SKIN')
ax3.plot(TIMESTEPS_SCM_WRF, TSK_SCM_WRF, color = 'C3',linestyle = 'dashed', linewidth = 1.5, label = 'SCM QC SKIN')     
ax3.plot(TIMESTEPS_SCM_WRF2, TSK_SCM_WRF2,color = 'C4', linestyle = 'dashed', linewidth = 1.5, label = 'SCM REF SKIN')

ax3.plot(ECMWF_TIMESTEPS, T2_TEMP_C_EC_TIDSSERIE, label = 'IFS')
plt.legend().get_frame().set_linewidth(1)
ax3.legend(fontsize = 'xx-large', loc='center left', bbox_to_anchor=(1, 1.1))
ax3.set_xlim(datum_START-0.05, datum_END)
#ax3.set_ylim(-30, 0)
ax3.set_xticks(TIMESTEPS_NC_OBS[datum_START:datum_END+1:1])
ax3.set_xticklabels(ALL_DATES[datum_START:datum_END+1:1], rotation=30, ha='right', fontsize=13)
#ax3.set_xlabel('Time (Hours)', fontsize=15)
ax3.set_ylabel('T $(C)$', fontsize=20)
ax3.tick_params(labelsize=14)

#ax3.set_title('2 meter and skin temperature', fontsize=12)
ax3.grid(linestyle="dotted")
ax3.text(7049.5, -10, 'c', fontsize = 18, fontweight='bold')
#fig1.subplots_adjust(hspace=0.4)



fig1.savefig('/home/sm_marha/FIGURES/FINAL/20180218/Final_comparison_3_subplot_180218_T_to_Td_1e-4_18h')       #ÄNDRA
plt.show()
