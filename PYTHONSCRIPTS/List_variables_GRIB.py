import numpy as N
import scipy.io.netcdf as S
import netCDF4 as N4
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import os
import pandas as pd
import datetime
from math import *

GRIBFIL_ua = 'wrf_mlth1ua2018012900_M_D.nc'
GRIBFIL_sfc = 'wrf_mlth1sfc2018012900_M_D.nc'


Filobjekt_GRIB = S.netcdf_file('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/' + GRIBFIL_ua, mode='r')
Filobjekt_GRIB_sfc = S.netcdf_file('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/'+ GRIBFIL_sfc, mode='r')


LAT=27    #SOD:27

LON=191   #SOD:191


'''Parametrar från surface filen'''



SURFACE_GEOPOTENTIAL = Filobjekt_GRIB_sfc.variables['z'][:]
SURFACE_GEOPOTENTIAL_SOD = SURFACE_GEOPOTENTIAL[0, LAT, LON]

SURFACE_PRESSURE = Filobjekt_GRIB_sfc.variables['sp'][:]
SURFACE_PRESSURE_SOD = SURFACE_PRESSURE[0, LAT, LON]




'''Parametrar från upper_air filen'''

MODEL_LEVELS = Filobjekt_GRIB.variables['level'][:]

SPECIFIC_HUMIDITY = Filobjekt_GRIB.variables['q'][:]
SPECIFIC_HUMIDITY_SOD = SPECIFIC_HUMIDITY[0, :, LAT, LON]


TEMPERATURE = Filobjekt_GRIB.variables['t'][:]
TEMPERATURE_SOD = TEMPERATURE[0, :, LAT, LON]
TEMPERATURE_C_SOD = TEMPERATURE_SOD-273.15

CLOUDCOVER = Filobjekt_GRIB.variables['cc'][:]
CLOUDCOVER_SOD = CLOUDCOVER[0, :, LAT, LON]

SNOWDEPTH = Filobjekt_GRIB_sfc.variables['sd'][:]
SNOWDEPTH_SOD = SNOWDEPTH[0, LAT, LON]                                   

CLOUDWATER = Filobjekt_GRIB.variables['clwc'][:]
CLOUDWATER_SOD = CLOUDWATER[0, :, LAT, LON]

CLOUDICE = Filobjekt_GRIB.variables['ciwc'][:]
CLOUDICE_SOD = CLOUDICE[0, :, LAT, LON]

CLOUD_TOTALWATER_SOD = CLOUDWATER_SOD + CLOUDICE_SOD

                                   
TIMESTEPS = Filobjekt_GRIB.variables['time'][:]
VERTICAL_LEVELS = Filobjekt_GRIB.variables['level'][:]


LN_SURFACE_PRESSURE = Filobjekt_GRIB.variables['lnsp'][:]
#SURFACE_PRESSURE = N.exp(LN_SURFACE_PRESSURE)    #Använder Surface pressure på surface i stället, då dessa p g a interpolering skiljer sig lite



'''Läser in a-. och b-koeffecienter'''

df = pd.read_csv(r'/home/sm_marha/TEXTFILER/A_B_Coefficient_ECMWF.csv', encoding='latin-1', delimiter=';', header = None)

#df_time_cloudbase = df.iloc[:,[4,19]]

df = df.iloc[:,:]

df_numpy_array = df.values

A= df_numpy_array[:, 1]
B= df_numpy_array[:, 2]
#TEMP = df_numpy_array[:136, 7]
#Density = df_numpy_array[:136, 8]
#PF = 100*df_numpy_array[:136, 4]
#FI= 9.80665*df_numpy_array[:136, 6]




'''Räknar ut virtuella temperaturen'''

T_VIRTUAL = (1.+0.609133*SPECIFIC_HUMIDITY)*TEMPERATURE

T_VIRTUAL_SOD = T_VIRTUAL[0, :, LAT, LON]


T_VIRTUAL_SOD = T_VIRTUAL_SOD[::-1] #Inverterar vektorn så nivå 137 kommer först)




'''Räknar ut de halva trycknivåerna'''

p_half = N.zeros(len(A))

for i in range(len(A)):

    p = A[i] + B[i] * SURFACE_PRESSURE[0,LAT,LON]

    p_half[i] = p


'''Räknar ut geopotentialen på de halva nivåerna



for i in range(len(p_half)-1):

Fi_half[i+1] = Fi_half[i] + Rd*(  T_VIRTUAL_SOD[i]  *  (  log(p_half[i])-log(p_half[i+1])  ) )'''
     


'''Räknar ut hela trycknivåerna'''
     
p_full = 0.5*(p_half[1:]+p_half[:-1])





'''Räknar ut geopotentialen på de fulla nivåerna.'''

Rd = 287.06
g0 = 9.80665

#TEMP_v = p_full/(Rd*Density) #används till att kontrollera algorithmen via ecmwf.int


Fi = N.zeros(len(p_full))
Fi[0]=10*g0

#FI = N.zeros(len(p_full))  #används till att kontrollera algorithmen via ecmwf.int
#FI[0]=10*g0

for i in range(len(p_full)-1):

    Fi[i+1] = Fi[i] + Rd*(  T_VIRTUAL_SOD[i]  *  (  log(p_full[i])-log(p_full[i+1])  ) )
                              
    #FI[i+1] = FI[i] + Rd*(  PF[i]/(Rd*Density[i])  *  (  log(PF[i])-log(PF[i+1])  )  )    #används till att kontrollera algorithmen via ecmwf.int



'''Plockar ut datat för Sodankylä'''

FULL_LEVELS_SOD = Fi/g0
FULL_LEVELS_SOD = FULL_LEVELS_SOD[::-1]

#FULL_LEVELS_SOD = FULL_LEVELS_SOD[::-1]+222   #lägger till modellhöjden för att få det på samma sätt som i de WRF-filker jag jämför med.

avr = lambda x: "%8.3e" % x

i=0
for CC, CW in zip(CLOUDCOVER_SOD, CLOUDWATER_SOD):
    if CC > 0 or CW > 0:
        print('HEIGHT(' + str(MODEL_LEVELS[i])+') ='  + ' ' + str(int(FULL_LEVELS_SOD[i])) + '      ' \
              + 'CLDFRA = ' + str(round(CLOUDCOVER_SOD[i], 2)) +  '     '  + 'CLW = ' +  str(avr(CLOUDWATER_SOD[i])) \
              +  '      ' + 'CLI = ' + str(avr(CLOUDICE_SOD[i])) +  '     ' + 'CLOUD_TOTALWATER = ' + \
              str(avr(CLOUD_TOTALWATER_SOD[i]))  +  '      ' + 'T =' + str(round(TEMPERATURE_C_SOD[i], 2)) + '     '  + \
              'qv = ' + str(avr(SPECIFIC_HUMIDITY_SOD[i])))
        print('\n')
    i+=1
        
print(SNOWDEPTH_SOD)
