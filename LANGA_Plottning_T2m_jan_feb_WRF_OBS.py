#!/usr/bin/env python3



'''Detta program plottar molnbasen från alla körningar i januari och februari utifrån befintliga textfiler Molnbas.txt och
Datum.txt. Dessutom plottas molnbasen för januari och februari enligt observationer.'''

import numpy as N
import scipy.io.netcdf as S
import netCDF4 as N4
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import sys 
import pandas as pd
import os
import datetime


print('Python version ' + sys.version)
print('Pandas version ' + pd.__version__)
print('Matplotlib version ' + matplotlib.__version__)



#####################################################
#PLOTTNING AV OBSERVATIONER FRÅN SODANKYLÄ I JANUARI#
#####################################################


df = pd.read_csv(r'/home/sm_marha/TEXTFILER/Automatic_Station_jan_feb_all_parameters.txt', encoding='latin-1', delimiter=',')

#df_time_cloudbase = df.iloc[:,[4,19]]

df = df.loc[:,['DATE TIME','T']]

df_numpy_array = df.values

TIDSVEKTOR = df_numpy_array[1:, 0]

TIDSVEKTOR_MINUTER = []

for datum in TIDSVEKTOR:
    datum = datum[0:16]
    datum = datetime.datetime.strptime(datum,'%Y-%m-%d %H:%M').strftime('%d/%m %H')
    TIDSVEKTOR_MINUTER.append(datum)

TEMP_C_TIDSSERIE = df_numpy_array[1:, 1].astype(N.float64)

TIDSVEKTOR_TIMMAR = TIDSVEKTOR_MINUTER[0::6]

TEMP_C_GLES = TEMP_C_TIDSSERIE[0::6]


TIDSVEKTOR = N.arange(len(TEMP_C_GLES))

print(len(TIDSVEKTOR_TIMMAR))
print(len(TEMP_C_GLES))

fig = plt.figure(1, figsize =(16, 12))

plt.legend(['WRF'])

#plt.ylim(0,200)
plt.plot(TIDSVEKTOR, TEMP_C_GLES)


#sys.exit()


############################################
#PLOTTNING AV ALLA KÖRNINGAR I JAN OCH FEB # 
############################################

'''Listar alla 00z-körningar'''
file_list = os.listdir('/home/sm_marha/WRF_JAN_FEB_KORNINGAR/ALLA_00_KORNINGAR/')

wrfout_list = []

for file in file_list:
    if file.startswith('wrfout'):
        wrfout_list.append(file)

wrfout_list.sort()

VARIABEL = 'T2'


TIDSVEKTOR_TIMMAR = N.arange(len(wrfout_list))

TEMP_C_2M_TIDSSERIE = []
TIDSVEKTOR_DATUM = []
wrfout_list_datum = []


for index, file in enumerate(wrfout_list):
        print('File:' + str(index))
        f="/home/sm_marha/WRF_JAN_FEB_KORNINGAR/ALLA_00_KORNINGAR/"+file
        Filobjekt = N4.Dataset(f, mode='r')
        TEMP_K_2M = Filobjekt.variables[VARIABEL][:]
        TEMPC_C_2M = TEMP_K_2M[0, 563, 352]-273.15

        wrfout_list_datum.append(file[11:])
        datum = datetime.datetime.strptime(str(file[11:]), '%Y-%m-%d_%H:%M:%S').strftime('%d/%m %H')
        TIDSVEKTOR_DATUM.append(datum)
        
        TEMP_C_2M_TIDSSERIE.append(TEMPC_C_2M)



print(len(TIDSVEKTOR_DATUM))
print(len(TEMP_C_2M_TIDSSERIE))        
        

fig = plt.figure(1, figsize=(16, 12))



plt.plot(TIDSVEKTOR_TIMMAR, TEMP_C_2M_TIDSSERIE)

#plt.locator_params(axis='x', nticks=59)


TIDSVEKTOR_TIMMAR_GLES = TIDSVEKTOR_TIMMAR[::96]
TIDSVEKTOR_DATUM_GLES = TIDSVEKTOR_DATUM[::96]



plt.xticks(TIDSVEKTOR_TIMMAR_GLES, TIDSVEKTOR_DATUM_GLES, rotation= 45, size=12)



plt.yticks(size=12)

#plt.minorticks_on()

#plt.tick_params(which = 'both', direction = 'in', pad = 10)

plt.legend(['OBS', 'WRF'])



#plt.xlim(0,24)
#plt.ylim(0,1000)
plt.ylabel('TYemperature (C)', fontsize = 16)
plt.xlabel('Time (UTC)', fontsize = 16)

plt.title('T2m WRF/OBS Sodankylä January and February', fontsize = 16)


plt.savefig('/home/sm_marha/FIGURES/T2m_WRF_vs_OBS_Sodankyla_jan_feb_1000m')
plt.show()
        

