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


TEMP_C_GLES = TEMP_C_TIDSSERIE[0::6]


TIDSVEKTOR = N.arange(len(TEMP_C_GLES))

data_dict_jan = {}


j=1

for i in range(0, 720, 24):
    data_dict_jan['Day_%s' % str(j).zfill(2)] = TEMP_C_GLES[i:i+24]
    j+=1


TIDSVEKTOR_TIMMAR = N.arange(len(data_dict_jan['Day_01'][:]))


fig_1, axs_1 = plt.subplots(5, 6, figsize=(12,10))
fig_1.subplots_adjust(hspace=0.3, wspace=0.3)


k=1

for i in range(5):
    for j in range(6):
        print('New day OBS')
        ax = axs_1[i, j]
        ax.set_title(str(k)+ '/1')
        if i==0 or i==1 or i==2 or i==3:
            ax.set_xticklabels([])
        if j==0:
            ax.set_ylabel('Temp (C)')
        if i==4:
            ax.set_xlabel('Hours')
            
        ax.plot(TIDSVEKTOR_TIMMAR, data_dict_jan['Day_%s' % str(k).zfill(2)], label = 'OBS')
        
        k+=1




data_dict_feb = {}

j=0

for i in range(720, len(TEMP_C_GLES), 24):
    data_dict_feb['Day_%s' % str(j).zfill(2)] = TEMP_C_GLES[i:i+24]  
    j+=1

fig_2, axs_2 = plt.subplots(5, 6, figsize=(12,10))
fig_2.subplots_adjust(hspace=0.3, wspace=0.3)

k=0


    
for i in range(5):
    for j in range(6):
        
        print('New day OBS')
        ax = axs_2[i, j]
        if k==0:
            ax.set_title('31/1')
        else:
            ax.set_title(str(k)+ '/2')
        if i==0 or i==1 or i==2 or i==3:
            ax.set_xticklabels([])
        if j==0:
            ax.set_ylabel('Temp (C)')
        if i==4:
                ax.set_xlabel('Hours')
            
        ax.plot(TIDSVEKTOR_TIMMAR, data_dict_feb['Day_%s' % str(k).zfill(2)], label = 'OBS')

        if k==28:
            break
        k+=1

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




data_dict_jan = {}

j=1

for i in range(0, 720, 24):
    data_dict_jan['Day_%s' % str(j).zfill(2)] = wrfout_list[i:i+24]
    j+=1
    




VARIABEL = 'T2'


TIDSVEKTOR_TIMMAR = N.arange(len(data_dict_jan['Day_01'][:]))
TIDSVEKTOR_DATUM = []
wrfout_list_datum = []


k=1

for i in range(5):
    for j in range(6):
        print('New day')
        ax = axs_1[i, j]
        TEMP_C_2M_TIDSSERIE = []
        for file in data_dict_jan['Day_%s' % str(k).zfill(2)]:
            f="/home/sm_marha/WRF_JAN_FEB_KORNINGAR/ALLA_00_KORNINGAR/"+file
            Filobjekt = N4.Dataset(f, mode='r')            
            TEMP_K_2M = Filobjekt.variables[VARIABEL][:]
            TEMPC_C_2M = TEMP_K_2M[0, 563, 352]-273.15
            TEMP_C_2M_TIDSSERIE.append(TEMPC_C_2M)
        ax.plot(TIDSVEKTOR_TIMMAR,  TEMP_C_2M_TIDSSERIE, label = 'WRF')
        
        k+=1

plt.legend()


fig_1.savefig('/home/sm_marha/FIGURES/T2m_WRF_vs_OBS_Sodankyla_jan_every_date_jan')


data_dict_feb = {}

j=0

for i in range(720, len(wrfout_list), 24):
    data_dict_feb['Day_%s' % str(j).zfill(2)] = wrfout_list[i:i+24]
    j+=1
    




VARIABEL = 'T2'


TIDSVEKTOR_TIMMAR = N.arange(len(data_dict_feb['Day_01'][:]))
TIDSVEKTOR_DATUM = []
wrfout_list_datum = []


k=0

for i in range(5):
    for j in range(6):
        print('New day')
        ax = axs_2[i, j]
        TEMP_C_2M_TIDSSERIE = []
        for file in data_dict_feb['Day_%s' % str(k).zfill(2)]:
            f="/home/sm_marha/WRF_JAN_FEB_KORNINGAR/ALLA_00_KORNINGAR/"+file
            Filobjekt = N4.Dataset(f, mode='r')            
            TEMP_K_2M = Filobjekt.variables[VARIABEL][:]
            TEMPC_C_2M = TEMP_K_2M[0, 563, 352]-273.15
            TEMP_C_2M_TIDSSERIE.append(TEMPC_C_2M)
        ax.plot(TIDSVEKTOR_TIMMAR,  TEMP_C_2M_TIDSSERIE, label = 'WRF')
        if k==28:
            break
        
        k+=1


plt.legend()


'''        

wrfout_list_datum = files[11:]
datum = datetime.datetime.strptime(str(file[11:]), '%Y-%m-%d_%H:%M:%S').strftime('%d/%m %H')
TIDSVEKTOR_DATUM.append(datum)

plt.xticks(TIDSVEKTOR_TIMMAR_GLES, TIDSVEKTOR_DATUM_GLES, rotation= 45, size=12)'''




fig_2.savefig('/home/sm_marha/FIGURES/T2m_WRF_vs_OBS_Sodankyla_jan_every_date_feb')
plt.show()
        

