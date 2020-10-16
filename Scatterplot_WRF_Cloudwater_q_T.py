import numpy as N
import scipy.io.netcdf as S
import netCDF4 as N4
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import os
import pandas as pd
from math import *
import sys


'''Detta skript beräknar molnklimatologin för WRF och ECMWF.'''





###########
# WRF OUT #
#######################################################################################


file_list = os.listdir('/nobackup/smhid12/sm_marha/2018-01-18_00z/run/')


wrfout_list = []

for file in file_list:
    if file.startswith('wrfout'):
        wrfout_list.append(file)

wrfout_list.sort()


wrfout_list = wrfout_list[1::3]   #Väljer ut vilka tidssteg jag vill ha, då alla tidssteg ligger i olika filer från WRF.


#For-loop som loopar över de olika tidssteg jag gjort ovan

data_dict = {}    #Initierar tom dictionary

i = 0

for file in wrfout_list:
    print('New file!')

    Filobjekt = S.netcdf_file('/nobackup/smhid12/sm_marha/2018-01-18_00z/run/'+ file, mode='r')
    

    PERT_THETA = Filobjekt.variables['T'][0, :, :, :]       #Plockar ut det nollte elementet från tidsdimensionen, vilket är det enda, eftersom det bara ligger 1 tid i varje fil.

    THETA = PERT_THETA+300

    P_PERT = Filobjekt.variables['P'][0, :, :, :]
    P_BASE = Filobjekt.variables['PB'][0, :, :, :]
    P= P_BASE+P_PERT

    TEMP = THETA/((1000/(P/100))**0.286)
    TEMP_C = TEMP-273.15

    QVAPOR = Filobjekt.variables['QVAPOR'][0, :, :, :]
    E_MATTNAD = N.exp(N.log(611.2)+(17.62*TEMP_C/(243.12+TEMP_C)))
    QVAPOR_MATTNAD = 0.622*E_MATTNAD/(P-E_MATTNAD)
    RH = QVAPOR/QVAPOR_MATTNAD

    
    QVAPOR = QVAPOR*1000
    QICE = Filobjekt.variables['QICE'][0, :, :, :]*1000
    QCLOUD = Filobjekt.variables['QCLOUD'][0, :, :, :]*1000

    data_dict['TEMP_C_%s' % str(i).zfill(2)] = TEMP_C   #Lagrar aktuell variabel i dictionary   ÄNDRA PÅ 2 STÄLLEN NÄR BYTE AV VARIABEL SKER
    data_dict['QCLOUD_%s' % str(i).zfill(2)] = QCLOUD   #ÄNDRA NÄR BYTE AV VARIABEL SKER

    i+=1


fig, axs = plt.subplots(3, 4, figsize=(12,8))   #Skapar tomma subplottar med 3 rader och 4 kolumner
fig.subplots_adjust(hspace=0.3, wspace=0.3)

k=0         # Iterationsvariabel för att få rätt överskrift på respektive subplot

for i in range(3):
    for j in range(4):
        ax = axs[i, j]
        ax.scatter(data_dict['TEMP_C_%s' % str(k).zfill(2)][1:20, :, :], data_dict['QCLOUD_%s' % str(k).zfill(2)][1:20, :, :] , s=0.5)   #ÄNDRA VID BYTE AV VARIABEL
        ax.set_ylim(0, 1.2)
        #ax.set_xlim(0.8, 1.02)
        t = 1+3*k
        if  i==2:
            ax.set_xlabel('Temperature (C)')        #ÄNDRA NÄR BYTE AV VARIABEL SKER
        if  j==0:  
            ax.set_ylabel('Cloud water (g/kg)')                  
        ax.set_title('Time step:%s' %t)
        k+=1

plt.suptitle('Cloudwater vs Temperature C')              #ÄNDRA NÄR BYTE AV VARIABEL SKER
fig.savefig('/home/sm_marha/FIGURES/FOKUS_18_20_JAN/18_JAN/WRF_Cloudwater_Temperature_C_clim.png')      #ÄNDRA NÄR BYTE AV VARIABEL SKER

plt.show()

