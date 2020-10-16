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
import matplotlib.ticker as ticker
import netCDF4 as N4


'''Detta skript plottar de första tidsstegen av körningen och man väljer olika variabler. Under överskriften
'skriver ut alla wrf-filer' väljer man hur många tidssteg man vill plotta. Längre ner väljer man variabel
och variabelnamn. Om det är en variabel som finns i utfilen väljer man längst ner automatiska, långa
xlabel och title. Vid självuträklnade variabler tar man de kort och får skriva in manuellt. Set.major.formatter
används när man vill ha i exponentform, men tas bort när man t ex vill plotta relativ fuiktighet eller
temperatur. I for-loopen finns också ax1.text som plottar en siffra på varje kurva. Tas bort om man vill!'''



file_list = os.listdir('/nobackup/smhid12/sm_marha/2018-02-26_00z_91_LEVELS_QC_QI_42h/run/WRFINPUT_W_TO_Wvmax_0/')
#file_list = os.listdir('/nobackup/smhid12/sm_marha/2018-02-17_18z_radt1/run/')

wrfout_list = []

for file in file_list:
    if file.startswith('wrfout'):
        wrfout_list.append(file)

wrfout_list.sort()




'''Skriver ut alla wrf-filer'''

wrfout_list = wrfout_list[5:10]

for file in wrfout_list:
    print (file)





for place, file in enumerate(wrfout_list):


    Filobjekt = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-26_00z_91_LEVELS_QC_QI_42h/run/WRFINPUT_W_TO_Wvmax_0/'+ file, mode='r')

    PERT_THETA = Filobjekt.variables['T'][0, :, 563, 352]

    THETA = PERT_THETA+300

    P_PERT = Filobjekt.variables['P'][0, :, 563, 352]
    P_BASE = Filobjekt.variables['PB'][0, :, 563, 352]
    P= P_BASE+P_PERT

    TEMP = THETA/((1000/(P/100))**0.286)
    TEMP_C = TEMP-273.15

    QVAPOR = Filobjekt.variables['QVAPOR'][0, :, 563, 352]
    E_MATTNAD = N.exp(N.log(611.2)+(17.62*TEMP_C/(243.12+TEMP_C)))
    QVAPOR_MATTNAD = 0.622*E_MATTNAD/(P-E_MATTNAD)
    RH = QVAPOR/QVAPOR_MATTNAD

    QICE = Filobjekt.variables['QICE'][0, :,  563, 352]
    QCLOUD = Filobjekt.variables['QCLOUD'][0, :,  563, 352]

    
    PH = Filobjekt.variables['PH'][0, :, 563, 352]
    PHB = Filobjekt.variables['PHB'][0, :, 563, 352]
    MODELLNIVAHOJD = (PH+PHB)/9.81
    MODELLNIVAHOJD_TER = MODELLNIVAHOJD - MODELLNIVAHOJD[0]


    MASSLEVELS = 0.5*(MODELLNIVAHOJD_TER[:-1] + MODELLNIVAHOJD_TER[1:])

    VARIABLE = RH
    VARIABLENAME = 'RH'

    fig = plt.figure(1)

    ax1 = plt.subplot(111)

    ax1.plot(VARIABLE[0:40], MASSLEVELS[0:40], label = 't=%d ' %(place))

    

    #ax1.text(VARIABLE[12], MASSLEVELS[12], str(i))


   

    
#plt.title(Filobjekt.variables[VARIABLENAME].description.decode("utf-8"))
plt.title(' Saturation mixing ratio')

#plt.xlabel(Filobjekt.variables[VARIABLENAME].description.decode("utf-8") + '    ' + Filobjekt.variables[VARIABLENAME].units.decode("utf-8"))
plt.xlabel('$W_vmax$ (g/kg)')          #('$\Theta$ (K)')

#ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%8.1e'))

ax1.grid(linestyle="dotted")

ax1.legend()

plt.ylabel('Height (m)')

#fig.savefig('/home/sm_marha/FIGURES/FOKUS_17_27_FEB/18_FEB/Initialisation_Thermodynamic_vertical_soundings_W_v_max')

plt.show()


