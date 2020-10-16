import numpy as N
import scipy.io.netcdf as S
import netCDF4 as N4
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import pandas as pd
import datetime
from math import *
from matplotlib.ticker import FormatStrFormatter

Filobjekt_GRIB = N4.Dataset('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1ua201802251800_M_D.nc', mode='r') #ÄNDRA DATUM
Filobjekt_GRIB_sfc = N4.Dataset('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1sfc201802251800_M_D.nc' , mode='r')   #ÄNDRA DATUM



VV_EC =30

LATITUD = 27 # Har värdet 0 uppe 
LONGITUD = 191 # Har värdet 0 på vänstra randen. ECMWF-området är längst i x-led till skillnad från WRF-området,som är längst i y-led.

nrows = 3
ncols = 4

fig1, axs = plt.subplots(nrows, ncols, figsize=(10,8))
plt.suptitle('IFS Cloud water and cloud ice 18022518')        #ÄNDRA DATUM
fig1.subplots_adjust(hspace=0.4, wspace=0.4)

fig2, axs2 = plt.subplots(nrows, ncols, figsize=(10,8))
plt.suptitle('IFS Relative humidity 18022518')                #ÄNDRA DATUM
fig2.subplots_adjust(hspace=0.4, wspace=0.4)


TIDSSTEG = 0

for ii in range(nrows):
    for jj in range (ncols):
        
        print(TIDSSTEG)
        QICE_EC = Filobjekt_GRIB.variables['ciwc'][:]*1000
        QICE_EC_SOD = QICE_EC[TIDSSTEG, :, LATITUD, LONGITUD]
        QICE_EC_SOD = QICE_EC_SOD[::-1]

        QCLOUD_EC = Filobjekt_GRIB.variables['clwc'][:]*1000
        QCLOUD_EC_SOD = QCLOUD_EC[TIDSSTEG, :, LATITUD, LONGITUD]
        QCLOUD_EC_SOD = QCLOUD_EC_SOD[::-1]

        TEMPERATURE_EC = Filobjekt_GRIB.variables['t'][:]
        TEMPERATURE_EC_SOD = TEMPERATURE_EC[TIDSSTEG, :, LATITUD, LONGITUD]
        TEMPERATURE_EC_SOD = TEMPERATURE_EC_SOD[::-1]
        TEMPERATURE_C_EC_SOD = TEMPERATURE_EC_SOD -273.15

        '''markgeopotentialen och marktrycket från surface-gribfilen'''

        SURFACE_GEOPOTENTIAL_EC = Filobjekt_GRIB_sfc.variables['z'][:]
        SURFACE_GEOPOTENTIAL_EC_SOD = SURFACE_GEOPOTENTIAL_EC[TIDSSTEG, LATITUD, LONGITUD]

        SURFACE_PRESSURE_EC = Filobjekt_GRIB_sfc.variables['sp'][:]
        SURFACE_PRESSURE_EC_SOD = SURFACE_PRESSURE_EC[TIDSSTEG, LATITUD, LONGITUD]

        SPECIFIC_HUMIDITY_EC = Filobjekt_GRIB.variables['q'][:]
        SPECIFIC_HUMIDITY_EC_SOD = SPECIFIC_HUMIDITY_EC[TIDSSTEG, :, LATITUD, LONGITUD]
        SPECIFIC_HUMIDITY_EC_SOD = SPECIFIC_HUMIDITY_EC_SOD[::-1]

        TIMESTEPS_EC = Filobjekt_GRIB.variables['time'][:]
        VERTICAL_LEVELS_EC = Filobjekt_GRIB.variables['level'][:]



        '''BERÄKNING AV DE VERTIKALA NIVÅERNA I ECMWF'''



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




        '''Räknar ut virtuella temperaturen samt plockar ut datat för Sodankylä'''

        T_VIRTUAL_EC = (1.+0.609133*SPECIFIC_HUMIDITY_EC)*TEMPERATURE_EC

        T_VIRTUAL_EC_SOD = T_VIRTUAL_EC[TIDSSTEG,:, LATITUD, LONGITUD]


        T_VIRTUAL_EC_SOD = T_VIRTUAL_EC_SOD[::-1] #Inverterar vektorn så nivå 137 kommer först)




        '''Räknar ut de halva trycknivåerna för Sodankylä'''

        p_half = N.zeros(len(A))

        for i in range(len(A)):

            p = A[i] + B[i] * SURFACE_PRESSURE_EC[TIDSSTEG,LATITUD,LONGITUD]

            p_half[i] = p


        '''Räknar ut geopotentialen på de halva nivåerna



        for i in range(len(p_half)-1):

            Fi_half[i+1] = Fi_half[i] + Rd*(  T_VIRTUAL_SOD[i]  *  (  log(p_half[i])-log(p_half[i+1])  ) )'''
         


        '''Räknar ut hela trycknivåerna för Sodankylä'''
                 
        p_full = 0.5*(p_half[1:]+p_half[:-1])



        '''Räknar ut geopotentialen på de fulla nivåerna över Sodankylä.'''

        Rd = 287.06
        g0 = 9.80665

        #TEMP_v = p_full/(Rd*Density) #används till att kontrollera algorithmen via ecmwf.int


        Fi = N.zeros(len(p_full))
        Fi[0]=10*g0

        #FI = N.zeros(len(p_full))  #används till att kontrollera algorithmen via ecmwf.int
        #FI[0]=10*g0

        for i in range(len(p_full)-1):

            Fi[i+1] = Fi[i] + Rd*(  T_VIRTUAL_EC_SOD[i]  *  (  log(p_full[i])-log(p_full[i+1])  ) )
                                  
            #FI[i+1] = FI[i] + Rd*(  PF[i]/(Rd*Density[i])  *  (  log(PF[i])-log(PF[i+1])  )  )    #används till att kontrollera algorithmen via ecmwf.int

        '''Fulla nivåer Sodankylä'''

        FULL_LEVELS_EC_SOD = Fi/g0


        '''Beräknar potentiella temperaturen. detta görs först här för att p_full behövs på varje nivå.'''

        THETA_EC_SOD = TEMPERATURE_EC_SOD *((1000/(p_full/100))**0.286)

        '''Beräknar relativa fuktigheten EC'''

        MIXING_RATIO_EC_SOD = SPECIFIC_HUMIDITY_EC_SOD/(1-SPECIFIC_HUMIDITY_EC_SOD)

        E_MATTNAD_EC_SOD = N.exp(N.log(611.2)+(17.62*TEMPERATURE_C_EC_SOD/(243.12+TEMPERATURE_C_EC_SOD)))

        MIXING_RATIO_MATTNAD_EC_SOD = 0.622*E_MATTNAD_EC_SOD/(p_full-E_MATTNAD_EC_SOD)

        RH_EC_SOD = MIXING_RATIO_EC_SOD/MIXING_RATIO_MATTNAD_EC_SOD


        axs[ii, jj].plot(QCLOUD_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[0:VV_EC], c = 'g', label = 'QC')
        axs[ii, jj].plot(QICE_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[0:VV_EC], c = 'b', label = 'QI')
        axs[ii, jj].set_title('Time step = %dh' %(TIDSSTEG*3))
        #axs[ii, jj].xaxis.set_major_formatter(FormatStrFormatter('%8.1e'))
        axs[ii, jj].set_xlim(0, 2e-1)
        
        if ii == 2:
            axs[ii, jj].set_xlabel('QC and QI (g/kg)')
        
        if jj == 0:
            axs[ii, jj].set_ylabel('Height (m)')
            
        axs[ii, jj].grid(linestyle="dotted")
        axs[ii, jj].legend(fontsize = 'x-small')
        


        axs2[ii, jj].plot(RH_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[0:VV_EC], c = 'g', label = 'RH')
        axs2[ii, jj].set_title('Time step = %dh' %(TIDSSTEG*3))
        axs2[ii, jj].set_xlim(0, 1)

        if ii == 2:
            axs2[ii, jj].set_xlabel('RH (%)')

        if jj == 0:
            axs2[ii, jj].set_ylabel('Height (m)')
            
    
        axs2[ii, jj].grid(linestyle="dotted")
        axs2[ii, jj].legend(fontsize = 'x-small')
        
        TIDSSTEG += 1

fig1.savefig('/home/sm_marha/FIGURES/FOKUS_17_27_FEB/26_FEB/ECMWF_SOUNDINGS_3x4_QC_QI_2518')   #ÄNDRA
fig2.savefig('/home/sm_marha/FIGURES/FOKUS_17_27_FEB/26_FEB/ECMWF_SOUNDINGS_3x4_RH_2518')     #ÄNDRA



plt.show()
