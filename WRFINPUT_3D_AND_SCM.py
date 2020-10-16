import numpy as N
import scipy.io.netcdf as S
import netCDF4 as N4
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import pandas as pd
import datetime
from math import *
import matplotlib.ticker as ticker
import Picking_values_Sounding_Dat_file_version_2



Filobjekt_WRF_3D_WRF_INPUT = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-18_00z_91_LEVELS_QC_QI_42h/run/wrfinput_d01')           #3D
Filobjekt_WRF_SCM_WRF_INPUT = N4.Dataset('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfinput_d01_0218_QC_QI_REFERENCE')     #SINGLE COLUMN
Filobjekt_GRIB = N4.Dataset('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1ua2018021800_M_D.nc', mode='r')         
Filobjekt_GRIB_sfc = N4.Dataset('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1sfc2018021800_M_D.nc' , mode='r')


##################
# WRF 3D wrfinput#
#######################################################################################


VV_WRF = 33



QICE_WRF_3D_WRFINPUT = Filobjekt_WRF_3D_WRF_INPUT.variables['QICE'][0, :, 563, 352]
QCLOUD_WRF_3D_WRFINPUT = Filobjekt_WRF_3D_WRF_INPUT.variables['QCLOUD'][0, :, 563, 352]

#print(QICE_WRF_3D_WRFOUT-QICE_WRF_3D_WRFINPUT)
print(QICE_WRF_3D_WRFINPUT)
#print(QCLOUD_WRF_3D_WRFOUT-QCLOUD_WRF_3D_WRFINPUT)

PERT_THETA_WRF_3D_WRFINPUT = Filobjekt_WRF_3D_WRF_INPUT.variables['T'][0, :, 563, 352]

P_PERT_WRF_3D_WRFINPUT = Filobjekt_WRF_3D_WRF_INPUT.variables['P'][0, :, 563, 352]
P_BASE_WRF_3D_WRFINPUT = Filobjekt_WRF_3D_WRF_INPUT.variables['PB'][0, :, 563, 352]
P_WRF_3D_WRFINPUT = P_BASE_WRF_3D_WRFINPUT+P_PERT_WRF_3D_WRFINPUT

THETA_WRF_3D_WRFINPUT = PERT_THETA_WRF_3D_WRFINPUT+300
TEMP_WRF_3D_WRFINPUT = THETA_WRF_3D_WRFINPUT/((1000/(P_WRF_3D_WRFINPUT/100))**0.286)
TEMP_C_WRF_3D_WRFINPUT = TEMP_WRF_3D_WRFINPUT-273.15
QVAPOR_WRF_3D_WRFINPUT = Filobjekt_WRF_3D_WRF_INPUT.variables['QVAPOR'][0, :, 563, 352]

E_MATTNAD_WRF_3D_WRFINPUT = N.exp(N.log(611.2)+(17.62*TEMP_C_WRF_3D_WRFINPUT/(243.12+TEMP_C_WRF_3D_WRFINPUT)))
QVAPOR_MATTNAD_WRF_3D_WRFINPUT = 0.622*E_MATTNAD_WRF_3D_WRFINPUT/(P_WRF_3D_WRFINPUT-E_MATTNAD_WRF_3D_WRFINPUT)
RH_WRF_3D_WRFINPUT = QVAPOR_WRF_3D_WRFINPUT/QVAPOR_MATTNAD_WRF_3D_WRFINPUT
T_DAGGPUNKT_WRF_3D_WRFINPUT = (243.12*N.log(611.2)-243.12*N.log(RH_WRF_3D_WRFINPUT*E_MATTNAD_WRF_3D_WRFINPUT)) / (N.log(RH_WRF_3D_WRFINPUT*E_MATTNAD_WRF_3D_WRFINPUT)-17.62-N.log(611.2))
#T_DAGGPUNKT_WRF_3D_WRFINPUT = TEMP_C_WRF_3D_WRFINPUT - ((100-RH_WRF_3D_WRFINPUT*100)/5.)    #Enklare formel för Td
          
PH_WRF_3D_WRFINPUT = Filobjekt_WRF_3D_WRF_INPUT.variables['PH'][0, :, 563, 352]
PHB_WRF_3D_WRFINPUT = Filobjekt_WRF_3D_WRF_INPUT.variables['PHB'][0, :, 563, 352]
MODELLNIVAHOJD_WRF_3D_WRFINPUT = (PH_WRF_3D_WRFINPUT+PHB_WRF_3D_WRFINPUT)/9.81
MODELLNIVAHOJD_TER_WRF_3D_WRFINPUT = MODELLNIVAHOJD_WRF_3D_WRFINPUT[:] -MODELLNIVAHOJD_WRF_3D_WRFINPUT[0]


MASSLEVELS_WRF_3D_WRFINPUT = 0.5*(MODELLNIVAHOJD_TER_WRF_3D_WRFINPUT[:-1] + MODELLNIVAHOJD_TER_WRF_3D_WRFINPUT[1:])



###################
# WRF SCM wrfinput#
#######################################################################################


QICE_WRF_SCM_WRFINPUT = Filobjekt_WRF_SCM_WRF_INPUT.variables['QICE'][0, :, 1, 1]*1000
QCLOUD_WRF_SCM_WRFINPUT = Filobjekt_WRF_SCM_WRF_INPUT.variables['QCLOUD'][0, :, 1, 1]*1000

#print(QICE_WRF_SCM_WRFOUT-QICE_WRF_SCM_WRFINPUT)
print(QICE_WRF_SCM_WRFINPUT)
#print(QCLOUD_WRF_SCM_WRFOUT-QCLOUD_WRF_SCM_WRFINPUT)

PERT_THETA_WRF_SCM_WRFINPUT = Filobjekt_WRF_SCM_WRF_INPUT.variables['T'][0, :, 1, 1]

P_PERT_WRF_SCM_WRFINPUT = Filobjekt_WRF_SCM_WRF_INPUT.variables['P'][0, :, 1, 1]
P_BASE_WRF_SCM_WRFINPUT = Filobjekt_WRF_SCM_WRF_INPUT.variables['PB'][0, :, 1, 1]
P_WRF_SCM_WRFINPUT = P_BASE_WRF_SCM_WRFINPUT+P_PERT_WRF_SCM_WRFINPUT

THETA_WRF_SCM_WRFINPUT = PERT_THETA_WRF_SCM_WRFINPUT+300
TEMP_WRF_SCM_WRFINPUT = THETA_WRF_SCM_WRFINPUT/((1000/(P_WRF_SCM_WRFINPUT/100))**0.286)
TEMP_C_WRF_SCM_WRFINPUT = TEMP_WRF_SCM_WRFINPUT-273.15
QVAPOR_WRF_SCM_WRFINPUT = Filobjekt_WRF_SCM_WRF_INPUT.variables['QVAPOR'][0, :, 1, 1]*1000

E_MATTNAD_WRF_SCM_WRFINPUT = N.exp(N.log(611.2)+(17.62*TEMP_C_WRF_SCM_WRFINPUT/(243.12+TEMP_C_WRF_SCM_WRFINPUT)))
QVAPOR_MATTNAD_WRF_SCM_WRFINPUT = 0.622*E_MATTNAD_WRF_SCM_WRFINPUT/(P_WRF_SCM_WRFINPUT-E_MATTNAD_WRF_SCM_WRFINPUT)*1000
RH_WRF_SCM_WRFINPUT = QVAPOR_WRF_SCM_WRFINPUT/QVAPOR_MATTNAD_WRF_SCM_WRFINPUT
T_DAGGPUNKT_WRF_SCM_WRFINPUT = (243.12*N.log(611.2)-243.12*N.log(RH_WRF_SCM_WRFINPUT*E_MATTNAD_WRF_SCM_WRFINPUT)) / (N.log(RH_WRF_SCM_WRFINPUT*E_MATTNAD_WRF_SCM_WRFINPUT)-17.62-N.log(611.2))
#T_DAGGPUNKT_WRF_3D_WRFINPUT = TEMP_C_WRF_3D_WRFINPUT - ((100-RH_WRF_3D_WRFINPUT*100)/5.)    #Enklare formel för Td
          
PH_WRF_SCM_WRFINPUT = Filobjekt_WRF_SCM_WRF_INPUT.variables['PH'][0, :, 1, 1]
PHB_WRF_SCM_WRFINPUT = Filobjekt_WRF_SCM_WRF_INPUT.variables['PHB'][0, :, 1, 1]
MODELLNIVAHOJD_WRF_SCM_WRFINPUT = (PH_WRF_SCM_WRFINPUT+PHB_WRF_SCM_WRFINPUT)/9.81
MODELLNIVAHOJD_TER_WRF_SCM_WRFINPUT = MODELLNIVAHOJD_WRF_SCM_WRFINPUT[:] -MODELLNIVAHOJD_WRF_SCM_WRFINPUT[0]


MASSLEVELS_WRF_SCM_WRFINPUT = 0.5*(MODELLNIVAHOJD_TER_WRF_SCM_WRFINPUT[:-1] + MODELLNIVAHOJD_TER_WRF_SCM_WRFINPUT[1:])



##########
#PLOTTING#
#######################################################################################




fig1 = plt.figure(1, figsize=(10,8)) 

ax1 = plt.subplot(111)

#ax1.plot(TEMP_C_WRF_3D_WRFINPUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFINPUT[0:VV_WRF], linewidth = 1, c = 'g', label = '3D T')
#ax1.plot(T_DAGGPUNKT_WRF_3D_WRFINPUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFINPUT[0:VV_WRF], linewidth = 1, c = 'g', linestyle = 'dashed', label = '3D DP')

ax1.plot(TEMP_C_WRF_SCM_WRFINPUT[0:VV_WRF], MASSLEVELS_WRF_SCM_WRFINPUT[0:VV_WRF], linewidth = 1, c = 'b', label = 'SCM T')
ax1.scatter(TEMP_C_WRF_SCM_WRFINPUT[0:VV_WRF], MASSLEVELS_WRF_SCM_WRFINPUT[0:VV_WRF], c = 'b', s = 3)
ax1.plot(T_DAGGPUNKT_WRF_SCM_WRFINPUT[0:VV_WRF], MASSLEVELS_WRF_SCM_WRFINPUT[0:VV_WRF], linewidth = 1, c = 'b', linestyle = 'dashed', label = 'SCM DP')

ax1.set_title('Temperature and Dew point WRF SCM', fontsize=18)      
ax1.set_ylabel('Height (m)',fontsize=15)
ax1.set_xlabel('Temperature (C)', fontsize=15)
ax1.set_xlim(-50, 0)
ax1.grid(linestyle="dotted")
ax1.legend()



fig2 = plt.figure(2, figsize=(10,8))

ax2 = plt.subplot(111)

#ax2.plot(THETA_WRF_3D_WRFINPUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFINPUT[0:VV_WRF], linewidth = 1, c = 'g', label = '3D Theta')


ax2.plot(THETA_WRF_SCM_WRFINPUT[0:VV_WRF], MASSLEVELS_WRF_SCM_WRFINPUT[0:VV_WRF], linewidth = 1, c = 'b', label = 'SCM Theta')



ax2.set_title('Potential temperature WRF SCM', fontsize=18)      
ax2.set_ylabel('Height (m)',fontsize=15)
ax2.set_xlabel('Potential temperature (C)', fontsize=15)
#ax2.set_xlim(-50, 0)
ax2.grid(linestyle="dotted")
ax2.legend(fontsize = 'x-large')




#######################
#HÄMTAR SONDERINGSDATA#
#######################################################################

j=48   #ÄNDRA (18 FEB)          #ÄNDRA





'''Hämtar sonderingsdata med hjälp av funktionen 'pick dates' i Reading_Dat_Files'''


MATRIX_2, datum = Picking_values_Sounding_Dat_file_version_2.pick_dates(j, '00')


MATRIX_2 = (N.array(MATRIX_2)).astype(float)




'''Beräknar relativ fuktighet och specifik fuktighet'''



VV_SOND = 400   #2000



TEMP_C_SONDERING = MATRIX_2[0:VV_SOND, 4]

TEMP_K_SONDERING = TEMP_C_SONDERING + 273.15

DAGGPUNKT_SONDERING = MATRIX_2[0:VV_SOND, 7]

P_TEMPERATURE_SONDERING = MATRIX_2[0:VV_SOND, 6] * 100  #Gör om från hPa till Pa

P_DEWPOINT_SONDERING = MATRIX_2[0:VV_SOND, 9] * 100  #Gör om från hPa till Pa

E_MATTNAD_SONDERING = N.exp(N.log(611.2)+(17.62*TEMP_C_SONDERING/(243.12+TEMP_C_SONDERING)))

QVAPOR_MATTNAD_SONDERING = 0.622*E_MATTNAD_SONDERING/(P_TEMPERATURE_SONDERING-E_MATTNAD_SONDERING)

E_AKTUELL_SONDERING = N.exp(N.log(611.2)+(17.62*DAGGPUNKT_SONDERING/(243.12+DAGGPUNKT_SONDERING)))

RH_SONDERING = E_AKTUELL_SONDERING/E_MATTNAD_SONDERING

QVAPOR_SONDERING = 0.622*E_AKTUELL_SONDERING/(P_TEMPERATURE_SONDERING-E_AKTUELL_SONDERING)

THETA_SONDERING = TEMP_K_SONDERING *((1000/(P_TEMPERATURE_SONDERING/100))**0.286)

HEIGHT_TEMP_TER_SONDERING = MATRIX_2[0:VV_SOND, 5]-180



##########
#PLOTTING#
#######################################################################################


ax1.plot(TEMP_C_SONDERING[0:VV_SOND], HEIGHT_TEMP_TER_SONDERING[0:VV_SOND], linewidth = 1, c = 'k', label = 'OBS T')
ax1.plot(DAGGPUNKT_SONDERING[0:VV_SOND], HEIGHT_TEMP_TER_SONDERING[0:VV_SOND], linewidth = 1, c = 'k', linestyle = 'dashed', label = 'OBS DP')

ax1.legend(fontsize = 'x-large')

ax1.set_xlim(-20, -5)





################
# ECMWF HYBRID #
#######################################################################################


VV_EC = 30
ECMWF_TIMESTEP = 0


LATITUD = 27 # Har värdet 0 uppe 
LONGITUD = 191 # Har värdet 0 på vänstra randen. ECMWF-området är längst i x-led till skillnad från WRF-området,som är längst i y-led.

TIMESTEPS = Filobjekt_GRIB.variables['time'][:]
time = Filobjekt_GRIB.variables['time']
dates = N4.num2date(TIMESTEPS, time.units, time.calendar)
ALLA_DATUM=[]

for date in dates:
    datum = date.strftime('%d/%m %Hz')
    ALLA_DATUM.append(datum)


LONGITUDE = Filobjekt_GRIB.variables['longitude'][:]
LATITUDE = Filobjekt_GRIB.variables['latitude'][:]

lon_matrix, lat_matrix = N.meshgrid(LONGITUDE, LATITUDE)

QICE_EC = Filobjekt_GRIB.variables['ciwc'][:]
QICE_EC_SOD = QICE_EC[ECMWF_TIMESTEP, :, LATITUD, LONGITUD]
QICE_EC_SOD = QICE_EC_SOD[::-1]

QCLOUD_EC = Filobjekt_GRIB.variables['clwc'][:]*1000
QCLOUD_EC_SOD = QCLOUD_EC[ECMWF_TIMESTEP, :, LATITUD, LONGITUD]
QCLOUD_EC_SOD = QCLOUD_EC_SOD[::-1]


CLOUD_FRACTION_EC = Filobjekt_GRIB.variables['cc'][:]
CLOUD_FRACTION_EC_SOD = CLOUD_FRACTION_EC[ECMWF_TIMESTEP, :, LATITUD, LONGITUD]
CLOUD_FRACTION_EC_SOD = CLOUD_FRACTION_EC_SOD[::-1]


TEMPERATURE_EC = Filobjekt_GRIB.variables['t'][:]
TEMPERATURE_EC_SOD = TEMPERATURE_EC[ECMWF_TIMESTEP, :, LATITUD, LONGITUD]
TEMPERATURE_EC_SOD = TEMPERATURE_EC_SOD[::-1]
TEMPERATURE_C_EC_SOD = TEMPERATURE_EC_SOD -273.15
E_MATTNAD_EC_SOD = N.exp(N.log(611.2)+(17.62*TEMPERATURE_C_EC_SOD/(243.12+TEMPERATURE_C_EC_SOD)))




    
TEMPERATURE_C_EC = TEMPERATURE_EC-273.15


SPECIFIC_HUMIDITY_EC = Filobjekt_GRIB.variables['q'][:]
SPECIFIC_HUMIDITY_EC_SOD = SPECIFIC_HUMIDITY_EC[ECMWF_TIMESTEP, :, LATITUD, LONGITUD]
SPECIFIC_HUMIDITY_EC_SOD = SPECIFIC_HUMIDITY_EC_SOD[::-1]


TIMESTEPS_EC = Filobjekt_GRIB.variables['time'][:]
VERTICAL_LEVELS_EC = Filobjekt_GRIB.variables['level'][:]




'''markgeopotentialen och marktrycket från surface-gribfilen'''

SURFACE_GEOPOTENTIAL_EC = Filobjekt_GRIB_sfc.variables['z'][:]
SURFACE_GEOPOTENTIAL_EC_SOD = SURFACE_GEOPOTENTIAL_EC[ECMWF_TIMESTEP, LATITUD, LONGITUD]

SURFACE_PRESSURE_EC = Filobjekt_GRIB_sfc.variables['sp'][:]
SURFACE_PRESSURE_EC_SOD = SURFACE_PRESSURE_EC[ECMWF_TIMESTEP, LATITUD, LONGITUD]







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

T_VIRTUAL_EC_SOD = T_VIRTUAL_EC[ECMWF_TIMESTEP,:, LATITUD, LONGITUD]


T_VIRTUAL_EC_SOD = T_VIRTUAL_EC_SOD[::-1] #Inverterar vektorn så nivå 137 kommer först)




'''Räknar ut de halva trycknivåerna för Sodankylä'''

p_half = N.zeros(len(A))

for i in range(len(A)):

    p = A[i] + B[i] * SURFACE_PRESSURE_EC[ECMWF_TIMESTEP,LATITUD,LONGITUD]

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

QVAPOR_MATTNAD_EC_SOD = 0.622*E_MATTNAD_EC_SOD/(p_full-E_MATTNAD_EC_SOD)

QVAPOR_EC_SOD = SPECIFIC_HUMIDITY_EC_SOD/(1- SPECIFIC_HUMIDITY_EC_SOD)

RH_EC_SOD =  QVAPOR_EC_SOD/QVAPOR_MATTNAD_EC_SOD

T_DAGGPUNKT_EC_SOD = (243.12*N.log(611.2)-243.12*N.log(RH_EC_SOD*E_MATTNAD_EC_SOD)) / (N.log(RH_EC_SOD*E_MATTNAD_EC_SOD)-17.62-N.log(611.2))



##########
#PLOTTING#
#######################################################################################

fig3 = plt.figure(3, figsize=(10,8)) 

ax3 = plt.subplot(111)


ax3.plot(QCLOUD_WRF_SCM_WRFINPUT[0:VV_WRF], MASSLEVELS_WRF_SCM_WRFINPUT[0:VV_WRF], linewidth = 1, c = 'b', label = 'SCM QV')
ax3.scatter(QCLOUD_WRF_SCM_WRFINPUT[0:VV_WRF], MASSLEVELS_WRF_SCM_WRFINPUT[0:VV_WRF], s = 3, c = 'b')
ax3.plot(QCLOUD_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[0:VV_EC], linewidth = 1, c = 'r', label = 'IFS QV')
ax3.scatter(QCLOUD_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[0:VV_EC], s= 3, c = 'r')




ax3.set_title('Cloud water IFS vs WRF SCM wrfinput', fontsize=18)      
ax3.set_ylabel('Height (m)',fontsize=15)
ax3.set_xlabel('Cloud water (g/kg)', fontsize=15)
#ax3.set_xlim(-50, 0)
ax3.grid(linestyle="dotted")
ax3.legend(fontsize = 'x-large')




plt.show()
