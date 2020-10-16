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
import Picking_values_Sounding_Dat_file

'''Detta skript jämför WRF SCM initialiserat med molnvatten med WRF 3D och ECMWF. För 3D-körningen ligger varje tidpunkt i olika wrfout-filer så
man får byta fil i sökvägen till filobjektet direkt här nedan för att byta tid. För SCM samt ECMWF ändrar man tiiden längre ned.

När man skiftar tidssteg för att ha samma tidssteg på alla körningar, leta efter #ÄNDRA. Om man ska ändra tid för sonderingsdata, ändra j-värdet
under sonderingsdata, men glöm heller inte att ändra i Reading_Dat_Files om Du vill ha 00-, eller 12z-sonderingen'''

Filobjekt_GRIB = N4.Dataset('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1ua2018021800_M_D.nc', mode='r')
Filobjekt_GRIB_sfc = N4.Dataset('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1sfc2018021800_M_D.nc' , mode='r')
Filobjekt_WRF_3D_WRFOUT = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-18_00z_91_LEVELS_QC_QI_42h/run/wrfout_d01_2018-02-18_06:00:00', mode='r')  #ÄNDRA 
Filobjekt_SCM_WRFOUT = N4.Dataset('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfout_d01_2018-02-18_00:00:00_QC_QI_T_to_Td_1e-5')



################
# ECMWF HYBRID #
#######################################################################################


TIDSSTEG = 2    #ÄNDRA
VV_EC =30

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
QICE_EC_SOD = QICE_EC[TIDSSTEG, :, LATITUD, LONGITUD]
QICE_EC_SOD = QICE_EC_SOD[::-1]

QCLOUD_EC = Filobjekt_GRIB.variables['clwc'][:]
QCLOUD_EC_SOD = QCLOUD_EC[TIDSSTEG, :, LATITUD, LONGITUD]
QCLOUD_EC_SOD = QCLOUD_EC_SOD[::-1]

QCLOUD_EC_2D = Filobjekt_GRIB.variables['clwc'][0:15, -VV_EC:, LATITUD, LONGITUD]
QCLOUD_EC_2D = N.transpose(QCLOUD_EC_2D)
QCLOUD_EC_2D = QCLOUD_EC_2D[::-1, :]

CLOUD_FRACTION_EC = Filobjekt_GRIB.variables['cc'][:]
CLOUD_FRACTION_EC_SOD = CLOUD_FRACTION_EC[TIDSSTEG, :, LATITUD, LONGITUD]
CLOUD_FRACTION_EC_SOD = CLOUD_FRACTION_EC_SOD[::-1]


TEMPERATURE_EC = Filobjekt_GRIB.variables['t'][:]
TEMPERATURE_C_EC = TEMPERATURE_EC-273.15
TEMPERATURE_EC_SOD = TEMPERATURE_EC[TIDSSTEG, :, LATITUD, LONGITUD]
TEMPERATURE_EC_SOD = TEMPERATURE_EC_SOD[::-1]
TEMPERATURE_C_EC_SOD = TEMPERATURE_EC_SOD -273.15
    
TEMPERATURE_C_EC = TEMPERATURE_EC-273.15


SPECIFIC_HUMIDITY_EC = Filobjekt_GRIB.variables['q'][:]
SPECIFIC_HUMIDITY_EC_SOD = SPECIFIC_HUMIDITY_EC[TIDSSTEG, :, LATITUD, LONGITUD]
SPECIFIC_HUMIDITY_EC_SOD = SPECIFIC_HUMIDITY_EC_SOD[::-1]

TIMESTEPS_EC = Filobjekt_GRIB.variables['time'][:]
VERTICAL_LEVELS_EC = Filobjekt_GRIB.variables['level'][:]




'''markgeopotentialen och marktrycket från surface-gribfilen'''

SURFACE_GEOPOTENTIAL_EC = Filobjekt_GRIB_sfc.variables['z'][:]
SURFACE_GEOPOTENTIAL_EC_SOD = SURFACE_GEOPOTENTIAL_EC[TIDSSTEG, LATITUD, LONGITUD]

SURFACE_PRESSURE_EC = Filobjekt_GRIB_sfc.variables['sp'][:]
SURFACE_PRESSURE_EC_SOD = SURFACE_PRESSURE_EC[TIDSSTEG, LATITUD, LONGITUD]




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





#################
# WRF 3D WRFOUT#
#######################################################################################

'''Ändra sökvägen för filobjektet i början av skriptet för att ändra tidpunkt, då varje steg i 3D ligger i olika NetCDF-filer.'''

VV_WRF = 40

PERT_THETA_WRF_3D_WRFOUT = Filobjekt_WRF_3D_WRFOUT.variables['T'][0, :, 563, 352]

P_PERT_WRF_3D_WRFOUT = Filobjekt_WRF_3D_WRFOUT.variables['P'][0, :, 563, 352]
P_BASE_WRF_3D_WRFOUT = Filobjekt_WRF_3D_WRFOUT.variables['PB'][0, :, 563, 352]
P_WRF_3D_WRFOUT = P_BASE_WRF_3D_WRFOUT+P_PERT_WRF_3D_WRFOUT

THETA_WRF_3D_WRFOUT = PERT_THETA_WRF_3D_WRFOUT+300

TEMP_WRF_3D_WRFOUT = THETA_WRF_3D_WRFOUT/((1000/(P_WRF_3D_WRFOUT/100))**0.286)
TEMP_C_WRF_3D_WRFOUT = TEMP_WRF_3D_WRFOUT-273.15
QVAPOR_WRF_3D_WRFOUT = Filobjekt_WRF_3D_WRFOUT.variables['QVAPOR'][0, :, 563, 352]
E_MATTNAD_WRF_3D_WRFOUT = N.exp(N.log(611.2)+(17.62*TEMP_C_WRF_3D_WRFOUT/(243.12+TEMP_C_WRF_3D_WRFOUT)))
QVAPOR_MATTNAD_WRF_3D_WRFOUT = 0.622*E_MATTNAD_WRF_3D_WRFOUT/(P_WRF_3D_WRFOUT-E_MATTNAD_WRF_3D_WRFOUT)
RH_WRF_3D_WRFOUT = QVAPOR_WRF_3D_WRFOUT/QVAPOR_MATTNAD_WRF_3D_WRFOUT
          
PH_WRF_3D_WRFOUT = Filobjekt_WRF_3D_WRFOUT.variables['PH'][0, :, 563, 352]
PHB_WRF_3D_WRFOUT = Filobjekt_WRF_3D_WRFOUT.variables['PHB'][0, :, 563, 352]
MODELLNIVAHOJD_WRF_3D_WRFOUT = (PH_WRF_3D_WRFOUT+PHB_WRF_3D_WRFOUT)/9.81
MODELLNIVAHOJD_TER_WRF_3D_WRFOUT = MODELLNIVAHOJD_WRF_3D_WRFOUT[:] -MODELLNIVAHOJD_WRF_3D_WRFOUT[0]


MASSLEVELS_WRF_3D_WRFOUT = 0.5*(MODELLNIVAHOJD_TER_WRF_3D_WRFOUT[:-1] + MODELLNIVAHOJD_TER_WRF_3D_WRFOUT[1:])


QCLOUD_SOD_WRF_3D_WRFOUT = Filobjekt_WRF_3D_WRFOUT.variables['QCLOUD'][0, :, 563, 352]*1000    # I g/kg
QICE_SOD_WRF_3D_WRFOUT = Filobjekt_WRF_3D_WRFOUT.variables['QICE'][0, :, 563, 352]*1000        # I g/kg




#################
# WRF SCM WRFOUT#
#######################################################################################

TIDPUNKT = 1081 #ÄNDRA   7561 Sätt hur mycket data som ska plockas ut. Kan även göras vid plottning genom att sätta variabeln t.

TIMESTEPS_SCM_WRFOUT = Filobjekt_SCM_WRFOUT.variables['XTIME'][TIDPUNKT]


PERT_THETA_SCM_WRFOUT = Filobjekt_SCM_WRFOUT.variables['T'][TIDPUNKT, :, 1, 1]

P_PERT_SCM_WRFOUT = Filobjekt_SCM_WRFOUT.variables['P'][TIDPUNKT, :, 1, 1]
P_BASE_SCM_WRFOUT = Filobjekt_SCM_WRFOUT.variables['PB'][TIDPUNKT, :, 1, 1]
P_SCM_WRFOUT = P_BASE_SCM_WRFOUT+P_PERT_SCM_WRFOUT

THETA_SCM_WRFOUT = PERT_THETA_SCM_WRFOUT+300

TEMP_SCM_WRFOUT = THETA_SCM_WRFOUT/((1000/(P_SCM_WRFOUT/100))**0.286)
TEMP_C_SCM_WRFOUT = TEMP_SCM_WRFOUT-273.15
QVAPOR_SCM_WRFOUT = Filobjekt_SCM_WRFOUT.variables['QVAPOR'][TIDPUNKT, :, 1, 1]
E_MATTNAD_SCM_WRFOUT = N.exp(N.log(611.2)+(17.62*TEMP_C_SCM_WRFOUT/(243.12+TEMP_C_SCM_WRFOUT)))
QVAPOR_MATTNAD_SCM_WRFOUT = 0.622*E_MATTNAD_SCM_WRFOUT/(P_SCM_WRFOUT-E_MATTNAD_SCM_WRFOUT)
          
PH_SCM_WRFOUT = Filobjekt_SCM_WRFOUT.variables['PH'][TIDPUNKT, :, 1, 1]
PHB_SCM_WRFOUT = Filobjekt_SCM_WRFOUT.variables['PHB'][0, :, 1, 1]
MODELLNIVAHOJD_SCM_WRFOUT = (PH_SCM_WRFOUT+PHB_SCM_WRFOUT)/9.81
MODELLNIVAHOJD_TER_SCM_WRFOUT = MODELLNIVAHOJD_SCM_WRFOUT[:] -MODELLNIVAHOJD_SCM_WRFOUT[0]


MASSLEVELS_SCM_WRFOUT = 0.5*(MODELLNIVAHOJD_TER_SCM_WRFOUT[:-1] + MODELLNIVAHOJD_TER_SCM_WRFOUT[1:])







#######################
#HÄMTAR SONDERINGSDATA#
#######################

j=48   #ÄNDRA
    
'''Hämtar sonderingsdata med hjälp av funktionen 'pick dates' i Reading_Dat_Files'''


MATRIX_2 = Picking_values_Sounding_Dat_file.pick_dates(j)


MATRIX_2 = (N.array(MATRIX_2)).astype(float)




'''Beräknar relativ fuktighet och specifik fuktighet'''



VV_SOND = 470   #2000



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

HEIGHT_TER_SONDERING = MATRIX_2[0:VV_SOND, 5]-180



###########
#PLOTTNING#
######################################################################



fig1 =plt.figure(1, figsize=(10,10))

#Rc-params måste sättas innan anrop till subplot, men gäller sedan för alla tills nytt anrop görs.
plt.rc('xtick', labelsize='x-large')
plt.rc('ytick', labelsize='x-large')

ax1 = plt.subplot(111)

ax1.grid(linestyle="dotted")

#ax1.plot(THETA_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[:VV_EC],linewidth = 1, c = 'b', label = 'IFS(Analysis)')
ax1.plot(THETA_WRF_3D_WRFOUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF], linewidth = 1, c = 'r', label = 'WRF 3D')
ax1.plot(THETA_SCM_WRFOUT[0:VV_WRF], MASSLEVELS_SCM_WRFOUT[0:VV_WRF], linewidth = 1, c = 'g', label = 'WRF SCM')
#ax1.plot(THETA_SONDERING[0:VV_SOND], HEIGHT_TER_SONDERING[0:VV_SOND], linewidth = 1, c = 'k', label = 'OBS(00z)')
ax1.scatter(THETA_SCM_WRFOUT[0:VV_WRF], MASSLEVELS_SCM_WRFOUT[0:VV_WRF], s=3, c = 'g')




ax1.set_title('Potential temperature 2018-02-18  00z', fontsize=15)      #ÄNDRA
ax1.set_ylabel('Height (m)',fontsize=15)
ax1.set_xlabel('Potential temperature (K)', fontsize=15)



ax1.legend(fontsize = 'xx-large')

#fig1.savefig('/home/sm_marha/FIGURES/FOKUS_17_27_FEB/18_FEB/In_Time_EC_WRF_OBS_00z')  #ÄNDRA




fig2 =plt.figure(2, figsize=(10,10))

ax1 = plt.subplot(111)

ax1.grid(linestyle="dotted")

ax1.plot(QCLOUD_SOD_WRF_3D_WRFOUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF], linewidth = 1, c = 'g', label = 'QC WRF 3D')
ax1.plot(QICE_SOD_WRF_3D_WRFOUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF], linewidth = 1, c = 'b', label = 'QI WRF 3D')


ax1.set_title('Qcloud and Qice 2018-02-18  00z', fontsize=15)      #ÄNDRA
ax1.set_ylabel('Height (m)',fontsize=15)
ax1.set_xlabel('Cloud water and Cloud ice (g/kg)', fontsize=15)       #ÄNDRA 
ax1.set_ylim(0, 3000)



fig3 =plt.figure(3, figsize=(10,10))

ax1 = plt.subplot(111)

ax1.grid(linestyle="dotted")

ax1.plot(RH_WRF_3D_WRFOUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF], linewidth = 1, c = 'g', label = 'RH WRF 3D')



ax1.set_title('Relative humidity 2018-02-18  00z', fontsize=15)      #ÄNDRA
ax1.set_ylabel('Height (m)',fontsize=15)
ax1.set_xlabel('Relative humidity (%)', fontsize=15)       #ÄNDRA 
ax1.set_ylim(0, 3000)

plt.show()









