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
import Picking_values_Sounding_Dat_file_version_2



'''Detta skript jämför WRF SCM initialiserat med molnvatten med WRF 3D och ECMWF. För 3D-körningen ligger varje tidpunkt i olika wrfout-filer så
man får byta fil i sökvägen till filobjektet direkt här nedan för att byta tid. För SCM samt ECMWF ändrar man tiden längre ned.

När man skiftar tidssteg för att ha samma tidssteg på alla körningar, leta efter #ÄNDRA. Om man ska ändra tid för sonderingsdata, ändra j-värdet
under sonderingsdata, men glöm heller inte att ändra i Reading_Dat_Files om Du vill ha 00-, eller 12z-sonderingen'''






Filobjekt_GRIB = N4.Dataset('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1ua2018021800_M_D.nc', mode='r')
Filobjekt_GRIB_sfc = N4.Dataset('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1sfc2018021800_M_D.nc' , mode='r')
Filobjekt_Coldstart_WRF = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-18_00z_91_LEVELS_QC_QI/run/wrfout_d01_2018-02-18_00:00:00', mode='r')
Filobjekt_SCM = N4.Dataset('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfinput_d01')
Filobjekt_SCM_WRFOUT = N4.Dataset('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfout_d01_2018-02-18_00:00:00')
#Filobjekt_Coldstart = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-18_00z_91_LEVELS_QC_QI_42h_T_to_Td_1e-4/z/run/wrfout_d01_2018-02-18_00:00:00', mode='r')


################
# ECMWF HYBRID #
#######################################################################################


TIDSSTEG = 0    #ÄNDRA
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
QICE_EC_SOD = QICE_EC_SOD[::-1]*1000

QCLOUD_EC = Filobjekt_GRIB.variables['clwc'][:]
QCLOUD_EC_SOD = QCLOUD_EC[TIDSSTEG, :, LATITUD, LONGITUD]
QCLOUD_EC_SOD = QCLOUD_EC_SOD[::-1]*1000

QCLOUD_EC_2D = Filobjekt_GRIB.variables['clwc'][0:15, -VV_EC:, LATITUD, LONGITUD]
QCLOUD_EC_2D = N.transpose(QCLOUD_EC_2D)
QCLOUD_EC_2D = QCLOUD_EC_2D[::-1, :]

CLOUD_FRACTION_EC = Filobjekt_GRIB.variables['cc'][:]
CLOUD_FRACTION_EC_SOD = CLOUD_FRACTION_EC[TIDSSTEG, :, LATITUD, LONGITUD]
CLOUD_FRACTION_EC_SOD = CLOUD_FRACTION_EC_SOD[::-1]


TEMPERATURE_EC = Filobjekt_GRIB.variables['t'][:]
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


'''Beräknar potentiella temperaturen. detta görs först här för att p_full behövs på varje nivå. Detsamma gäller Mättnadsblandningsförhållandet'''

THETA_EC_SOD = TEMPERATURE_EC_SOD *((1000/(p_full/100))**0.286)


MIXING_RATIO_EC_SOD = SPECIFIC_HUMIDITY_EC_SOD/(1-SPECIFIC_HUMIDITY_EC_SOD)
E_MATTNAD_Water_EC_SOD = N.exp(N.log(611.2)+(17.62*TEMPERATURE_C_EC_SOD/(243.12+TEMPERATURE_C_EC_SOD)))
MIXING_RATIO_MATTNAD_Water_EC_SOD = 0.622*E_MATTNAD_Water_EC_SOD/(p_full-E_MATTNAD_Water_EC_SOD)

RH_EC_Water_SOD = MIXING_RATIO_EC_SOD/MIXING_RATIO_MATTNAD_Water_EC_SOD






####################
# WRF 3D COLDSTART #
#######################################################################################
LAT = 563
LON = 352


ETANIVAER_MASS_CS = Filobjekt_Coldstart_WRF.variables['ZNU'][:]
ETANIVAER_FULL_CS = Filobjekt_Coldstart_WRF.variables['ZNW'][:]
ZETATOP_CS = Filobjekt_Coldstart_WRF.variables['ZETATOP'][:]
HGT_CS = Filobjekt_Coldstart_WRF.variables['HGT'][:]
PH_CS = Filobjekt_Coldstart_WRF.variables['PH'][:]
PHB_CS = Filobjekt_Coldstart_WRF.variables['PHB'][:]
GEOPOTENTIAL_CS = PH_CS + PHB_CS
QCLOUD_CS = Filobjekt_Coldstart_WRF.variables['QCLOUD'][:]
QVAPOR_CS = Filobjekt_Coldstart_WRF.variables['QVAPOR'][:]
MODELLNIVAHOJD_CS = (PH_CS+PHB_CS)/9.81
MODELLNIVAHOJD_TER = MODELLNIVAHOJD_CS[:] -MODELLNIVAHOJD_CS[0, 0, :, :]
MASSLEVELS_TER_SOD = 0.5*(MODELLNIVAHOJD_TER[0,:-1, LAT,  LON] + MODELLNIVAHOJD_TER[0, 1:, LAT,  LON])
W_CS = Filobjekt_Coldstart_WRF.variables['W'][:]
'''PERT_THETA är den variabeln som i detta fall är vald längst upp. I många andra fall behövs inga omräkningar av variabeln utan
VARIABELVARDE är den viktiga variabeln genom hela skriptet.'''
PERT_THETA_CS = Filobjekt_Coldstart_WRF.variables['T'][:]
THETA_CS = PERT_THETA_CS+300
THETA_CS_SOD = THETA_CS[0, :, LAT,  LON]
P_PERT_CS = Filobjekt_Coldstart_WRF.variables['P'][:]
P_BASE_CS = Filobjekt_Coldstart_WRF.variables['PB'][:]
P_CS = P_BASE_CS+P_PERT_CS
P_CS_SOD = P_CS[0, :, LAT,  LON]
T2_CS = Filobjekt_Coldstart_WRF.variables['T2'][:]
TSK_CS = Filobjekt_Coldstart_WRF.variables['TSK'][:]
SNOWH_CS = Filobjekt_Coldstart_WRF.variables['SNOWH'][:]
TEMP_CS = THETA_CS/((1000/(P_CS/100))**0.286)
TEMP_C_CS = TEMP_CS-273.15
TEMP_C_SOD_CS = TEMP_C_CS[0, :, LAT,  LON]
TEMP_2M_C_SOD_CS = T2_CS[0, LAT,  LON]-273.15
TSK_C_SOD_CS = TSK_CS[0, LAT,  LON]-273.15
QVAPOR_SOD_CS = Filobjekt_Coldstart_WRF.variables['QVAPOR'][0, :, LAT,  LON]     
SNOW_DEPTH_SOD_CS = SNOWH_CS[0, LAT,  LON]

E_MATTNAD_Water_CS_SOD = N.exp(N.log(611.2)+(17.62*TEMP_C_SOD_CS/(243.12+TEMP_C_SOD_CS)))
QVAPOR_MATTNAD_Water_CS_SOD = 0.622*E_MATTNAD_Water_CS_SOD/(P_CS_SOD-E_MATTNAD_Water_CS_SOD)

RH_CS_Water_SOD = QVAPOR_SOD_CS/QVAPOR_MATTNAD_Water_CS_SOD

QCLOUD_SOD_CS = Filobjekt_Coldstart_WRF.variables['QCLOUD'][0, :, LAT,  LON]*1000    # I g/kg
QICE_SOD_CS = Filobjekt_Coldstart_WRF.variables['QICE'][0, :, LAT, LON]*1000        # I g/kg


####################
# WRF SCM WRFINPUT #
#######################################################################################

#TIMESTEPS_SCM = Filobjekt_SCM.variables['XTIME'][0:TIDSLANGD]

VV_WRF = 34


PERT_THETA_SCM = Filobjekt_SCM.variables['T'][0, :, 1, 1]

P_PERT_SCM = Filobjekt_SCM.variables['P'][0, :, 1, 1]
P_BASE_SCM = Filobjekt_SCM.variables['PB'][0, :, 1, 1]
P_SCM = P_BASE_SCM+P_PERT_SCM

THETA_SCM = PERT_THETA_SCM+300

TEMP_SCM = THETA_SCM/((1000/(P_SCM/100))**0.286)
TEMP_C_SCM = TEMP_SCM-273.15
QVAPOR_SCM = Filobjekt_SCM.variables['QVAPOR'][0, :, 1, 1]
E_MATTNAD_SCM = N.exp(N.log(611.2)+(17.62*TEMP_C_SCM/(243.12+TEMP_C_SCM)))
QVAPOR_MATTNAD_SCM = 0.622*E_MATTNAD_SCM/(P_SCM-E_MATTNAD_SCM)
QCLOUD_SCM = Filobjekt_SCM.variables['QCLOUD'][1, :, 1,  1]*1000
          
PH_SCM = Filobjekt_SCM.variables['PH'][0, :, 1, 1]
PHB_SCM = Filobjekt_SCM.variables['PHB'][0, :, 1, 1]
MODELLNIVAHOJD_SCM = (PH_SCM+PHB_SCM)/9.81
MODELLNIVAHOJD_TER_SCM = MODELLNIVAHOJD_SCM[:] -MODELLNIVAHOJD_SCM[0]


MASSLEVELS_SCM = 0.5*(MODELLNIVAHOJD_TER_SCM[:-1] + MODELLNIVAHOJD_TER_SCM[1:])

'''
#################
# WRF SCM WRFOUT#
#######################################################################################

TIDPUNKT = 7561 #7561 Sätt hur mycket data som ska plockas ut. Kan även göras vid plottning genom att sätta variabeln t.

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

'''





#######################
#HÄMTAR SONDERINGSDATA#
#######################

j=57   #ÄNDRA
    
'''Hämtar sonderingsdata med hjälp av funktionen 'pick dates' i Reading_Dat_Files'''


MATRIX_2, datum = Picking_values_Sounding_Dat_file_version_2.pick_dates(j, '00')


MATRIX_2 = (N.array(MATRIX_2)).astype(float)




'''Beräknar relativ fuktighet och specifik fuktighet'''



VV_SOND = 420   #2000




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


fig1 = plt.figure(1, figsize=(15,10))

#rc-parametrs måste anropas före anrop till subplot och gäller då för alla subplottar

plt.rc('xtick', labelsize='x-large')   #When using rc, xtick is lika xticklabel, i.e setting xticklabels bigger
plt.rc('ytick', labelsize='x-large')

ax1 = plt.subplot(121)

ax1.grid(linestyle="dotted")

ax1.plot(THETA_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[:VV_EC],linewidth = 1, c = 'b', label = 'EC(00z)')
ax1.plot(THETA_CS_SOD[0:VV_WRF], MASSLEVELS_TER_SOD[0:VV_WRF], linewidth = 1, c = 'r', label = 'WRF 3D(00z)')
ax1.plot(THETA_SONDERING[0:VV_SOND], HEIGHT_TER_SONDERING[0:VV_SOND], linewidth = 1, c = 'k', label = 'OBS(00z)')
ax1.plot(THETA_SCM[0:VV_WRF], MASSLEVELS_SCM[0:VV_WRF], linewidth = 1, c = 'g', label = 'WRF SCM(00z)')
ax1.scatter(THETA_SCM[0:VV_WRF], MASSLEVELS_SCM[0:VV_WRF], s=3, c = 'g')

ax1.set_title('Potential temperature 20180118 00z', fontsize=15)   #ÄNDRA
ax1.set_ylabel('Height (m)', fontsize=15)
ax1.set_xlabel('Potential temperature (K)', fontsize=15)
ax1.legend(fontsize = 'large')



ax2 = plt.subplot(122)

ax2.grid(linestyle="dotted")

ax2.plot(TEMPERATURE_C_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[:VV_EC],linewidth = 1, c = 'b', label = 'EC(00z)')
ax2.plot(TEMP_C_SOD_CS[0:VV_WRF], MASSLEVELS_TER_SOD[0:VV_WRF], linewidth = 1, c = 'r', label = 'WRF 3D(00z)')
ax2.plot(TEMP_C_SONDERING[0:VV_SOND], HEIGHT_TER_SONDERING[0:VV_SOND], linewidth = 1, c = 'k', label = 'OBS(00z)')
ax2.plot(TEMP_C_SCM[0:VV_WRF], MASSLEVELS_SCM[0:VV_WRF], linewidth = 1, c = 'g', label = 'WRF SCM(00z)')
ax2.scatter(TEMP_C_SCM[0:VV_WRF], MASSLEVELS_SCM[0:VV_WRF], s=3, c = 'g')


ax2.set_title('Temperature 20180118 00z', fontsize=15)  #ÄNDRA
ax2.set_ylabel('Height (m)', fontsize=15)
ax2.set_xlabel('Temperature (C)', fontsize=15)
ax2.legend(fontsize = 'large')


fig2 = plt.figure(2, figsize=(15,10))

ax1 = plt.subplot(121)

ax1.grid(linestyle="dotted")

ax1.plot(QVAPOR_SOD_CS[0:VV_WRF], MASSLEVELS_TER_SOD[0:VV_WRF], linewidth = 1, c = 'r', label = 'WRF 3D(00z)')
ax1.set_title('Mixing ratio and  20180227 00z', fontsize=15)   #ÄNDRA
ax1.set_ylabel('Height (m)', fontsize=15)
ax1.set_xlabel('Mixing ratio (g/kg)', fontsize=15)
ax1.legend(fontsize = 'large')


ax2 = plt.subplot(122)

ax2.grid(linestyle="dotted")

ax2.plot(RH_CS_Water_SOD[0:VV_WRF], MASSLEVELS_TER_SOD[0:VV_WRF], linewidth = 1, c = 'g', label = 'WRF 3D(00z)')
ax2.set_title('Relative humidity  20180118 00z', fontsize=15)   #ÄNDRA
ax2.set_ylabel('Height (m)', fontsize=15)
ax2.set_xlabel('Relative humidity (%)', fontsize=15)
ax2.legend(fontsize = 'large')



fig3 = plt.figure(3, figsize=(7, 9))

ax1 = plt.subplot(111)

ax1.grid(linestyle="dotted")

#ax1.plot(QCLOUD_SCM[0:VV_WRF], MASSLEVELS_SCM[0:VV_WRF], linewidth = 1, c = 'g', label = 'QC(00z)')
ax1.plot(QCLOUD_SOD_CS[0:VV_WRF], MASSLEVELS_TER_SOD[0:VV_WRF], linewidth = 1, c = 'g', label = 'QC(00z)')
ax1.scatter(QCLOUD_SOD_CS[0:VV_WRF], MASSLEVELS_TER_SOD[0:VV_WRF], c = 'g')
ax1.plot(QICE_SOD_CS[0:VV_WRF], MASSLEVELS_TER_SOD[0:VV_WRF], linewidth = 1, c = 'b', label = 'QI(00z)')
ax1.plot(QCLOUD_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[:VV_EC],linewidth = 1, c = 'r', label = 'EC(00z)')
ax1.set_title('WRF 3D Qcloud and Qice  20180118 00z', fontsize=15)   #ÄNDRA
ax1.set_ylabel('Height (m)', fontsize=15)
ax1.set_xlabel('Cloud water and cloud ice (g/kg)', fontsize=15)
ax1.legend(fontsize = 'large')



fig4 = plt.figure(4, figsize=(7, 9))

ax1 = plt.subplot(111)

ax1.plot(P_CS_SOD[0:VV_WRF], MASSLEVELS_TER_SOD[0:VV_WRF], linewidth = 1, c = 'k', label = 'P(00z)')
ax1.grid(linestyle="dotted")

#fig1.savefig('/home/sm_marha/FIGURES/FOKUS_17_27_FEB/18_FEB/Initialisation_EC_WRF_OBS_00z_plus_00')

plt.show()









