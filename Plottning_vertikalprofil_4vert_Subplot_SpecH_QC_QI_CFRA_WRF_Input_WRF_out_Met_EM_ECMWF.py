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


'''Detta skript plottar vertikala profiler av specifik fuktighet, molnvatten, molnis samt molnmängd.
Man kan skifta mellan lågupplöst och högupplöst GRIB för att påvisa vikten av horisontell interpolering.'''

D = 20


'''Öppnar filobjekten'''

Filobjekt_Coldstart_WRF = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-'+str(D)+'_00z_91_LEVELS_QC_QI_42h/run/wrfinput_d01', mode='r')
Filobjekt_GRIB = N4.Dataset('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1ua201802'+str(D)+'00_M_D.nc', mode='r')
Filobjekt_GRIB_sfc = N4.Dataset('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1sfc201802'+str(D)+'00_M_D.nc' , mode='r')
Filobjekt_wrf_out_WRF = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-'+str(D)+'_00z_91_LEVELS_QC_QI_42h/run/wrfout_d01_2018-02-'+str(D)+'_00:00:00', mode='r')
Filobjekt_Met_EM_WRF = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-'+str(D)+'_00z_91_LEVELS_QC_QI_42h/run/met_em.d01.2018-02-'+str(D)+'_00:00:00.nc' , mode='r')


'''EC-DATA'''

LATITUD = 17 # Har värdet 0 uppe 
LONGITUD = 170 # Har värdet 0 på vänstra randen. ECMWF-området är längst i x-led till skillnad från WRF-området,som är längst i y-led.



'''EC-DATA HI-RES'''

#LATITUD = 53 # Har värdet 0 uppe 
#LONGITUD = 381 # Har värdet 0 på vänstra randen. ECMWF-området är längst i x-led till skillnad från WRF-området,som är längst i y-led.


'''WRF DATA'''

LATITUD_WRF = 563   # Har värdet 0 på nedre randen
LONGITUD_WRF = 352  # Har värdet 0 på vänstra randen



##############
# COLD START #
####################################################################################


ETANIVAER_MASS_CS = Filobjekt_Coldstart_WRF.variables['ZNU'][:]
ETANIVAER_FULL_CS = Filobjekt_Coldstart_WRF.variables['ZNW'][:]
ZETATOP_CS = Filobjekt_Coldstart_WRF.variables['ZETATOP'][:]
HGT_CS = Filobjekt_Coldstart_WRF.variables['HGT'][:]
PH_CS = Filobjekt_Coldstart_WRF.variables['PH'][:]
PHB_CS = Filobjekt_Coldstart_WRF.variables['PHB'][:]
GEOPOTENTIAL_CS = PH_CS + PHB_CS
MODELLNIVAHOJD_CS = (PH_CS+PHB_CS)/9.81
W_CS = Filobjekt_Coldstart_WRF.variables['W'][:]
'''PERT_THETA är den variabeln som i detta fall är vald längst upp. I många andra fall behövs inga omräkningar av variabeln utan
VARIABELVARDE är den viktiga variabeln genom hela skriptet.'''
PERT_THETA_CS = Filobjekt_Coldstart_WRF.variables['T'][:]
THETA_CS = PERT_THETA_CS+300
P_PERT_CS = Filobjekt_Coldstart_WRF.variables['P'][:]
P_BASE_CS = Filobjekt_Coldstart_WRF.variables['PB'][:]
P_CS = P_BASE_CS+P_PERT_CS
T2_CS = Filobjekt_Coldstart_WRF.variables['T2'][:]
TSK_CS = Filobjekt_Coldstart_WRF.variables['TSK'][:]
SNOWH_CS = Filobjekt_Coldstart_WRF.variables['SNOWH'][:]
TEMP_CS = THETA_CS/((1000/(P_CS/100))**0.286)
TEMP_C_CS = TEMP_CS-273.15
TEMP_C_SOD_CS = TEMP_C_CS[0, :, LATITUD_WRF, LONGITUD_WRF]
TEMP_2M_C_SOD_CS = T2_CS[0, LATITUD_WRF, LONGITUD_WRF]-273.15
TSK_C_SOD_CS = TSK_CS[0, LATITUD_WRF, LONGITUD_WRF]-273.15
SNOW_DEPTH_CS = SNOWH_CS[0, LATITUD_WRF, LONGITUD_WRF]
QVAPOR_SOD_CS = Filobjekt_Coldstart_WRF.variables['QVAPOR'][0, :, LATITUD_WRF, LONGITUD_WRF]
QICE_SOD_CS = Filobjekt_Coldstart_WRF.variables['QICE'][0, :, LATITUD_WRF, LONGITUD_WRF]
QCLOUD_SOD_CS = Filobjekt_Coldstart_WRF.variables['QCLOUD'][0, :, LATITUD_WRF, LONGITUD_WRF]

'''Beräknar MASSLEVELS Cold_start'''

MASSLEVELS_1D_CS = 0.5*(MODELLNIVAHOJD_CS[0,:-1, LATITUD_WRF, LONGITUD_WRF] + MODELLNIVAHOJD_CS[0, 1:, LATITUD_WRF, LONGITUD_WRF])
MASSLEVELS_1D_MINUS_TER_CS = MASSLEVELS_1D_CS - MODELLNIVAHOJD_CS[0, 0, LATITUD_WRF, LONGITUD_WRF]

lon = Filobjekt_Coldstart_WRF.variables['XLONG'][0, :, :]
lon_units = Filobjekt_Coldstart_WRF.variables['XLONG'].units
lat = Filobjekt_Coldstart_WRF.variables['XLAT'][0, :, :]
lat_units = Filobjekt_Coldstart_WRF.variables['XLAT'].units



#E_MATTNAD_CS = N.exp(N.log(611.2)+(22.46*TEMP_C_CS/(272.62+TEMP_C_CS)))  #Magnus formula over ice more accurate than Clausius_Clapeyron)

E_MATTNAD_CS = N.exp(N.log(611.2)+(17.62*TEMP_C_CS/(243.12+TEMP_C_CS)))  #Magnus formula over water more accurate than Clausius_Clapeyron)

#E_MATTNAD_CS = 1000*0.61078*N.exp(17.27*TEMP_C_CS/(TEMP_C_CS+237.3))  #Tetens Eq för temp>0, tryck i Pascal

#E_MATTNAD_CS = 1000*0.61078*N.exp(21.875*TEMP_C_CS/(TEMP_C_CS+265.5))  #Tetens Eq för temp<0, tryck i Pascal

QVAPOR_MATTNAD_CS = 0.622*E_MATTNAD_CS/(P_CS-E_MATTNAD_CS)

QVAPOR_MATTNAD_SOD_CS = QVAPOR_MATTNAD_CS[0, :, LATITUD_WRF, LONGITUD_WRF]

RH_SOD_CS = QVAPOR_SOD_CS/QVAPOR_MATTNAD_SOD_CS





###########
# WRF OUT #
#######################################################################################


ETANIVAER_MASS_WO = Filobjekt_wrf_out_WRF.variables['ZNU'][:]
ETANIVAER_FULL_WO = Filobjekt_wrf_out_WRF.variables['ZNW'][:]
ZETATOP_WO = Filobjekt_wrf_out_WRF.variables['ZETATOP'][:]
HGT_WO = Filobjekt_wrf_out_WRF.variables['HGT'][:]
PH_WO = Filobjekt_wrf_out_WRF.variables['PH'][:]
PHB_WO = Filobjekt_wrf_out_WRF.variables['PHB'][:]
GEOPOTENTIAL_WO = PH_WO + PHB_WO
MODELLNIVAHOJD_WO = (PH_WO+PHB_WO)/9.81
W_WO = Filobjekt_wrf_out_WRF.variables['W'][:]
'''PERT_THETA är den variabeln som i detta fall är vald längst upp. I många andra fall behövs inga omräkningar av variabeln utan
VARIABELVARDE är den viktiga variabeln genom hela skriptet.'''
PERT_THETA_WO = Filobjekt_wrf_out_WRF.variables['T'][:]
THETA_WO = PERT_THETA_WO+300
P_PERT_WO = Filobjekt_wrf_out_WRF.variables['P'][:]
P_BASE_WO = Filobjekt_wrf_out_WRF.variables['PB'][:]
P_WO = P_BASE_CS+P_PERT_CS
T2_WO = Filobjekt_wrf_out_WRF.variables['T2'][:]
TSK_WO = Filobjekt_wrf_out_WRF.variables['TSK'][:]
SNOWH_WO = Filobjekt_wrf_out_WRF.variables['SNOWH'][:]
TEMP_WO = THETA_WO/((1000/(P_WO/100))**0.286)
TEMP_C_WO = TEMP_WO-273.15
TEMP_C_SOD_WO = TEMP_C_WO[0, :, LATITUD_WRF, LONGITUD_WRF]
TEMP_2M_C_SOD_WO = T2_WO[0, LATITUD_WRF, LONGITUD_WRF]-273.15
TSK_C_SOD_WO = TSK_WO[0, LATITUD_WRF, LONGITUD_WRF]-273.15
SNOW_DEPTH_WO = SNOWH_WO[0, LATITUD_WRF, LONGITUD_WRF]
QVAPOR_SOD_WO = Filobjekt_wrf_out_WRF.variables['QVAPOR'][0, :, LATITUD_WRF, LONGITUD_WRF]
QICE_SOD_WO = Filobjekt_Coldstart_WRF.variables['QICE'][0, :, LATITUD_WRF, LONGITUD_WRF]
QCLOUD_SOD_WO = Filobjekt_Coldstart_WRF.variables['QCLOUD'][0, :, LATITUD_WRF, LONGITUD_WRF]




'''Beräknar MASSLEVELS WRF_out'''

MASSLEVELS_1D_WO = 0.5*(MODELLNIVAHOJD_WO[0,:-1, LATITUD_WRF, LONGITUD_WRF] + MODELLNIVAHOJD_WO[0, 1:, LATITUD_WRF, LONGITUD_WRF])
MASSLEVELS_1D_MINUS_TER_WO = MASSLEVELS_1D_WO - MODELLNIVAHOJD_WO[0, 0, LATITUD_WRF, LONGITUD_WRF]

lon = Filobjekt_wrf_out_WRF.variables['XLONG'][0, :, :]
lon_units = Filobjekt_wrf_out_WRF.variables['XLONG'].units
lat = Filobjekt_wrf_out_WRF.variables['XLAT'][0, :, :]
lat_units = Filobjekt_wrf_out_WRF.variables['XLAT'].units






##############
# WRF MET_EM #
#######################################################################################


T_C_EM_SOD = Filobjekt_Met_EM_WRF.variables['TT'][0, :, LATITUD_WRF, LONGITUD_WRF] - 273.15

SKINTEMP_EM_SOD = Filobjekt_Met_EM_WRF.variables['SKINTEMP'][0, LATITUD_WRF, LONGITUD_WRF]

SPECHUMD_EM_SOD = Filobjekt_Met_EM_WRF.variables['SPECHUMD'][0, :, LATITUD_WRF, LONGITUD_WRF]

RH_EM_SOD= Filobjekt_Met_EM_WRF.variables['RH'][0, :, LATITUD_WRF, LONGITUD_WRF]/100

PRESSURE_EM_SOD= Filobjekt_Met_EM_WRF.variables['PRESSURE'][0, :, LATITUD_WRF, LONGITUD_WRF]

SOILHEIGHT_EM_SOD = Filobjekt_Met_EM_WRF.variables['SOILHGT'][0, LATITUD_WRF, LONGITUD_WRF]

HEIGHT_EM_SOD = Filobjekt_Met_EM_WRF.variables['GHT'][0, :, LATITUD_WRF, LONGITUD_WRF]

H_MINUS_SH_EM_SOD = HEIGHT_EM_SOD-SOILHEIGHT_EM_SOD

HGT_M_EM_SOD = Filobjekt_Met_EM_WRF.variables['HGT_M'][0, LATITUD_WRF, LONGITUD_WRF]

H_MINUS_HGTM_EM_SOD = HEIGHT_EM_SOD-HGT_M_EM_SOD

QICE_EM_SOD = Filobjekt_Met_EM_WRF.variables['QI'][0, :, LATITUD_WRF, LONGITUD_WRF]

QCLOUD_EM_SOD = Filobjekt_Met_EM_WRF.variables['QC'][0, :, LATITUD_WRF, LONGITUD_WRF]





################
# ECMWF HYBRID #
#######################################################################################


#VARIABEL_GRIB = Filobjekt_GRIB.variables[VARIABEL_GRIB][:]
LONGITUDE = Filobjekt_GRIB.variables['longitude'][:]
LATITUDE = Filobjekt_GRIB.variables['latitude'][:]

lon_matrix, lat_matrix = N.meshgrid(LONGITUDE, LATITUDE)

QICE_EC = Filobjekt_GRIB.variables['ciwc'][:]
QICE_EC_SOD = QICE_EC[0, :, LATITUD, LONGITUD]
QICE_EC_SOD = QICE_EC_SOD[::-1]

QCLOUD_EC = Filobjekt_GRIB.variables['clwc'][:]
QCLOUD_EC_SOD = QCLOUD_EC[0, :, LATITUD, LONGITUD]
QCLOUD_EC_SOD = QCLOUD_EC_SOD[::-1]


CLOUD_FRACTION_EC = Filobjekt_GRIB.variables['cc'][:]
CLOUD_FRACTION_EC_SOD = CLOUD_FRACTION_EC[0, :, LATITUD, LONGITUD]
CLOUD_FRACTION_EC_SOD = CLOUD_FRACTION_EC_SOD[::-1]


TEMPERATURE_EC = Filobjekt_GRIB.variables['t'][:]
TEMPERATURE_C_EC = TEMPERATURE_EC-273.15
TEMPERATURE_C_EC_SOD = TEMPERATURE_C_EC[0, :, LATITUD, LONGITUD]
TEMPERATURE_C_EC_SOD = TEMPERATURE_C_EC_SOD[::-1]




'''markgeopotentialen och marktrycket från surface-gribfilen'''

SURFACE_GEOPOTENTIAL_EC = Filobjekt_GRIB_sfc.variables['z'][:]
SURFACE_GEOPOTENTIAL_EC_SOD = SURFACE_GEOPOTENTIAL_EC[0, LATITUD, LONGITUD]

SURFACE_PRESSURE_EC = Filobjekt_GRIB_sfc.variables['sp'][:]
SURFACE_PRESSURE_EC_SOD = SURFACE_PRESSURE_EC[0, LATITUD, LONGITUD]




'''Parametrar från upper_air filen''' 

SPECIFIC_HUMIDITY_EC = Filobjekt_GRIB.variables['q'][:]
SPECIFIC_HUMIDITY_EC_SOD = SPECIFIC_HUMIDITY_EC[0, :, LATITUD, LONGITUD]
SPECIFIC_HUMIDITY_EC_SOD = SPECIFIC_HUMIDITY_EC_SOD[::-1]

TIMESTEPS_EC = Filobjekt_GRIB.variables['time'][:]
VERTICAL_LEVELS_EC = Filobjekt_GRIB.variables['level'][:]



'''Använder Surface pressure på surface i stället, då dessa p g a interpolering skiljer sig lite'''

#LN_SURFACE_PRESSURE = Filobjekt_GRIB.variables['lnsp'][:]

#SURFACE_PRESSURE = N.exp(LN_SURFACE_PRESSURE)    
#SURFACE_PRESSURE_SOD = SURFACE_PRESSURE[0,2,LATITUD, LONGITUD]






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




'''Räknar ut virtuella temperaturen'''

T_VIRTUAL_EC = (1.+0.609133*SPECIFIC_HUMIDITY_EC)*TEMPERATURE_EC

T_VIRTUAL_EC_SOD = T_VIRTUAL_EC[0,:, LATITUD, LONGITUD]


T_VIRTUAL_EC_SOD = T_VIRTUAL_EC_SOD[::-1] #Inverterar vektorn så nivå 137 kommer först)




'''Räknar ut de halva trycknivåerna'''

p_half = N.zeros(len(A))

for i in range(len(A)):

    p = A[i] + B[i] * SURFACE_PRESSURE_EC[0,LATITUD,LONGITUD]

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

    Fi[i+1] = Fi[i] + Rd*(  T_VIRTUAL_EC_SOD[i]  *  (  log(p_full[i])-log(p_full[i+1])  ) )
                              
    #FI[i+1] = FI[i] + Rd*(  PF[i]/(Rd*Density[i])  *  (  log(PF[i])-log(PF[i+1])  )  )    #används till att kontrollera algorithmen via ecmwf.int



'''Plockar ut datat för Sodankylä'''

FULL_LEVELS_EC_SOD = Fi/g0




VV = 34   #20
VV_SOND = 2400   #20000
VV_EC = 60     #30






############
# PLOTTING #
#################################################################################################################

    
fig1 = plt.figure(1, figsize=(15,10))

ax = plt.subplot(141)

ax.plot(QVAPOR_SOD_CS[0:VV], MASSLEVELS_1D_MINUS_TER_CS[0:VV], linewidth = 1, label = 'WRF_Input(00z)')

ax.plot(QVAPOR_SOD_WO[0:VV], MASSLEVELS_1D_MINUS_TER_WO[0:VV],linewidth = 1, label = 'WRF_OUT(00z)' )

ax.plot(SPECHUMD_EM_SOD[1:VV_EC], H_MINUS_SH_EM_SOD[1:VV_EC],linewidth = 1, label = 'WRF_MET_EM(00z)')

ax.plot(SPECIFIC_HUMIDITY_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[0:VV_EC],linewidth = 1, label = 'ECMWF(00z)' )



#ax.plot(QVAPOR_SONDERING, MATRIX_2[0:VV_SOND, 5]-180, 'k',linewidth = 1, label = 'SOUNDING(00z)')    

ax.grid(linestyle="dotted")

ax.legend()

#ax.set_ylim(0,500)

plt.title('Specific humidity')

plt.suptitle('Coldstart vs ECMWF 201801'+str(D))

ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%8.1e'))

ax.tick_params(axis="x", labelsize=10)


plt.xlabel('(kg/kg)')

plt.ylabel('Height [m]')






ax2 = plt.subplot(142)

ax2.plot(QCLOUD_SOD_CS[0:VV], MASSLEVELS_1D_MINUS_TER_CS[0:VV],linewidth = 1, label = 'WRF_Input(00z)')

ax2.plot(QCLOUD_SOD_WO[0:VV], MASSLEVELS_1D_MINUS_TER_WO[0:VV],linewidth = 1, label = 'WRF_OUT(00z)' )

ax2.plot(QCLOUD_EM_SOD[0:VV_EC], H_MINUS_SH_EM_SOD[0:VV_EC],linewidth = 1, label = 'WRF_MET_EM(00z)')

ax2.plot(QCLOUD_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[0:VV_EC],linewidth = 1, label = 'ECMWF(00z)')

#ax2.plot(RH_SONDERING, MATRIX_2[0:VV_SOND, 5]-180, 'k',linewidth = 1, label = 'SOUNDING(00z)')

ax2.grid(linestyle="dotted")

ax2.legend(loc=2)

plt.xlabel('(kg/kg)')

ax2.set_yticklabels([])

ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('%8.1e'))

ax2.tick_params(axis="x", labelsize=9)

#ax2.set_ylim(0,500)



plt.title('Cloud water')



ax3 = plt.subplot(143)

ax3.plot(QICE_SOD_CS[0:VV], MASSLEVELS_1D_MINUS_TER_CS[0:VV], linewidth = 1, label = 'WRF_Input(00z)')

ax3.plot(QICE_SOD_WO[0:VV], MASSLEVELS_1D_MINUS_TER_WO[0:VV], linewidth = 1, label = 'WRF_OUT(00z)' )

ax3.plot(QICE_EM_SOD[0:VV_EC], H_MINUS_SH_EM_SOD[0:VV_EC], linewidth = 1, label = 'WRF_MET_EM(00z)')

ax3.plot(QICE_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[0:VV_EC], linewidth = 1, label = 'ECMWF(00z)' )

#ax3.plot(TEMP_C_SONDERING, MATRIX_2[0:VV_SOND, 5]-180, 'k',linewidth = 1, label = 'SOUNDING(00z)')
#ax3.plot(THETA_FG[0, 0:VV, LATITUD_WRF, LONGITUD_WRF], MASSLEVELS_1D_MINUS_TER_CS[0:VV],'r',linewidth = 1, label = 'WRF_CS(18z)')

#ax3.set_ylim(0,500)
 

ax3.grid(linestyle="dotted")

ax3.legend()

plt.title('Cloud ice')

ax3.set_yticklabels([])

#labels = ["{8.1e}".format(x) for x in QICE_EC_SOD[0:VV_EC]]
#ax3.set_xticklabels(labels)

ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter('%8.1e'))

ax3.tick_params(axis="x", labelsize=9)

plt.xlabel('(kg/kg)')

   


ax4 = plt.subplot(144)

#ax4.plot([0:11], MASSLEVELS_1D_MINUS_TER_CS[0:11], 'r', label = 'WRF_CS')

ax4.plot(CLOUD_FRACTION_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[0:VV_EC],linewidth = 1, label = 'ECMWF(00z)' )



#ax4.scatter(0, MOLNBAS[j], c = 'k', s=12)

ax4.grid(linestyle="dotted")

ax4.legend()

plt.title('Cloud fraction')
    
ax4.set_yticklabels([])

ax4.tick_params(axis="x", labelsize=9)

#ax4.set_ylim(0,500)

plt.xlabel('(%)')




'''Sparar figuren och visar plotten'''


plt.show()

#fig1.savefig('/home/sm_marha/FIGURES/FOKUS_18_20_JAN/20180118_CW_CI_CLDFRA_HIGHRES_GRIB')
    


   












    
