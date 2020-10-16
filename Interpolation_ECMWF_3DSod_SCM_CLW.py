import numpy as N
import scipy.io.netcdf as S
import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4 as N4
import pandas as pd
from math import *


######
#GRIB#
######################################################



Filobjekt_GRIB = S.netcdf_file('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1ua2018022000_M_D.nc', mode='r')
Filobjekt_GRIB_sfc = S.netcdf_file('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1sfc2018022000_M_D.nc', mode='r')


LAT=27    #SOD:53

LON=191   #SOD:381

VARIABLE_EC = 'clwc'


'''Parametrar från surface filen'''

SURFACE_GEOPOTENTIAL_GRIB = Filobjekt_GRIB_sfc.variables['z'][:]
SURFACE_GEOPOTENTIAL_SOD_GRIB = SURFACE_GEOPOTENTIAL_GRIB[0, LAT, LON]

SURFACE_PRESSURE_GRIB = Filobjekt_GRIB_sfc.variables['sp'][:]
SURFACE_PRESSURE_SOD_GRIB = SURFACE_PRESSURE_GRIB[0, LAT, LON]


'''Parametrar från upper_air filen'''

CLOUD_VARIABLE_GRIB = Filobjekt_GRIB.variables[VARIABLE_EC][:]
CLOUD_VARIABLE_SOD_GRIB = CLOUD_VARIABLE_GRIB[0, ::-1, LAT, LON]*1000

SURFACE_GEOPOTENTIAL_GRIB = Filobjekt_GRIB_sfc.variables['z'][:]
SURFACE_GEOPOTENTIAL_SOD_GRIB = SURFACE_GEOPOTENTIAL_GRIB[0, LAT, LON]

SURFACE_PRESSURE_GRIB = Filobjekt_GRIB_sfc.variables['sp'][:]
SURFACE_PRESSURE_SOD_GRIB = SURFACE_PRESSURE_GRIB[0, LAT, LON]

TEMPERATURE_GRIB = Filobjekt_GRIB.variables['t'][:]
TEMPERATURE_SOD_GRIB = TEMPERATURE_GRIB[0, :, LAT, LON]

SPECIFIC_HUMIDITY_GRIB = Filobjekt_GRIB.variables['q'][:]
SPECIFIC_HUMIDITY_SOD_GRIB = SPECIFIC_HUMIDITY_GRIB[0, :, LAT, LON]

TIMESTEPS_GRIB = Filobjekt_GRIB.variables['time'][:]
VERTICAL_LEVELS_GRIB = Filobjekt_GRIB.variables['level'][:]



'''Läser in a-. och b-koeffecienter'''

df = pd.read_csv(r'/home/sm_marha/TEXTFILER/A_B_Coefficient_ECMWF.csv', encoding='latin-1', delimiter=';', header = None)

df = df.iloc[:,:]

df_numpy_array = df.values

A= df_numpy_array[:, 1]
B= df_numpy_array[:, 2]




'''Räknar ut virtuella temperaturen'''

T_VIRTUAL_GRIB = (1.+0.609133*SPECIFIC_HUMIDITY_GRIB)*TEMPERATURE_GRIB

T_VIRTUAL_SOD_GRIB = T_VIRTUAL_GRIB[0, :, LAT, LON]


T_VIRTUAL_SOD_GRIB = T_VIRTUAL_SOD_GRIB[::-1] #Inverterar vektorn så nivå 137 kommer först)




'''Räknar ut de halva trycknivåerna'''

p_half = N.zeros(len(A))

for i in range(len(A)):

    p = A[i] + B[i] * SURFACE_PRESSURE_GRIB[0,LAT,LON]

    p_half[i] = p
    


'''Räknar ut hela trycknivåerna'''
     
p_full = 0.5*(p_half[1:]+p_half[:-1])



'''Räknar ut geopotentialen på de fulla nivåerna.'''

Rd = 287.06
g0 = 9.80665


Fi = N.zeros(len(p_full))
Fi[0]=10*g0


for i in range(len(p_full)-1):

    Fi[i+1] = Fi[i] + Rd*(  T_VIRTUAL_SOD_GRIB[i]  *  (  log(p_full[i])-log(p_full[i+1])  ) )
    


'''Plockar ut datat för Sodankylä'''

'''Välj med, eller utan terräng'''

FULL_LEVELS_SOD_GRIB = Fi/g0

#FULL_LEVELS_SOD_GRIB = FULL_LEVELS_SOD_GRIB+222   #lägger till modellhöjden för att få det på samma sätt som i de WRF-filer jag jämför med.






########
#WRF 3D#
######################################################

LATITUDE = 563    #Sodankylä = 563

LONGITUDE = 352    #Sodankylä = 352

VARIABLE_WRF = 'QCLOUD'

'''Välj 46 eller 91 nivåer'''


#Filobjekt_3D_wrfinput = N4.Dataset('/nobackup/smhid12/sm_marha/2018-01-18_00z_91_LEVELS_QC_QI/run/wrfinput_d01', mode='r')
Filobjekt_3D_wrfinput = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-20_00z_91_LEVELS_QC_QI_42h/run/wrfinput_d01', mode='r')
Filobjekt_Met_EM_WRF = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-20_00z_91_LEVELS_QC_QI_42h/run/met_em.d01.2018-02-20_00:00:00.nc', mode='r')

SOILHEIGHT = Filobjekt_Met_EM_WRF.variables['SOILHGT'][0, 563, 352]
HEIGHT = Filobjekt_Met_EM_WRF.variables['GHT'][0, :, 563, 352]
HALF_LEVELS_SOD = HEIGHT-SOILHEIGHT

CLOUD_VARIABLE_SOD_MET_EM = Filobjekt_Met_EM_WRF.variables['QC'][0, :, 563, 352]*1000




CLOUD_VARIABLE_SOD_3D = Filobjekt_3D_wrfinput.variables[VARIABLE_WRF][0, :, LATITUDE, LONGITUDE]*1000

'''Beräknar parametrar Cold_start'''

ETANIVAER_MASS_3D = Filobjekt_3D_wrfinput.variables['ZNU'][:]
ETANIVAER_FULL_3D = Filobjekt_3D_wrfinput.variables['ZNW'][:]
ZETATOP_3D = Filobjekt_3D_wrfinput.variables['ZETATOP'][:]
HGT_3D = Filobjekt_3D_wrfinput.variables['HGT'][:]
PH_3D = Filobjekt_3D_wrfinput.variables['PH'][:]
PHB_3D = Filobjekt_3D_wrfinput.variables['PHB'][:]
GEOPOTENTIAL_3D = PH_3D + PHB_3D
MODELLNIVAHOJD_3D = (PH_3D+PHB_3D)/9.81
W_3D = Filobjekt_3D_wrfinput.variables['W'][:]
'''PERT_THETA är den variabeln som i detta fall är vald längst upp. I många andra fall behövs inga omräkningar av variabeln utan
VARIABELVARDE är den viktiga variabeln genom hela skriptet.'''
PERT_THETA_3D = Filobjekt_3D_wrfinput.variables['T'][:]
THETA_3D = PERT_THETA_3D+300
P_PERT_3D = Filobjekt_3D_wrfinput.variables['P'][:]
P_BASE_3D = Filobjekt_3D_wrfinput.variables['PB'][:]
P_3D = P_BASE_3D+P_PERT_3D
TEMP_3D = THETA_3D/((1000/(P_3D/100))**0.286)
TEMP_C_3D = TEMP_3D-273.15
TEMP_C_SOD_3D = TEMP_C_3D[0, :, LATITUDE, LONGITUDE]



'''Beräknar MASSLEVELS WRF 3D'''

'''välj terräng, eller utan terräng'''

MASSLEVELS_1D_3D = 0.5*(MODELLNIVAHOJD_3D[0,:-1, LATITUDE, LONGITUDE] + MODELLNIVAHOJD_3D[0, 1:, LATITUDE, LONGITUDE])
MASSLEVELS_1D_MINUS_TER_3D = MASSLEVELS_1D_3D - MODELLNIVAHOJD_3D[0, 0, LATITUDE, LONGITUDE]



lon = Filobjekt_3D_wrfinput.variables['XLONG'][0, :, :]
lon_units = Filobjekt_3D_wrfinput.variables['XLONG'].units
lat = Filobjekt_3D_wrfinput.variables['XLAT'][0, :, :]
lat_units = Filobjekt_3D_wrfinput.variables['XLAT'].units







#########
#WRF SCM#
######################################################

'''Välj 46 eller 91 nivåer'''

Filobjekt_SCM_wrfinput = N4.Dataset('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfinput_d01', mode='r')
#Filobjekt_SCM_wrfinput = N4.Dataset('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY/test/em_scm_xy/wrfinput_d01', mode='r')



'''Parametrar'''

CLOUD_VARIABLE_SCM = Filobjekt_SCM_wrfinput.variables[VARIABLE_WRF][0, :, 1, 1]*1000



'''Beräknar vertikala nivåer'''

PH_SCM = Filobjekt_SCM_wrfinput.variables['PH'][0, :, 1, 1]
PHB_SCM = Filobjekt_SCM_wrfinput.variables['PHB'][0, :, 1, 1]
MODELLNIVAHOJD_SCM = (PH_SCM+PHB_SCM)/9.81
MODELLNIVAHOJD_MINUS_TER_SCM = MODELLNIVAHOJD_SCM -MODELLNIVAHOJD_SCM[0]


'''Välj med, eller utan terränghöjd'''

#MASSLEVELS_SCM = 0.5*(MODELLNIVAHOJD_SCM[:-1] + MODELLNIVAHOJD_SCM[1:])
MASSLEVELS_SCM = 0.5*(MODELLNIVAHOJD_MINUS_TER_SCM[:-1] + MODELLNIVAHOJD_MINUS_TER_SCM[1:])




#################################### 
#            PLOTTNING             # 
#################################### 


##########
#FIGURE 1#
##########


fig1 = plt.figure(1)

ax1 = plt.subplot(111)

ax1.plot(CLOUD_VARIABLE_SOD_GRIB[:37], FULL_LEVELS_SOD_GRIB[:37], linewidth = 1, label = 'ECMWF')
ax1.scatter(CLOUD_VARIABLE_SOD_GRIB[:37], FULL_LEVELS_SOD_GRIB[:37], c='b', s=5)
ax1.plot(CLOUD_VARIABLE_SOD_3D[:], MASSLEVELS_1D_MINUS_TER_3D[:], linewidth = 1, label = 'WRF 3D')
ax1.scatter(CLOUD_VARIABLE_SOD_3D[:], MASSLEVELS_1D_MINUS_TER_3D[:], color='orange', s=5)
ax1.plot(CLOUD_VARIABLE_SCM[:], MASSLEVELS_SCM[:], linewidth = 1, label = 'WRF_SCM')
ax1.scatter(CLOUD_VARIABLE_SCM[:], MASSLEVELS_SCM[:], c='g', s=5)

ax1.plot(CLOUD_VARIABLE_SOD_MET_EM[:], HALF_LEVELS_SOD[:], linewidth = 1, label = 'MET_EM')
ax1.scatter(CLOUD_VARIABLE_SOD_MET_EM[:], HALF_LEVELS_SOD[:], color='black', s=5)


#ax1.set_xlim(0, 0.0001)
ax1.set_ylim(0, 3000)
ax1.legend()
ax1.grid(linestyle="dotted")
plt.ylabel('Height (m)')
plt.xlabel(Filobjekt_3D_wrfinput[VARIABLE_WRF].description + ' (g/kg)')
plt.title('Comparing ' + Filobjekt_3D_wrfinput[VARIABLE_WRF].description + ' ECMWF, WRF 3D')


'''Välj 46 eller 91 nivåer'''

#fig1.savefig('/home/sm_marha/FIGURES/SINGLE_COLUMN/18011800Z/91_LEVELS/COMPARING_INITIAL_PROFILES_EC_3D_SCM/QC_QI_3_LEVELS')
#fig1.savefig('/home/sm_marha/FIGURES/SINGLE_COLUMN/18011800Z/COMPARING_INITIAL_PROFILES_EC_3D_SCM/QC_QI_2_LEVELS')'''

plt.show()






