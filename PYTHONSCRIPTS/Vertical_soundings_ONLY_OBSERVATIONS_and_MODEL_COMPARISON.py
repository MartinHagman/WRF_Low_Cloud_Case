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
from pylab import *

rc('axes', linewidth=1.5)

fig1, axs = plt.subplots(2, 4, figsize=(10,8))
fig1.subplots_adjust(hspace=0.4)
#fig.subplots_adjust(right=0.50)



ax1 = axs[0, 0]
ax1w = axs[0, 1]
ax2 = axs[0, 2]
ax2w = axs[0, 3]
ax3 = axs[1, 0]
ax3w = axs[1, 1]
ax4 = axs[1, 2]
ax4w = axs[1, 3]



'''Lägger till vind diagrammen på rätt position, genom att flytta temperatureaxlarna åt vänster och sen malla in vindaxlarna'''

print(ax1.get_position())
box_ax1 =  ax1.get_position()

ax1.set_position([box_ax1.x0, box_ax1.y0, box_ax1.width*1.6, box_ax1.height])
ax1.set_ylabel('Height ($m$)', fontsize='18')
ax1.xaxis.set_visible(True)
ax1.set_xticklabels([])


print("\n")


print(ax1w.get_position())
box_ax1w =  ax1w.get_position()

ax1w.set_position([0.407, box_ax1w.y0, box_ax1w.width*0.5, box_ax1w.height])
ax1w.xaxis.set_visible(True)
ax1w.set_yticklabels([])
ax1w.set_frame_on(True)
ax1w.set_xticklabels([])
print("\n")


print(ax2.get_position())
box_ax2 =  ax2.get_position()

ax2.set_position([box_ax2.x0, box_ax2.y0, box_ax2.width*1.6, box_ax2.height])
ax2.set_yticklabels([])
ax2.set_xticklabels([])
ax2.xaxis.set_visible(True)
print("\n")


print(ax2w.get_position())
box_ax2w =  ax2w.get_position()

ax2w.set_position([0.811, box_ax2w.y0, box_ax2w.width*0.5, box_ax2w.height])
ax2w.xaxis.set_visible(True)
ax2w.set_yticklabels([])
ax2w.set_frame_on(True)
ax2w.set_xticklabels([])
print("\n")


print(ax3.get_position())
box_ax3 =  ax3.get_position()
ax3.set_position([box_ax3.x0, box_ax3.y0, box_ax3.width*1.6, box_ax3.height])
ax3.set_xlabel('T ($^oC$)', fontsize='18')
ax3.set_ylabel('Height ($m$)', fontsize='18')
ax3.xaxis.set_visible(True)
#ax3.set_xticklabels([])

print("\n")


print(ax3w.get_position())
box_ax3w =  ax3w.get_position()

ax3w.set_position([0.407, box_ax3w.y0, box_ax3w.width*0.5, box_ax3w.height])
ax3w.xaxis.set_visible(True)
ax3w.set_yticklabels([])
ax3w.set_frame_on(True)
ax3w.set_xlabel('Wind ($ms^{-1}$)', fontsize='17')
print("\n")


print(ax4.get_position())
box_ax4 =  ax4.get_position()

ax4.set_position([box_ax4.x0, box_ax4.y0, box_ax4.width*1.6, box_ax4.height])
ax4.set_yticklabels([])
ax4.set_xlabel('T ($^oC$)', fontsize='18')
ax4.xaxis.set_visible(True)
#ax4.set_xticklabels([])
print("\n")


print(ax4w.get_position())
box_ax4w =  ax4w.get_position()
ax4w.xaxis.set_visible(True)
ax4w.set_yticklabels([])
ax4w.set_frame_on(True)
ax4w.set_xlabel('Wind ($ms^{-1}$)', fontsize='17')
ax4w.set_position([0.811, box_ax4w.y0, box_ax4w.width*0.5, box_ax4w.height])


print("\n")




'''Vind-diagrammen'''


AXS  = [ax1, ax2, ax3, ax4]

AX_WINDS = [ax1w, ax2w, ax3w, ax4w]


ax11 = ax1.twiny()
ax22 = ax2.twiny()
ax33 = ax3.twiny()
ax44 = ax4.twiny()




'''Den osynliga andra x-axeln.'''


AX_THETAS = [ax11, ax22, ax33, ax44]




#ax11.spines['right'].set_visible(False)
#ax11.set_frame_on(False)
ax11.set_position([box_ax1.x0, box_ax1.y0, box_ax1.width*1.6, box_ax1.height])
ax11.set_position([box_ax1.x0, box_ax1.y0, box_ax1.width*1.6, box_ax1.height])
ax11.xaxis.set_visible(False)
ax11.yaxis.set_visible(False)
ax11.set_frame_on(True)


ax22.set_position([box_ax2.x0, box_ax2.y0, box_ax2.width*1.6, box_ax2.height])
ax22.xaxis.set_visible(False)
ax22.yaxis.set_visible(False)
ax22.set_frame_on(False)

ax33.set_position([box_ax3.x0, box_ax3.y0, box_ax3.width*1.6, box_ax3.height])
ax33.xaxis.set_visible(False)
ax33.yaxis.set_visible(False)
ax33.set_frame_on(False)

ax44.set_position([box_ax4.x0, box_ax4.y0, box_ax4.width*1.6, box_ax4.height])
ax44.xaxis.set_visible(False)
ax44.yaxis.set_visible(False)
ax44.set_frame_on(False)


WRF_TIMES = N.array(['00', '06', '12', '18'])   #18z finns tyvärr ej 

IFS_TIMESTEPS = N.array([0, 2, 4, 6])

SOUNDING_TIMES = N.array(['00', '6', '12', '18'])

#SCM_TIMESTEPS = N.array([1080, 2160, 3240, 4320])
SCM_TIMESTEPS = N.array([0, 1080, 2160, 3240])

r=0

for AX, AX_WIND, AX_THETA, SOUNDING_TIME, WRF_TIME, IFS_TIMESTEP, SCM_TIMESTEP in zip(AXS, AX_WINDS, AX_THETAS, SOUNDING_TIMES, WRF_TIMES, IFS_TIMESTEPS, SCM_TIMESTEPS):



##############
#OBSERVATIONS#
##################################################

    
  
    #######################
    #HÄMTAR SONDERINGSDATA#
    #######################

    j=49    #ÄNDRA (18 FEB)          #ÄNDRA

    



    '''Hämtar sonderingsdata med hjälp av funktionen 'pick dates' i Reading_Dat_Files'''


    MATRIX_2, datum = Picking_values_Sounding_Dat_file_version_2.pick_dates(j, SOUNDING_TIME)


    MATRIX_2 = (N.array(MATRIX_2)).astype(float)




    '''Beräknar relativ fuktighet och specifik fuktighet'''



    VV_SOND = 700  # 470   #2000



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

    WIND_VELOCITY_SONDERING = MATRIX_2[0:VV_SOND, 0]

    WIND_VELOCITY_SONDERING = MATRIX_2[0:VV_SOND, 0]

    WIND_DIRECTION_SONDERING = MATRIX_2[0:VV_SOND, 2]

    U_WIND_SONDERING = -WIND_VELOCITY_SONDERING*N.sin(pi/180*WIND_DIRECTION_SONDERING)

    V_WIND_SONDERING = -WIND_VELOCITY_SONDERING*N.cos(pi/180*WIND_DIRECTION_SONDERING)

    ALTITUDE_WIND_SONDERING = MATRIX_2[0:VV_SOND, 1]



    #####
    #WRF#
    #####################################################################

    #Filobjekt_WRF_3D_WRFOUT = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-26_00z_91_LEVELS_QC_QI_24h_T_to_Td_1e-5/run/wrfout_d01_2018-02-26_%s:00:00' %WRF_TIME, mode='r')  #ÄNDRA

    #Filobjekt_WRF_3D_WRFOUT = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-20_00z_91_LEVELS_QC_QI_Icloud_2/run/wrfout_d01_2018-02-20_%s:00:00' %WRF_TIME, mode='r')  #ÄNDRA

    Filobjekt_WRF_3D_WRFOUT = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-18_00z_91_LEVELS_QC_QI_42h/run/wrfout_d01_2018-02-18_%s:00:00' %WRF_TIME, mode='r')  #ÄNDRA

    

    VV_WRF = 40  #40

    LAT = 563 #563

    LON = 352 #352

    PERT_THETA_WRF_3D_WRFOUT = Filobjekt_WRF_3D_WRFOUT.variables['T'][0, :, LAT, LON]

    P_PERT_WRF_3D_WRFOUT = Filobjekt_WRF_3D_WRFOUT.variables['P'][0, :, LAT, LON]
    P_BASE_WRF_3D_WRFOUT = Filobjekt_WRF_3D_WRFOUT.variables['PB'][0, :, LAT, LON]
    P_WRF_3D_WRFOUT = P_BASE_WRF_3D_WRFOUT+P_PERT_WRF_3D_WRFOUT

    THETA_WRF_3D_WRFOUT = PERT_THETA_WRF_3D_WRFOUT+300

    

    TEMP_WRF_3D_WRFOUT = THETA_WRF_3D_WRFOUT/((1000/(P_WRF_3D_WRFOUT/100))**0.286)
    TEMP_C_WRF_3D_WRFOUT = TEMP_WRF_3D_WRFOUT-273.15
    QVAPOR_WRF_3D_WRFOUT = Filobjekt_WRF_3D_WRFOUT.variables['QVAPOR'][0, :, LAT, LON]
    E_MATTNAD_WRF_3D_WRFOUT = N.exp(N.log(611.2)+(17.62*TEMP_C_WRF_3D_WRFOUT/(243.12+TEMP_C_WRF_3D_WRFOUT)))
    QVAPOR_MATTNAD_WRF_3D_WRFOUT = 0.622*E_MATTNAD_WRF_3D_WRFOUT/(P_WRF_3D_WRFOUT-E_MATTNAD_WRF_3D_WRFOUT)
    CLOUD_FRACTION_WRF_3D_WRFOUT = Filobjekt_WRF_3D_WRFOUT.variables['CLDFRA'][0, :, LAT, LON]
    
    RH_WRF_3D_WRFOUT = QVAPOR_WRF_3D_WRFOUT/QVAPOR_MATTNAD_WRF_3D_WRFOUT

    T_DAGGPUNKT_3D_WRFOUT = (243.12*N.log(611.2)-243.12*N.log(RH_WRF_3D_WRFOUT*E_MATTNAD_WRF_3D_WRFOUT)) / (N.log(RH_WRF_3D_WRFOUT*E_MATTNAD_WRF_3D_WRFOUT)-17.62-N.log(611.2))
              
    PH_WRF_3D_WRFOUT = Filobjekt_WRF_3D_WRFOUT.variables['PH'][0, :, LAT, LON]
    PHB_WRF_3D_WRFOUT = Filobjekt_WRF_3D_WRFOUT.variables['PHB'][0, :, LAT, LON]
    MODELLNIVAHOJD_WRF_3D_WRFOUT = (PH_WRF_3D_WRFOUT+PHB_WRF_3D_WRFOUT)/9.81
    MODELLNIVAHOJD_TER_WRF_3D_WRFOUT = MODELLNIVAHOJD_WRF_3D_WRFOUT[:] -MODELLNIVAHOJD_WRF_3D_WRFOUT[0]

    


    MASSLEVELS_WRF_3D_WRFOUT = 0.5*(MODELLNIVAHOJD_TER_WRF_3D_WRFOUT[:-1] + MODELLNIVAHOJD_TER_WRF_3D_WRFOUT[1:])


    QCLOUD_SOD_WRF_3D_WRFOUT = Filobjekt_WRF_3D_WRFOUT.variables['QCLOUD'][0, :, LAT, LON]*1000    # I g/kg
    QICE_SOD_WRF_3D_WRFOUT = Filobjekt_WRF_3D_WRFOUT.variables['QICE'][0, :, LAT, LON]*1000        # I g/kg

    THETA_LIQUID_SOD_3D_WRFOUT = THETA_WRF_3D_WRFOUT - (2260*THETA_WRF_3D_WRFOUT)/(1004*TEMP_WRF_3D_WRFOUT)* QCLOUD_SOD_WRF_3D_WRFOUT/1000

    U = Filobjekt_WRF_3D_WRFOUT.variables['U'][:]
    V = Filobjekt_WRF_3D_WRFOUT.variables['V'][:]
    U_10 = Filobjekt_WRF_3D_WRFOUT.variables['U10'][:]
    V_10 = Filobjekt_WRF_3D_WRFOUT.variables['V10'][:]
    XLAT = Filobjekt_WRF_3D_WRFOUT.variables['XLAT'][:]
    XLONG = Filobjekt_WRF_3D_WRFOUT.variables['XLONG'][:]
    X_LAT_U = Filobjekt_WRF_3D_WRFOUT.variables['XLAT_U'][:]
    X_LAT_U = Filobjekt_WRF_3D_WRFOUT.variables['XLAT_U'][:]
    X_LON_U = Filobjekt_WRF_3D_WRFOUT.variables['XLONG_U'][:]
    X_LAT_V = Filobjekt_WRF_3D_WRFOUT.variables['XLAT_V'][:]
    X_LONG_U = Filobjekt_WRF_3D_WRFOUT.variables['XLONG_V'][:]
    SINALPHA = Filobjekt_WRF_3D_WRFOUT.variables['SINALPHA'][:]
    COSALPHA = Filobjekt_WRF_3D_WRFOUT.variables['COSALPHA'][:]

   

    U_3D = U[0, :, :, :] 
    V_3D = V[0, :, :, :]

    U_3D_UNSTAG = 0.5*(U_3D[:,:,:-1] + U_3D[:,:,1:])
    V_3D_UNSTAG = 0.5*(V_3D[:,:-1,:] + V_3D[:,1:,:])

    '''print(N.shape(U_3D_UNSTAG))
    print(N.shape(V_3D_UNSTAG))'''

    U_3D_UNSTAG_EARTH = U_3D_UNSTAG*COSALPHA - V_3D_UNSTAG*SINALPHA
    V_3D_UNSTAG_EARTH = V_3D_UNSTAG*COSALPHA + U_3D_UNSTAG*SINALPHA

    U_SOD = U_3D_UNSTAG_EARTH[:, 563, 352]
    V_SOD = V_3D_UNSTAG_EARTH[:, 563, 352]

    VELOCITY_SOD = N.sqrt(U_SOD**2+V_SOD**2)

    '''for place, element in  enumerate(QCLOUD_SOD_WRF_3D_WRFOUT):
        if element < 0:
            QCLOUD_SOD_WRF_3D_WRFOUT[place] = 0'''
            

    ################
    # ECMWF HYBRID #
    #######################################################################################


    Filobjekt_GRIB = N4.Dataset('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1ua2018021800_M_D.nc', mode='r')
    Filobjekt_GRIB_sfc = N4.Dataset('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1sfc2018021800_M_D.nc' , mode='r')



    TIDSSTEG = IFS_TIMESTEP    #ÄNDRA
    VV_EC =35

    LATITUD = 27 # Har värdet 0 uppe 
    LONGITUD = 191 # Har värdet 0 på vänstra randen. ECMWF-området är längst i x-led till skillnad från WRF-området,som är längst i y-led.

    '''TIMESTEPS = Filobjekt_GRIB.variables['time'][:]
    time = Filobjekt_GRIB.variables['time']
    dates = N4.num2date(TIMESTEPS, time.units, time.calendar)
    ALLA_DATUM=[]

    for date in dates:
        datum = date.strftime('%d/%m %Hz')
        ALLA_DATUM.append(datum)'''


    LONGITUDE = Filobjekt_GRIB.variables['longitude'][:]
    LATITUDE = Filobjekt_GRIB.variables['latitude'][:]

    lon_matrix, lat_matrix = N.meshgrid(LONGITUDE, LATITUDE)

    U_WIND= Filobjekt_GRIB.variables['u'][TIDSSTEG, :, LATITUD, LONGITUD]
    V_WIND = Filobjekt_GRIB.variables['v'][TIDSSTEG, :, LATITUD, LONGITUD]

    U_WIND = U_WIND[::-1]
    V_WIND = V_WIND[::-1]

    VELOCITY_EC = N.sqrt(U_WIND**2 + V_WIND**2)

    QICE_EC = Filobjekt_GRIB.variables['ciwc'][:]*1000
    QICE_EC_SOD = QICE_EC[TIDSSTEG, :, LATITUD, LONGITUD]
    QICE_EC_SOD = QICE_EC_SOD[::-1]

    QCLOUD_EC = Filobjekt_GRIB.variables['clwc'][:]*1000
    QCLOUD_EC_SOD = QCLOUD_EC[TIDSSTEG, :, LATITUD, LONGITUD]    
    QCLOUD_EC_SOD = QCLOUD_EC_SOD[::-1]
    print(QCLOUD_EC_SOD)

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

    E_MATTNAD_EC_SOD = N.exp(N.log(611.2)+(17.62*TEMPERATURE_C_EC_SOD/(243.12+TEMPERATURE_C_EC_SOD)))
    


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
   

    THETA_LIQUID_EC_SOD = THETA_EC_SOD - ((2.26e6*THETA_EC_SOD)/(1004*TEMPERATURE_EC_SOD)* QCLOUD_EC_SOD/1000)

    #print('THETA-THETA_liq')
    #print(THETA_EC_SOD - THETA_LIQUID_EC_SOD)

    
    '''Beräknar QVAPOR mättnad. Detta görs också först här för att p_full behövs på varje nivå.'''

    QVAPOR_MATTNAD_EC_SOD = 0.622*E_MATTNAD_EC_SOD/(p_full-E_MATTNAD_EC_SOD)

    SPECIFIC_HUMIDITY_MATTNAD_EC_SOD = QVAPOR_MATTNAD_EC_SOD/(1 + QVAPOR_MATTNAD_EC_SOD)

    RH_EC_SOD = SPECIFIC_HUMIDITY_EC_SOD/SPECIFIC_HUMIDITY_MATTNAD_EC_SOD

    T_DAGGPUNKT_EC_SOD = (243.12*N.log(611.2)-243.12*N.log(RH_EC_SOD*E_MATTNAD_EC_SOD)) / (N.log(RH_EC_SOD*E_MATTNAD_EC_SOD)-17.62-N.log(611.2))
   
    



    #################
    # WRF SCM WRFOUT#
    #######################################################################################

    Filobjekt_SCM_WRFOUT = N4.Dataset('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfout_d01_2018-02-20_00:00:00_QC_QI_REFERENCE')


    #TIDPUNKT = 0 #ÄNDRA   7561 Sätt hur mycket data som ska plockas ut. Kan även göras vid plottning genom att sätta variabeln t.

    TIMESTEPS_SCM_WRFOUT = Filobjekt_SCM_WRFOUT.variables['XTIME'][SCM_TIMESTEP]


    PERT_THETA_SCM_WRFOUT = Filobjekt_SCM_WRFOUT.variables['T'][SCM_TIMESTEP, :, 1, 1]

    P_PERT_SCM_WRFOUT = Filobjekt_SCM_WRFOUT.variables['P'][SCM_TIMESTEP, :, 1, 1]
    P_BASE_SCM_WRFOUT = Filobjekt_SCM_WRFOUT.variables['PB'][SCM_TIMESTEP, :, 1, 1]
    P_SCM_WRFOUT = P_BASE_SCM_WRFOUT+P_PERT_SCM_WRFOUT

    THETA_SCM_WRFOUT = PERT_THETA_SCM_WRFOUT+300

    TEMP_SCM_WRFOUT = THETA_SCM_WRFOUT/((1000/(P_SCM_WRFOUT/100))**0.286)
    TEMP_C_SCM_WRFOUT = TEMP_SCM_WRFOUT-273.15
    QVAPOR_SCM_WRFOUT = Filobjekt_SCM_WRFOUT.variables['QVAPOR'][SCM_TIMESTEP, :, 1, 1]
    E_MATTNAD_SCM_WRFOUT = N.exp(N.log(611.2)+(17.62*TEMP_C_SCM_WRFOUT/(243.12+TEMP_C_SCM_WRFOUT)))
    QVAPOR_MATTNAD_SCM_WRFOUT = 0.622*E_MATTNAD_SCM_WRFOUT/(P_SCM_WRFOUT-E_MATTNAD_SCM_WRFOUT)
              
    PH_SCM_WRFOUT = Filobjekt_SCM_WRFOUT.variables['PH'][SCM_TIMESTEP, :, 1, 1]
    PHB_SCM_WRFOUT = Filobjekt_SCM_WRFOUT.variables['PHB'][SCM_TIMESTEP, :, 1, 1]
    MODELLNIVAHOJD_SCM_WRFOUT = (PH_SCM_WRFOUT+PHB_SCM_WRFOUT)/9.81
    MODELLNIVAHOJD_TER_SCM_WRFOUT = MODELLNIVAHOJD_SCM_WRFOUT[:] -MODELLNIVAHOJD_SCM_WRFOUT[0]


    MASSLEVELS_SCM_WRFOUT = 0.5*(MODELLNIVAHOJD_TER_SCM_WRFOUT[:-1] + MODELLNIVAHOJD_TER_SCM_WRFOUT[1:])
    
    
    ###########
    #PLOTTNING#
    ######################################################################

    rc('axes', linewidth=1.5)

    #Rc-params måste sättas innan anrop till subplot, men gäller sedan för alla tills nytt anrop görs.
    plt.rc('xtick', labelsize='x-large')
    plt.rc('ytick', labelsize='x-large')


    AX.grid(linestyle="dotted")
    AX_WIND.grid(linestyle="dotted")
    '''Lns1 osv för att få legenderna i samma lager'''
    lns1 = AX.plot(TEMP_C_SONDERING[0:VV_SOND], HEIGHT_TER_SONDERING[0:VV_SOND], linewidth = 1.4, c = 'k', label = 'T/Wind OBS')
    lns2 = AX.plot(DAGGPUNKT_SONDERING[0:VV_SOND], HEIGHT_TER_SONDERING[0:VV_SOND], linewidth = 1.4, c = 'k', linestyle = 'dashed', label = 'DP OBS')
    #AX_WIND.barbs(0.5*N.ones(VV_SOND-15)[0::18], HEIGHT_TER_SONDERING[0:VV_SOND-15:18], U_WIND_SONDERING[0:VV_SOND-15:18], V_WIND_SONDERING[0:VV_SOND-15:18], barbcolor = 'k')
    AX_WIND.plot(WIND_VELOCITY_SONDERING[0:VV_SOND], HEIGHT_TER_SONDERING[0:VV_SOND],linewidth = 1.4, c = 'k', label = 'OBS')
    AX_WIND.plot(VELOCITY_SOD[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF], linewidth = 1, c = 'b', label = 'WRF')
    AX_WIND.plot(VELOCITY_EC[0:VV_EC], FULL_LEVELS_EC_SOD[:VV_EC], linewidth = 1, c = 'r', label = 'IFS')

    
    #lns3 = AX_THETA.plot(THETA_SONDERING[0:VV_SOND], HEIGHT_TER_SONDERING[0:VV_SOND], linewidth = 1.5, c = 'k', label = 'PT OBS')

    #lns4 = AX_THETA.plot(THETA_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[:VV_EC],linewidth = 1, c = 'r', label = 'PT IFS')

    #if WRF_TIME != '09':
        #lns5 = AX_THETA.plot(THETA_WRF_3D_WRFOUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF], linewidth = 1, c = 'b', label = 'PT WRF')
        #lns6 = AX_THETA.plot(THETA_LIQUID_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[:VV_EC], linewidth = 1, c = 'r', label = 'IFS_TH_L')

    #lns7 = AX_THETA.plot(THETA_SCM_WRFOUT[0:VV_WRF], MASSLEVELS_SCM_WRFOUT[0:VV_WRF], linewidth = 1,linestyle = 'dashed', c = 'b', label = 'PT SCM')

    lns8 = AX.plot(TEMP_C_WRF_3D_WRFOUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF], linewidth = 0.9, c = 'b', label = 'T/Wind WRF')
    lns9 = AX.plot(T_DAGGPUNKT_3D_WRFOUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF], linewidth = 0.9, c = 'b', linestyle = 'dashed', label = 'DP WRF')

    lns10 = AX.plot(TEMPERATURE_C_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[:VV_EC], linewidth = 0.9, c = 'r', label = 'T/Wind IFS')
    lns11 = AX.plot(T_DAGGPUNKT_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[:VV_EC], linewidth = 0.9, c = 'r', linestyle = 'dashed', label = 'DP IFS')

    
    
    '''För att få legenderna i olika lager i samma legend'''

    #lns = lns1 + lns2 + lns3 + lns4 + lns5 + lns6 + lns7 + lns8 + lns9 + 1ns10 + lns11       # ALLA

    lns = lns1 + lns2 +lns8 + lns9 + lns10 + lns11         # T och Td

    #lns = lns3 + lns4 +lns5 + lns7       # THETA
    
    labs = [l.get_label() for l in lns]

    
    AX.text(-9, 1800,'%s' %(datum[5:8]), fontsize=18, fontweight='bold')      
    #AX.set_ylabel('Height (m)',fontsize=18)
    #AX.set_xlabel('Temperature (C)', fontsize=18)
    AX.set_xlim(-25, -5)      #(-25, -5)
    AX.set_ylim(-60, 2000)
    AX_WIND.set_ylim(-60, 2000)
    AX_WIND.set_xticks([0, 5, 10, 15, 20])
    AX.tick_params(which='major', direction='in', labelsize = '17')
    
    AX_WIND.set_xlim(0,11)
    AX_WIND.tick_params(which='major', direction='in', labelsize= '15')
    AX_THETA.set_xlim(260, 270)
    AX_THETA.tick_params(axis='x',direction="in", labelsize = '17')

    if WRF_TIME == '09':
        AX_THETA.legend(lns, labs, loc='upper left', fontsize = 'x-small')
    elif WRF_TIME == '00':
        AX_THETA.legend(lns[0:6], labs, loc='upper left', fontsize = 'x-large', ncol = 3,  bbox_to_anchor=(0.30, -0.03) )

    #AX_THETA.legend(loc = 2)

    #fig1.suptitle('Temperature, Dewpoint and Potential temperature', fontsize=20)

    #fig1.savefig('/home/sm_marha/FIGURES/FINAL/Sounding_OBSERVATION_WRF_T_Td_180218')  #ÄNDRA


    if r == 0:
        fig2,axs = plt.subplots(2, 4, figsize=(16,8))
        fig2.subplots_adjust(bottom=0.2)
        fig2.subplots_adjust(hspace=0.05)
        fig2.subplots_adjust(wspace=0.05)
        fig2.subplots_adjust(right=0.85)
        
        marker = '.'

    if r == 0:
        cr = 'C0'
    elif r == 1:
        cr = 'C1'
    elif r == 2:
        cr = 'C2'
    else:
        cr = 'C3'
        
       
            

    axs[0,0].plot(QCLOUD_SOD_WRF_3D_WRFOUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF], linestyle = 'dashed',c= cr)
    if r == 0:
        lns00 = axs[0,0].plot(QCLOUD_SOD_WRF_3D_WRFOUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF], linestyle = 'dashed', c= cr, label ='t = 00h')
    if r == 1:
        lns111 = axs[0,0].plot(QCLOUD_SOD_WRF_3D_WRFOUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF], linestyle = 'dashed',c= cr, label ='t = 06h')
    if r == 2:
        lns22 = axs[0,0].plot(QCLOUD_SOD_WRF_3D_WRFOUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF], linestyle = 'dashed',c= cr, label ='t = 12h')
    if r == 3:
        lns33 = axs[0,0].plot(QCLOUD_SOD_WRF_3D_WRFOUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF], linestyle = 'dashed',c= cr, label ='t = 18h')
    
    axs[0,0].set_ylabel('Height ($m$)', fontsize = '18')
    axs[0,0].grid(linestyle="dotted")
    axs[0,0].set_ylim(0, 2000)
    axs[0,0].set_xlim(0, 0.4)
    axs[0,0].set_xticklabels([])
    axs[0,0].text(0.025, 1800,'a', fontsize=18, fontweight='bold')
    axs[0,0].tick_params(which='major', direction='in', labelsize = '17')
    axs[0,0].set_yticks([500, 1000, 1500, 2000])
    axs[0,0].set_yticklabels([500, 1000, 1500, 2000])

    axs[0,1].plot(QICE_SOD_WRF_3D_WRFOUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF],linestyle = 'dashed', c= cr)    
    axs[0,1].grid(linestyle="dotted")
    axs[0,1].set_ylim(0, 2000)
    axs[0,1].set_xlim(0, 1e-1)
    axs[0,1].set_yticklabels([])
    axs[0,1].set_xticklabels([])
    axs[0,1].text(0.006, 1800,'b', fontsize=18, fontweight='bold')
    axs[0,1].tick_params(which='major', direction='in', labelsize = '17')

    axs[0,2].plot(CLOUD_FRACTION_WRF_3D_WRFOUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF],linestyle = 'dashed', c= cr)
    axs[0,2].grid(linestyle="dotted")
    axs[0,2].set_ylim(0, 2000)
    axs[0,2].set_xlim(-0.05, 1.05)
    axs[0,2].set_yticklabels([])
    axs[0,2].set_xticklabels([])
    axs[0,2].text(0.05, 1800,'c', fontsize=18, fontweight='bold')
    axs[0,2].tick_params(which='major', direction='in', labelsize = '17')
    
    axs[0,3].plot(RH_WRF_3D_WRFOUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF], linestyle = 'dashed',c= cr)
    axs[0,3].grid(linestyle="dotted")
    axs[0,3].set_ylim(0, 2000)
    axs[0,3].set_xlim(0, 1.05)
    axs[0,3].set_yticklabels([])
    axs[0,3].set_xticklabels([])
    axs[0,3].text(0.03, 1800,'d', fontsize=18, fontweight='bold')
    axs[0,3].tick_params(which='major', direction='in', labelsize = '17')

    axs[1,0].plot(QCLOUD_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[0:VV_EC], linestyle = 'dashed', c= cr)
    axs[1,0].set_ylabel('Height ($m$)', fontsize = '18')
    axs[1,0].grid(linestyle="dotted")
    axs[1,0].set_ylim(0, 2000)
    axs[1,0].set_yticks([500, 1000, 1500, 2000])
    axs[1,0].set_yticklabels([500, 1000, 1500, 2000])
    axs[1,0].set_xlabel('Cloud water ($gkg^{-1}$)', fontsize = '17')
    axs[1,0].set_xlim(0, 0.4)
    axs[1,0].text(0.025, 1800,'e', fontsize=18, fontweight='bold')
    axs[1,0].tick_params(which='major', direction='in', labelsize = '17')
    
    axs[1,1].plot(QICE_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[0:VV_EC], linestyle = 'dashed', c= cr)
    axs[1,1].grid(linestyle="dotted")
    axs[1,1].set_ylim(0, 2000)
    axs[1,1].set_xlim(0, 1e-1)
    axs[1,1].set_yticklabels([])
    axs[1,1].set_xticks([ 2.5e-2, 5e-2, 7.5e-2, 1e-1])
    axs[1,1].set_xticklabels([ 2.5e-2, 5e-2, 7.5e-2, 1e-1])
    #axs[1,1].xaxis.set_major_formatter(FormatStrFormatter('%8.1e'))
    axs[1,1].set_xlabel('Cloud ice ($gkg^{-1}$)', fontsize = '17')
    axs[1,1].text(0.006, 1800,'f', fontsize=18, fontweight='bold')
    axs[1,1].tick_params(which='major', direction='in', labelsize = '17')

    axs[1,2].plot(CLOUD_FRACTION_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[0:VV_EC], linestyle = 'dashed', c=cr)
    axs[1,2].grid(linestyle="dotted")
    axs[1,2].set_xlim(-0.05, 1.05)
    axs[1,2].set_ylim(0, 2000)
    axs[1,2].set_yticklabels([])
    axs[1,2].set_xticks([0.5, 1.0])
    axs[1,2].set_xticklabels([0.5, 1.0])
    axs[1,2].set_xlabel('Cloud fraction (%)', fontsize = '17')
    axs[1,2].text(0.05, 1800,'g', fontsize=18, fontweight='bold')
    axs[1,2].tick_params(which='major', direction='in', labelsize = '17')

    axs[1,3].plot(RH_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[0:VV_EC], linestyle = 'dashed', c= cr)
    axs[1,3].grid(linestyle="dotted")
    axs[1,3].set_ylim(0, 2000)
    axs[1,3].set_xlim(0, 1.05)
    axs[1,3].set_yticklabels([])
    axs[1,3].set_xticks([0.5, 1.0])
    axs[1,3].set_xticklabels([0.5, 1.0])
    axs[1,3].set_xlabel('RH (%)', fontsize = '17')
    axs[1,3].text(0.03, 1800,'h', fontsize=18, fontweight='bold')
    axs[1,3].tick_params(which='major', direction='in', labelsize = '17')
   


    if r == 1:
        lns44 = axs[0,3].plot(RH_SONDERING[0:VV_SOND], HEIGHT_TER_SONDERING[0:VV_SOND], linewidth = 1.7, c = 'k', label = 'OBS 06h')
        axs[1,3].plot(RH_SONDERING[0:VV_SOND], HEIGHT_TER_SONDERING[0:VV_SOND], linewidth = 1.7, c = 'k', label = 'RH_OBS')

    if r == 2:
        lns55 = axs[0,3].plot(RH_SONDERING[0:VV_SOND], HEIGHT_TER_SONDERING[0:VV_SOND], linewidth = 1.7, linestyle = 'dashed', c = 'k', label = 'OBS 12h')
        axs[1,3].plot(RH_SONDERING[0:VV_SOND], HEIGHT_TER_SONDERING[0:VV_SOND], linewidth = 1.7, c = 'k', linestyle = 'dashed', label = 'RH_OBS')

    '''

    if r == 0:
    
        fig2 = plt.figure(2, figsize=(16,4))

        ax1 = plt.subplot(111)

        ax1.plot(QCLOUD_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[0:VV_EC], c= 'r')
        ax1.scatter(QCLOUD_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[0:VV_EC], c= 'r')
        ax1.plot(QCLOUD_SOD_WRF_3D_WRFOUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF], c= 'b')
        ax1.scatter(QCLOUD_SOD_WRF_3D_WRFOUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF], c= 'b')
        ax1.set_title('Cloud water IFS WRF')
        ax1.set_xlabel('Qcloud (g/kg)')
        ax1.set_ylabel('Height (m)')


        fig3 = plt.figure(3, figsize=(10,8))

        ax1 = plt.subplot(111)

        ax1.plot(QICE_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[0:VV_EC], c= 'r')
        ax1.scatter(QICE_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[0:VV_EC], c= 'r')
        ax1.plot(QICE_SOD_WRF_3D_WRFOUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF], c= 'b')
        ax1.scatter(QICE_SOD_WRF_3D_WRFOUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF], c= 'b')
        ax1.set_title('Cloud ice IFS WRF')
        ax1.set_xlabel('Qice (g/kg)')
        ax1.set_ylabel('Height (m)')


        fig4 = plt.figure(4, figsize=(10,8))

        ax1 = plt.subplot(111)

        ax1.plot(CLOUD_FRACTION_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[0:VV_EC], c= 'r')
        ax1.scatter(CLOUD_FRACTION_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[0:VV_EC], c= 'r')
        ax1.plot(CLOUD_FRACTION_WRF_3D_WRFOUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF], c= 'b')
        ax1.scatter(CLOUD_FRACTION_WRF_3D_WRFOUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF], c= 'b')
        ax1.set_title('Cloud fraction IFS WRF')
        ax1.set_xlabel('Cloud fraction (g/kg)')
        ax1.set_ylabel('Height (m)')


        fig5 = plt.figure(5, figsize=(10,8))

        ax1 = plt.subplot(111)

        ax1.plot(RH_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[0:VV_EC], c= 'r')
        ax1.scatter(RH_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[0:VV_EC], c= 'r')
        ax1.plot(RH_WRF_3D_WRFOUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF], c= 'b')
        ax1.scatter(RH_WRF_3D_WRFOUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF], c= 'b')
        ax1.set_title('RH IFS WRF')
        ax1.set_xlabel('RH (g/kg)')
        ax1.set_ylabel('Height (m)')
    '''
        

        

    r+=1

lns_tva = lns00 + lns111 +lns22 + lns33    
labs_tva = [l.get_label() for l in lns_tva]

#axs[0,0].legend(lns_tva, labs_tva, loc='upper right', fontsize = 'small')
#axs[0,1].legend(lns_tva, labs_tva, loc='upper right', fontsize = 'small')
#axs[0,2].legend(lns_tva, labs_tva, loc='upper right', fontsize = 'small')
lns_tva = lns00 + lns111 +lns22 + lns33 + lns44 + lns55
labs_tva = [l.get_label() for l in lns_tva]
#axs[0,3].legend(lns_tva, labs_tva, loc='lower left', fontsize = 'small')


lns_tva = lns00 + lns111 +lns22 + lns33    
labs_tva = [l.get_label() for l in lns_tva]

#axs[1,0].legend(lns_tva, labs_tva, loc='upper right', fontsize = 'small')
#axs[1,1].legend(lns_tva, labs_tva, loc='upper right', fontsize = 'small')
#axs[1,2].legend(lns_tva, labs_tva, loc='upper right', fontsize = 'small')
lns_tva = lns00 + lns111 +lns22 + lns33 + lns44 + lns55
labs_tva = [l.get_label() for l in lns_tva]
axs[0,3].legend(lns_tva, labs_tva, loc='lower right', fontsize = 'xx-large',  bbox_to_anchor=(1.15,1), ncol = 6)

    
fig1.savefig('/home/sm_marha/FIGURES/FINAL/20180218/Vertical_soundings_IFS_WRF3D_OBS_NEW')
fig2.savefig('/home/sm_marha/FIGURES/FINAL/20180218/Vertical_soundings_CLW_CI_CF_RH_NEW')
    
plt.show()
