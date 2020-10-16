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






fig, axs = plt.subplots(2, 2, figsize=(10,8))
fig.subplots_adjust(hspace=0.4, wspace=0.3)

ax1 = axs[0, 0]
ax2 = axs[0, 1]
ax3 = axs[1, 0]
ax4 = axs[1, 1]

AXS  = [ax1, ax2, ax3, ax4]



WRF_OUTPUT_TIMES = N.array(['2018-02-19_00:00:00', '2018-02-19_06:00:00', '2018-02-19_12:00:00', '2018-02-19_18:00:00'])   #ÄNDRA DATUM

ECMWF_TIMESTEPS = N.array([0, 2, 4, 6])

SOUNDING_TIMES = N.array(['00', '6', '12', '18'])



for AX, WRF_OUTPUT_TIME, ECMWF_TIMESTEP, SOUNDING_TIME in zip(AXS, WRF_OUTPUT_TIMES, ECMWF_TIMESTEPS, SOUNDING_TIMES):



    Filobjekt_GRIB = N4.Dataset('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1ua2018021800_M_D.nc', mode='r')         #ÄNDRA DATUM
    Filobjekt_GRIB_sfc = N4.Dataset('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1sfc2018021800_M_D.nc' , mode='r')   #ÄNDRA DATUM
    Filobjekt_WRF_3D_WRFOUT = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-19_00z_91_LEVELS_QC_QI_42h_T_to_Td_1e-4_and_gt_1e-4_Icloud_2/run/wrfout_d01_%s' %(WRF_OUTPUT_TIME), mode='r')   #ÄNDRA SÖKVÄG TILL RÄTT WRF-KÖRNING
    Filobjekt_WRF_3D_WRF_INPUT = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-19_00z_91_LEVELS_QC_QI_42h_T_to_Td_1e-4_and_gt_1e-4_Icloud_2/run/wrfinput_d01')




    ################
    # ECMWF HYBRID #
    #######################################################################################


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
    QICE_EC_SOD = QICE_EC[ECMWF_TIMESTEP, :, LATITUD, LONGITUD]
    QICE_EC_SOD = QICE_EC_SOD[::-1]

    QCLOUD_EC = Filobjekt_GRIB.variables['clwc'][:]
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

    #################
    # WRF 3D WRFOUT#
    #######################################################################################

    '''Ändra sökvägen för filobjektet i början av skriptet för att ändra tidpunkt, då varje steg i 3D ligger i olika NetCDF-filer.'''

    VV_WRF = 33


    QICE_WRF_3D_WRFOUT = Filobjekt_WRF_3D_WRFOUT.variables['QICE'][0, :, 563, 352]
    QCLOUD_WRF_3D_WRFOUT = Filobjekt_WRF_3D_WRFOUT.variables['QCLOUD'][0, :, 563, 35]

    

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
    T_DAGGPUNKT_WRF_3D_WRFOUT = (243.12*N.log(611.2)-243.12*N.log(RH_WRF_3D_WRFOUT*E_MATTNAD_WRF_3D_WRFOUT)) / (N.log(RH_WRF_3D_WRFOUT*E_MATTNAD_WRF_3D_WRFOUT)-17.62-N.log(611.2))
              
    PH_WRF_3D_WRFOUT = Filobjekt_WRF_3D_WRFOUT.variables['PH'][0, :, 563, 352]
    PHB_WRF_3D_WRFOUT = Filobjekt_WRF_3D_WRFOUT.variables['PHB'][0, :, 563, 352]
    MODELLNIVAHOJD_WRF_3D_WRFOUT = (PH_WRF_3D_WRFOUT+PHB_WRF_3D_WRFOUT)/9.81
    MODELLNIVAHOJD_TER_WRF_3D_WRFOUT = MODELLNIVAHOJD_WRF_3D_WRFOUT[:] -MODELLNIVAHOJD_WRF_3D_WRFOUT[0]


    MASSLEVELS_WRF_3D_WRFOUT = 0.5*(MODELLNIVAHOJD_TER_WRF_3D_WRFOUT[:-1] + MODELLNIVAHOJD_TER_WRF_3D_WRFOUT[1:])




    ##################
    # WRF 3D wrfinput#
    #######################################################################################

    '''Ändra sökvägen för filobjektet i början av skriptet för att ändra tidpunkt, då varje steg i 3D ligger i olika NetCDF-filer.'''

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
              
    PH_WRF_3D_WRFINPUT = Filobjekt_WRF_3D_WRFOUT.variables['PH'][0, :, 563, 352]
    PHB_WRF_3D_WRFINPUT = Filobjekt_WRF_3D_WRFOUT.variables['PHB'][0, :, 563, 352]
    MODELLNIVAHOJD_WRF_3D_WRFINPUT = (PH_WRF_3D_WRFINPUT+PHB_WRF_3D_WRFINPUT)/9.81
    MODELLNIVAHOJD_TER_WRF_3D_WRFINPUT = MODELLNIVAHOJD_WRF_3D_WRFINPUT[:] -MODELLNIVAHOJD_WRF_3D_WRFINPUT[0]


    MASSLEVELS_WRF_3D_WRFINPUT = 0.5*(MODELLNIVAHOJD_TER_WRF_3D_WRFINPUT[:-1] + MODELLNIVAHOJD_TER_WRF_3D_WRFINPUT[1:])





    #######################
    #HÄMTAR SONDERINGSDATA#
    #######################

    j=56    #ÄNDRA (26 FEB)          #ÄNDRA

    



    '''Hämtar sonderingsdata med hjälp av funktionen 'pick dates' i Reading_Dat_Files'''


    MATRIX_2, datum = Picking_values_Sounding_Dat_file_version_2.pick_dates(j, SOUNDING_TIME)


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



    #Rc-params måste sättas innan anrop till subplot, men gäller sedan för alla tills nytt anrop görs.
    plt.rc('xtick', labelsize='x-large')
    plt.rc('ytick', labelsize='x-large')


    AX.grid(linestyle="dotted")

    
    #AX.plot(TEMPERATURE_C_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[:VV_EC],linewidth = 1, c = 'b', label = 'IFS T')
    #AX.plot(T_DAGGPUNKT_EC_SOD[0:VV_EC], FULL_LEVELS_EC_SOD[:VV_EC],linewidth = 1, c = 'b', linestyle = 'dashed', label = 'IFS DP')
    AX.plot(THETA_WRF_3D_WRFOUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF], linewidth = 1, c = 'r', label = 'WRF Theta')
    #AX.plot(TEMP_C_WRF_3D_WRFOUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF], linewidth = 1, c = 'r', label = 'WRF T')
    #AX.plot(T_DAGGPUNKT_WRF_3D_WRFOUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFOUT[0:VV_WRF], linewidth = 1, c = 'r', linestyle = 'dashed', label = 'WRF DP')
    #AX.plot(TEMP_C_SONDERING[0:VV_SOND], HEIGHT_TER_SONDERING[0:VV_SOND], linewidth = 1, c = 'k', label = 'OBS PT')
    #AX.plot(DAGGPUNKT_SONDERING[0:VV_SOND], HEIGHT_TER_SONDERING[0:VV_SOND], linewidth = 1, c = 'k', linestyle = 'dashed', label = 'OBS DP')
    
    AX.plot(THETA_WRF_3D_WRFINPUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFINPUT[0:VV_WRF], linewidth = 1, c = 'g', label = 'WRFIN Theta')
    #AX.plot(T_DAGGPUNKT_WRF_3D_WRFINPUT[0:VV_WRF], MASSLEVELS_WRF_3D_WRFINPUT[0:VV_WRF], linewidth = 1, c = 'g', linestyle = 'dashed', label = 'WRFIN DP')


    AX.set_title('%sz' %(WRF_OUTPUT_TIME[11:13]), fontsize=18)      
    AX.set_ylabel('Height (m)',fontsize=15)
    AX.set_xlabel('Temperature (C)', fontsize=15)
    #AX.set_xlim(-30, -5)


    #AX.legend()

fig.suptitle('Temperature and Dewpoint %s' %(WRF_OUTPUT_TIME[0:10]), fontsize=20)

#fig.savefig('/home/sm_marha/FIGURES/FOKUS_17_27_FEB/19_FEB/Sounding_Intercomparison_EC_WRF_OBS_TC_DP')  #ÄNDRA


plt.show()
