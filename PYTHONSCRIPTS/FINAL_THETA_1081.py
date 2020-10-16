import numpy as N
import scipy.io.netcdf as S
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import os
import datetime
import matplotlib.colors as mplc
import netCDF4 as N4
import Picking_values_Sounding_Dat_file_version_2
import pandas as pd
from math import *






#############
# SCM MODEL #
###########################################################


fig1 = plt.figure(1, figsize=(6, 12))

ax = plt.subplot(111)


fig2 = plt.figure(2, figsize=(6, 12))

ax2 = plt.subplot(111)



KORNING = ['YSU', 'MYNN', 'MYNN+SODTEMP']

i = 0

FILOBJEKTLISTA = ['/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfout_d01_2018-02-27_00:00:00_025T_LIQ_AND_VARIATION_YSU_TD0',
             '/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfout_d01_2018-02-27_00:00:00_025T_LIQ_AND_VARIATION_MYNN',
             '/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfout_d01_2018-02-27_00:00:00_T_LIQ_MYNN_SOUNDING_SODANKYLA_QVAPOR_100']




for FILOBJEKT in FILOBJEKTLISTA:


    Filobjekt = N4.Dataset(FILOBJEKT)

    TIDSLANGD = 2701   #7561 totalt....      540, 1080, 1620, 2160, 2701, 3241, 3781, 4321, 4861.......osv
    TIDPUNKT = 1081

       

    TIMESTEPS = Filobjekt.variables['XTIME'][0:TIDSLANGD]
    print(N.shape(TIMESTEPS))
    time = Filobjekt.variables['XTIME']
    dates = N4.num2date(TIMESTEPS, time.units)


    ALLA_DATUM=[]

    for date in dates:
        datum = date.strftime('%d/%-m %Hz')
        ALLA_DATUM.append(datum)

    ALLA_DATUM_3h = ALLA_DATUM[0:len(ALLA_DATUM):540]




    PERT_THETA_2D = Filobjekt.variables['T'][0:TIDSLANGD, :, 1, 1]

    P_PERT_2D = Filobjekt.variables['P'][0:TIDSLANGD, :, 1, 1]
    P_BASE_2D = Filobjekt.variables['PB'][0:TIDSLANGD, :, 1, 1]
    P_2D= P_BASE_2D+P_PERT_2D

    THETA_2D = PERT_THETA_2D+300

    TEMP_2D = THETA_2D/((1000/(P_2D/100))**0.286)
    TEMP_C_2D = TEMP_2D-273.15
    TEMP_C = TEMP_C_2D[:, 0]
    T2M = Filobjekt.variables['T2'][0:TIDSLANGD, 1, 1]-273.15
    T2M[0]= T2M[1]
    QVAPOR_2D = Filobjekt.variables['QVAPOR'][0:TIDSLANGD, :, 1, 1]
    E_MATTNAD_2D = N.exp(N.log(611.2)+(17.62*TEMP_C_2D/(243.12+TEMP_C_2D)))
    QVAPOR_MATTNAD_2D = 0.622*E_MATTNAD_2D/(P_2D-E_MATTNAD_2D)
          


    RH_2D = QVAPOR_2D/QVAPOR_MATTNAD_2D
    RH_2D = N.transpose(RH_2D)

    TSK = Filobjekt.variables['TSK'][0:TIDSLANGD, 1, 1]-273.15
    EMISS = Filobjekt.variables['EMISS'][0:TIDSLANGD, 1, 1]

    CLDFRA = Filobjekt.variables['CLDFRA'][0:TIDSLANGD, 0, 1, 1]
    CLDFRA_2D = Filobjekt.variables['CLDFRA'][0:TIDSLANGD, :, 1, 1]


    QCLOUD = Filobjekt.variables['QCLOUD'][0:TIDSLANGD, 0, 1, 1]*1000
    QCLOUD_2D = Filobjekt.variables['QCLOUD'][0:TIDSLANGD, :, 1, 1]*1000


    QRAIN = Filobjekt.variables['QRAIN'][0:TIDSLANGD, 0, 1, 1]*1000
    QRAIN_2D = Filobjekt.variables['QRAIN'][0:TIDSLANGD, :, 1, 1]*1000


    QSNOW = Filobjekt.variables['QSNOW'][0:TIDSLANGD, 0, 1, 1]*1000
    QSNOW_2D = Filobjekt.variables['QSNOW'][0:TIDSLANGD, :, 1, 1]*1000


    QICE = Filobjekt.variables['QICE'][0:TIDSLANGD, 0, 1, 1]*1000
    QICE_2D = Filobjekt.variables['QICE'][0:TIDSLANGD, :, 1, 1]*1000


    QGRAUP = Filobjekt.variables['QGRAUP'][0:TIDSLANGD, 0, 1, 1]
    QGRAUP_2D = Filobjekt.variables['QGRAUP'][0:TIDSLANGD, :, 1, 1]

                                                                                                                 


    RTHRATEN = Filobjekt.variables['RTHRATEN'][0:TIDSLANGD, 0, 1, 1]
    RTHRATEN_2D = Filobjekt.variables['RTHRATEN'][0:TIDSLANGD, :, 1, 1]
    RTHRATEN_2D = N.transpose(RTHRATEN_2D)

    RTHBLTEN = Filobjekt.variables['RTHBLTEN'][0:TIDSLANGD, 0, 1, 1]
    RTHBLTEN_2D = Filobjekt.variables['RTHBLTEN'][0:TIDSLANGD, :, 1, 1]
    RTHBLTEN_2D = N.transpose(RTHBLTEN_2D)

    EXCH_H = Filobjekt.variables['EXCH_H'][0:TIDSLANGD, 0, 1, 1]
    EXCH_H_2D = Filobjekt.variables['EXCH_H'][0:TIDSLANGD, :, 1, 1]
    EXCH_H_2D = N.transpose(EXCH_H_2D)

    EXCH_M = Filobjekt.variables['EXCH_M'][0:TIDSLANGD, 0, 1, 1]
    EXCH_M_2D = Filobjekt.variables['EXCH_M'][0:TIDSLANGD, :, 1, 1]
    EXCH_M_2D = N.transpose(EXCH_M_2D)
    '''
    TKE_PBL = Filobjekt.variables['TKE_PBL'][0:TIDSLANGD, 0, 1, 1]
    TKE_PBL_2D = Filobjekt.variables['TKE_PBL'][0:TIDSLANGD, :, 1, 1]
    TKE_PBL_2D = N.transpose(TKE_PBL_2D)

    DTKE = Filobjekt.variables['DTKE'][0:TIDSLANGD, 0, 1, 1]                   #ÄNDRA TKE
    QSHEAR = Filobjekt.variables['QSHEAR'][0:TIDSLANGD, 0, 1, 1]
    QBUOY = Filobjekt.variables['QBUOY'][0:TIDSLANGD, 0, 1, 1]
    QDISS = Filobjekt.variables['QDISS'][0:TIDSLANGD, 0, 1, 1]
    QWT = Filobjekt.variables['QWT'][0:TIDSLANGD, 0, 1, 1]
    '''



    QTOTAL_2D = QVAPOR_2D + QCLOUD_2D + QICE_2D + QRAIN_2D + QSNOW_2D
    POISSON_C = 0.2854*(1-0.24*QVAPOR_2D)
    GAMMA = QVAPOR_2D*461/(1005.7+QTOTAL_2D*1.996e3)


    THETA_L_2D = THETA_2D - ((2.26e6*THETA_2D)/(1004*TEMP_2D)* QCLOUD_2D)

        
    '''Väntar med att transponera dessa tills Theta_L har beräknats'''
    #THETA_2D = N.transpose(THETA_2D)
    THETA_L_2D = N.transpose(THETA_L_2D)
    CLDFRA_2D = N.transpose(CLDFRA_2D)  #Transponerar för att det ska passa TIME och LEVELS nedan
    QCLOUD_2D = N.transpose(QCLOUD_2D)
    QRAIN_2D = N.transpose(QRAIN_2D)
    QSNOW_2D = N.transpose(QSNOW_2D)
    QICE_2D = N.transpose(QICE_2D)
    QGRAUP_2D = N.transpose(QGRAUP_2D)



                                                                                                                    

    RAINNC = Filobjekt.variables['RAINNC'][0:TIDSLANGD, 1, 1]
    SNOWNC = Filobjekt.variables['SNOWNC'][0:TIDSLANGD, 1, 1]



    GLW = Filobjekt.variables['GLW'][0:TIDSLANGD, 1, 1]
    SWDOWN = Filobjekt.variables['SWDOWN'][0:TIDSLANGD, 1, 1]
    GRDFLX = Filobjekt.variables['GRDFLX'][0:TIDSLANGD, 1, 1]
    HFX = Filobjekt.variables['HFX'][0:TIDSLANGD, 1, 1]
    LH = Filobjekt.variables['LH'][0:TIDSLANGD, 1, 1]
    EMISS = Filobjekt.variables['EMISS'][0:TIDSLANGD, 1, 1]
    ALBEDO = Filobjekt.variables['ALBEDO'][0:TIDSLANGD, 1, 1]
    NOAHRES = Filobjekt.variables['NOAHRES'][0:TIDSLANGD, 1, 1]
    FLX1 = Filobjekt.variables['FLX1'][0:TIDSLANGD, 1, 1]
    FLX2 = Filobjekt.variables['FLX2'][0:TIDSLANGD, 1, 1]
    FLX3 = Filobjekt.variables['FLX3'][0:TIDSLANGD, 1, 1]
        
    SWUP = ALBEDO*SWDOWN
    GUPLW = EMISS*5.67e-8*(TSK+273.15)**4
    GUPLW[0] = 0
    GLW[0] = GLW[1]

    LONGW_RESIDUAL = EMISS*GLW - GUPLW
    HFX = -HFX
    SWUP = -SWUP



    U = Filobjekt.variables['U'][:]
    V = Filobjekt.variables['V'][:]
    U_10 = Filobjekt.variables['U10'][:]
    V_10 = Filobjekt.variables['V10'][:]
    XLAT = Filobjekt.variables['XLAT'][:]
    XLONG = Filobjekt.variables['XLONG'][:]
    X_LAT_U = Filobjekt.variables['XLAT_U'][:]
    X_LAT_U = Filobjekt.variables['XLAT_U'][:]
    X_LON_U = Filobjekt.variables['XLONG_U'][:]
    X_LAT_V = Filobjekt.variables['XLAT_V'][:]
    X_LONG_U = Filobjekt.variables['XLONG_V'][:]
    SINALPHA = Filobjekt.variables['SINALPHA'][:]
    COSALPHA = Filobjekt.variables['COSALPHA'][:]

    SINALPHA = SINALPHA[TIDPUNKT, :, :]      #Detta är ett tillägg i SCM jämfört med 3D
    COSALPHA = COSALPHA[TIDPUNKT, :, :]      #Detta är ett tillägg i SCM jämfört med 3D'''

    U_3D = U[2160, :, :, :] 
    V_3D = V[2160, :, :, :]

    U_3D_UNSTAG = 0.5*(U_3D[:,:,:-1] + U_3D[:,:,1:])
    V_3D_UNSTAG = 0.5*(V_3D[:,:-1,:] + V_3D[:,1:,:])

    '''print(N.shape(U_3D_UNSTAG))
    print(N.shape(V_3D_UNSTAG))'''

    U_3D_UNSTAG_EARTH = U_3D_UNSTAG*COSALPHA - V_3D_UNSTAG*SINALPHA
    V_3D_UNSTAG_EARTH = V_3D_UNSTAG*COSALPHA + U_3D_UNSTAG*SINALPHA

    U_SOD = U_3D_UNSTAG_EARTH[:, 1, 1]
    V_SOD = V_3D_UNSTAG_EARTH[:, 1, 1]

    VELOCITY_SOD = N.sqrt(U_SOD**2+V_SOD**2)

    RESIDUAL = EMISS*GLW + SWDOWN + SWUP - GUPLW + HFX - LH + GRDFLX - FLX1 - FLX2 - FLX3



    ZNW = Filobjekt.variables['ZNW'][0:TIDSLANGD, :]

    P = Filobjekt.variables['P'][0:TIDSLANGD, :, 1, 1]

    PB = Filobjekt.variables['PB'][0:TIDSLANGD, :, 1, 1]



    PH = Filobjekt.variables['PH'][:, :, 1, 1]
    PHB = Filobjekt.variables['PHB'][:, :, 1, 1]
    MODELLNIVAHOJD = (PH+PHB)/9.81
    MODELLNIVAHOJD_TER = MODELLNIVAHOJD[0:TIDSLANGD, :] -MODELLNIVAHOJD[0, 0]


    MASSLEVELS = 0.5*(MODELLNIVAHOJD_TER[0:TIDSLANGD, :-1] + MODELLNIVAHOJD_TER[0:TIDSLANGD, 1:])  #Olika masslevels för varje tidssteg, om än mycket små skillnader!

    ax.plot( THETA_2D[1081, 0:30],MASSLEVELS[1081, 0:30], label = KORNING[i])
    ax.legend()

    ax2.plot(VELOCITY_SOD[0:30], MASSLEVELS[0, 0:30], label = KORNING[i])
    ax2.legend()

    

    i+=1


#############
# SOUNDINGS #
###########################################################


j=58    #ÄNDRA (25 FEB)          #ÄNDRA

    
SOUNDING_TIME = '00'


'''Hämtar sonderingsdata med hjälp av funktionen 'pick dates' i Reading_Dat_Files'''


MATRIX_2, datum = Picking_values_Sounding_Dat_file_version_2.pick_dates(j, SOUNDING_TIME)


MATRIX_2 = (N.array(MATRIX_2)).astype(float)




'''Beräknar relativ fuktighet och specifik fuktighet'''



VV_SOND = 530  # 470   #2000



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

ax.plot(THETA_SONDERING[0:330], HEIGHT_TER_SONDERING[0:330],linewidth = 1.2, c = 'k', label = 'OBS' )
ax2.plot(WIND_VELOCITY_SONDERING[0:330], HEIGHT_TER_SONDERING[0:330],linewidth = 1.2, c = 'k', label = 'OBS' )

ax.grid(linestyle="dotted")
ax.set_xlabel('Potential temperature (K)', fontsize = 14)
ax.set_ylabel('Height (m)', fontsize = 14)
ax.set_xlim(263, 266)
ax.set_title('Potential temperature', fontsize = 14)
ax.legend()

ax2.grid(linestyle="dotted")
ax2.set_xlabel('Wind (m/s)', fontsize = 14)
ax2.set_ylabel('Height (m)', fontsize = 14)
#ax2.set_xlim(263, 266)
ax2.set_title('Wind velocity', fontsize = 14)
ax2.legend()
ax2.plot


plt.show()
