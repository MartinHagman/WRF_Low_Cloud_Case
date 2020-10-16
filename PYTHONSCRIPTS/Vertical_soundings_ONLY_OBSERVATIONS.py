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





print(ax1.get_position())
box_ax1 =  ax1.get_position()

ax1.set_position([box_ax1.x0, box_ax1.y0, box_ax1.width*1.6, box_ax1.height])
print("\n")


print(ax1w.get_position())
box_ax1w =  ax1w.get_position()

ax1w.set_position([0.407, box_ax1w.y0, box_ax1w.width*0.5, box_ax1w.height])
ax1w.xaxis.set_visible(False)
ax1w.yaxis.set_visible(False)
ax1w.set_frame_on(False)
print("\n")


print(ax2.get_position())
box_ax2 =  ax2.get_position()

ax2.set_position([box_ax2.x0, box_ax2.y0, box_ax2.width*1.6, box_ax2.height])
print("\n")


print(ax2w.get_position())
box_ax2w =  ax2w.get_position()

ax2w.set_position([0.811, box_ax2w.y0, box_ax2w.width*0.5, box_ax2w.height])
ax2w.xaxis.set_visible(False)
ax2w.yaxis.set_visible(False)
ax2w.set_frame_on(False)
print("\n")


print(ax3.get_position())
box_ax3 =  ax3.get_position()

ax3.set_position([box_ax3.x0, box_ax3.y0, box_ax3.width*1.6, box_ax3.height])
print("\n")


print(ax3w.get_position())
box_ax3w =  ax3w.get_position()

ax3w.set_position([0.407, box_ax3w.y0, box_ax3w.width*0.5, box_ax3w.height])
ax3w.xaxis.set_visible(False)
ax3w.yaxis.set_visible(False)
ax3w.set_frame_on(False)
print("\n")


print(ax4.get_position())
box_ax4 =  ax4.get_position()

ax4.set_position([box_ax4.x0, box_ax4.y0, box_ax4.width*1.6, box_ax4.height])
print("\n")


print(ax4w.get_position())
box_ax4w =  ax4w.get_position()
ax4w.xaxis.set_visible(False)
ax4w.yaxis.set_visible(False)
ax4w.set_frame_on(False)

ax4w.set_position([0.811, box_ax4w.y0, box_ax4w.width*0.5, box_ax4w.height])


print("\n")







AXS  = [ax1, ax2, ax3, ax4]

AX_WINDS = [ax1w, ax2w, ax3w, ax4w]


ax11 = ax1.twiny()
ax22 = ax2.twiny()
ax33 = ax3.twiny()
ax44 = ax4.twiny()


AX_THETAS = [ax11, ax22, ax33, ax44]




#ax11.spines['right'].set_visible(False)
#ax11.set_frame_on(False)
ax11.set_position([box_ax1.x0, box_ax1.y0, box_ax1.width*1.6, box_ax1.height])
ax11.set_position([box_ax1.x0, box_ax1.y0, box_ax1.width*1.6, box_ax1.height])
ax11.xaxis.set_visible(False)
ax11.yaxis.set_visible(False)
ax11.set_frame_on(False)


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




SOUNDING_TIMES = N.array(['00', '6', '12', '18'])
#SOUNDING_TIMES = N.array(['00', '12'])



for AX, AX_WIND, AX_THETA, SOUNDING_TIME in zip(AXS, AX_WINDS, AX_THETAS, SOUNDING_TIMES):




    
  
    #######################
    #HÄMTAR SONDERINGSDATA#
    #######################

    j=51    #ÄNDRA (25 FEB)          #ÄNDRA

    



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

    WIND_VELOCITY_SONDERING = MATRIX_2[0:VV_SOND, 0]

    WIND_VELOCITY_SONDERING = MATRIX_2[0:VV_SOND, 0]

    WIND_DIRECTION_SONDERING = MATRIX_2[0:VV_SOND, 2]

    U_WIND_SONDERING = -WIND_VELOCITY_SONDERING*N.sin(pi/180*WIND_DIRECTION_SONDERING)

    V_WIND_SONDERING = -WIND_VELOCITY_SONDERING*N.cos(pi/180*WIND_DIRECTION_SONDERING)

    ALTITUDE_WIND_SONDERING = MATRIX_2[0:VV_SOND, 1]


###########
#PLOTTNING#
######################################################################



    #Rc-params måste sättas innan anrop till subplot, men gäller sedan för alla tills nytt anrop görs.
    plt.rc('xtick', labelsize='x-large')
    plt.rc('ytick', labelsize='x-large')


    AX.grid(linestyle="dotted")
    '''Lns1 osv för att få legenderna i samma lager'''
    lns1 = AX.plot(TEMP_C_SONDERING[0:VV_SOND], HEIGHT_TER_SONDERING[0:VV_SOND], linewidth = 1, c = 'b', label = 'T')
    lns2 = AX.plot(DAGGPUNKT_SONDERING[0:VV_SOND], HEIGHT_TER_SONDERING[0:VV_SOND], linewidth = 1, c = 'b', linestyle = 'dashed', label = 'DP')
    AX_WIND.barbs(0.5*N.ones(VV_SOND-15)[0::18], HEIGHT_TER_SONDERING[0:VV_SOND-15:18], U_WIND_SONDERING[0:VV_SOND-15:18], V_WIND_SONDERING[0:VV_SOND-15:18], barbcolor = 'k')
    lns3 = AX_THETA.plot(THETA_SONDERING[0:VV_SOND], HEIGHT_TER_SONDERING[0:VV_SOND], linewidth = 1, c = 'r', label = 'PT')


    '''För att få legenderna i olika lager i samma legend'''
    lns = lns1 + lns2 + lns3
    labs = [l.get_label() for l in lns]
    
    AX.set_title('%s' %(datum), fontsize=18)      
    #AX.set_ylabel('Height (m)',fontsize=15)
    #AX.set_xlabel('Temperature (C)', fontsize=15)
    AX.set_xlim(-35, 0)
    AX.set_ylim(-60, 2800)
    AX_WIND.set_ylim(-60, 2800)
    AX_THETA.set_xlim(245, 285)
    AX_THETA.tick_params(axis='x',direction="in")
    AX_THETA.legend(lns, labs, loc='upper left', fontsize = 'x-small')
    #AX_THETA.legend(loc = 2)

fig1.suptitle('Temperature, Dewpoint and Potential temperature', fontsize=20)

#fig1.savefig('/home/sm_marha/FIGURES/FOKUS_17_27_FEB/Sounding_OBSERVATION_180227')  #ÄNDRA


plt.show()
