import numpy as N
import netCDF4 as N4
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import os
import pandas as pd
import datetime
from math import *
import Picking_values_Sounding_Dat_file_version_2






#######################
#HÄMTAR SONDERINGSDATA#
##########################################################


VV_SOND = 700   #2000

#Definierar arrays med nollor. Fyra kolumner och VV_SOND rader

THETA_SONDERING_MATRIX = N.zeros((VV_SOND,4))

TEMP_C_MATRIX = N.zeros((VV_SOND,4))

ALTITUDE_TEMPERATURE_MATRIX = N.zeros((VV_SOND, 4))

TIME_MATRIX = N.ones((VV_SOND, 4))

DATE_LIST = []


#Endast för plottning

VARIABLE = 'THETA'              #ÄNDRA
DATE = '20180227-28'            #ÄNDRA


Sounding_times = N.array(['00', '6', '12', '18'])

counter = 0
counter_time = 0    #Till för att ge olika avstånd mellan xticks om en sondering saknas. På dessa sätts sedan xticklabels ut från DATE_LIST

for j in range(57, 58):#ÄNDRA   Loopar över dagar. VIKTIGT: Börja på datumet FÖRE det Du vill ha, då
                               #Picking_values_Sounding_Dat_file_version_2.pick_dates börjar på 20180101 
    for i in Sounding_times:    #Loopar över sounding times ovan
        
        '''Hämtar sonderingsdata med hjälp av funktionen 'pick dates' i Picking_values_Sounding_Dat_file_version_2'''


        MATRIX_2, datum = Picking_values_Sounding_Dat_file_version_2.pick_dates(j, i)

        MATRIX_2 = (N.array(MATRIX_2)).astype(float)




        '''Beräknar relativ fuktighet och specifik fuktighet mm mm'''       
        
        

        TEMP_C_SONDERING = MATRIX_2[0:VV_SOND, 4]



        '''Tar bort kolumner ur matrisen när sonderingar saknas.'''
        
        if len(TEMP_C_SONDERING) == 0:
            
            THETA_SONDERING_MATRIX = N.delete(THETA_SONDERING_MATRIX, [counter], axis=1)
            TEMP_C_MATRIX = N.delete(TEMP_C_MATRIX, [counter], axis=1)
            ALTITUDE_TEMPERATURE_MATRIX = N.delete(ALTITUDE_TEMPERATURE_MATRIX, [counter], axis=1)
            TIME_MATRIX = N.delete(TIME_MATRIX, [counter], axis=1)
            counter_time+=1
            print(N.shape(THETA_SONDERING_MATRIX))
                  
            continue   #Hoppar till nästa iterering i for-loopen

        
       
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

        
        '''Spar undan data i olika kolumner i matrisen. Varje kolumn, ett tidssteg.'''
        
        THETA_SONDERING_MATRIX[:, counter] = THETA_SONDERING

        TEMP_C_MATRIX[:, counter] = TEMP_C_SONDERING

        ALTITUDE_TEMPERATURE_MATRIX[:, counter] = HEIGHT_TER_SONDERING

        TIME_MATRIX[:, counter] =  TIME_MATRIX[:, counter] * counter_time  #Fixar så att kolumnerna har värdet 0, 1, 2, 3... alternativ 0, 2 , 3 om t ex 1 saknas.)

        DATE_LIST.append(datum)  


        counter+=1
        counter_time+=1






########################
#LIDAR BACKSCATTER DATA#
##########################################################



#Filobjekt_LIDAR = N4.Dataset('/home/sm_marha/201801-201803_sodankylä_ct25k')





        
        



############
# PLOTTING #
#################################################################################################################

fig1 = plt.figure(1, figsize=(14,6))

ax1 = plt.subplot(111)


contour_levels = N.arange(0, 330, 0.5)      #ÄNDRA
#contour_levels = N.arange(-80, 10, 0.5)    #ÄNDRA


cp = ax1.contour(TIME_MATRIX, ALTITUDE_TEMPERATURE_MATRIX,  THETA_SONDERING_MATRIX, contour_levels, colors='red',  linestyles='dotted', linewidth='0.2' ) #ÄNDRA

ax1.set_xticks(TIME_MATRIX[0, :])
ax1.set_ylim(0,2000)            #ÄNDRA
plt.tick_params(labelsize=15)          #Ändra tick-storlek utan att ändra labels
ax1.set_xticklabels(DATE_LIST, rotation=45, ha='right', fontsize=15)
ax1.set_title('Sounding potential temperature (K)', fontsize=18)
ax1.set_ylabel('Height (m)', fontsize=18)
plt.clabel(cp, inline=True, fmt='%1.1f', fontsize=10)

fig1.tight_layout()
fig1.subplots_adjust(right=0.9)

#plt.savefig('/home/sm_marha/FIGURES/FOKUS_17_27_FEB/18_FEB/SOUNDING_Crossection_2000_%s_%s' % (VARIABLE, DATE))     #ÄNDRA

plt.show()
