import numpy as N
import scipy.io.netcdf as S
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import os
import datetime
import netCDF4 as N4
from matplotlib.colors import LogNorm



file_list = os.listdir('/nobackup/smhid12/sm_marha/2018-02-27_00z_91_LEVELS_Clouds_Where_Clouds_In_Sounding_YSU_TOPD1_100_Reduced_Stability_in_clouds/run/')    #ÄNDRA


wrfout_list = []

for file in file_list:
    if file.startswith('wrfout'):
        wrfout_list.append(file)

wrfout_list.sort()


wrfout_list_datum = []
EXCH_H_MATRIX = []
EXCH_M_MATRIX = []
THETA_MATRIX = []
RTHRATEN_MATRIX =[]
RTHBLTEN_MATRIX = []
RQCBLTEN_MATRIX = []
RQIBLTEN_MATRIX = []

MODELLNIVAHOJD_TER_MATRIX = []
MASSLEVELS_1D_MINUS_TER_MATRIX = []

TIDSVEKTOR_MINUTER = range(len(wrfout_list[:]))
TIDSVEKTOR_DATUM = []

LATITUDE = 563

LONGITUDE = 352

TIDSLANGD= 1

for file in wrfout_list:          #Ha med [0:37] om det finns fler körningar än 6h
    print (file)


    Filobjekt = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-27_00z_91_LEVELS_Clouds_Where_Clouds_In_Sounding_YSU_TOPD1_100_Reduced_Stability_in_clouds/run/'+file, mode='r')  #ÄNDRA


    PH_1D = Filobjekt.variables['PH'][0, :, LATITUDE, LONGITUDE]
    PHB_1D = Filobjekt.variables['PHB'][0, :, LATITUDE, LONGITUDE]
    MODELLNIVAHOJD_1D = (PH_1D+PHB_1D)/9.81
    MASSLEVELS = 0.5*(MODELLNIVAHOJD_1D[:-1] + MODELLNIVAHOJD_1D[1:])
    TERRAIN_HEIGHT = Filobjekt.variables['HGT'][0, LATITUDE, LONGITUDE]
    MODELLNIVAHOJD_TER_1D = MODELLNIVAHOJD_1D -MODELLNIVAHOJD_1D[0]   
    MASSLEVELS_1D_MINUS_TER = MASSLEVELS - TERRAIN_HEIGHT
    MODELLNIVAHOJD_TER_MATRIX.append(MODELLNIVAHOJD_TER_1D)
    MASSLEVELS_1D_MINUS_TER_MATRIX.append(MASSLEVELS_1D_MINUS_TER)

    THETA_PERT_1D = Filobjekt.variables['T'][0, :, LATITUDE, LONGITUDE]
    THETA_1D = THETA_PERT_1D+300
    THETA_MATRIX.append(THETA_1D)
    
    RTHRATEN_1D = Filobjekt.variables['RTHRATEN'][0, :, LATITUDE, LONGITUDE] *3600
    RTHRATEN_MATRIX.append(RTHRATEN_1D)

    
    RTHBLTEN_1D = Filobjekt.variables['RTHBLTEN'][0, :, LATITUDE, LONGITUDE] *3600
    RTHBLTEN_MATRIX.append(RTHBLTEN_1D)

    RQCBLTEN_1D = Filobjekt.variables['RQCBLTEN'][0, :, LATITUDE, LONGITUDE] 
    RQCBLTEN_MATRIX.append(RQCBLTEN_1D)

    RQIBLTEN_1D = Filobjekt.variables['RQIBLTEN'][0, :, LATITUDE, LONGITUDE]
    RQIBLTEN_MATRIX.append(RQIBLTEN_1D)


    
    EXCH_H_1D = Filobjekt.variables['EXCH_H'][0, :, LATITUDE, LONGITUDE]
    EXCH_H_MATRIX.append(EXCH_H_1D)

    EXCH_M_1D = Filobjekt.variables['EXCH_M'][0, :, LATITUDE, LONGITUDE]
    EXCH_M_MATRIX.append(EXCH_M_1D)
    
    

    wrfout_list_datum.append(file[11:])
    datum = datetime.datetime.strptime(str(file[11:]), '%Y-%m-%d_%H:%M:%S').strftime('%d/%-m %Hz')
    TIDSVEKTOR_DATUM.append(datum)


EXCH_M_MATRIX = N.transpose(EXCH_M_MATRIX)
EXCH_H_MATRIX = N.transpose(EXCH_H_MATRIX)
THETA_MATRIX = N.transpose(THETA_MATRIX)
RTHRATEN_MATRIX = N.transpose(RTHRATEN_MATRIX)
RTHBLTEN_MATRIX = N.transpose(RTHBLTEN_MATRIX)
RQCBLTEN_MATRIX = N.transpose(RQCBLTEN_MATRIX)
RQIBLTEN_MATRIX = N.transpose(RQIBLTEN_MATRIX)

TOTTHTEN_MATRIX = RTHRATEN_MATRIX + RTHBLTEN_MATRIX


##########
#PLOTTING#
##########################################################################
VARIABLE = 'EXCH_H'

fig1 = plt.figure(1, figsize=(12,5))
fig1.subplots_adjust(bottom=0.3)
fig1.subplots_adjust(right=0.90)

ax1 = plt.subplot(111)



'''TURBULENSVARIABLER'''
    
TIME_FULL, LEVELS_FULL = N.meshgrid(TIDSVEKTOR_MINUTER, MODELLNIVAHOJD_TER_1D)   #Behövs för K-värdena på halva nivåerna. Där finns turbulensparametrarna.

TIME, LEVELS = N.meshgrid(TIDSVEKTOR_MINUTER, MASSLEVELS_1D_MINUS_TER)     #För tendenserna på massnivåerna


'''Välj interpolering i fälten med contourf, eller plotta med scatter för att se exakt på vilka modellnivåer det finns moln'''
myplot = ax1.contourf(TIME_FULL, LEVELS_FULL, EXCH_H_MATRIX,[0.01, 0.1, 1, 10, 100], cmap ='binary', norm = LogNorm())                  #Diffusion coefficients          
#myplot = ax1.contourf(TIME, LEVELS, RTHBLTEN_MATRIX ,[-3, -2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2, 3], extend = 'both', cmap ='seismic') #cmap ='binary', norm = LogNorm())        #Tendenser
#myplot = ax1.contourf(TIME, LEVELS, RQCBLTEN_MATRIX ,[-4e-7, -3e-7, -2e-7, -1e-7, 0, 1e-7, 2e-7, 3e-7, 4e-7], extend = 'both', cmap ='seismic') #cmap ='binary', norm = LogNorm())        #Tendenser nederbörd
#myplot = ax3.scatter(TIME_FULL, LEVELS_FULL, c=VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]], cmap='binary'
#myplot = ax3.pcolormesh(TIME_FULL, LEVELS_FULL, VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]], vmin=-0.5, vmax=0.5, cmap = 'seismic')


'''Plottar contours av Theta tillsammans med contourf'''

contour_levels = N.arange(250, 290, 0.5)
cp = plt.contour(TIME, LEVELS, THETA_MATRIX, contour_levels, colors='red', linestyles='dotted', linewidth='0.2' )
plt.clabel(cp, inline=True, fmt='%1.1f', fontsize=8)

    


    
'''Plottar masslevels vid t=0 längs y-axeln''' 
'''X = N.zeros_like(MASSLEVELS[0, :])+3
ax1.scatter(X, MASSLEVELS[0, :], marker= '>', color='black', s=20)
ax1.set_xlim(0, N.max(TIME))'''



ax1.set_ylim(0,2000)

plt.tick_params(labelsize=12)
plt.xticks(TIDSVEKTOR_MINUTER[0:len(TIDSVEKTOR_MINUTER):3], TIDSVEKTOR_DATUM[0:len(TIDSVEKTOR_DATUM):3], rotation= 45, size=6, ha='right', fontsize=12)
ax1.set_ylabel('Height (m)', fontsize=12)
ax1.set_xlabel('Time (UTC)', fontsize=12)
plt.ylim(0,2500)            #ÄNDRA
#plt.xlim(0,360)
#ax1.legend(loc=1,prop={'size': 12})
ax1.grid(linestyle="dotted")
#ax1.set_title(Filobjekt.variables[VARIABLE].description + ' (m2/s)', fontsize=16)

ax1.set_title('Km MYNN ')

'''OBS, detta måste komma efter man använder plt., eftersom inställningarna görs på axes för colorbaren'''
#cbar_ax = fig1.add_axes([0.92, 0.30, 0.01, 0.58])
#cbar = fig1.colorbar(myplot, cax=cbar_ax, ticks =[0.01, 0.1, 1, 10, 100, 100])#, format = "%8.1e")
cbar = fig1.colorbar(myplot)#, format = "%8.1e")


plt.show()
