import numpy as N
import netCDF4 as N4
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
import pandas as pd
import datetime
from math import *
import os
#from metpy.plots import ctables


file_list = os.listdir('/home/sm_marha/201801-201803_sodankyla_ct25k')

nc_list =[]

for file in file_list:
    if file.endswith('.nc'):
        nc_list.append(file)

nc_list.sort()

nc_list = nc_list[0:59]




fig1 = plt.figure(1, figsize=(20,10))
fig1.subplots_adjust(bottom=0.3)
ax1 = plt.subplot(111)

ALLA_DATUM_ACCUMULATED= []

for place, file in enumerate(nc_list):
    print(file)
    

    Filobjekt_LIDAR = N4.Dataset('/home/sm_marha/201801-201803_sodankyla_ct25k/'+file, mode = 'r')


    BACKSCATTER_COEFF = Filobjekt_LIDAR.variables['beta'][:]
    BC_scale_factor =  Filobjekt_LIDAR.variables['beta'].scale_factor
    #BACKSCATTER_COEFF = BACKSCATTER_COEFF*BC_scale_factor
    BACKSCATTER_COEFF = N.transpose(BACKSCATTER_COEFF)

    RANGE = Filobjekt_LIDAR.variables['range'][:]

    TIME_ACCUMULATED = place*24 + Filobjekt_LIDAR.variables['time'][:]
    TIME = Filobjekt_LIDAR.variables['time'][:]
    time = Filobjekt_LIDAR.variables['time']

    dates = N4.num2date(TIME, time.units)

    ALLA_DATUM = []

    for date in dates:
        datum = date.strftime('%d/%-m')
        ALLA_DATUM.append(datum)
    print(N.shape(ALLA_DATUM))

    ALLA_DATUM = ALLA_DATUM[0:len(ALLA_DATUM):5760]
    ALLA_DATUM_ACCUMULATED.append(ALLA_DATUM[0])
    

    #Tar bort dubletter
    '''for place, DATUM in enumerate(ALLA_DATUM):
        print(place, DATUM)
        print(ALLA_DATUM[place-1])
        if DATUM == ALLA_DATUM[place-1]:
            ALLA_DATUM.remove(DATUM)'''
    #print(ALLA_DATUM)


    TIME_MATRIX, RANGE_MATRIX = N.meshgrid(TIME_ACCUMULATED, RANGE)


    
    

    #cmap = ctables.registry.get_colortable('NWSReflectivity')

    myplot = ax1.contourf(TIME_MATRIX, RANGE_MATRIX, BACKSCATTER_COEFF,  norm=mplc.LogNorm(vmin=BACKSCATTER_COEFF.min(), vmax=BACKSCATTER_COEFF.max()), levels=N.logspace(start=-7,stop=-3, num=100),cmap='jet')

    #fig1.tight_layout()
    fig1.subplots_adjust(right=0.90)
    cbar_ax = fig1.add_axes([0.92, 0.30, 0.01, 0.6])
    cbar = fig1.colorbar(myplot, cax=cbar_ax)#, format = "%8.1e")
    cbar.set_ticks([10**-7, 10**-6, 10**-5, 10**-4, 10**-3])
    ax1.set_ylim(0,2)
    ax1.set_title('Attenuated backscatter coefficient (s/mr)', fontsize=18)
    ax1.set_ylabel('Height (km)', fontsize=18)
    ax1.set_xlabel('Time (Hours)', fontsize=18)
    
ALLA_DATUM_ACCUMULATED = ALLA_DATUM_ACCUMULATED[::2]     

ax1.set_xticks(N.arange(0,(TIME_ACCUMULATED[-1]),48))
ax1.set_xticklabels(ALLA_DATUM_ACCUMULATED,rotation=45, fontsize=10)#, ha='right',)
#plt.tick_params(labelsize=15)
plt.show()

fig1.savefig('/home/sm_marha/FIGURES/FOKUS_17_27_FEB/18_FEB/BACKSCATTER_WHOLE_PERIOD')
