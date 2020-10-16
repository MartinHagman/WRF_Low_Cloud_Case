import numpy as N
import netCDF4 as N4
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
import pandas as pd
import datetime
from math import *
import os





file_list = os.listdir('/home/sm_marha/201801-201803_sodankyla_ct25k')

nc_list =[]

for file in file_list:
    if file.endswith('.nc'):
        nc_list.append(file)

nc_list.sort()




BACKSCATTER_COEFF_TIDSSERIE = N.zeros((256, len(nc_list)*5760))
RANGE_TIDSSERIE = N.zeros(len(nc_list)*256)
TIME_TIDSSERIE = N.zeros(len(nc_list)*5760)
ALLA_DATUM_TIDSSERIE = []


for place, file in enumerate(nc_list):
    print(file)
    print(place)

    Filobjekt_LIDAR = N4.Dataset('/home/sm_marha/201801-201803_sodankyla_ct25k/'+file, mode = 'r')


    BACKSCATTER_COEFF = Filobjekt_LIDAR.variables['beta'][:]
    BC_scale_factor =  Filobjekt_LIDAR.variables['beta'].scale_factor
    #BACKSCATTER_COEFF = BACKSCATTER_COEFF*BC_scale_factor
    BACKSCATTER_COEFF = N.transpose(BACKSCATTER_COEFF)

    RANGE = Filobjekt_LIDAR.variables['range'][:]

    TIME = Filobjekt_LIDAR.variables['time'][:]

    '''for iterator, TIME_STEP in enumerate(TIME):
        if TIME_STEP == TIME[iterator-1]:
            TIME.remove(TIME_STEP)'''

    time = Filobjekt_LIDAR.variables['time']
    dates = N4.num2date(TIME, time.units)

    

    ALLA_DATUM = []

    for date in dates:
        datum = date.strftime('%d/%-m %H:%M:%S')
        ALLA_DATUM.append(datum)
    print(len(ALLA_DATUM))

    
    BACKSCATTER_COEFF_TIDSSERIE[:, (place*len(TIME)):(place*5760+5760)] = BACKSCATTER_COEFF
    RANGE_TIDSSERIE[(place*256):(place*256+256)] = RANGE
    TIME_TIDSSERIE[(place*5760):(place*5760+5760)] = TIME+place*24
    ALLA_DATUM_TIDSSERIE.append(ALLA_DATUM)

'''
#Tar bort dubletter
for iterator, DATUM in enumerate(ALLA_DATUM):
    if DATUM == ALLA_DATUM[iterator-1]:
            ALLA_DATUM.remove(DATUM)'''
    

TIME_MATRIX, RANGE_MATRIX = N.meshgrid(TIME_TIDSSERIE, RANGE_TIDSSERIE)


fig1 = plt.figure(1, figsize=(14,6))
fig1.subplots_adjust(bottom=0.3)

ax1 = plt.subplot(111)

#cmap = ctables.registry.get_colortable('NWSReflectivity')

myplot = ax1.contourf(TIME_MATRIX, RANGE_MATRIX, BACKSCATTER_COEFF_TIDSSERIE,  norm=mplc.LogNorm(vmin=BACKSCATTER_COEFF.min(), vmax=BACKSCATTER_COEFF.max()), levels=N.logspace(start=-7,stop=-3, num=100),cmap='jet')

#fig1.tight_layout()
fig1.subplots_adjust(right=0.90)
cbar_ax = fig1.add_axes([0.92, 0.30, 0.01, 0.6])
cbar = fig1.colorbar(myplot, cax=cbar_ax)#, format = "%8.1e")
cbar.set_ticks([10**-7, 10**-6, 10**-5, 10**-4, 10**-3])
ax1.set_ylim(0,2)
ax1.set_title('Attenuated backscatter coefficient (s/mr)', fontsize=18)
ax1.set_ylabel('Height (km)', fontsize=18)
ax1.set_xlabel('Time (Hours)', fontsize=18)
ax1.set_xticks(N.arange(0, 25))
ax1.set_xticklabels(ALLA_DATUM[0:len(ALLA_DATUM)],rotation=45, ha='right', fontsize=15)
plt.tick_params(labelsize=15)
plt.show()

#fig1.savefig('/home/sm_marha/FIGURES/FOKUS_17_27_FEB/18_FEB/20180218_BACKSCATTER')
