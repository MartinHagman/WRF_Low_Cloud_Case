import numpy as N
import netCDF4 as N4
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
import pandas as pd
import datetime
from math import *
from pylab import *
#from metpy.plots import ctables




Filobjekt_LIDAR = N4.Dataset('/home/sm_marha/201801-201803_sodankyla_ct25k/20180227_sodankyla_ct25k.nc')


BACKSCATTER_COEFF = Filobjekt_LIDAR.variables['beta'][:]
BC_scale_factor =  Filobjekt_LIDAR.variables['beta'].scale_factor
#BACKSCATTER_COEFF = BACKSCATTER_COEFF*BC_scale_factor
BACKSCATTER_COEFF = N.transpose(BACKSCATTER_COEFF)

RANGE = Filobjekt_LIDAR.variables['range'][:]

TIME = Filobjekt_LIDAR.variables['time'][:]
time = Filobjekt_LIDAR.variables['time']
dates = N4.num2date(TIME, time.units)

ALLA_DATUM = []

for date in dates:
    datum = date.strftime('%d/%-m %H:%M:%S')
    ALLA_DATUM.append(datum)
print(N.shape(ALLA_DATUM))

LAST_DATE = ALLA_DATUM[-1]

ALLA_DATUM = ALLA_DATUM[0:len(ALLA_DATUM):720]

ALLA_DATUM.append(LAST_DATE)

for datum in ALLA_DATUM:
    datum = datum[0:10]



#Tar bort dubletter
'''for place, DATUM in enumerate(ALLA_DATUM):
    print(place, DATUM)
    print(ALLA_DATUM[place-1])
    if DATUM == ALLA_DATUM[place-1]:
        ALLA_DATUM.remove(DATUM)'''
print(ALLA_DATUM)


TIME_MATRIX, RANGE_MATRIX = N.meshgrid(TIME, RANGE)


rc('axes', linewidth=1.5)


fig1 = plt.figure(1, figsize=(10,4))

fig1.subplots_adjust(bottom=0.3)

ax1 = plt.subplot(111)

ax1.tick_params(which = 'major', direction ='in', width = 1   )

#cmap = ctables.registry.get_colortable('NWSReflectivity')

myplot = ax1.contourf(TIME_MATRIX, RANGE_MATRIX, BACKSCATTER_COEFF,  norm=mplc.LogNorm(vmin=BACKSCATTER_COEFF.min(), vmax=BACKSCATTER_COEFF.max()), levels=N.logspace(start=-7,stop=-3, num=100),cmap='jet')

#fig1.tight_layout()
fig1.subplots_adjust(right=0.85)
cbar_ax = fig1.add_axes([0.87, 0.30, 0.01, 0.58])
cbar = fig1.colorbar(myplot, cax=cbar_ax)#, format = "%8.1e")
cbar.set_ticks([10**-7, 10**-6, 10**-5, 10**-4, 10**-3])
#ax1.set_title('Attenuated backscatter coefficient (s/mr)', fontsize=18)
ax1.set_ylabel('Height ($km$)', fontsize=16)
#ax1.set_xlabel('Time (Hours)', fontsize=18)
ax1.set_yticks([0, 0.5, 1.0, 1.5, 2.0])
ax1.set_yticklabels([0, 0.5, 1.0, 1.5, 2.0], fontsize=14)
ax1.set_xticks(N.arange(0, 26, 3))
#ax1.set_xticklabels(ALLA_DATUM[0:len(ALLA_DATUM)],rotation=45, fontsize=12, ha='right')
ax1.set_xticklabels(['27/2 00z', '27/2 03z', '27/2 06z', '27/2 09z', '27/2 12z', '27/2 15z', '27/2 18z', '27/2 21z', '28/2 00z'],rotation=45, fontsize=13)
plt.tick_params(labelsize=15)
ax1.set_ylim(0,1.5)
ax1.grid(linestyle="dotted")

fig1.savefig('/home/sm_marha/FIGURES/FINAL/20180227/Backscatter_20180227')

plt.show()
