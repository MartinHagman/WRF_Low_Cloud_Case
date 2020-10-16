import numpy as N
import scipy.io.netcdf as S
import netCDF4 as N4
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import os
import pandas as pd
import datetime
from math import *
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize

#import matplotlib.ticker as ticker

'''Detta skript plottar ECMWF tvärsnitt 2 variabler i taget. Två saker ändras, variabel i plottningskommandot samt ändring
av variabel som string för titeln.'''

Filobjekt_GRIB = N4.Dataset('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1ua2018022000_M_D.nc', mode='r')    #ÄNDRA
Filobjekt_GRIB_sfc = N4.Dataset('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1sfc2018022000_M_D.nc' , mode='r')   #ÄNDRA


VV_EC = 35

ANTAL_TIDSSTEG = 14

FULL_LEVELS_EC_SOD_alla_tidssteg = N.zeros((VV_EC, ANTAL_TIDSSTEG))
THETA_EC_2D_SOD = N.zeros((VV_EC, ANTAL_TIDSSTEG ))               #ÄNDRA VID 2518 vs 2600 (ANTAL_TIDSSTEG-2)        




LATITUD = 47 # Har värdet 0 uppe 
LONGITUD = 186 # Har värdet 0 på vänstra randen. ECMWF-området är längst i x-led till skillnad från WRF-området,som är längst i y-led.

TIMESTEPS = Filobjekt_GRIB.variables['time'][:]
time = Filobjekt_GRIB.variables['time']     #OBS, ingen array
dates = N4.num2date(TIMESTEPS, time.units, time.calendar)
ALLA_DATUM=[]


'''Gör om tiden till ett visst datumformat.'''

for date in dates:
    datum = date.strftime('%d/%-m %Hz')
    ALLA_DATUM.append(datum)


LONGITUDE = Filobjekt_GRIB.variables['longitude'][:]
LATITUDE = Filobjekt_GRIB.variables['latitude'][:]

lon_matrix, lat_matrix = N.meshgrid(LONGITUDE, LATITUDE)



'''Variabler som finns på varje nivå och vid alla tidpunkter direkt från början'''

QICE_EC_2D = Filobjekt_GRIB.variables['ciwc'][0:ANTAL_TIDSSTEG, -VV_EC:, LATITUD, LONGITUD]*1000  #ÄNDRA VID 2518 vs 2600, (2:ANTAL_TIDSSTEG)      #Shape=(15, VV_EC)). Plockar de VV_EC sista vertikala nivåerna bakifrån, fortfarande felvända.
QICE_EC_2D = N.transpose(QICE_EC_2D)    #shape=(VV_EC, 15)
QICE_EC_2D = QICE_EC_2D[::-1, :]     #Vänder på alla kolumner, som är de utvalda VV_EC-värdena. På så sätt blir värdena vid marken också de första värdena i varje kolumn.
                                     #Shapen blir då fortsatt (VV_EC, Antal tidpunkter)

QCLOUD_EC_2D = Filobjekt_GRIB.variables['clwc'][0:ANTAL_TIDSSTEG, -VV_EC:, LATITUD, LONGITUD]*1000      #ÄNDRA VID 2518 vs 2600, (2:ANTAL_TIDSSTEG)
QCLOUD_EC_2D = N.transpose(QCLOUD_EC_2D)
QCLOUD_EC_2D = QCLOUD_EC_2D[::-1, :]


CLOUD_FRACTION_EC_2D = Filobjekt_GRIB.variables['cc'][0:ANTAL_TIDSSTEG, -VV_EC:, LATITUD, LONGITUD]     #ÄNDRA VID 2518 vs 2600, (2:ANTAL_TIDSSTEG)
CLOUD_FRACTION_EC_2D = N.transpose(CLOUD_FRACTION_EC_2D)
CLOUD_FRACTION_EC_2D = CLOUD_FRACTION_EC_2D[::-1, :]





'''Variabler som ska loopas för varje tidssteg för att de ska användas vid varje tidsssteg för att beräkna nya variabler,
som i sin tur behövs när de vertikala nivåerna vid varje tidpunkt skall beräknas.'''



for TIDSSTEG in range(ANTAL_TIDSSTEG):         #ÄNDRA VID 2518 vs 2600, (2,ANTAL_TIDSSTEG)

    print('TIDSSTEG= ' + str(TIDSSTEG))
    QICE_EC = Filobjekt_GRIB.variables['ciwc'][:]*1000
    QICE_EC_SOD = QICE_EC[TIDSSTEG, -VV_EC:, LATITUD, LONGITUD]
    QICE_EC_SOD = QICE_EC_SOD[::-1]         #Allt inverteras så att lägsta nivån hamnar längst ner i Arrayen.


    QCLOUD_EC = Filobjekt_GRIB.variables['clwc'][:]*1000
    QCLOUD_EC_SOD = QCLOUD_EC[TIDSSTEG, -VV_EC:, LATITUD, LONGITUD]
    QCLOUD_EC_SOD = QCLOUD_EC_SOD[::-1]


    CLOUD_FRACTION_EC = Filobjekt_GRIB.variables['cc'][:]
    CLOUD_FRACTION_EC_SOD = CLOUD_FRACTION_EC[TIDSSTEG, -VV_EC::, LATITUD, LONGITUD]
    CLOUD_FRACTION_EC_SOD = CLOUD_FRACTION_EC_SOD[::-1]


    TEMPERATURE_EC = Filobjekt_GRIB.variables['t'][:]
    TEMPERATURE_EC_SOD = TEMPERATURE_EC[TIDSSTEG, -VV_EC:, LATITUD, LONGITUD]
    TEMPERATURE_EC_SOD = TEMPERATURE_EC_SOD[::-1]
    TEMPERATURE_C_EC_SOD = TEMPERATURE_EC_SOD -273.15
        
    TEMPERATURE_C_EC = TEMPERATURE_EC-273.15     #Denna är inte inverterad i höjdled här och fortfarande 4-dimensionell!!


    SPECIFIC_HUMIDITY_EC = Filobjekt_GRIB.variables['q'][:]
    SPECIFIC_HUMIDITY_EC_SOD = SPECIFIC_HUMIDITY_EC[TIDSSTEG, -VV_EC:, LATITUD, LONGITUD]
    SPECIFIC_HUMIDITY_EC_SOD = SPECIFIC_HUMIDITY_EC_SOD[::-1]

    TIMESTEPS_EC = Filobjekt_GRIB.variables['time'][:]
    VERTICAL_LEVELS_EC = Filobjekt_GRIB.variables['level'][:]  #Ger vektor med värden 1-138 på plats 0 till 137




    '''markgeopotentialen och marktrycket från surface-gribfilen'''

    SURFACE_GEOPOTENTIAL_EC = Filobjekt_GRIB_sfc.variables['z'][:]
    SURFACE_GEOPOTENTIAL_EC_SOD = SURFACE_GEOPOTENTIAL_EC[TIDSSTEG, LATITUD, LONGITUD]  #Endimensionell, MEN EJ inverterad

    SURFACE_PRESSURE_EC = Filobjekt_GRIB_sfc.variables['sp'][:]
    SURFACE_PRESSURE_EC_SOD = SURFACE_PRESSURE_EC[TIDSSTEG, LATITUD, LONGITUD]  #Endimensionell, MEN EJ inverterad




    '''BERÄKNING AV DE VERTIKALA NIVÅERNA I ECMWF för varje tidssteg i GRIB-filen'''
    

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

    T_VIRTUAL_EC = (1.+0.609133*SPECIFIC_HUMIDITY_EC)*TEMPERATURE_EC       # Beräknas med 4-dimensionella vektorer.....

    T_VIRTUAL_EC_SOD = T_VIRTUAL_EC[TIDSSTEG,:, LATITUD, LONGITUD]      # Först här görs nya variabeln endimensionell.


    T_VIRTUAL_EC_SOD = T_VIRTUAL_EC_SOD[::-1] #Inverterar vektorn så nivå 137 kommer först)




    '''Räknar ut de halva trycknivåerna för Sodankylä'''

    p_half = N.zeros(len(A))

    for i in range(len(A)):

        p = A[i] + B[i] * SURFACE_PRESSURE_EC[TIDSSTEG,LATITUD,LONGITUD]  #Ett tryck för Sodankylä i marknivå, men då vi itererar över i blir det en kolumn tryck på halva nivåer(138 stycken)

        p_half[i] = p


    '''Räknar ut geopotentialen på de halva nivåerna



    for i in range(len(p_half)-1):

        Fi_half[i+1] = Fi_half[i] + Rd*(  T_VIRTUAL_SOD[i]  *  (  log(p_half[i])-log(p_half[i+1])  ) )'''
     


    '''Räknar ut hela trycknivåerna för Sodankylä'''
             
    p_full = 0.5*(p_half[1:]+p_half[:-1])   #Fulla tryck beräknas för kolumnen över Sodankylä(137 st)



    '''Räknar ut geopotentialen på de fulla nivåerna över Sodankylä.'''

    Rd = 287.06
    g0 = 9.80665

    #TEMP_v = p_full/(Rd*Density) #Används till att kontrollera algorithmen via ecmwf.int


    Fi = N.zeros(len(p_full))
    Fi[0]=10*g0

    #FI = N.zeros(len(p_full))  #Används till att kontrollera algorithmen via ecmwf.int
    #FI[0]=10*g0

    for i in range(len(p_full)-1):

        Fi[i+1] = Fi[i] + Rd*(  T_VIRTUAL_EC_SOD[i]  *  (  log(p_full[i])-log(p_full[i+1])  ) )
                              
        #FI[i+1] = FI[i] + Rd*(  PF[i]/(Rd*Density[i])  *  (  log(PF[i])-log(PF[i+1])  )  )    #Används till att kontrollera algorithmen via ecmwf.int

    '''Fulla nivåer Sodankylä'''

    FULL_LEVELS_EC_SOD = Fi/g0   #Denna är alla 137, eftersom vi började med 138 stycken A och B
    FULL_LEVELS_EC_SOD = FULL_LEVELS_EC_SOD[:VV_EC]

    FULL_LEVELS_EC_SOD_alla_tidssteg[:, TIDSSTEG] =  FULL_LEVELS_EC_SOD

    


    '''Beräknar potentiella temperaturen. Detta görs först här för att p_full behövs på varje nivå.'''


    THETA_EC_SOD = TEMPERATURE_EC_SOD *((1000/(p_full[:VV_EC]/100))**0.286)   #p_full är alla 137 så vi tar de VV_EC första av dessa 137 eftersom TEMPERATURE_EC_SOD är VV_EC lång

    THETA_EC_2D_SOD[:, TIDSSTEG] = THETA_EC_SOD   #ÄNDRA VID 2518 vs 2600, THETA_EC_2D_SOD[:, TIDSSTEG-2]





############
# PLOTTING #
#################################################################################################################
DATE = '180129'     #Endast för title, ylabel osv           #ÄNDRA

VARIABLE_1 = 'clwc' #Endast för title, ylabel osv           #ÄNDRA

VARIABLE_2 = 'cc'   #Endast för title, ylabel osv           #ÄNDRA

VARIABLE_3 = 'ciwc' #Endast för title, ylabel osv           #ÄNDRA



fig1 = plt.figure(1, figsize=(12,5))   


ax1 = plt.subplot(111)

TIME, LEVELS = N.meshgrid(N.array(N.arange(ANTAL_TIDSSTEG)), FULL_LEVELS_EC_SOD)  #Välj samma antal tidssteg som är gjort för 2D-variabeln som plottas här nedanför.


c_levels = [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1]


#myplot = ax1.contourf(TIME, LEVELS, QCLOUD_EC_2D, c_levels, colors=('snow','whitesmoke', 'lightgray', 'silver','darkgrey','gray', 'dimgrey', 'black', 'magenta'))#, norm = LogNorm())
myplot = ax1.contourf(TIME, LEVELS, QCLOUD_EC_2D,cmap ='binary')#, norm = LogNorm()) 
#myplot = ax1.scatter(TIME, LEVELS, QCLOUD_EC_2D, c_levels, colors=('snow','whitesmoke', 'lightgray', 'silver','darkgrey','gray', 'dimgrey', 'black', 'magenta'))
#myplot = ax1.pcolormesh(TIME, LEVELS, QCLOUD_EC_2D,  c_levels, colors=('snow','whitesmoke', 'lightgray', 'silver','darkgrey','gray', 'dimgrey', 'black', 'magenta'))


'''Plottar contours för Theta i SAMMA figur som sonderingarna'''

TIME, LEVELS = N.meshgrid(N.array(N.arange(ANTAL_TIDSSTEG)), FULL_LEVELS_EC_SOD)  #Tidssteg range(0,15) ger 0-14 ovan. När man kör arange på detta IGEN tappar vi ett tidssteg. Vi får 0-13.
                                            #ÄNDRA VID 2518 VS 2600                       #Dessutom.....när TIDSSTEG=0 blir vektorn tom. Därav TIDSSTEG+1.
                                            #N.arange(2,ANTAL_TIDSSTEG)
contour_levels = N.arange(0, 330, 0.5)
cp = ax1.contour(TIME, LEVELS, THETA_EC_2D_SOD,contour_levels, colors='red', linestyles='dotted', linewidth='0.2' )
plt.clabel(cp, inline=True, fmt='%1.1f', fontsize=10)





'''Plottar vertikala nivåer längs y-axeln'''


X = N.zeros_like(FULL_LEVELS_EC_SOD)[0:VV_EC]
ax1.scatter(X, FULL_LEVELS_EC_SOD[0:VV_EC], marker= '>', color='black', s=10)

                                                


ax1.grid(linestyle="dotted")

#ax.legend()

ax1.set_ylabel('Height (m)', fontsize=18)

#ax1.set_xlabel('Vertical model soundings every 3rd hour')

ax1.set_xticks(N.arange(0,15))
plt.tick_params(labelsize=15)
ax1.set_xticklabels(ALLA_DATUM, rotation=45, ha='right', fontsize=15)
ax1.set_xlim(0, 14)
ax1.set_ylim(0,2500)

ax1.set_title(Filobjekt_GRIB.variables[VARIABLE_1].long_name, fontsize=18)

fig1.tight_layout()
fig1.subplots_adjust(right=0.90)
fig1.subplots_adjust(bottom=0.3)
cbar_ax = fig1.add_axes([0.92, 0.30, 0.01, 0.62])
#fig1.colorbar(myplot, cax=cbar_ax,ticks=[1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1], format = "%8.1e")
fig1.colorbar(myplot, cax=cbar_ax, format = "%8.1e")
#plt.savefig('/home/sm_marha/FIGURES/FOKUS_17_27_FEB/18_FEB/ECMWF_Crossection_2500_%s_%s' % (VARIABLE_1, DATE))   #ÄNDRA










fig2 = plt.figure(2, figsize=(12,5))   


ax1 = plt.subplot(111)

TIME, LEVELS = N.meshgrid(N.array(N.arange(ANTAL_TIDSSTEG)), FULL_LEVELS_EC_SOD[0:VV_EC])   #ÄNDRA VID 2518 VS 2600, (N.arange(2,ANTAL_TIDSSTEG))
                                    


myplot2 = ax1.contourf(TIME, LEVELS, CLOUD_FRACTION_EC_2D, extend='both', cmap ='binary')
#myplot2 = ax1.scatter(TIME, LEVELS, c=QCLOUD_EC_2D, cmap ='binary')
#myplot2 = ax1.pcolormesh(TIME, LEVELS, QCLOUD_EC_2D, cmap ='binary')




'''Plottar masslevels längs y-axeln'''

X = N.zeros_like(FULL_LEVELS_EC_SOD)[0:VV_EC]
ax1.scatter(X, FULL_LEVELS_EC_SOD[0:VV_EC], marker= '>', color='black', s=10)


ax1.grid(linestyle="dotted")

#ax.legend()

ax1.set_ylabel('Height (m)', fontsize=18)

#ax1.set_xlabel('Vertical model soundings every 3rd hour')
ax1.set_xticks(N.arange(0,15))
ax1.set_xticklabels(ALLA_DATUM, rotation=45, ha='right', fontsize=15)
ax1.set_xlim(0, 14)
ax1.set_ylim(0,2500)
ax1.set_title(Filobjekt_GRIB.variables[VARIABLE_2].long_name + '(g/kg)', fontsize=18)
fig2.tight_layout()
fig2.subplots_adjust(right=0.90)
fig2.subplots_adjust(bottom=0.3)
cbar_ax = fig2.add_axes([0.92, 0.30, 0.01, 0.62])
fig2.colorbar(myplot2, cax=cbar_ax, format = "%8.1e")


#plt.savefig('/home/sm_marha/FIGURES/FOKUS_17_27_FEB/18_FEB/ECMWF_Crossection_2500_%s_%s' % (VARIABLE_2, DATE))







fig3 = plt.figure(3, figsize=(12,5))   


ax1 = plt.subplot(111)

TIME, LEVELS = N.meshgrid(N.array(N.arange(ANTAL_TIDSSTEG)), FULL_LEVELS_EC_SOD)     #ÄNDRA VID 2518 VS 2600, (N.arange(2,ANTAL_TIDSSTEG))  #Välj samma antal tidssteg som är gjort för 2D-variabeln som plottas här nedanför.

myplot3 = ax1.contourf(TIME, LEVELS, QICE_EC_2D,cmap ='binary')#, norm = LogNorm()) 
#myplot3 = ax1.scatter(TIME, LEVELS, c=QICE_EC_2D, cmap ='binary')
#myplot3 = ax1.pcolormesh(TIME, LEVELS, QICE_EC_2D, cmap ='binary')


'''Plottar contours för Theta i SAMMA figur som sonderingarna'''

TIME, LEVELS = N.meshgrid(N.array(N.arange(ANTAL_TIDSSTEG)), FULL_LEVELS_EC_SOD)  #Tidssteg range(0,15) ger 0-14 ovan. När man kör arange på detta IGEN tappar vi ett tidssteg. Vi får 0-13.
                                                                                         #Dessutom.....när TIDSSTEG=0 blir vektorn tom. Därav TIDSSTEG+1.

contour_levels = N.arange(0, 330, 0.5)
cp = ax1.contour(TIME, LEVELS, THETA_EC_2D_SOD,contour_levels, colors='red', linestyles='dotted', linewidth='0.2' )
plt.clabel(cp, inline=True, fmt='%1.1f', fontsize=10)




'''Plottar vertikala nivåer längs y-axeln'''


X = N.zeros_like(FULL_LEVELS_EC_SOD)[0:VV_EC]
ax1.scatter(X, FULL_LEVELS_EC_SOD[0:VV_EC], marker= '>', color='black', s=10)

                                                


ax1.grid(linestyle="dotted")

#ax.legend()

ax1.set_ylabel('Height (m)', fontsize=18)

#ax1.set_xlabel('Vertical model soundings every 3rd hour')

ax1.set_xticks(N.arange(0,15))
plt.tick_params(labelsize=15)
ax1.set_xticklabels(ALLA_DATUM, rotation=45, ha='right', fontsize=15)
ax1.set_xlim(0, 14)
ax1.set_ylim(0,2500)

ax1.set_title(Filobjekt_GRIB.variables[VARIABLE_3].long_name, fontsize=18)

fig3.tight_layout()
fig3.subplots_adjust(right=0.90)
fig3.subplots_adjust(bottom=0.3)
cbar_ax = fig3.add_axes([0.92, 0.30, 0.01, 0.62])
fig3.colorbar(myplot3, cax=cbar_ax, format = "%8.1e")

#plt.savefig('/home/sm_marha/FIGURES/FOKUS_17_27_FEB/18_FEB/ECMWF_Crossection_2500_%s_%s' % (VARIABLE_3, DATE))   #ÄNDRA



plt.show()
