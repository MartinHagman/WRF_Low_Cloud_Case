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
import matplotlib.ticker as ticker




################
# ECMWF HYBRID #
#######################################################################################

'''Detta skript plottar ECMWF-sonderingar, QCLOUD isoarea(valbar parameter), samt isolinjer Theta. Välj endast GRIB-filer överst
i skriptet, en ua och en sfc från samma datum.'''



Filobjekt_GRIB = N4.Dataset('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1ua2018012900_M_D.nc', mode='r')
Filobjekt_GRIB_sfc = N4.Dataset('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1sfc2018012900_M_D.nc' , mode='r')



VV_EC = 32

THETA_EC_2D_SOD = N.zeros((137, 19))




#I GRIB-filerna finns 19 stycken tidssteg. Här väljs index från 0 till 15, 0-42h. For-loopen rullar så många man väljer!


'''Här skapar vi 'TIDSSTEG' antal vertikala sonderingar. Även plottkommandot loopas, men figuren plottas inte förrns på slutet. '''

for TIDSSTEG in range(0,15):

    print(TIDSSTEG)

    '''EC-DATA'
       Upper air'''

    LATITUD = 27 # Har värdet 0 uppe 
    LONGITUD = 191 # Har värdet 0 på vänstra randen. ECMWF-området är längst i x-led till skillnad från WRF-området,som är längst i y-led.

    TIMESTEPS = Filobjekt_GRIB.variables['time'][:]
    time = Filobjekt_GRIB.variables['time']
    dates = N4.num2date(TIMESTEPS, time.units, time.calendar)
    ALLA_DATUM=[]


    '''Gör om tiden till ett visst datumformat.'''

    for date in dates:
        datum = date.strftime('%d/%m %Hz')
        ALLA_DATUM.append(datum)

    
    LONGITUDE = Filobjekt_GRIB.variables['longitude'][:]
    LATITUDE = Filobjekt_GRIB.variables['latitude'][:]

    lon_matrix, lat_matrix = N.meshgrid(LONGITUDE, LATITUDE)

    QICE_EC = Filobjekt_GRIB.variables['ciwc'][:]
    QICE_EC_SOD = QICE_EC[TIDSSTEG, :, LATITUD, LONGITUD]
    QICE_EC_SOD = QICE_EC_SOD[::-1]         #Allt inverteras så att lägsta nivån hamnar längst ner i Arrayen.

    QCLOUD_EC = Filobjekt_GRIB.variables['clwc'][:]
    QCLOUD_EC_SOD = QCLOUD_EC[TIDSSTEG, :, LATITUD, LONGITUD]
    QCLOUD_EC_SOD = QCLOUD_EC_SOD[::-1]

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
        
    TEMPERATURE_C_EC = TEMPERATURE_EC-273.15     #Denna är inte inverterad i höjdled här och fortfarande 4-dimensionell!!
 

    SPECIFIC_HUMIDITY_EC = Filobjekt_GRIB.variables['q'][:]
    SPECIFIC_HUMIDITY_EC_SOD = SPECIFIC_HUMIDITY_EC[TIDSSTEG, :, LATITUD, LONGITUD]
    SPECIFIC_HUMIDITY_EC_SOD = SPECIFIC_HUMIDITY_EC_SOD[::-1]

    TIMESTEPS_EC = Filobjekt_GRIB.variables['time'][:]
    VERTICAL_LEVELS_EC = Filobjekt_GRIB.variables['level'][:]  #Ger vektor med värden 1-138 på plats 0 till 137




    '''markgeopotentialen och marktrycket från surface-gribfilen'''

    SURFACE_GEOPOTENTIAL_EC = Filobjekt_GRIB_sfc.variables['z'][:]
    SURFACE_GEOPOTENTIAL_EC_SOD = SURFACE_GEOPOTENTIAL_EC[TIDSSTEG, LATITUD, LONGITUD]  #Endimensionell, MEN EJ inverterad

    SURFACE_PRESSURE_EC = Filobjekt_GRIB_sfc.variables['sp'][:]
    SURFACE_PRESSURE_EC_SOD = SURFACE_PRESSURE_EC[TIDSSTEG, LATITUD, LONGITUD]  #Endimensionell, MEN EJ inverterad




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

    FULL_LEVELS_EC_SOD = Fi/g0

    
    '''Beräknar potentiella temperaturen. Detta görs först här för att p_full behövs på varje nivå.'''

    THETA_EC_SOD = TEMPERATURE_EC_SOD *((1000/(p_full/100))**0.286)

    THETA_EC_2D_SOD[:, TIDSSTEG] = THETA_EC_SOD

    


    ############
    # PLOTTING #
    #################################################################################################################
    
   
    

    fig1 = plt.figure(1, figsize=(16,6))

    
    '''Om tidssteget är lika med 0, så plottas contourf för QCLOUD_EC_2D för hela perioden'''

    if TIDSSTEG == 0:
        
        ax1 = plt.subplot(111)
        
        TIME, LEVELS = N.meshgrid(N.array(N.arange(0,15)), FULL_LEVELS_EC_SOD[0:VV_EC])

        myplot = ax1.contourf(TIME, LEVELS, QCLOUD_EC_2D, extend='both', cmap ='binary')
        #myplot = ax1.scatter(TIME, LEVELS, c=QCLOUD_EC_2D, cmap ='binary')
        #myplot = ax1.pcolormesh(TIME, LEVELS, QCLOUD_EC_2D, cmap ='binary')

        
        
        ax2 = ax1.twiny()   #Gemensam y-axel, men olika x-axlar.

    ax2.plot(THETA_EC_SOD[0:VV_EC]+5+TIDSSTEG*30, FULL_LEVELS_EC_SOD[0:VV_EC],linewidth = 1, c = 'b', label = 'ECMWF tidssteg ' + str(TIDSSTEG) )    #Förskjutning 30 för varje loop för att passa med tiden i comtourf. +5 är finlir med förskjutningen. Detta får göras om, när man väljer färre sonderingar.
    

    if TIDSSTEG%2 == 0:    #Olika färger på tidsangivelserna ovanför varje tempkurva!
        
        ax2.text(THETA_EC_SOD[VV_EC-1]+5+TIDSSTEG*30, FULL_LEVELS_EC_SOD[VV_EC-1]+30, ALLA_DATUM[TIDSSTEG], horizontalalignment='center', color = 'g', fontsize=7)
        
    else:
        ax2.text(THETA_EC_SOD[VV_EC-1]+5+TIDSSTEG*30, FULL_LEVELS_EC_SOD[VV_EC-1]+30, ALLA_DATUM[TIDSSTEG], horizontalalignment='center', color = 'r', fontsize=7)
        


'''Plottar contours för Theta i SAMMA figur som sonderingarna'''

TIME, LEVELS = N.meshgrid(N.array(N.arange(0,TIDSSTEG+1)), FULL_LEVELS_EC_SOD[0:VV_EC])  #Tidssteg range(0,15) ger 0-14 ovan. När man kör arange på detta IGEN tappar vi ett tidssteg. Vi får 0-13.
                                                                                         #Dessutom.....när TIDSSTEG=0 blir vektorn tom. Därav TIDSSTEG+1.

contour_levels = N.arange(0, 330, 0.5)
cp = ax1.contour(TIME, LEVELS, THETA_EC_2D_SOD[0:VV_EC, 0:TIDSSTEG+1],contour_levels, colors='red', linestyles='dotted', linewidth='0.2' )
plt.clabel(cp, inline=True, fmt='%1.1f', fontsize=8)



'''Plottar masslevels längs y-axeln''' 
X = N.zeros_like(FULL_LEVELS_EC_SOD)[0:VV_EC]-0.2
ax1.scatter(X, FULL_LEVELS_EC_SOD[0:VV_EC], marker= '>', color='black', s=10)


ax1.grid(linestyle="dotted")

#ax.legend()

ax1.set_ylabel('Height (m)')

#ax1.set_xlabel('Vertical model soundings every 3rd hour')

ax1.set_xticks(N.arange(0,15))
ax1.set_xticklabels(ALLA_DATUM, rotation=45)
ax1.set_xlim(-2, 15)
#cbar = fig1.colorbar(myplot, ax=ax2, orientation='horizontal', format = "%8.1e")
ax2.set_xticklabels([])
ax2.set_xticks([])
ax1.set_xlim(-2, 15)


ax2.set_xlim(210, 720)

plt.suptitle('Model soundings PT (K), isolines PT (K)  and  isoarea cloud water (g/kg)')





'''Ny enskild förgur av 1 sondering i taget, för att se detaljerna bättre.'''

fig2 = plt.figure(2)

ax1 = plt.subplot(111)

ax1.grid(linestyle="dotted")

ax1.plot(THETA_EC_2D_SOD[0:VV_EC, 0], FULL_LEVELS_EC_SOD[:VV_EC],linewidth = 1, c = 'b', label = 'ECMWF(00z)' + str(TIDSSTEG) )

ax1.set_title('Potential temperature')
ax1.set_ylabel('Height (m)')
ax1.set_xlabel('Potential temperature (K)')
ax1.legend(fontsize = 'large')

plt.show()






