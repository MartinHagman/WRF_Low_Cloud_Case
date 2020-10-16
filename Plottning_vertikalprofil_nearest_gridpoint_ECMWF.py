import numpy as N
import scipy.io.netcdf as S
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import pandas as pd
from math import *
import matplotlib.ticker as ticker
from matplotlib.legend import Legend

'''Igår, 190818, har jag korrigerat förskjutningen i höjdled mellan ECMWF och Met_em-filerna, genom att börja plotta andra elementet i GHT, dvs 9-10 meter.'''


'''I detta skript kan man välja om man vill plotta data från en högupplöst körning med varje gridpunkt från ECMWF, eller från
lågupplöst med varannan. Välj mellan filobjekten här nedanför. Man måste också välja uppsättningen gridpunktsnummer nedanför
filobjekten, beroende på vilken körning man valt.Välj även filobjekt för WRF-data längre ned Dessutom väljer man variabel och
variabelnamn under plottning.

I båda if-satserna i plottkommandot, måste man också ändra indexnummer till den högupplösta varianten när man ändrar.'''

'''Öppnar filobjektet'''

#Filobjekt coarse resolution

#Filobjekt_GRIB = S.netcdf_file('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1ua2018021800_M_D.nc', mode='r')
#Filobjekt_GRIB_sfc = S.netcdf_file('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1sfc2018021800_M_D.nc', mode='r')


#Filobjekt high resolution

Filobjekt_GRIB = S.netcdf_file('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1ua2018011800_MH_D.nc', mode='r')
Filobjekt_GRIB_sfc = S.netcdf_file('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1sfc2018011800_MH_D.nc', mode='r')




#ECMWF coarse resolution (27, 191) = Sodankylä


#LATITUDE_LIST = [27, 27, 27, 26, 26, 26, 28, 28, 28]

#LONGITUDE_LIST = [190, 191, 192, 190, 191, 192, 190, 191, 192]



#ECMWF high resolution (53,381) = Sodankylä

LATITUDE_LIST = [53, 53, 53, 52, 52, 52, 54, 54, 54]

LONGITUDE_LIST = [380, 381, 382, 380, 381, 382, 380, 381, 382]








for LATITUD, LONGITUD in zip(LATITUDE_LIST, LONGITUDE_LIST): 

    ################
    # ECMWF HYBRID #
    #######################################################################################


    #VARIABEL_GRIB = Filobjekt_GRIB.variables[VARIABEL_GRIB][:]
    LONGITUDE = Filobjekt_GRIB.variables['longitude'][:]
    LATITUDE = Filobjekt_GRIB.variables['latitude'][:]

    lon_matrix, lat_matrix = N.meshgrid(LONGITUDE, LATITUDE)
    QICE_EC = Filobjekt_GRIB.variables['ciwc'][:]
    QICE_EC_SOD = QICE_EC[0, :, LATITUD, LONGITUD]
    QICE_EC_SOD = QICE_EC_SOD[::-1]

    QCLOUD_EC = Filobjekt_GRIB.variables['clwc'][:]
    QCLOUD_EC_SOD = QCLOUD_EC[0, :, LATITUD, LONGITUD]
    QCLOUD_EC_SOD = QCLOUD_EC_SOD[::-1]


    CLOUD_FRACTION_EC = Filobjekt_GRIB.variables['cc'][:]
    CLOUD_FRACTION_EC_SOD = CLOUD_FRACTION_EC[0, :, LATITUD, LONGITUD]
    CLOUD_FRACTION_EC_SOD = CLOUD_FRACTION_EC_SOD[::-1]


    TEMPERATURE_EC = Filobjekt_GRIB.variables['t'][:]
    TEMPERATURE_C_EC = TEMPERATURE_EC-273.15
    TEMPERATURE_C_EC_SOD = TEMPERATURE_C_EC[0, :, LATITUD, LONGITUD]
    TEMPERATURE_C_EC_SOD = TEMPERATURE_C_EC_SOD[::-1]

    GEOPOTENTIAL = Filobjekt_GRIB_sfc.variables['z'][:]
    GEOPOTENTIAL_1D_SOD = GEOPOTENTIAL[0, 27, 192]




    '''markgeopotentialen och marktrycket från surface-gribfilen'''

    SURFACE_GEOPOTENTIAL_EC = Filobjekt_GRIB_sfc.variables['z'][:]
    SURFACE_GEOPOTENTIAL_EC_SOD = SURFACE_GEOPOTENTIAL_EC[0, LATITUD, LONGITUD]

    SURFACE_PRESSURE_EC = Filobjekt_GRIB_sfc.variables['sp'][:]
    SURFACE_PRESSURE_EC_SOD = SURFACE_PRESSURE_EC[0, LATITUD, LONGITUD]




    '''Plottningarametrar från upper_air filen''' 

    SPECIFIC_HUMIDITY_EC = Filobjekt_GRIB.variables['q'][:]
    SPECIFIC_HUMIDITY_EC_SOD = SPECIFIC_HUMIDITY_EC[0, :, LATITUD, LONGITUD]
    SPECIFIC_HUMIDITY_EC_SOD = SPECIFIC_HUMIDITY_EC_SOD[::-1]

    TIMESTEPS_EC = Filobjekt_GRIB.variables['time'][:]
    VERTICAL_LEVELS_EC = Filobjekt_GRIB.variables['level'][:]



    '''Använder Surface pressure på surface i stället, då dessa p g a interpolering skiljer sig lite'''

    #LN_SURFACE_PRESSURE = Filobjekt_GRIB.variables['lnsp'][:]

    #SURFACE_PRESSURE = N.exp(LN_SURFACE_PRESSURE)    
    #SURFACE_PRESSURE_SOD = SURFACE_PRESSURE[0,2,LATITUD, LONGITUD]






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




    '''Räknar ut virtuella temperaturen'''

    T_VIRTUAL_EC = (1.+0.609133*SPECIFIC_HUMIDITY_EC)*TEMPERATURE_EC

    T_VIRTUAL_EC_SOD = T_VIRTUAL_EC[0,:, LATITUD, LONGITUD]


    T_VIRTUAL_EC_SOD = T_VIRTUAL_EC_SOD[::-1] #Inverterar vektorn så nivå 137 kommer först)




    '''Räknar ut de halva trycknivåerna'''

    p_half = N.zeros(len(A))

    for i in range(len(A)):

        p = A[i] + B[i] * SURFACE_PRESSURE_EC[0,LATITUD,LONGITUD]

        p_half[i] = p


    '''Räknar ut geopotentialen på de halva nivåerna



    for i in range(len(p_half)-1):

        Fi_half[i+1] = Fi_half[i] + Rd*(  T_VIRTUAL_SOD[i]  *  (  log(p_half[i])-log(p_half[i+1])  ) )'''
     


    '''Räknar ut hela trycknivåerna'''
     
    p_full = 0.5*(p_half[1:]+p_half[:-1])



    '''Räknar ut geopotentialen på de fulla nivåerna.'''

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





    FULL_LEVELS_EC_SOD = Fi/g0





    ################
    # PLOTTNING    #
    #######################################################################################

    '''Välj variabel och variabelnamn innan plottning. Ändra också filnamn på det som sparas.'''


    VARIABLE = QICE_EC_SOD  #Ändra variabel här
    VARIABLE_NAME = 'ciwc'  #Ändra variabelnamn här


   
    VV_EC = 25    #60

    

    if LATITUD == int(53) and LONGITUD == int(381):       #Välj index för Sodankylä beroende på upplösning.(27, 191), (53, 381)
        
        c = 'r'
        

    else:

        c='g'
        

    

    fig = plt.figure(1)
    ax = plt.subplot(111)
    
    if LATITUD == int(53) and LONGITUD == int(381):                                                                  #Välj index för Sodankylä beroende på upplösning.(27,191),(53, 381)
        ax.plot(VARIABLE[0:VV_EC], FULL_LEVELS_EC_SOD[0:VV_EC], linewidth = 1, color = c, label = 'Sodankylä GRIB'  )    
    elif LATITUD == int(53) and LONGITUD == int(380):                                                                #Välj index för första gridpunkten i vektorn ovan beroende på upplösning.(27, 190),(53, 380) 
        ax.plot(VARIABLE[0:VV_EC], FULL_LEVELS_EC_SOD[0:VV_EC], linewidth = 1, color = c, label = 'Nearest points GRIB'  )       
    else:
        ax.plot(VARIABLE[0:VV_EC], FULL_LEVELS_EC_SOD[0:VV_EC], linewidth = 1, color = c)
                
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%8.1e'))
    ax.tick_params(axis="x", labelsize=9)

    plt.title('Nearest gridpoints Sodankylä in ECMWF data')
    plt.ylabel('Height (m)')
    plt.xlabel(Filobjekt_GRIB.variables[VARIABLE_NAME].long_name.decode("utf-8") + '    (' + Filobjekt_GRIB.variables[VARIABLE_NAME].units.decode("utf-8") + ')')



#################################################################
# WRF DATA MET_EM + PLOTTNING EFTER LOOPEN AV ECMWF DATA OVAN   #
#######################################################################################

'''Välj rätt WRF körning, highres_GRIB eller den vanliga. Välj även variabel'''



'''WRF från coarse resolution'''

#Filobjekt = S.netcdf_file('/nobackup/smhid12/sm_marha/2018-02-18_00z_QC_QI/run/met_em.d01.2018-02-18_00:00:00.nc', mode='r')



'''WRF från high resolution'''

Filobjekt = S.netcdf_file('/nobackup/smhid12/sm_marha/2018-01-18_00z_QC_QI_HIGHRES_GRIB/run/met_em.d01.2018-01-18_00:00:00.nc', mode='r')





VARIABEL = 'QI'   #Ändra variabel här!


VARIABELVARDE = Filobjekt.variables[VARIABEL][:]
VARIABELVARDE_1D_SOD = VARIABELVARDE[0, :, 563, 352]
P = Filobjekt.variables['PRES'][:]
TEMP = Filobjekt.variables['TT'][:]
GHT = Filobjekt.variables['GHT'][:]
TEMP_C = TEMP-273.15
MODELLNIVAHOJD = GHT
MODELLNIVAHOJD_1D_SOD = MODELLNIVAHOJD[0, :, 563, 352]
TERRAIN_HEIGHT_SOD = Filobjekt.variables['SOILHGT'][0, 563, 352]
MODELLNIVAHOJD_1D_MINUS_TER = MODELLNIVAHOJD_1D_SOD - TERRAIN_HEIGHT_SOD


'''Höjdfilerna ser lite olika ut. GHT, efter avdrag av terränghöjden, är 0 meter, d v s marken. I FULL_LEVELS_EC är första elementet
10 meter, första fulla nivån. Detta måste jag korrigera för enligt nedan, där jag börjar med 2a elementet i GHT, då detta blir första
fulla nivån ovan marken. Det är den första nivån som har molnvatten. På marken finns t ex inget molnvatten, så den parametern kommer
alltid att vara 0 där. Eftersom variabelvärdet också kommer från met_em-filerna(även i fallet med ECs vertikala nivåer), så hoppar
jag över även det första elementet där, då ju detta är markvärdet av den aktuella variabeln. I och med detta måste jag dra av ett element,
enligt FULL_LEVELS_EC_SOD[0:VV_EC-1], för att vektorerna ska bli lika långa.'''

ax.plot(VARIABELVARDE_1D_SOD[1:VV_EC], MODELLNIVAHOJD_1D_MINUS_TER[1:VV_EC], color = 'b', label = 'Met_em levels')
ax.plot(VARIABELVARDE_1D_SOD[1:VV_EC], FULL_LEVELS_EC_SOD[0:VV_EC-1], color = 'y', label = 'GRIB levels')
ax.legend()


#fig.savefig('/home/sm_marha/FIGURES/FOKUS_17_27_FEB/18_FEB/20180218_nearest_gridpoint_highres_QI')

plt.show()
