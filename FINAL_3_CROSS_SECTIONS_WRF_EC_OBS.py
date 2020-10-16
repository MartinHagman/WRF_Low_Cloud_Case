import numpy as N
import netCDF4 as N4
import matplotlib.pyplot as plt
import os
import pandas as pd
import datetime
import sys
import matplotlib.colors as mplc
from matplotlib.colors import LogNorm
from pylab import *
import Picking_values_Sounding_Dat_file_version_2
import matplotlib.gridspec as gridspec


'''Viktigt i detta skript är att ändra i xlim samtidigt som man ändrar i antalet tidssteg! Så om man sätter t ex [0:19] för WRF ska xlim vara
ax1.set_xlim(0, 18). samma me ECMWF, om ANTAL TIDSSTEG = 7 blir det xlim(0,6).

Näe man har en figur utan colorbar används gridspec för att dra in den figuren åt höger. Gridspec behövs inte om man u stället ersätter
alla "barsen" med en gemensam, för då skapar man en axel först som drar in alla figurer åt vänster.'''


#######
# WRF #
################################################################################################################
print('WRF')



file_list = os.listdir('/nobackup/smhid12/sm_marha/2018-02-27_00z_91_LEVELS_Clouds_Where_Clouds_In_Sounding_YSU_1_100/run/')    #ÄNDRA

#file_list = os.listdir('/nobackup/smhid12/sm_marha/2018-02-18_00z_91_LEVELS_QC_QI_42h_T_to_Td_1e-4/run/')


wrfout_list = []

for file in file_list:
    if file.startswith('wrfout'):
        wrfout_list.append(file)

wrfout_list.sort()

#wrfout_list = wrfout_list[0:6]     #Om man vill ha färre filer än "ALLA"

wrfout_list = wrfout_list[0:19]


'''Skriver ut alla wrf-filer'''


for file in wrfout_list[0:19]:         # Ex vis 0 till 18h = [0:19]
    print (file)

    


'''Listar diverse variabler'''


wrfout_list_datum = []

TSK_TIDSSERIE = []

T2_TIDSSERIE = []

TEMP_C_MATRIX = []

CLOUD_WATER_MATRIX = []

CLOUD_ICE_MATRIX = []

CLOUD_FRACTION_MATRIX =[]

CLOUD_GRAUPEL_MATRIX = []

SNOW_MATRIX = []

RAIN_MATRIX = []

THETA_MATRIX = []

TIDSVEKTOR_MINUTER = range(len(wrfout_list[:]))  #Definierar vektor från 0 till 36 med längd 37. Plottar sedan denna på x-axeln,
                                              # men gör om xticklabels till aktuella datum och aktuell tid. Att den heter
                                              # minuter har ingen betydelse. Skulle lika gärna kunnat heta TIDSVEKTOR.


                                              
#TIDSVEKTOR_MINUTER = TIDSVEKTOR_MINUTER   #Ha med denna om det finns fler körningar än till än till 6h med tidssteget 10 min


TIDSVEKTOR_DATUM = []


#Sätter vilken punkt som ska plottas
#Sodankylä = 563, 352

LATITUDE = 563

LONGITUDE = 352


for file in wrfout_list:     #Ha med [0:37] om det finns fler körningar än 6h

        print(file)
        
        Filobjekt = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-27_00z_91_LEVELS_Clouds_Where_Clouds_In_Sounding_YSU_1_100/run/' + file, mode='r')
        #Filobjekt = Filobjekt.Dataset('/nobackup/smhid12/sm_marha/2018-02-27_00z_91_LEVELS_Clouds_Where_Clouds_In_Sounding_YSU_TOPD1/run/'+file, mode='r')  #ÄNDRA
        #Filobjekt = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-18_00z_91_LEVELS_QC_QI_42h_T_to_Td_1e-4/run/' + file, mode='r')
        
        TSK = Filobjekt.variables['TSK'][0, LATITUDE, LONGITUDE]-273.15
        T2 = Filobjekt.variables['T2'][0, LATITUDE, LONGITUDE]-273.15      
        
        THETA_PERT_1D = Filobjekt.variables['T'][0, :, LATITUDE, LONGITUDE]
        P_PERT_1D = Filobjekt.variables['P'][0, :, LATITUDE, LONGITUDE]
        P_BASE_1D = Filobjekt.variables['PB'][0, :, LATITUDE, LONGITUDE]   
        
        P_TOT_1D = P_BASE_1D + P_PERT_1D
        THETA_1D = THETA_PERT_1D+300
        TEMP_K_1D = THETA_1D/((1000/(P_TOT_1D/100))**0.286)
        TEMP_C_1D = TEMP_K_1D-273.15     


        PH = Filobjekt.variables['PH'][:]
        PHB = Filobjekt.variables['PHB'][:]
        MODELLNIVAHOJD = (PH+PHB)/9.81
        MASSLEVELS = 0.5*(MODELLNIVAHOJD[0,:-1, LATITUDE, LONGITUDE] + MODELLNIVAHOJD[0, 1:, LATITUDE, LONGITUDE])
        TERRAIN_HEIGHT = Filobjekt.variables['HGT'][0, LATITUDE, LONGITUDE]
        MASSLEVELS_1D_MINUS_TER = MASSLEVELS - TERRAIN_HEIGHT

        CLOUD_ICE_1D= Filobjekt.variables['QICE'][0, :, LATITUDE, LONGITUDE]*1000
        CLOUD_WATER_1D = Filobjekt.variables['QCLOUD'][0, :, LATITUDE, LONGITUDE]*1000
        CLOUD_GRAUPEL_1D = Filobjekt.variables['QGRAUP'][0, :, LATITUDE, LONGITUDE]*1000
        SNOW_1D = Filobjekt.variables['QSNOW'][0, :, LATITUDE, LONGITUDE]*1000
        RAIN_1D = Filobjekt.variables['QRAIN'][0, :, LATITUDE, LONGITUDE]*1000
        CLDFRA_1D = Filobjekt.variables['CLDFRA'][0, :, LATITUDE, LONGITUDE]*1000
        


        # Lägger till 1-D vektorer till listorna och får listor som innehåller numpy-arrayer.....Nybörjar-röra...
        CLOUD_WATER_MATRIX.append(CLOUD_WATER_1D)
        CLOUD_ICE_MATRIX.append(CLOUD_ICE_1D)
        CLOUD_FRACTION_MATRIX.append(CLDFRA_1D)
        CLOUD_GRAUPEL_MATRIX.append(CLOUD_GRAUPEL_1D)
        SNOW_MATRIX.append(SNOW_1D)
        RAIN_MATRIX.append(RAIN_1D)
        THETA_MATRIX.append(THETA_1D)
        
        TSK_TIDSSERIE.append(TSK)
        T2_TIDSSERIE.append(T2)
        TEMP_C_MATRIX.append(TEMP_C_1D)
        
        
        
        '''Fixar vektor med datum som ska tjäna som xticklabels'''
        
        wrfout_list_datum.append(file[11:])
        datum = datetime.datetime.strptime(str(file[11:]), '%Y-%m-%d_%H:%M:%S').strftime('%d/%-m %Hz')
        TIDSVEKTOR_DATUM.append(datum)


'''Vektorerna görs om till 2D, där espektive vektor får samma shape som DATA, så att de kan plottas med contourf'''

TIME, LEVELS = N.meshgrid(TIDSVEKTOR_MINUTER, MASSLEVELS_1D_MINUS_TER)

'''Transponerar för att få överensstämmelse med TIME och LEVELS ovan'''

CLOUD_WATER_MATRIX = N.transpose(CLOUD_WATER_MATRIX)
CLOUD_ICE_MATRIX = N.transpose(CLOUD_ICE_MATRIX)
CLOUD_FRACTION_MATRIX = N.transpose(CLOUD_FRACTION_MATRIX)
CLOUD_GRAUPEL_MATRIX = N.transpose(CLOUD_GRAUPEL_MATRIX)
SNOW_MATRIX = N.transpose(SNOW_MATRIX)
RAIN_MATRIX = N.transpose(RAIN_MATRIX)
THETA_MATRIX = N.transpose(THETA_MATRIX)



#########
# ECMWF #
#################################################################################################################

print('ECMWF')


Filobjekt_GRIB = N4.Dataset('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1ua2018022700_M_D.nc', mode='r')    #ÄNDRA
Filobjekt_GRIB_sfc = N4.Dataset('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1sfc2018022700_M_D.nc' , mode='r')   #ÄNDRA


VV_EC = 35

ANTAL_TIDSSTEG = 7

FULL_LEVELS_EC_SOD_alla_tidssteg = N.zeros((VV_EC, ANTAL_TIDSSTEG))
THETA_EC_2D_SOD = N.zeros((VV_EC, ANTAL_TIDSSTEG ))               #ÄNDRA VID 2518 vs 2600 (ANTAL_TIDSSTEG-2)        




LATITUD = 27 # Har värdet 0 uppe 
LONGITUD = 192 # Har värdet 0 på vänstra randen. ECMWF-området är längst i x-led till skillnad från WRF-området,som är längst i y-led.

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





#######################
#HÄMTAR SONDERINGSDATA#
#####################################################################################
    
print('SOUNDINGS')

VV_SOND = 700   #2000

#Definierar arrays med nollor. Fyra kolumner och VV_SOND rader

THETA_SONDERING_MATRIX = N.zeros((VV_SOND,4))

TEMP_C_MATRIX = N.zeros((VV_SOND,8))

ALTITUDE_TEMPERATURE_MATRIX = N.zeros((VV_SOND, 4))

TIME_MATRIX = N.ones((VV_SOND, 4))

DATE_LIST = []


#Endast för plottning

VARIABLE = 'THETA'              #ÄNDRA
DATE = '20180218-19'            #ÄNDRA


Sounding_times = N.array(['00', '6', '12', '18'])

counter = 0
counter_time = 0    #Till för att ge olika avstånd mellan xticks om en sondering saknas. På dessa sätts sedan xticklabels ut från DATE_LIST

for j in range(58,59):#ÄNDRA   Loopar över dagar. VIKTIGT: Börja på aktuellt datum, då
                               #Picking_values_Sounding_Dat_file_version_2.pick_dates börjar på 20180100 
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













############
# PLOTTING #
#################################################################################################################




print('PLOTTING CLOUD WATER')


gs = gridspec.GridSpec(ncols=20, nrows=3)        




#cbar_ax = fig1.add_axes([0.87, 0.53, 0.01, 0.35])

rc('axes', linewidth=1.5)

'''Här kan man välja om man vill ha WRF_SCM eller WRF_3D i första subploten. Vill man byta till WRF2 eller WRF_SCM2 byter man inne i plotkommandona'''


'''Plottning av sex olika figurer i form av tvärsnitt'''


################ PLOTTNING WRF ###############################

i=1

#DATA_LIST = [CLOUD_WATER_MATRIX, CLOUD_ICE_MATRIX, CLOUD_FRACTION_MATRIX, CLOUD_GRAUPEL_MATRIX, SNOW_MATRIX, RAIN_MATRIX]    #ÄNDRA

#VARIABLE_LIST =['QCLOUD', 'QICE', 'CLDFRA', 'QGRAUP', 'QSNOW', 'QRAIN']            #ÄNDRA

DATA_LIST = [CLOUD_WATER_MATRIX]        #ÄNDRA

VARIABLE_LIST =['QCLOUD']       #ÄNDRA


for DATA, VARIABLE in zip(DATA_LIST, VARIABLE_LIST):


    fig1 = plt.figure(str(i), figsize=(12,8))                       #Här skapas figur objektet, som bara är 1 DATA_LIST har ett element, t ex QCLOUD
    fig1.subplots_adjust(right=0.85)
    #fig1.subplots_adjust(hspace=0.2)
    rc('axes', linewidth=1.5)

    ax1 = fig1.add_subplot(gs[0, :])
    ax1 = plt.subplot(311)
    ax1.tick_params(which='major', direction='in', width = 2, labelsize = '17')
    


    if VARIABLE == 'QCLOUD':
        norm = mplc.Normalize(0, 5e-1)
        Antal = 5e-1/20
        myplot = ax1.contourf(TIME, LEVELS, DATA, N.arange(0, 5e-1, Antal), norm = norm, extend = 'both', cmap ='binary')
        #myplot = ax1.contourf(TIME, LEVELS, DATA, cmap ='binary', extend = 'both')   #decode("utf-8") Denna behövs med scipy, men ej med Netcdf4
        contour_levels = N.arange(0, 330, 1.5)
        myplot_2 = ax1.contour(TIME, LEVELS, THETA_MATRIX, contour_levels, colors='red', linestyles='solid', linewidths= 0.5)
        ax1.clabel(myplot_2, contour_levels[::2], inline=True, fmt='%1.0f', fontsize=10)
    else:
        myplot = ax1.contourf(TIME, LEVELS, DATA, cmap ='binary')   #decode("utf-8") Denna behövs med scipy, men ej med Netcdf4





    plt.tick_params(labelsize=15)
    plt.xticks(TIDSVEKTOR_MINUTER[0:len(TIDSVEKTOR_MINUTER):3], TIDSVEKTOR_DATUM[0:len(TIDSVEKTOR_DATUM):3], rotation= 45, size=6, ha='right', fontsize=15)
    #ax1.set_ylabel('Height (m)', fontsize=18)
    #ax1.set_xlabel('Time (UTC)', fontsize=18)
    ax1.set_ylim(0,2000)            #ÄNDRA
    ax1.set_xlim(0, 18)
    #ax1.legend(loc=1,prop={'size': 12})
    ax1.grid(linestyle="dotted")
    ax1.set_xticklabels([])
    
    # Plottar modellnivåer
    X_W = N.zeros_like(MASSLEVELS_1D_MINUS_TER)
    ax1.scatter(X_W, MASSLEVELS_1D_MINUS_TER, marker= '>', color='black', s=8)
    ax1.text(17.5, 1700, 'a', fontsize = 18, fontweight='bold')

    '''OBS, detta måste komma efter man använder plt., eftersom inställningarna görs på axes för colorbaren'''

  
    
    
    #fig1.colorbar(myplot)#, format = "%8.1e")               # Enskild color bar
    






################ PLOTTNING ECMWF ###############################

ax2 = fig1.add_subplot(gs[1, :])
ax2.tick_params(which='major', direction='in', width =2, labelsize=15)
#ax2 = plt.subplot(312)

TIME, LEVELS = N.meshgrid(N.array(N.arange(ANTAL_TIDSSTEG)), FULL_LEVELS_EC_SOD[0:VV_EC])   #ÄNDRA VID 2518 VS 2600, (N.arange(2,ANTAL_TIDSSTEG))
                                    



norm = mplc.Normalize(0, 5e-1)
Antal = 5e-1/20
myplot3 = ax2.contourf(TIME, LEVELS, QCLOUD_EC_2D, N.arange(0, 5e-1, Antal), norm = norm, extend = 'both', cmap ='binary')
#myplot3 = ax2.contourf(TIME, LEVELS, QCLOUD_EC_2D, extend='both', cmap ='binary')
#myplot3 = ax2.scatter(TIME, LEVELS, c=QCLOUD_EC_2D, cmap ='binary')
#myplot3 = ax2.pcolormesh(TIME, LEVELS, QCLOUD_EC_2D, cmap ='binary')

contour_levels = N.arange(0, 330, 1.5)
cp = ax2.contour(TIME, LEVELS, THETA_EC_2D_SOD,contour_levels, colors='red', linestyles='solid', linewidths= 0.5 )

ax2.clabel(cp,contour_levels[::2], inline=True, fmt='%1.0f', fontsize=10)





'''Plottar masslevels längs y-axeln'''

X = N.zeros_like(FULL_LEVELS_EC_SOD)[0:VV_EC]
ax2.scatter(X, FULL_LEVELS_EC_SOD[0:VV_EC], marker= '>', color='black', s=10)


ax2.grid(linestyle="dotted")

#ax.legend()

ax2.set_ylabel('Height (m)', fontsize=15)

#ax2.set_xlabel('Vertical model soundings every 3rd hour')
ax2.set_xticks(N.arange(0,15))
ax2.set_yticklabels([0, 500, 1000, 1500, 2000], fontsize=15)
#ax2.set_xticklabels(ALLA_DATUM, rotation=45, ha='right', fontsize=14)
ax2.set_xticklabels([])
ax2.set_xlim(0, 6)
ax2.set_ylim(0,2000)
ax2.text(5.8, 1700, 'b', fontsize = 18, fontweight='bold')
#cbar_ax = fig2.add_axes([0.92, 0.30, 0.01, 0.62])


cbar_ax = fig1.add_axes([0.87, 0.38, 0.01, 0.50])       #Gemensam bar
cbar = fig1.colorbar(myplot, cax=cbar_ax)#, format = "%8.1e")
cbar_ax.set_title(' ($gkg^{-1})$', fontsize = 16)
cbar.ax.tick_params(labelsize=16)

#fig1.colorbar(myplot3)#, format = "%8.1e")              # Enskild color bar





################ PLOTTNING SOUNDINGS CROSS SECTIONS ###############################


ax3 = fig1.add_subplot(gs[2, :])             #Om enskild colorbar för plot 1 och 2 sätts ax3 = fig1.add_subplot(gs[2, :16])
#ax3 = plt.subplot(313)
ax3.tick_params(which='major', direction='in', width =2, labelsize = '17')

contour_levels = N.arange(0, 330, 1.5)      #ÄNDRA
#contour_levels = N.arange(-80, 10, 0.5)    #ÄNDRA


cp = ax3.contour(TIME_MATRIX, ALTITUDE_TEMPERATURE_MATRIX,  THETA_SONDERING_MATRIX, contour_levels, colors='red', linestyles='solid', linewidths= 0.5 ) #ÄNDRA

ax3.set_xticks(TIME_MATRIX[0, :])
ax3.set_ylim(0,2000)            #ÄNDRA
plt.tick_params(labelsize=15)          #Ändra tick-storlek utan att ändra labels
ax3.set_xticklabels(DATE_LIST, rotation=45, fontsize=15)
#ax3.set_title('Sounding potential temperature (K)', fontsize=18)
#ax3.set_ylabel('Height (m)', fontsize=12)
ax3.clabel(cp, contour_levels[::2], inline=True, fmt='%1.0f', fontsize=10)
ax3.grid(linestyle="dotted")
ax3.text(2.9, 1700, 'c', fontsize = 18, fontweight='bold')


plt.savefig('/home/sm_marha/FIGURES/FINAL/20180227/Contours_WRF_THETA_LIQ_ECMWF_SOUNDINGS_20180227_YSU_1')    #ÄNDRA
#plt.savefig('/home/sm_marha/FIGURES/FINAL/20180218/Contours_WRF_ECMWF_SOUNDINGS_20180218_NEW')



plt.show()
