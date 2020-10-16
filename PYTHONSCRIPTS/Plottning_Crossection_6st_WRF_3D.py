import numpy as N
import scipy.io.netcdf as S
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import os
import datetime
import netCDF4 as N4
from matplotlib.colors import LogNorm



'''Detta skript plottar WRF contourf för 6 variabler samtidigt
OBS!!! Ändra även sökvägen för filobjektet längre ner, inte bara sökvägen för file_list direkt här nedanför.
Behövs även ändra i filnamnet när man sparar längst ner. T ex om det är med, eller utan, QC_QI'''




file_list = os.listdir('/nobackup/smhid12/sm_marha/2018-02-27_00z_91_LEVELS_Clouds_Where_Clouds_In_Sounding_YSU_TOPD1/run/')    #ÄNDRA


wrfout_list = []

for file in file_list:
    if file.startswith('wrfout'):
        wrfout_list.append(file)

wrfout_list.sort()

#wrfout_list = wrfout_list[0:6]     #Om man vill ha färre filer än "ALLA"

wrfout_list = wrfout_list[0:19]


'''Skriver ut alla wrf-filer'''


for file in wrfout_list[0:19]:          #Ha med [0:37] om det finns fler körningar än 6h
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

LATITUDE = 400

LONGITUDE = 400


for file in wrfout_list:     #Ha med [0:37] om det finns fler körningar än 6h

        print(file)
        
        Filobjekt = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-27_00z_91_LEVELS_Clouds_Where_Clouds_In_Sounding_YSU_TOPD1/run/'+file, mode='r')  #ÄNDRA

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

############
# Plotting #   Loopar över två listor för att göra 6 figurer samtidigt.
###################################################################################################

'''Plottning av sex olika figurer i form av tvärsnitt'''


i=1

#DATA_LIST = [CLOUD_WATER_MATRIX, CLOUD_ICE_MATRIX, CLOUD_FRACTION_MATRIX, CLOUD_GRAUPEL_MATRIX, SNOW_MATRIX, RAIN_MATRIX]    #ÄNDRA

#VARIABLE_LIST =['QCLOUD', 'QICE', 'CLDFRA', 'QGRAUP', 'QSNOW', 'QRAIN']            #ÄNDRA

DATA_LIST = [CLOUD_WATER_MATRIX]        #ÄNDRA

VARIABLE_LIST =['QCLOUD']       #ÄNDRA


for DATA, VARIABLE in zip(DATA_LIST, VARIABLE_LIST):


    fig = plt.figure(str(i), figsize=(12,5))

    ax1 = plt.subplot(111)

    if VARIABLE == 'QCLOUD':
        myplot = ax1.contourf(TIME, LEVELS, DATA, cmap ='binary')   #decode("utf-8") Denna behövs med scipy, men ej med Netcdf4
        contour_levels = N.arange(0, 330, 0.5)
        myplot_2 = ax1.contour(TIME, LEVELS, THETA_MATRIX, contour_levels, colors='red', linestyles='dotted', linewidth='0.2')
        plt.clabel(myplot_2, inline=True, fmt='%1.1f', fontsize=8)
    else:
        myplot = ax1.contourf(TIME, LEVELS, DATA, cmap ='binary')   #decode("utf-8") Denna behövs med scipy, men ej med Netcdf4


    if VARIABLE=='QCLOUD':
        plt.title(Filobjekt.variables[VARIABLE].description + ' ' + '(g/kg)' + '  QC QI 42h T_to_Td Whole MYJ ', fontsize=18)# + ' and Potential temperature (K)', fontsize=18)           #ÄNDRA
    elif VARIABLE=='CLDFRA':
        plt.title(Filobjekt.variables[VARIABLE].description, fontsize=18)
    else:
        plt.title(Filobjekt.variables[VARIABLE].description + ' ' + '(g/kg)', fontsize=18)





    plt.tick_params(labelsize=15)
    plt.xticks(TIDSVEKTOR_MINUTER[0:len(TIDSVEKTOR_MINUTER):3], TIDSVEKTOR_DATUM[0:len(TIDSVEKTOR_DATUM):3], rotation= 45, size=6, ha='right', fontsize=15)
    ax1.set_ylabel('Height (m)', fontsize=18)
    ax1.set_xlabel('Time (UTC)', fontsize=18)
    plt.ylim(0,2500)            #ÄNDRA
    #plt.xlim(0,360)
    #ax1.legend(loc=1,prop={'size': 12})
    ax1.grid(linestyle="dotted")


    # Plottar modellnivåer
    X = N.zeros_like(MASSLEVELS_1D_MINUS_TER)
    ax1.scatter(X+0.2, MASSLEVELS_1D_MINUS_TER, marker= '>', color='black', s=8)
    

    '''OBS, detta måste komma efter man använder plt., eftersom inställningarna görs på axes för colorbaren'''

    fig.tight_layout()
    fig.subplots_adjust(right=0.90)
    cbar_ax = fig.add_axes([0.92, 0.29, 0.01, 0.63])
    fig.colorbar(myplot, cax=cbar_ax, format = "%8.1e")
    

    ''' JAKOBS TIPS!!!!!!!

    filename = '/home/sm_marha/FIGURES/FOKUS_18_20_JAN/18_JAN/'
    filename += wrfout_list[1][11:21]
    filename += '_%s_00z' % (VARIABLE)
    filename += '.%s' % outformat
    fig.savefig(filename)''' 
    
    #fig.savefig('/home/sm_marha/FIGURES/FOKUS_17_27_FEB/26_FEB/' + wrfout_list[1][11:21] + '_%s_00z_2500m_91_LEVELS_QC_QI_1e-5_42h_Icloud_2' % (VARIABLE)) #ÄNDRA
    #fig.savefig('/home/sm_marha/FIGURES/FOKUS_18_20_JAN/29_JAN/' + wrfout_list[1][11:21] + '_%s_00z_2000m' % (VARIABLE))
    i+=1

    




plt.show()
        
