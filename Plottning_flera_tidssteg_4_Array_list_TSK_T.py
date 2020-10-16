import numpy as N
import scipy.io.netcdf as S
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import os
import datetime
import netCDF4 as N4




'''Detta skript plottar fyra olika figurer för den variabel och för de directoryn som anges i skriptet'''




FILE_LIST_ARRAY = [os.listdir('/nobackup/smhid12/sm_marha/2018-02-18_00z_91_LEVELS/run')]#),
            #os.listdir('/nobackup/smhid12/sm_marha/2018-01-18_00z_QC_QI/run'),
            #os.listdir('/nobackup/smhid12/sm_marha/2018-02-18_00z/run'),
            #os.listdir('/nobackup/smhid12/sm_marha/2018-02-18_00z_QC_QI/run')]


FILE_DIR_ARRAY = '/nobackup/smhid12/sm_marha/2018-02-18_00z_91_LEVELS/run/',   #Vid bara EN sträng och ingen lista, inga hakparanteser.
            #'/nobackup/smhid12/sm_marha/2018-01-18_00z_QC_QI/run/',
            #'/nobackup/smhid12/sm_marha/2018-02-18_00z/run/',
            #'/nobackup/smhid12/sm_marha/2018-02-18_00z_QC_QI/run/']



LAT = 563
LON = 352


wrfout_list_datum = []

TSK_TIDSSERIE = []

T2_TIDSSERIE = []

TEMP_C_TIDSSERIE = []




TIDSVEKTOR_DATUM = []



i=1  #För plottningen

for file_list, file_dir in zip(FILE_LIST_ARRAY, FILE_DIR_ARRAY): 


    wrfout_list = []

    for file in file_list:
       if file.startswith('wrfout'):
            wrfout_list.append(file)

    wrfout_list.sort()

    wrfout_list = wrfout_list[0:24]



    '''Definierar tomma vektorer inuti loopen för att de ska nollställas i varje loop'''

    wrfout_list_datum = []

    TSK_TIDSSERIE = []

    T2_TIDSSERIE = []

    TEMP_C_TIDSSERIE = []

    TIDSVEKTOR_DATUM = []

    TIDSVEKTOR_MINUTER = range(len(wrfout_list))  #Definierar vektor från 0 till 36 med längd 37. Plottar sedan denna på x-axeln,
                                                  # men gör om xticklabels till aktuella datum och aktuell tid. att den heter                                              # minuter har ingen betydelse. Skulle lika gärna kunnat heta TIDSVEKTOR.  
                                                  # minuter har ingen betydelse. Skulle lika gärna kunnat heta TIDSVEKTOR.


    '''Skriver ut alla wrf-filer'''


    for file in wrfout_list:
        print (file)




    
    for file in wrfout_list:
        
        print(file)
        
        Filobjekt = N4.Dataset(file_dir + file, mode='r')

        TSK = Filobjekt.variables['TSK'][0, LAT, LON]-273.15
        
        T2 = Filobjekt.variables['T2'][0, LAT, LON]-273.15      
        
        THETA_PERT = Filobjekt.variables['T'][0, 0, LAT, LON]
        P_PERT = Filobjekt.variables['P'][0, 0, LAT, LON]
        P_BASE = Filobjekt.variables['PB'][0, 0, LAT, LON]   
        
        P_TOT = P_BASE + P_PERT
        THETA = THETA_PERT+300
        TEMP_K = THETA/((1000/(P_TOT/100))**0.286)
        TEMP_C = TEMP_K-273.15     


        PH = Filobjekt.variables['PH'][:]
        PHB = Filobjekt.variables['PHB'][:]
        MODELLNIVAHOJD = (PH+PHB)/9.81
        MASSLEVELS = 0.5*(MODELLNIVAHOJD[0,:-1, LAT, LON] + MODELLNIVAHOJD[0, 1:, LAT, LON])
        TERRAIN_HEIGHT = Filobjekt.variables['HGT'][0, LAT, LON]
        MASSLEVELS_1D_MINUS_TER = MASSLEVELS - TERRAIN_HEIGHT


        TSK_TIDSSERIE.append(TSK)
        T2_TIDSSERIE.append(T2)
        TEMP_C_TIDSSERIE.append(TEMP_C)

        
        
        
        '''Fixar vektor med datum som ska tjäna som xticklabels'''
        
        wrfout_list_datum.append(file[11:])
        datum = datetime.datetime.strptime(str(file[11:]), '%Y-%m-%d_%H:%M:%S').strftime('%d/%m %H:%M')
        TIDSVEKTOR_DATUM.append(datum)

        




############
# Plotting #   
###################################################################################################


    fig = plt.figure(str(i), figsize=(12,5))

    ax1 = plt.subplot(111)

    ax1.plot(TIDSVEKTOR_MINUTER, TSK_TIDSSERIE, linewidth = 2)

    plt.tick_params(labelsize=15)
    plt.xticks(TIDSVEKTOR_MINUTER[0:len(TIDSVEKTOR_MINUTER):3], TIDSVEKTOR_DATUM[0:len(TIDSVEKTOR_DATUM):3], rotation= 45, size=6, ha='right', fontsize=15)
    plt.title('SKIN TEMPERATURE 2018-02-26 REFERENCE', fontsize=18)
    ax1.set_ylabel('Skin temperature (C)' , fontsize=16)
    ax1.set_xlabel('Time (UTC)', fontsize=18)
    #plt.ylim(0,1000)            #ÄNDRA
    #plt.xlim(0,36)
    #ax1.legend(loc=1,prop={'size': 12})
    ax1.grid(linestyle="dotted")


   
    fig.tight_layout()
    fig.subplots_adjust(right=0.90)












    '''VARIABLE = 'TSK'  #Endast för y-label

 

    fig = plt.figure(str(i), figsize=(12,5))

    ax1 = plt.subplot(111)

    fig .subplots_adjust(bottom=0.3)


    ax1.plot(TIDSVEKTOR_MINUTER, TSK_TIDSSERIE, linewidth = 2)


    plt.title('SKIN TEMPERATURE 2018-02-18 REFERENCE', fontsize=18)

    plt.xticks(TIDSVEKTOR_MINUTER, TIDSVEKTOR_DATUM  , rotation= 45, ha='right')

    plt.xlabel('Time', fontsize=16)

    plt.ylabel('Skin temperature (C)' , fontsize=16)#(Filobjekt.variables[VARIABLE].description + '    (' + Filobjekt.variables[VARIABLE].units + ')')

    plt.tick_params(labelsize=12)

    ax1.grid(linestyle="dotted")

    i+=1
'''

plt.show()



        
