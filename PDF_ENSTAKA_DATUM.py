import numpy as N
import scipy.io.netcdf as S
import netCDF4 as N4
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import sys 
import pandas as pd
import os
import datetime
from matplotlib.colors import LogNorm
import matplotlib.colors as mplc
from pylab import *


########################################
#OBSERVATIONER FRÅN SODANKYLÄ I JANUARI#
########################################


df = pd.read_csv(r'/home/sm_marha/TEXTFILER/Automatic_Station_jan_feb_all_parameters.txt', encoding='latin-1', delimiter=',')

#df_time_cloudbase = df.iloc[:,[4,19]]

df = df.loc[:,['DATE TIME','CH1']]

df_numpy_array = df.values



TIDSVEKTOR_OBS = df_numpy_array[1:, 0]





TIDSVEKTOR_MINUTER_OBS = []

for datum in TIDSVEKTOR_OBS:
    datum = datum[0:16]
    datum = datetime.datetime.strptime(datum,'%Y-%m-%d %H:%M').strftime('%d/%m %H:%M')
    TIDSVEKTOR_MINUTER_OBS.append(datum)

MOLNBAS_TIDSSERIE_OBS = df_numpy_array[1:, 1]

TIDSVEKTOR_TIMMAR_OBS = TIDSVEKTOR_MINUTER_OBS[0::6]


#Plockar ut varje hel timme
MOLNBAS_GLES_OBS = MOLNBAS_TIDSSERIE_OBS[0::6]
#Gör om från sträng till ifloat
MOLNBAS_GLES_OBS = pd.to_numeric(MOLNBAS_GLES_OBS, errors='coerce')

#Sätter värdet 100 när molnbas är NaN, vilker detsamma som klart i detta fall
for place, element in enumerate(MOLNBAS_GLES_OBS):
    if element != element:
        MOLNBAS_GLES_OBS[place] = -200

   

#Gör en tidsvektor från 0 till 1415
TIDSVEKTOR_OBS = range(len(MOLNBAS_GLES_OBS))   #Värde 0:1416


#Plockar ut det datum jag vill ha.I detta fall dygn 57 0ch 49, dvs ANDRA SIFFRAN
MOLNBAS_GLES_OBS_1 = MOLNBAS_GLES_OBS[56*24:57*24]                                      #ÄNDRA
#MOLNBAS_GLES_OBS_2 = MOLNBAS_GLES_OBS[48*24:49*24]                                     #ÄNDRA    

#Slår ihop dessa till en lång array
#MOLNBAS_GLES_OBS_TOT = N.concatenate((MOLNBAS_GLES_OBS_1, MOLNBAS_GLES_OBS_2))         #ÄNDRA
MOLNBAS_GLES_OBS_TOT = MOLNBAS_GLES_OBS_1#ÄNDRA
MOLNBAS_GLES_OBS_TOT = MOLNBAS_GLES_OBS_TOT[0:]
print(MOLNBAS_GLES_OBS_TOT)


print(len(TIDSVEKTOR_TIMMAR_OBS))
print(len(MOLNBAS_GLES_OBS))
print(len(MOLNBAS_GLES_OBS_TOT))



###############################
#AKTUELL WRF KÖRNING/KÖRNINGAR#
###############################

#Initialiserar lista där jag sedan spar undan molnbaserna från ALLA patherna jag loopar över nedan




PATH_LIST = ['/nobackup/smhid12/sm_marha/2018-02-18_00z_91_LEVELS_42h/run/',
             '/nobackup/smhid12/sm_marha/2018-02-18_00z_91_LEVELS_QC_QI_T_to_Td_1e-4_Icloud_2/run/']

'''
PATH_LIST = ['/nobackup/smhid12/sm_marha/2018-02-26_00z_91_LEVELS_Icloud1/run/',
             '/nobackup/smhid12/sm_marha/2018-02-26_00z_91_LEVELS_QC_QI_24h_T_to_Td_Whole_domain_0_icloud2/run/']'''

'''
PATH_LIST = ['/nobackup/smhid12/sm_marha/2018-02-20_00z_91_LEVELS_15h/run/',
             '/nobackup/smhid12/sm_marha/2018-02-20_00z_91_LEVELS_QC_QI_18h_T-to_Td_0_Icloud2/run/']'''


Figure_letter_list= ['a', 'b']

#Loopar över patherna ovan

i=0
for PATH in PATH_LIST:

    CLOUDBASE_TIDSSERIE_WRF = []

    print(PATH)

    file_list = os.listdir(PATH)

    



    wrfout_list = []

    for file in file_list:
        if file.startswith('wrfout'):
            wrfout_list.append(file)

    wrfout_list.sort()
    wrfout_list = wrfout_list[0:24]         #Välj hur många filer Du vill ha om det finns fler i mappen.

    VARIABEL = 'CLDFRA'


#Plockar ut första nivån där det finns CLDFRA >0.8
    

    for place, file in enumerate(wrfout_list):
        print('file '+ str(place))
        f = PATH+file
        Filobjekt = N4.Dataset(f, mode='r')
        VARIABELVARDE = Filobjekt.variables[VARIABEL][:]
        VARIABELVARDE_1D = VARIABELVARDE[0, :, 563, 352]
        P_PERT = Filobjekt.variables['P'][:]
        P_BASE = Filobjekt.variables['PB'][:]
        P = P_BASE+P_PERT
        PH = Filobjekt.variables['PH'][:]
        PHB = Filobjekt.variables['PHB'][:]
        MODELLNIVAHOJD = (PH+PHB)/9.81
        MASSLEVELS = 0.5*(MODELLNIVAHOJD[0,:-1, 563, 352] + MODELLNIVAHOJD[0, 1:, 563, 352])


        for k in range(len(VARIABELVARDE_1D)):
            if VARIABELVARDE_1D[k] >=0.8:
                CLOUDBASE = MASSLEVELS[k]-MODELLNIVAHOJD[0,0,563,352]
                break
            else:
                CLOUDBASE = -200

    
        CLOUDBASE_TIDSSERIE_WRF.append(CLOUDBASE)

       
    CLOUDBASE_TIDSSERIE_WRF = N.array(CLOUDBASE_TIDSSERIE_WRF)
    print(CLOUDBASE_TIDSSERIE_WRF)
    

#####################
#PLOTTNING HISTOGRAM#
#####################

    if i == 0:
        rc('axes', linewidth=1.5)
        fig1, axs = plt.subplots(1, 2, figsize = (9,4))
        fig1.subplots_adjust(bottom = 0.2)
        fig1.subplots_adjust(hspace = 0.1)
        fig1.subplots_adjust(right = 0.85)



    '''Deklarerar variabler för de räta linjerna i figuren'''

    X = N.arange(0,10000)
    Y = N.arange(0,10000)
    X_0 = N.zeros(10000)
    Y_0 = N.zeros(10000)


    ax = axs[0]

    #Plottar räta linjer i figuren
    axs[i].plot(X, Y , 'r', linewidth = 2)
    axs[i].plot(X_0, Y , 'k', linewidth = 2)
    axs[i].plot(X, Y_0 , 'k', linewidth = 2)


    norm = mplc.Normalize(0, 10)  # Linjär normerin
    norm = mplc.LogNorm(1, 10)  # Logaritmisk normering
    
    hsist = axs[i].hist2d(MOLNBAS_GLES_OBS_TOT, CLOUDBASE_TIDSSERIE_WRF, bins=(22,22),range=[[-200, 2000], [-200, 2000]], norm = norm, cmap = 'PuRd') #0 ger första dygnet osv
    #cbar = fig1.colorbar(hsist[3], ax = axs[i], format = '%4d')    # Egen bar
    axs[i].set_ylim(-200, 2000)
    axs[i].set_xlim(-200, 2000)
    #ax.set_title('180226 OBS vs QC QI Reference Icloud2,bins(28, 32)', fontsize = 16)
    #ax.set_title('180226 OBS vs reference, Icloud2, bins(28, 32)', fontsize = 16)
    axs[i].set_xlabel('Cloudbase Observation (m)', fontsize = 15)
    axs[0].set_ylabel('Cloudbase WRF (m)', fontsize = 15)
    axs[i].tick_params(which='major', direction='in', labelsize = 15)
    axs[i].text(-180, 1850, Figure_letter_list[i], fontsize=15, fontweight='bold')
    #axs[i].tick_params(labelsize=16)
    #cbar.set_label('Number of hits', fontsize = 16)     # Egen bar
    #cbar.ax.tick_params(labelsize=13)  #Egen bar

    i+=1
    
axs[1].set_yticklabels([])
cbar_ax = fig1.add_axes([0.87, 0.20, 0.01, 0.67])       #Gemensam bar
cbar = fig1.colorbar(hsist[3], cax=cbar_ax, format = '%4d')
cbar.ax.tick_params(labelsize=14)
cbar.set_label('Number of hits', fontsize = 17)

plt.savefig('/home/sm_marha/FIGURES/FINAL/20180218/Molnbas_WRF_vs_OBS_Sodankyla__QC_QI_24h_T_to_Td_Whole_domain_0_icloud2_PDF_20180218_ENSTAKA')




plt.show()
