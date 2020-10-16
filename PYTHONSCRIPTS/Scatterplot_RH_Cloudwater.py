import numpy as N
import scipy.io.netcdf as S
import netCDF4 as N4
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import os
import pandas as pd
from math import *
import sys



'''Detta skript beräknar molnklimatologin för WRF och ECMWF.'''



'''
###########
# WRF OUT #
#######################################################################################


file_list = os.listdir('/nobackup/smhid12/sm_marha/2018-01-18_00z/run/')


wrfout_list = []

for file in file_list:
    if file.startswith('wrfout'):
        wrfout_list.append(file)

wrfout_list.sort()


wrfout_list = wrfout_list[6:]   #Väljer ut vilka tidssteg jag vill ha, då alla tidssteg ligger i olika filer från WRF.


#For-loop som loopar över de olika tidssteg jag gjort ovan

for file in wrfout_list:
    print('New file!')

    Filobjekt = S.netcdf_file('/nobackup/smhid12/sm_marha/2018-01-18_00z/run/'+ file, mode='r')
    

    PERT_THETA = Filobjekt.variables['T'][0, :, :, :]       #Plockar ut det nollte elementet från tidsdimensionen, vilket är det enda, eftersom det bara ligger 1 tid i varje fil.

    THETA = PERT_THETA+300

    P_PERT = Filobjekt.variables['P'][0, :, :, :]
    P_BASE = Filobjekt.variables['PB'][0, :, :, :]
    P= P_BASE+P_PERT

    TEMP = THETA/((1000/(P/100))**0.286)
    TEMP_C = TEMP-273.15

    QVAPOR = Filobjekt.variables['QVAPOR'][0, :, :, :]
    E_MATTNAD = N.exp(N.log(611.2)+(17.62*TEMP_C/(243.12+TEMP_C)))
    QVAPOR_MATTNAD = 0.622*E_MATTNAD/(P-E_MATTNAD)
    RH = QVAPOR/QVAPOR_MATTNAD


    QICE = Filobjekt.variables['QICE'][0, :, :, :]*1000
    QCLOUD = Filobjekt.variables['QCLOUD'][0, :, :, :]*1000

    fig = plt.figure(1)                                     #Eftersom plottningen ligger i loopen och vi plottar i ax hela tiden hamnar alla tidssteg i samma plot, även om den inte visas i varje steg.
    ax = plt.subplot(111)
    ax.scatter(RH[1:20, :, :], QICE[1:20, :, :], s=0.5)         #Här väljer jag hur många vertikala nivåer jag vill plotta. Är ju mest intresserad av St/sc-moln.

plt.title('WRF cloud climatology level 1-20, timestep 6 to 36' ) #All timesteps här är från utfil sex, efter spin-upen är klar. Alltså INTE från början.
plt.xlabel('RH (%)')
plt.ylabel('Cloud ice (g/kg)')

#ax.set_xlim(0.8, 1.0)
#ax.set_ylim(-0.1, 1)        #Vid interpoleringen uppstår negativa cloud water-värden. Därför minus på y-axeln.

fig.savefig('/home/sm_marha/FIGURES/FOKUS_18_20_JAN/18_JAN/Cloudice_clim_WRF_timestep_ALL_timestep.png')


sys.exit()'''


   

################
# ECMWF HYBRID #
#######################################################################################

'''Plockar ut variabler från GRIB-filerna'''


Filobjekt_GRIB = S.netcdf_file('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1ua2018011800_M_D.nc', mode='r')
Filobjekt_GRIB_sfc = S.netcdf_file('/nobackup/smhid12/sm_marha/EGNA_GRIBFILER/wrf_mlth1sfc2018011800_M_D.nc', mode='r')

#Här ligger alla tidssteg i samma fil, så därför får vi loopa över tidsdimensionen, i, nedan.


for i in range(0, 1):             #range(1:5) =1, 2, 3, 4 osv....

    print('New ECMWF file!')
    QICE_EC = Filobjekt_GRIB.variables['ciwc'][i, ::-1, :, :]*1000

    QCLOUD_EC = Filobjekt_GRIB.variables['clwc'][i, ::-1, :, :]*1000
    QCLOUD_FRACTION_EC = Filobjekt_GRIB.variables['cc'][i, ::-1, :, :]
    
    TEMPERATURE_EC = Filobjekt_GRIB.variables['t'][i, ::-1, :, :]
    TEMPERATURE_C_EC = TEMPERATURE_EC-273.15
    
    SPECIFIC_HUMIDITY_EC = Filobjekt_GRIB.variables['q'][i,::-1, :, :]

    MIXING_RATIO_EC = SPECIFIC_HUMIDITY_EC/(1-SPECIFIC_HUMIDITY_EC)

    E_MATTNAD_EC = N.exp(N.log(611.2)+(17.62*TEMPERATURE_C_EC/(243.12+TEMPERATURE_C_EC)))


    SURFACE_PRESSURE_EC = Filobjekt_GRIB_sfc.variables['sp'][i, :, :]



    '''BERÄKNING AV DE VERTIKALA NIVÅERNA I ECMWF'''


    
    '''Läser in a-. och b-koeffecienter'''

    df = pd.read_csv(r'/home/sm_marha/TEXTFILER/A_B_Coefficient_ECMWF.csv', encoding='latin-1', delimiter=';', header = None)



    df = df.iloc[:,:]

    df_numpy_array = df.values

    A= df_numpy_array[:, 1]
    B= df_numpy_array[:, 2]




    '''Räknar ut virtuella temperaturen'''

    T_VIRTUAL_EC = (1.+0.609133*SPECIFIC_HUMIDITY_EC)*TEMPERATURE_EC



    '''Räknar ut de halva trycknivåerna(edges of the layer. U, V, T, q is valid in the middle of each layer. https://rda.ucar.edu/datasets/ds115.4/docs/levels.hybrid.html'''

    p_half = N.empty([len(A), 113, 269])

    for i in range(len(A)):

        p = A[i] + B[i] * SURFACE_PRESSURE_EC

        p_half[i, :, :] = p


     


    '''Räknar ut hela trycknivåerna'''
     
    p_full = 0.5*(p_half[1:, :, :]+p_half[:-1, :, :])   #Väljer först ut alla element utom det första i z-led olch sedan alla utom det sista i z-led, adderar och delar på 2.



    '''Räknar ut geopotentialen på de fulla nivåerna.'''

    Rd = 287.06
    g0 = 9.80665



    Fi = N.zeros_like(p_full)
    Fi[0, :, :]=10*g0        #Sätter geopotentialhöjden till 10m ggr g på lägsta fulla nivån. Detta behövs för att stara loopen nedan.


    for i in range(len(p_full[:, 112, 268])-1):     #Tar bara en random horisontell punkt, då det enda jag är intresserad av är längden vertikalt. Denna är samma i alla horisontella punkter.

        Fi[i+1, :, :] = Fi[i, :, :] + Rd*(  T_VIRTUAL_EC[i, :, :]  *  (  N.log(p_full[i, :, :])-N.log(p_full[i+1, :, :])  ) ) #Beräknar geopotentalen nerifrån och upp poå alla vertikala nivåer.
                              
   

    

    FULL_LEVELS = Fi/g0


    MIXING_RATIO_MATTNAD_EC = 0.622*E_MATTNAD_EC/(p_full-E_MATTNAD_EC)

    RH_EC = MIXING_RATIO_EC/MIXING_RATIO_MATTNAD_EC

    
    fig = plt.figure(1)                                             #På samma sätt som i WRF ovan plottar vi i for-loopen i samma ax i figur2, men figuren visas ej förrns efter loopen.
    ax = plt.subplot(111)
    ax.scatter(RH_EC[0:30, :, :], QICE_EC[0:30, :, :], s=0.5)



plt.title('ECMWF cloud climatology analysis, vertical level 1-29'  )
plt.xlabel('RH (%)')
plt.ylabel('Cloud ice (g/kg)')

#ax.set_xlim(0.5, 1.1)
#ax.set_ylim(-0.1, 1)
fig.savefig('/home/sm_marha/FIGURES/FOKUS_18_20_JAN/18_JAN/Cloudice_clim_ECMWF_Analysis.png')

plt.show()
            

   

    





    
