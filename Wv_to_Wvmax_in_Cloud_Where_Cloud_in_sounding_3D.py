import numpy as N
import scipy.io.netcdf as S
import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4 as N4



###########
#VARIABLES#
###########


#Plockar wrfinput_d01 från den trdimensionella modellen

Filobjekt = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-27_00z_91_LEVELS_Clouds_Where_Clouds_In_Sounding_YSU_TOPD0_100_KOPIA/run/wrfinput_d01', mode='r+')


QCLOUD = Filobjekt.variables['QCLOUD'][0, :, :, :]
QICE = Filobjekt.variables['QICE'][0, :, :, :]


P = Filobjekt.variables['P'][0,:, :, :]
PB = Filobjekt.variables['PB'][0, :, :, :]
P_TOT = P+PB

PERT_THETA = Filobjekt.variables['T'][0, :, :, :]
THETA = PERT_THETA + 300
T_K = THETA/((1000/(P_TOT/100))**0.286)
T_C = T_K - 273.15

QVAPOR = Filobjekt.variables['QVAPOR'][0, :, :, :]
E_MATTNAD = N.exp(N.log(611.2)+(17.62*T_C/(243.12+T_C)))
QVAPOR_MATTNAD = 0.622*E_MATTNAD/(P_TOT-E_MATTNAD)
RH = QVAPOR/QVAPOR_MATTNAD
T_DAGGPUNKT = (243.12*N.log(611.2)-243.12*N.log(RH*E_MATTNAD)) / (N.log(RH*E_MATTNAD)-17.62-N.log(611.2))

PH = Filobjekt.variables['PH'][0, :, :, :]
PHB = Filobjekt.variables['PHB'][0, :, :, :]
MODELLNIVAHOJD = (PH+PHB)/9.81
MODELLNIVAHOJD_TER = MODELLNIVAHOJD[:, :, :] -MODELLNIVAHOJD[0, :, :]

MASSLEVELS = 0.5*(MODELLNIVAHOJD_TER[:-1, :, :] + MODELLNIVAHOJD_TER[1:, :, :])



Filobjekt.close()



#Om QCLOUD-värdet inte är lika med 0 så ska blandningsförhållandet höjas till mättnadsblandningsförhållande för att få RH=100%
QVAPOR_NEW = QVAPOR.copy()


for j in range(500, QCLOUD.shape[1]):
    for i in range(250, QCLOUD.shape[2]):
        for k, value in enumerate(MASSLEVELS[:, j, i]):
            #if value >= 1e-5:
            if value >1000 and value <1400:
                #print(k, j, i,  QCLOUD[k, j, i])
                #print(T_C[k, j, i])
                #print(T_DAGGPUNKT[k, j, i])
                QVAPOR_NEW[k, j, i] = QVAPOR_MATTNAD[k, j, i]
                #print(T_C[k, j, i])

                

LIQ_WATER = N.zeros_like(THETA)
THETA_LIQ = THETA.copy()
THETA_BELOW_CLOUD = N.zeros_like(THETA)




#Plockar ut första nivån under moln och ser vad Theta är där. Här är Theta = Theta_liquid och denna vill jag ha konstant i molnet sedan.


for j in range(500, QCLOUD.shape[1]):
    for i in range(250, QCLOUD.shape[2]):
        for k, value in enumerate(MASSLEVELS[:, j, i]):
    
            #if value >= 1e-5 and QCLOUD[k-1, j, i] < 1e-5:
            if value >1000 and value <1400:
                THETA_BELOW_CLOUD = THETA[k-1, j, i]
                LIQ_WATER[k, j, i] = (THETA[k, j, i] - THETA_BELOW_CLOUD)/((2.260e6*THETA[k, j, i])/(1004*T_K[k, j, i]))
                #THETA[k, j, i] = THETA_BELOW_CLOUD      #Denna rad är till för att ta bort stabiliteten där molnet har initialiserats, dvs sätta rheta konstant med höjden.
                THETA_LIQ[k, j, i] = THETA_BELOW_CLOUD

                
            '''elif value >= 1e-5 and QCLOUD[k-1, j, i] > 1e-5:
                LIQ_WATER[k, j, i] = (THETA[k, j, i] - THETA_BELOW_CLOUD)/((2.260e6*THETA[k, j, i])/(1004*T_K[k, j, i]))    #om det finns moln på nivån under definieras ingeet nytt Theta_Below Cloud
                THETA_LIQ[k, j, i] = THETA_BELOW_CLOUD'''




            
#Det är Pert Theta som skrivs in nedan om vi hoppar över nedanstående avsnitt
            

PERT_THETA = THETA-300
'''
################
#OBS!!!!!!!!!! #
######################################################################################################
 
# Efter att ha ändrat Theta till att vara mer omblandat måste jag ändra daggpunkten igen så att inte temperaturen blir lägre än daggpunkten.
# Hela detta avsnitt kan tas bort om inte molnet ska göras mera omblandat

Filobjekt = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-27_00z_91_LEVELS_Clouds_Where_Clouds_In_Sounding_YSU_TOPD1_100_Reduced_Stability_in_clouds/run/wrfinput_d01', mode='r+')

#PERT_THETA = Filobjekt.variables['T'][0, :, :, :]
THETA = PERT_THETA + 300
T_K = THETA/((1000/(P_TOT/100))**0.286)
T_C = T_K - 273.15

QVAPOR = Filobjekt.variables['QVAPOR'][0, :, :, :]
E_MATTNAD = N.exp(N.log(611.2)+(17.62*T_C/(243.12+T_C)))
QVAPOR_MATTNAD = 0.622*E_MATTNAD/(P_TOT-E_MATTNAD)
RH = QVAPOR/QVAPOR_MATTNAD
T_DAGGPUNKT = (243.12*N.log(611.2)-243.12*N.log(RH*E_MATTNAD)) / (N.log(RH*E_MATTNAD)-17.62-N.log(611.2))

for j in range(500, QCLOUD.shape[1]):
    for i in range(250, QCLOUD.shape[2]):
        for k, value in enumerate(MASSLEVELS[:, j, i]):
            #if value >= 1e-5:
            if value >1000 and value <1400:
                #print(k, j, i,  QCLOUD[k, j, i])
                #print(T_C[k, j, i])
                #print(T_DAGGPUNKT[k, j, i])
                QVAPOR_NEW[k, j, i] = QVAPOR_MATTNAD[k, j, i]
                #print(T_C[k, j, i])

Filobjekt.close()

#######################################################################################################'''



#Öppnar netCDF-filen på nytt i "append"-läge

dataset = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-27_00z_91_LEVELS_Clouds_Where_Clouds_In_Sounding_YSU_TOPD0_100_KOPIA/run/wrfinput_d01', 'r+')


#Skriver in nya världen för blandningsförhållandet i de gridpunkter där det höjts till mättnadsblandningsförhållandet.
#För QCLOUD skrivs LIQ_WATER-värdena in.


for j in range(500, QCLOUD.shape[1]):
    for i in range(250, QCLOUD.shape[2]):
        for k, value in enumerate(MASSLEVELS[:, j, i]):
            #if value > 1e-5:
            if value >1000 and value <1400:
                dataset['QVAPOR'][0, k, j, i] = QVAPOR_NEW[k, j, i]                
            #if value > 1e-5:
            if value >1000 and value <1400:
                dataset['QCLOUD'][0, k, j, i] = LIQ_WATER[k, j, i]
                dataset['T'][0, k, j, i] = PERT_THETA[k, j, i]
            #else:
            #    dataset['QCLOUD'][0, k, j, i] = 0
                

dataset.close()         #VIKTIGT att stänga filen för annars kan man inte öppna upp den med nya värden nedan. Den går att öppna, men med gamla QCLOUD-värden.
    
                 





dataset_new = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-27_00z_91_LEVELS_Clouds_Where_Clouds_In_Sounding_YSU_TOPD0_100_KOPIA/run/wrfinput_d01', 'r+')

#dataset_old = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-19_00z_91_LEVELS_QC_QI_42h/run/wrfinput_d01_original', 'r+')


P = dataset_new.variables['P'][0,:, :, :]
PB = dataset_new.variables['PB'][0, :, :, :]
P_TOT = P+PB

PERT_THETA = dataset_new.variables['T'][0, :, :, :]
THETA = PERT_THETA + 300
T_K = THETA/((1000/(P_TOT/100))**0.286)
T_C = T_K - 273.15

QVAPOR = dataset_new.variables['QVAPOR'][0, :, :, :]
E_MATTNAD = N.exp(N.log(611.2)+(17.62*T_C/(243.12+T_C)))
QVAPOR_MATTNAD = 0.622*E_MATTNAD/(P_TOT-E_MATTNAD)
RH = QVAPOR/QVAPOR_MATTNAD
T_DAGGPUNKT = (243.12*N.log(611.2)-243.12*N.log(RH*E_MATTNAD)) / (N.log(RH*E_MATTNAD)-17.62-N.log(611.2))


dataset_new.close()


