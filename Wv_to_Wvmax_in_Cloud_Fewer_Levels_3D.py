import numpy as N
import scipy.io.netcdf as S
import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4 as N4



###########
#VARIABLES#
###########


#Plockar wrfinput_d01 från den trdimensionella modellen

Filobjekt = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-26_00z_Wv_to_Wvmax_0_Icloud2/run/wrfinput_d01', mode='r+')


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

Filobjekt.close()



#Om QCLOUD-värdet inte är lika med 0 så ska blandningsförhållandet höjas till mättnadsblandningsförhållande för att få RH=100%
QVAPOR_NEW = QVAPOR


for j in range(QCLOUD.shape[1]): #500
    for i in range(QCLOUD.shape[2]): #250
        for k, value in enumerate(QCLOUD[:, j, i]):
            if value >= 0:
                #print(k, j, i,  QCLOUD[k, j, i])
                #print(T_C[k, j, i])
                #print(T_DAGGPUNKT[k, j, i])
                QVAPOR_NEW[k, j, i] = QVAPOR_MATTNAD[k, j, i]
                #print(T_C[k, j, i])




#Öppnar netCDF-filen på nytt i "append"-läge

dataset = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-26_00z_Wv_to_Wvmax_0_Icloud2/run/wrfinput_d01', 'r+')


#Skriver in nya världen för temperaturen i de gridpunkter där den sänkts till daggpunktstemperaturen.
#Skilladen mot SCM-fallet är att här finns redan QCLOUD-värdena så de behöver inte skrivas in, men de som är mindre än ett visst värde plockas bort.
#vi väljer att bara sänka temperaturen på de nivåer där det finns mycket molnvatten.

for j in range(QCLOUD.shape[1]):
    for i in range(QCLOUD.shape[2]):
        for k, value in enumerate(QCLOUD[:, j, i]):
            if value > 0:
                dataset['QVAPOR'][0, k, j, i] = QVAPOR_NEW[k, j, i]
            if value < 0:
                dataset['QCLOUD'][0, k, j, i] = 0

dataset.close()         #VIKTIGT att stänga filen för annars kan man inte öppna upp den med nya värden nedan. Den går att öppna, men med gamla QCLOUD-värden.
    
                 





dataset_new = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-26_00z_Wv_to_Wvmax_0_Icloud2/run/wrfinput_d01', 'r+')

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


