import numpy as N
import scipy.io.netcdf as S
import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4 as N4


'''Detta skript ändrar blandningsförhållandet till mättnadsblandningsförhållandet på alla nivåer där det finns molnvatten.
Molnvattnet stoppas också in i wrfinput_d01-filen'''

###########
#VARIABLES#
###########

#Välj varifrån wrfinput_d01 ska tas ifrån. Från 3-dim, eller från SCM.

#Plockar wrfinput_d01 från den trdimensionella modellen


Filobjekt = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-27_00z_91_LEVELS/run/wrfinput_d01', mode='r')

#QCLOUD och QICE måste plockas från 3D-filen, för de finns ej i 1D/SCM

QCLOUD = Filobjekt.variables['QCLOUD'][0, :, 563, 352]
QICE = Filobjekt.variables['QICE'][0, :, 563, 352]


P = Filobjekt.variables['P'][0,:, 563, 352]
PB = Filobjekt.variables['PB'][0, :, 563, 352]
P_TOT = P+PB

PERT_THETA = Filobjekt.variables['T'][0, :, 563, 352]
THETA = PERT_THETA + 300
T_K = THETA/((1000/(P_TOT/100))**0.286)
T_C = T_K - 273.15

QVAPOR = Filobjekt.variables['QVAPOR'][0, :, 563, 352]
E_MATTNAD = N.exp(N.log(611.2)+(17.62*T_C/(243.12+T_C)))
QVAPOR_MATTNAD = 0.622*E_MATTNAD/(P_TOT-E_MATTNAD)
RH = QVAPOR/QVAPOR_MATTNAD
T_DAGGPUNKT = (243.12*N.log(611.2)-243.12*N.log(RH*E_MATTNAD)) / (N.log(RH*E_MATTNAD)-17.62-N.log(611.2))


Filobjekt.close()



#Plockar wrfinput_d01 från SCM, efter att ideal.exe har körts.

Filobjekt = N4.Dataset('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfinput_d01', mode='r')

P = Filobjekt.variables['P'][0,:, 1, 1]
PB = Filobjekt.variables['PB'][0,:, 1, 1]
P_TOT = P+PB

PERT_THETA = Filobjekt.variables['T'][0,:, 1, 1]
THETA = PERT_THETA + 300

T_K = THETA/((1000/(P_TOT/100))**0.286)
T_C = T_K - 273.15

QVAPOR = Filobjekt.variables['QVAPOR'][0,:, 1, 1]
E_MATTNAD = N.exp(N.log(611.2)+(17.62*T_C/(243.12+T_C)))
QVAPOR_MATTNAD = 0.622*E_MATTNAD/(P_TOT-E_MATTNAD)
RH = QVAPOR/QVAPOR_MATTNAD
T_DAGGPUNKT = (243.12*N.log(611.2)-243.12*N.log(RH*E_MATTNAD)) / (N.log(RH*E_MATTNAD)-17.62-N.log(611.2))

PH = Filobjekt.variables['PH'][0, :, 1, 1]
PHB = Filobjekt.variables['PHB'][0, :, 1, 1]
MODELLNIVAHOJD = (PH+PHB)/9.81
MODELLNIVAHOJD_TER = MODELLNIVAHOJD[:]# -MODELLNIVAHOJD[0]

MASSLEVELS = 0.5*(MODELLNIVAHOJD_TER[:-1] + MODELLNIVAHOJD_TER[1:])


Filobjekt.close()



#Skriver ut blandningsförhållande och respektive värdes plats i arrayen

'''for place, value in enumerate(QVAPOR):
    print(place, value)'''


#Om QCLOUD-värdet inte är lika med 0 så ska blandningsförhållandet höjas till mättnadsblandningsförhållande för att få RH=100%
QVAPOR_NEW = QVAPOR

for place, value in enumerate(MASSLEVELS):
    
    if value >1260 and value <1500:            #Om man tittar i en sondering direkt
    #if QCLOUD[place] >1e-5:                 #Om man ska använda nivåer i Centret, men initialisera med konstant Theta_liquid
        print(place,'     ', value, 1)
        QVAPOR_NEW[place] = QVAPOR_MATTNAD[place]
        



LIQ_WATER = N.zeros_like(THETA)
THETA_LIQ = THETA.copy()            #Om jag här skriver THETA_LIQ = THETA, så kommer även THETA ändras när jag senare ändrar THETA_LIQ........


#Plockar ut första nivån under moln och ser vad Theta är där. Här är Theta = Theta_liquid och denna vill jag ha konstant i molnet sedan.

for place, value in enumerate(MASSLEVELS):    
    #if QCLOUD[place] >1e-5:
    if value >1260 and value <1500:     # Lagt på 180 för terrängen 1080  1320
        print(place,'     ',value,2)
        THETA_BELOW_CLOUD = THETA[place-1]
        print(THETA_BELOW_CLOUD)
        break



#Genom att sedan ta skillnaden på denna konstanta liquid-temperatur och Theta i modellen får jag ut hur mycket molnvatten jag ska  1080 1320
#pytsa in på varje nivå, enligt formeln i Stull kapitel 13.


''' 100% adiabatisk kondensation'''  
'''
for place, value in enumerate(MASSLEVELS):
    #if QCLOUD[place] >1e-5:
    if value >1260 and value <1500:        
        LIQ_WATER[place] = (THETA[place] - THETA_BELOW_CLOUD)/((2.260e6*THETA[place])/(1004*T_K[place]))*0.25        
        if  LIQ_WATER[place] < 0:
            LIQ_WATER[place] = -LIQ_WATER[place]
        THETA_LIQ[place] = THETA_BELOW_CLOUD
        print(THETA_LIQ[place])                                                                               
'''


''' x% adiabatisk kondensation''' 


for place, value in enumerate(MASSLEVELS):
    #if QCLOUD[place] >1e-5:
    if value >1260 and value <1500:        
        LIQ_WATER[place] = (THETA[place] - THETA_BELOW_CLOUD)/((2.260e6*THETA[place])/(1004*T_K[place]))
        LIQ_WATER[place] = 0.25*LIQ_WATER[place]
        if  LIQ_WATER[place] < 0:
            LIQ_WATER[place] = -LIQ_WATER[place]
        THETA_LIQ[place] = THETA[place] - (THETA[place]/T_K[place])*(2.260e6/1004)*LIQ_WATER[place]
        THETA_BELOW_CLOUD = THETA_LIQ[place]
        print(THETA_LIQ[place]) 
        print(LIQ_WATER[place])
       



                               
#Öppnar netCDF-filen i "append"-läge

dataset = N4.Dataset('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfinput_d01', 'r+')



#Skriver in värden på QCLOUD och Wv-max på respektive plats i arrayen i netCDF-filen wrfinput_d01 som man får EFTER att ha kört
#ideal.exe. Detta motsvarar att köra fix_snow.ncl-skriptet.Vi skriver på varje vertikalnivå, men samma värde på alla
#x- och y-värden, eftersom det är "single column".

for place, value in enumerate(MASSLEVELS):
        #if QCLOUD[place] >1e-5:
        if value >1260 and value <1500:
            dataset['QCLOUD'][:, place, :, :] = LIQ_WATER[place]



for place, value in enumerate(MASSLEVELS):
        #if QCLOUD[place] >1e-5:
        if value >1260 and value <1500:
            dataset['QVAPOR'][:, place, :, :] = QVAPOR_NEW[place]




dataset.close()                                 #VIKTIGT att stänga filen för annars kan man inte öppna upp den med nya värden nedan. Den går att öppna, men med gamla QCLOUD-värden.



#Öppnar igen för att plocka fram nya värden på QCLOUD i Python.

dataset = N4.Dataset('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfinput_d01', 'r')


#Ominitialiserar QCLOUD för att kunna titta på den i Python
                               
QCLOUD = dataset['QCLOUD'][:]


