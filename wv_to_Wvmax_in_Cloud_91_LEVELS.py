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


Filobjekt = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-18_00z_91_LEVELS_QC_QI/run/wrfinput_d01', mode='r')

#QCLOUD och QICE måste plockas från 3D-filen, för de finns ej i 1D/SCM

QCLOUD = Filobjekt.variables['QCLOUD'][0, :, 563, 352]
print(QCLOUD)
QICE = Filobjekt.variables['QICE'][0, :, 563, 352]

'''
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
'''

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



Filobjekt.close()



#Skriver ut blandningsförhållande och respektive värdes plats i arrayen
print('QVAPOR')
for place, value in enumerate(QVAPOR):
    print(place, value)


#Om QCLOUD-värdet inte är lika med 0 så ska blandningsförhållandet höjas till mättnadsblandningsförhållande för att få RH=100%
QVAPOR_NEW = QVAPOR

print('QCLOUD')
for place, value in enumerate(QCLOUD):
    if value > 1e-5:
        print(place, value)
        QVAPOR_NEW[place] = QVAPOR_MATTNAD[place]
        

#Skriver ut blandningsförhållande och respektive värdes plats i arrayen igen för att se att värden ändrats
print('QVAPOR NEW')
for place, value in enumerate(QVAPOR_NEW):
    print(place, value)

                               


#Öppnar netCDF-filen i "append"-läge

dataset = N4.Dataset('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfinput_d01', 'r+')



#Skriver in värden på QCLOUD och QVAAPOR på respektive plats i arrayen i netCDF-filen wrfinput_d01 som man får EFTER att ha kört
#ideal.exe. Detta motsvarar att köra fix_snow.ncl-skriptet.Vi skriver på varje vertikalnivå, men samma värde på alla
#x- och y-värden, eftersom det är "single column".

for place, element in enumerate(QCLOUD):
    if element > 1e-5:
        dataset['QCLOUD'][:, place, :, :] = element


for place, element in enumerate(QVAPOR_NEW):
    dataset['QVAPOR'][:, place, :, :] = element

dataset.close()                                 #VIKTIGT att stänga filen för annars kan man inte öppna upp den med nya värden nedan. Den går att öppna, men med gamla QCLOUD-värden.



#Öppnar igen för att plocka fram nya värden på QCLOUD i Python.

dataset = N4.Dataset('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfinput_d01', 'r')


#Ominitialiserar QCLOUD för att kunna titta på den i Python
                               
QCLOUD = dataset['QCLOUD'][:] 
