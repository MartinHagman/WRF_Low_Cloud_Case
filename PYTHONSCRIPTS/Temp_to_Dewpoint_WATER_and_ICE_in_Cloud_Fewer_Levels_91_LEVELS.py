import numpy as N
import scipy.io.netcdf as S
import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4 as N4

'''91 LEVELS'''

'''Detta skript ändrar temperaturen till daggpunktstemperaturen på alla nivåer där det finns molnvatten. Molnvattnet stoppas också in
i wrfinput_d01-filen'''

###########
#VARIABLES#
###########

#Välj varifrån wrfinput_d01 ska tas ifrån. Från 3-dim, eller från SCM.

#Plockar wrfinput_d01 från den trdimensionella modellen


Filobjekt = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-18_00z_91_LEVELS_QC_QI/run/wrfinput_d01', mode='r')

#QCLOUD och QICE måste plockas från 3D-filen, för de finns ej i 1D/SCM

PH = Filobjekt.variables['PH'][0, :, 563, 352]
PHB = Filobjekt.variables['PHB'][0, :, 563, 352]
MODELLNIVAHOJD = (PH+PHB)/9.81
MODELLNIVAHOJD_TER = MODELLNIVAHOJD -MODELLNIVAHOJD[0]


MASSLEVELS = 0.5*(MODELLNIVAHOJD_TER[:-1] + MODELLNIVAHOJD_TER[1:])

QCLOUD = Filobjekt.variables['QCLOUD'][0, :, 563, 352]
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


'''Mättnad över vatten'''

QVAPOR = Filobjekt.variables['QVAPOR'][0,:, 1, 1]
E_MATTNAD = N.exp(N.log(611.2)+(17.62*T_C/(243.12+T_C)))
QVAPOR_MATTNAD = 0.622*E_MATTNAD/(P_TOT-E_MATTNAD)
RH = QVAPOR/QVAPOR_MATTNAD
T_DAGGPUNKT = (243.12*N.log(611.2)-243.12*N.log(RH*E_MATTNAD)) / (N.log(RH*E_MATTNAD)-17.62-N.log(611.2))


'''Mättnad över is'''

QVAPOR = Filobjekt.variables['QVAPOR'][0,:, 1, 1]
E_MATTNAD_ICE = N.exp(N.log(611.2)+(22.46*T_C/(272.62+T_C)))
QVAPOR_MATTNAD_ICE = 0.622*E_MATTNAD_ICE/(P_TOT-E_MATTNAD_ICE)
RH_ICE = QVAPOR/QVAPOR_MATTNAD_ICE
T_DAGGPUNKT_ICE = (272.62*N.log(611.2)-272.62*N.log(RH_ICE*E_MATTNAD_ICE)) / (N.log(RH_ICE*E_MATTNAD_ICE)-22.46-N.log(611.2))




Filobjekt.close()



#Skriver ut temperatur och respektive värdes plats i arrayen

print('TEMPERATURE')
for place, value in enumerate(T_C):
    print(place, value)


#Om QCLOUD-värdet inte är lika med 0 så ska temperaturen sänkas till daggpunkten för att få RH=100%
#Om QICE-värdet inte är lika med 0 samtidigt som QCLOUD-värdet är lika med 0 så ska temperaturen sänkas till daggpunkten för is för att få RH=100%

print('CLOUD WATER and CLOUD ICE')
for place, (value_water, value_ice) in enumerate(zip(QCLOUD, QICE)):
    if value_water > 1e-4:
        print('WATER', place, value_water, MASSLEVELS[place])
        T_C[place] = T_DAGGPUNKT[place]

    '''if value_water == 0 and value_ice > 0:
        print('ICE', place, value_ice, MASSLEVELS[place])
        T_C[place] = T_DAGGPUNKT_ICE[place]'''



#Skriver ut temperatur och respektive värdes plats i arrayen igen för att se att värden ändrats

print('CHANGED TEMPERATURE')
for place, value in enumerate(T_C):
    print(place, value)


#Återgår till potentiell temperatur, då detta ska skrivas in i wrfinput_d01 igen    


T_K = T_C + 273.15
THETA = T_K*((1000/(P_TOT/100))**0.286)
PERT_THETA = THETA-300                           

        
    

#Skriver ut även daggpunkts-arrayen.

print('DEW POINT')
for place, value in enumerate(T_DAGGPUNKT):
    print(place, value)
                               


#Öppnar netCDF-filen i "append"-läge

dataset = N4.Dataset('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfinput_d01', 'r+')



#Skriver in värden på QCLOUD och THETA på respektive plats i arrayen i netCDF-filen wrfinput_d01 som man får EFTER att ha kört
#ideal.exe. Detta motsvarar att köra fix_snow.ncl-skriptet.Vi skriver på varje vertikalnivå, men samma värde på alla
#x- och y-värden, eftersom det är "single column".

for place, value in enumerate(QCLOUD):

    if value > 0:
        dataset['QCLOUD'][:, place, :, :] = value


for place, value in enumerate(QICE):

    if value > 0:
        dataset['QICE'][:, place, :, :] = value
        


for place, value  in enumerate(PERT_THETA):
    dataset['T'][:, place, :, :] = value

dataset.close()                                 #VIKTIGT att stänga filen för annars kan man inte öppna upp den med nya värden nedan. Den går att öppna, men med gamla QCLOUD-värden.



#Öppnar igen för att plocka framnya värden på QCLOUD i Python.

dataset = N4.Dataset('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfinput_d01', 'r')


#Ominitialiserar QCLOUD för att kunna titta på den i Python
                               
QCLOUD_CHANGED = dataset['QCLOUD'][:]
QICE_CHANGED = dataset['QICE'][:]
