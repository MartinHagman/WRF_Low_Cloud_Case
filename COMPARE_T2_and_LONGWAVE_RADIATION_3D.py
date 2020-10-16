import numpy as N
import scipy.io.netcdf as S
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import os
import datetime  
import netCDF4 as N4
import pandas as pd



###############
# WRF FILES   #
#################################################################################################################
print('WRF FILES')




'''Ordnar lista med med de olika WRF-körningarna som senare loopas över. '''


wrfout_list = []

file_list = os.listdir('/nobackup/smhid12/sm_marha/2018-02-18_00z_91_LEVELS_42h/run')         #ÄNDRA  Detta är bara en lista av körningar. Dessa är lika när datumet är lika och behöver ej ändras mellan REF och ASS

for file in file_list:
    if file.startswith('wrfout_d01'):
        wrfout_list.append(file)

wrfout_list.sort()

wrfout_list = wrfout_list[0:24]        #Om körningarna innehåller fler än 24 tidssteg


'''Initialisering av tomma arrays'''

T2_C_TIME_SERIES_WRF_REF = N.empty_like(wrfout_list)
T2_C_TIME_SERIES_WRF_ASS = N.empty_like(wrfout_list)

GLW_TIME_SERIES_WRF_REF = N.empty_like(wrfout_list)
GLW_TIME_SERIES_WRF_ASS = N.empty_like(wrfout_list)




'''Loopar över respektive wrfout-fil i 2 olika directoryn och lagrar de aktuella variabelvärdena i de tomma arrayerna.'''


for place, file in enumerate(wrfout_list):

        print(file)


        '''FILOBJEKTEN samt variabler från både REF och ASDS som ej behöver beräknas'''
        
        Filobjekt_WRF_REF = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-18_00z_91_LEVELS_42h/run/' + file, mode='r')       #ÄNDRA
        Filobjekt_WRF_ASS = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-18_00z_91_LEVELS_QC_QI_T_to_Td_0_Whole_Dom/run/' + file, mode='r')       #ÄNDRA

        LONGWAVE_DOWN_SURFACE_WRF_REF = Filobjekt_WRF_REF.variables['GLW'][0, 563, 352]
        LONGWAVE_DOWN_SURFACE_WRF_ASS = Filobjekt_WRF_ASS.variables['GLW'][0, 563, 352]

        T2_K_WRF_REF = Filobjekt_WRF_REF.variables['T2'][0, 563, 352]
        T2_K_WRF_ASS = Filobjekt_WRF_ASS.variables['T2'][0, 563, 352]

        T2_C_WRF_REF = T2_K_WRF_REF-273.15
        T2_C_WRF_ASS = T2_K_WRF_ASS-273.15

        
        '''WRF_REF'''
        
        THETA_PERT_WRF_REF = Filobjekt_WRF_REF.variables['T'][0, :, 563, 352]
        P_PERT_WRF_REF = Filobjekt_WRF_REF.variables['P'][0, :, 563, 352]
        P_BASE_WRF_REF = Filobjekt_WRF_REF.variables['PB'][0, :, 563, 352]   
        
        P_TOT_WRF_REF= P_BASE_WRF_REF + P_PERT_WRF_REF
        THETA_WRF_REF = THETA_PERT_WRF_REF+300
        TEMP_K_WRF_REF = THETA_WRF_REF/((1000/(P_TOT_WRF_REF/100))**0.286)
        TEMP_C_WRF_REF = TEMP_K_WRF_REF-273.15     


        PH_WRF_REF = Filobjekt_WRF_REF.variables['PH'][:]
        PHB_WRF_REF = Filobjekt_WRF_REF.variables['PHB'][:]
        MODELLNIVAHOJD_WRF_REF = (PH_WRF_REF+PHB_WRF_REF)/9.81
        MASSLEVELS_WRF_REF = 0.5*(MODELLNIVAHOJD_WRF_REF[0,:-1, 563, 352] + MODELLNIVAHOJD_WRF_REF[0, 1:, 563, 352])
        TERRAIN_HEIGHT_WRF_REF = Filobjekt_WRF_REF.variables['HGT'][0, 563, 352]
        MASSLEVELS_1D_MINUS_TER_WRF_REF = MASSLEVELS_WRF_REF - TERRAIN_HEIGHT_WRF_REF


        '''WRF_ASS'''

        THETA_PERT_WRF_ASS = Filobjekt_WRF_ASS.variables['T'][0, :, 563, 352]
        P_PERT_WRF_ASS = Filobjekt_WRF_ASS.variables['P'][0, :, 563, 352]
        P_BASE_WRF_ASS = Filobjekt_WRF_ASS.variables['PB'][0, :, 563, 352]   
        
        P_TOT_WRF_ASS= P_BASE_WRF_ASS + P_PERT_WRF_ASS
        THETA_WRF_ASS = THETA_PERT_WRF_ASS+300
        TEMP_K_WRF_ASS = THETA_WRF_ASS/((1000/(P_TOT_WRF_ASS/100))**0.286)
        TEMP_C_WRF_ASS = TEMP_K_WRF_ASS-273.15     


        PH_WRF_ASS = Filobjekt_WRF_ASS.variables['PH'][:]
        PHB_WRF_ASS = Filobjekt_WRF_ASS.variables['PHB'][:]
        MODELLNIVAHOJD_WRF_ASS= (PH_WRF_ASS+PHB_WRF_ASS)/9.81
        MASSLEVELS_WRF_ASS = 0.5*(MODELLNIVAHOJD_WRF_ASS[0,:-1, 563, 352] + MODELLNIVAHOJD_WRF_ASS[0, 1:, 563, 352])
        TERRAIN_HEIGHT_WRF_ASS = Filobjekt_WRF_ASS.variables['HGT'][0, 563, 352]
        MASSLEVELS_1D_MINUS_TER_WRF_ASS = MASSLEVELS_WRF_ASS - TERRAIN_HEIGHT_WRF_ASS



        GLW_TIME_SERIES_WRF_REF[place] =  LONGWAVE_DOWN_SURFACE_WRF_REF
        GLW_TIME_SERIES_WRF_ASS[place] =  LONGWAVE_DOWN_SURFACE_WRF_ASS

        T2_C_TIME_SERIES_WRF_REF[place] = T2_C_WRF_REF
        T2_C_TIME_SERIES_WRF_ASS[place] = T2_C_WRF_ASS



GLW_TIME_SERIES_WRF_REF[0] = GLW_TIME_SERIES_WRF_REF[1]
GLW_TIME_SERIES_WRF_ASS[0] = GLW_TIME_SERIES_WRF_ASS[1]

###############
# OBS NC-FILE #
#################################################################################################################
print('Nc-obs')


Filobjekt_OBS = N4.Dataset('/home/sm_marha/sodankyla_yopp.nc', mode='r')

LONGWAVE_DOWN_SURFACE_OBS = Filobjekt_OBS.variables['rlds'][:]

TIMESTEPS_NC_OBS = Filobjekt_OBS.variables['time'][:]
time = Filobjekt_OBS.variables['time']     #OBS, ingen array
dates = N4.num2date(TIMESTEPS_NC_OBS, time.units, time.calendar)

ALL_DATES = []

'''Gör om tiden till ett visst datumformat.'''

for date in dates:
    datum = date.strftime('%y%m%d %Hz')
    ALL_DATES.append(datum)
    
for place, date in enumerate(ALL_DATES):
    if date == '180227 00z':                                #ÄNDRA
        print(place)
        datum_START = place
    if date == '180228 00z':                                #ÄNDRA
        print(place)
        datum_END = place




#####################
# AUTOMATIC STATION #
#################################################################################################################
print('Automatic station')


df = pd.read_csv(r'/home/sm_marha/TEXTFILER/Automatic_Station_jan_feb_all_parameters.txt', encoding='latin-1', delimiter=',')

df = df.loc[:,['DATE TIME','T']]

df_numpy_array = df.values

TIDSVEKTOR_OBS = df_numpy_array[1:, 0]

TIDSVEKTOR_MINUTER_OBS = []

'''Letar upp platserna i arrayen för ändpunkterna på de datum jag är intresserad av.'''


for place, date in enumerate(TIDSVEKTOR_OBS):
    if date == '2018-02-18 00:00:00+00':                #ÄNDRA
        print(place)
        date_START = place
    if date == '2018-02-19 00:00:00+00':                #ÄNDRA
        print(place)
        date_END = place


'''Gör om formatet på datumen.'''

for date in TIDSVEKTOR_OBS[date_START:date_END+1]: 
    date = date[0:16]
    date = datetime.datetime.strptime(date,'%Y-%m-%d %H:%M').strftime('%d/%m %H')
    TIDSVEKTOR_MINUTER_OBS.append(date)

T2_TIDSSERIE_OBS = df_numpy_array[date_START:date_END+1, 1]

TIDSVEKTOR_TIMMAR_OBS = TIDSVEKTOR_MINUTER_OBS[0::6]

T2_GLES_OBS = (T2_TIDSSERIE_OBS[0::6])


'''Strängarna görs om till flytal.'''

for place, temp in enumerate(T2_GLES_OBS):
    T2_GLES_OBS[place] = float(temp)
    



############
# PLOTTING #
#################################################################################################################

fig1 = plt.figure(1, figsize=(14,8))



ax1 = plt.subplot(211)

ax1.plot(TIMESTEPS_NC_OBS[datum_START:datum_END], LONGWAVE_DOWN_SURFACE_OBS[datum_START:datum_END], label = 'OBS')
ax1.plot(TIMESTEPS_NC_OBS[datum_START:datum_END], GLW_TIME_SERIES_WRF_REF, label = 'WRF 3D REF')   #Använder TIMESTEPS från obs nc-file, då dessa har värden 7032-7074
ax1.plot(TIMESTEPS_NC_OBS[datum_START:datum_END], GLW_TIME_SERIES_WRF_ASS, label = 'WRF 3D QC QI') #De andra(TIDSVEKTOR_TIMMAR_OBS och TIME_SERIES) startar på 0 och går till 43)
                                                                                                   #TIMESTEPS_SCM_WRF_T2 kan lika gärna vara TIMESTEPS_SCM_WRF_LONG_WAVE, då det bara var ett sätt att namnge tidssteget för 2-dimensionella tidsserierna.
ax1.set_xticks(TIMESTEPS_NC_OBS[datum_START:datum_END+1:3])
#ax1.set_xticklabels(ALL_DATES[datum_START:datum_END+1:3], rotation=45, ha='right', fontsize=12)
ax1.set_xticklabels([])
ax1.set_xlim(datum_START, datum_END)
ax1.set_ylabel('Radiation  (W/m2)', fontsize=15) 
#ax1.set_xlabel('Time (Hours)')
ax1.set_title('Longwave radiation down', fontsize=18)
ax1.grid(linestyle="dotted")
ax1.legend(loc=4)




ax2 = plt.subplot(212)

ax2.plot(TIMESTEPS_NC_OBS[datum_START:datum_END], T2_GLES_OBS[0:24], label = 'OBS')
ax2.plot(TIMESTEPS_NC_OBS[datum_START:datum_END], T2_C_TIME_SERIES_WRF_REF, label = 'WRF 3D REF')
ax2.plot(TIMESTEPS_NC_OBS[datum_START:datum_END], T2_C_TIME_SERIES_WRF_ASS, label = 'WRF 3D QC QI')

        
ax2.set_xticks(TIMESTEPS_NC_OBS[datum_START:datum_END+1:3])
ax2.set_xticklabels(ALL_DATES[datum_START:datum_END+1:3], rotation=45, ha='right', fontsize=12)
#ax2.set_xticklabels([])
ax2.set_xlim(datum_START, datum_END)
ax2.set_xlabel('Time (Hours)', fontsize=15)
ax2.set_ylabel('Temperature  (C)', fontsize=15)
ax2.set_title('2 meter temperature', fontsize=18)
ax2.grid(linestyle="dotted")
ax2.legend(loc=4)

fig1.subplots_adjust(hspace=0.4)
fig1.tight_layout()


#plt.savefig('/home/sm_marha/FIGURES/FOKUS_17_27_FEB/27_FEB/COMPARE_T2_LW_OBS_WRF')

plt.show()
