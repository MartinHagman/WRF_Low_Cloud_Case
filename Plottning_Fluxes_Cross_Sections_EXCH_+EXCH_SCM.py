import numpy as N
import scipy.io.netcdf as S
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import os
import datetime
import netCDF4 as N4
import matplotlib.colors as mplc
from matplotlib.colors import LogNorm



fig1, AXS = plt.subplots(2, 1, figsize=(10,6))
fig1.subplots_adjust(bottom=0.2)
fig1.subplots_adjust(right=0.90)

TIDSLANGD = 7561 #Sätt hur mycket data som ska plockas ut. Kan även göras vid plottning genom att sätta variabeln t.



DIRECTORIES = ['/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfout_d01_2018-02-27_00:00:00_025T_LIQ_AND_VARIATION_YSU_TD1',
                '/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfout_d01_2018-02-27_00:00:00_025T_LIQ_AND_VARIATION_MYNN']

i = 0

for DIRECTORY, ax in zip(DIRECTORIES, AXS):
    Filobjekt = N4.Dataset(DIRECTORY)






    TIMESTEPS = Filobjekt.variables['XTIME'][0:TIDSLANGD]
    time = Filobjekt.variables['XTIME']
    dates = N4.num2date(TIMESTEPS, time.units)


    ALLA_DATUM=[]

    for date in dates:
        datum = date.strftime('%d/%-m %Hz')
        ALLA_DATUM.append(datum)

    ALLA_DATUM_3h = ALLA_DATUM[0:len(ALLA_DATUM):540]





    PERT_THETA_2D = Filobjekt.variables['T'][0:TIDSLANGD, :, 1, 1]
    P_PERT_2D = Filobjekt.variables['P'][0:TIDSLANGD, :, 1, 1]
    P_BASE_2D = Filobjekt.variables['PB'][0:TIDSLANGD, :, 1, 1]
    P_2D= P_BASE_2D+P_PERT_2D

    THETA_2D = PERT_THETA_2D+300

    TEMP_2D = THETA_2D/((1000/(P_2D/100))**0.286)
    TEMP_C_2D = TEMP_2D-273.15
    TEMP_C = TEMP_C_2D[:, 0]
    QVAPOR_2D = Filobjekt.variables['QVAPOR'][0:TIDSLANGD, :, 1, 1]
    E_MATTNAD_2D = N.exp(N.log(611.2)+(17.62*TEMP_C_2D/(243.12+TEMP_C_2D)))
    QVAPOR_MATTNAD_2D = 0.622*E_MATTNAD_2D/(P_2D-E_MATTNAD_2D)

    THETA_2D = N.transpose(THETA_2D)

    RH_2D = QVAPOR_2D/QVAPOR_MATTNAD_2D
    RH_2D = N.transpose(RH_2D)

    TSK = Filobjekt.variables['TSK'][0:TIDSLANGD, 1, 1]-273.15
    EMISS = Filobjekt.variables['EMISS'][0:TIDSLANGD, 1, 1]

    CLDFRA = Filobjekt.variables['CLDFRA'][0:TIDSLANGD, 0, 1, 1]
    CLDFRA_2D = Filobjekt.variables['CLDFRA'][0:TIDSLANGD, :, 1, 1]
    CLDFRA_2D = N.transpose(CLDFRA_2D)  #Transponerar för att det ska passa TIME och LEVELS nedan

    QCLOUD = Filobjekt.variables['QCLOUD'][0:TIDSLANGD, 0, 1, 1]
    QCLOUD_2D = Filobjekt.variables['QCLOUD'][0:TIDSLANGD, :, 1, 1]
    QCLOUD_2D = N.transpose(QCLOUD_2D)

    QRAIN = Filobjekt.variables['QRAIN'][0:TIDSLANGD, 0, 1, 1]
    QRAIN_2D = Filobjekt.variables['QRAIN'][0:TIDSLANGD, :, 1, 1]
    QRAIN_2D = N.transpose(QRAIN_2D)


    QSNOW = Filobjekt.variables['QSNOW'][0:TIDSLANGD, 0, 1, 1]
    QSNOW_2D = Filobjekt.variables['QSNOW'][0:TIDSLANGD, :, 1, 1]
    QSNOW_2D = N.transpose(QSNOW_2D)


    QICE = Filobjekt.variables['QICE'][0:TIDSLANGD, 0, 1, 1]
    QICE_2D = Filobjekt.variables['QICE'][0:TIDSLANGD, :, 1, 1]
    QICE_2D = N.transpose(QICE_2D)

    QGRAUP = Filobjekt.variables['QGRAUP'][0:TIDSLANGD, 0, 1, 1]
    QGRAUP_2D = Filobjekt.variables['QGRAUP'][0:TIDSLANGD, :, 1, 1]
    QGRAUP_2D = N.transpose(QGRAUP_2D)

    RTHRATEN = Filobjekt.variables['RTHRATEN'][0:TIDSLANGD, 0, 1, 1]
    RTHRATEN_2D = Filobjekt.variables['RTHRATEN'][0:TIDSLANGD, :, 1, 1]
    RTHRATEN_2D = N.transpose(RTHRATEN_2D)

    RTHBLTEN = Filobjekt.variables['RTHBLTEN'][0:TIDSLANGD, 0, 1, 1]
    RTHBLTEN_2D = Filobjekt.variables['RTHBLTEN'][0:TIDSLANGD, :, 1, 1]
    RTHBLTEN_2D = N.transpose(RTHBLTEN_2D)

    EXCH_H = Filobjekt.variables['EXCH_H'][0:TIDSLANGD, 0, 1, 1]
    EXCH_H_2D = Filobjekt.variables['EXCH_H'][0:TIDSLANGD, :, 1, 1]
    EXCH_H_2D = N.transpose(EXCH_H_2D)

    EXCH_M = Filobjekt.variables['EXCH_M'][0:TIDSLANGD, 0, 1, 1]
    EXCH_M_2D = Filobjekt.variables['EXCH_M'][0:TIDSLANGD, :, 1, 1]
    EXCH_M_2D = N.transpose(EXCH_M_2D)

    RAINNC = Filobjekt.variables['RAINNC'][0:TIDSLANGD, 1, 1]
    SNOWNC = Filobjekt.variables['SNOWNC'][0:TIDSLANGD, 1, 1]



    GLW = Filobjekt.variables['GLW'][0:TIDSLANGD, 1, 1]
    SWDOWN = Filobjekt.variables['SWDOWN'][0:TIDSLANGD, 1, 1]
    GRDFLX = Filobjekt.variables['GRDFLX'][0:TIDSLANGD, 1, 1]
    HFX = Filobjekt.variables['HFX'][0:TIDSLANGD, 1, 1]
    LH = Filobjekt.variables['LH'][0:TIDSLANGD, 1, 1]
    EMISS = Filobjekt.variables['EMISS'][0:TIDSLANGD, 1, 1]
    ALBEDO = Filobjekt.variables['ALBEDO'][0:TIDSLANGD, 1, 1]
    NOAHRES = Filobjekt.variables['NOAHRES'][0:TIDSLANGD, 1, 1]
    FLX1 = Filobjekt.variables['FLX1'][0:TIDSLANGD, 1, 1]
    FLX2 = Filobjekt.variables['FLX2'][0:TIDSLANGD, 1, 1]
    FLX3 = Filobjekt.variables['FLX3'][0:TIDSLANGD, 1, 1]
            
    SWUP = ALBEDO*SWDOWN
    GUPLW = EMISS*5.67e-8*(TSK+273.15)**4
    GUPLW[0] = 0

    LONGW_RESIDUAL = EMISS*GLW - GUPLW
    HFX = -HFX
    SWUP = -SWUP

    RESIDUAL = EMISS*GLW + SWDOWN + SWUP - GUPLW + HFX - LH + GRDFLX - FLX1 - FLX2 - FLX3



    ZNW = Filobjekt.variables['ZNW'][0:TIDSLANGD, :]

    P = Filobjekt.variables['P'][0:TIDSLANGD, :, 1, 1]

    PB = Filobjekt.variables['PB'][0:TIDSLANGD, :, 1, 1]




    PH = Filobjekt.variables['PH'][:, :, 1, 1]
    PHB = Filobjekt.variables['PHB'][:, :, 1, 1]
    MODELLNIVAHOJD = (PH+PHB)/9.81
    MODELLNIVAHOJD_TER = MODELLNIVAHOJD[0:TIDSLANGD, :] -MODELLNIVAHOJD[0, 0]    #Olika modellnivåhöjder för varje tidssteg, men mycket små skillnader.


    MASSLEVELS = 0.5*(MODELLNIVAHOJD_TER[0:TIDSLANGD:, :-1] + MODELLNIVAHOJD_TER[0:TIDSLANGD, 1:])




    #################################### 
    #            PLOTTNING             # 
    ###################################


    if i==0:
        VARIABLE_DICT = {'EXCH_H_2D': EXCH_H_2D}
    elif i==1:
        VARIABLE_DICT = {'EXCH_M_2D': EXCH_M_2D}

    list(VARIABLE_DICT.keys())[0]




    '''MASSVARIABLER'''

    '''Använder massnivåerna vid tidpunkt 0 för att skapa matriser TIME och LEVELS med rätt dimension'''
    TIME, LEVELS = N.meshgrid(TIMESTEPS, MASSLEVELS[0, :])   #Behövs för contour Theta på massnivåerna(fulla nivåerna)

    '''Nollställer LEVELS med med fortsatt samma dimension'''
    LEVELS = N.zeros_like(LEVELS)

    '''Skriver in massnivåer för de 7561 olika tidpunkterna. Antal rows är 90 stycken. 0-89 om man indexerar enligt Python'''
    for row in range(N.shape(LEVELS)[0]):
        LEVELS[row, :] = (N.transpose(MASSLEVELS))[row, :]



    '''TURBULENSVARIABLER'''
        
    '''Använder massnivåerna vid tidpunkt 0 för att skapa matriser TIME och LEVELS med rätt dimension'''
    TIME_FULL, LEVELS_FULL = N.meshgrid(TIMESTEPS, MODELLNIVAHOJD_TER[0, :])   #Behövs för K-värdena på halva nivåerna. Där finns turbulensparametrarna.

    '''Nollställer LEVELS med med fortsatt samma dimension'''
    LEVELS_FULL = N.zeros_like(LEVELS_FULL)

    '''Skriver in massnivåer för de 7561 olika tidpunkterna. Antal rows är 90 stycken. 0-89 om man indexerar enligt Python'''
    for row in range(N.shape(LEVELS_FULL)[0]):
        LEVELS_FULL[row, :] = (N.transpose(MODELLNIVAHOJD_TER))[row, :]
        
        


    '''Denna tas ggr 1000 för alla utom RH och CLDFRA för att göra om till g/kg i stället för kg/kg'''

    '''if list(VARIABLE_DICT.keys())[0] != 'CLDFRA_2D' and list(VARIABLE_DICT.keys())[0] != 'RH_2D':
        VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]] = VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]]*1000'''


    #print(N.max(VARIABLE_DICT['EXCH_H_2D']))

    #norm = mplc.Normalize(0, 3e-1)
    norm = mplc.Normalize(0, 70)
    
    #A = 3e-1/20
    A = 70/20
    #myplot = ax.contourf(TIME_FULL, LEVELS_FULL,VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]], N.arange(0, 70, A), norm=norm, extend = 'both', cmap ='binary')
    

    '''Välj interpolering i fälten med contourf, eller plotta med scatter för att se exakt på vilka modellnivåer det finns moln'''
    #myplot = ax.contourf(TIME_FULL, LEVELS_FULL, VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]], [10, 20, 30, 40, 50, 100, 1000, 10000, 100000, 1000000], cmap ='binary', norm = LogNorm())
    myplot = ax.contourf(TIME_FULL, LEVELS_FULL,VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]],[0.001, 0.01, 0.1, 1, 10, 100], cmap ='binary', norm = LogNorm())
    #myplot = ax.scatter(TIME_FULL, LEVELS_FULL, c=VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]], cmap='binary')
    #myplot = ax.pcolormesh(TIME_FULL, LEVELS_FULL, VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]], vmin=-0.5, vmax=0.5, cmap = 'seismic')
    

    '''Plottar contours av Theta tillsammans med contourf'''
    
    contour_levels = N.arange(250, 290, 0.5)
    cp = ax.contour(TIME, LEVELS, THETA_2D, contour_levels, colors='red', linestyles='dotted', linewidth='0.2' )
    ax.clabel(cp,contour_levels[::2], inline=True, fmt='%1.0f', fontsize=8)
    
        
    '''Plottar masslevels vid t=0 längs y-axeln''' 
    X = N.zeros_like(MASSLEVELS[0, :])+3
    ax.scatter(X, MASSLEVELS[0, :], marker= '>', color='black', s=20)
    ax.set_xlim(0, N.max(TIME))
    ax.set_ylim(1000,1500)     #Finns ylim längre ner i skriptet också
    
    
    
        
    cbar_ax = fig1.add_axes([0.92, 0.20, 0.01, 0.68])
    #cbar = fig1.colorbar(myplot, ticks=, cax=cbar_ax)#, format = "%8.1e")
    cbar = fig1.colorbar(myplot, cax=cbar_ax)#, ticks =[0.01, 0.1, 1, 10, 100])#, format = "%8.1e")


        


    ax.set_xticks(N.arange(0, len(ALLA_DATUM)/3, 180 ))  #X-axeln har enheten minuter, medan alla datum är var 20 sekund. Utgå därför från minuter,
                                                          # då ett "hopp" längs x-axeln är en minut, när den ska splittas upp i 3h-intervall för x-ticksen

    #ax.tick_params(labelsize=15)

    if i==0:
        ax.set_xticklabels([])
        ax.text(7, 1750,'a', fontsize = 18, fontweight='bold')
        
    elif i==1:
        ax.set_xticklabels(ALLA_DATUM_3h, rotation=45, ha='right', fontsize=13)
        #ax.set_xlabel('Time', fontsize=15)
        ax.text(7, 1750,'b', fontsize = 18, fontweight='bold')
        
    ax.grid(linestyle="dotted")
    
    ax.set_ylabel('Height (m)', fontsize =15)

    #ax.set_xlim(0, 360) 

    ax.set_ylim(0, 8000)

    

    i+=1

#fig1.savefig('/home/sm_marha/FIGURES/FINAL/20180227/EXCH_+_EXCH_YSU_TD1_025')
plt.show()






