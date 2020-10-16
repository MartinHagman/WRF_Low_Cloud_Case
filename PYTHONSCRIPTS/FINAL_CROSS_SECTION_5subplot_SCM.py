import numpy as N
import scipy.io.netcdf as S
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import pandas as pd
import os
import datetime
import matplotlib.colors as mplc
import netCDF4 as N4
from matplotlib.colors import LogNorm
from pylab import *
import matplotlib.gridspec as gridspec

'''Detta skript plottar tvärsnitt samt tidsserier av  variabler. Det fungerar för alla parametrar som har sina värden på massnivåerna. För turbulensparametrar
finns ett annat skript som plottar dessa'''


TIDSLANGD = 2161 #7561 #Sätt hur mycket data som ska plockas ut. Kan även göras vid plottning genom att sätta variabeln t.

rc('axes', linewidth=1.5)



fig1 = plt.figure(1,figsize=(11,8))


gs = gridspec.GridSpec(ncols=20, nrows=5)



AXS = [ fig1.add_subplot(gs[0, :]), fig1.add_subplot(gs[1, :]), fig1.add_subplot(gs[2, :]), fig1.add_subplot(gs[3, :]), fig1.add_subplot(gs[4, :])]



fig1.subplots_adjust(right=0.8)
fig1.subplots_adjust(hspace=0.1)

#cbar_ax = fig1.add_axes([0.87, 0.53, 0.01, 0.35])


rc('axes', linewidth=1.5)



TYPE_OF_RUN = 'Wv_TO_Wvmax'   #Kan vara 'REFERENCE', 'T_TO_Td', 'Wv_TO_Wvmax'
DATE = '18021800'
INTERPOLATION = '_SCM'   #Kan vara '_SCM', '_3D', ''   (Tomt!!)

Figure_letter_list = ['a', 'b', 'c', 'd', 'e']


'''Välj 46 eller 91 nivåer! Se till att ändra även var figurerna skall sparas utifrån det.'''


Fileobjects = ['/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfout_d01_2018-02-18_00:00:00_REFERENCE',
               '/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfout_d01_2018-02-18_00:00:00_QC_QI_REFERENCE',
               '/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfout_d01_2018-02-18_00:00:00',
               '/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfout_d01_2018-02-18_00:00:00_QC_QI_T_to_Td_1e-5']


'''Fileobjects = ['/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfout_d01_2018-02-26_00:00:00_REFERENCE',
               '/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfout_d01_2018-02-26_00:00:00_QC_QI_REFERENCE',
               '/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfout_d01_2018-02-26_00:00:00_QC_QI_Wv_to_Wvmax_1e-5',
               '/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfout_d01_2018-02-26_00:00:00_QC_QI_T_to_Td_1e-5']'''


i=0

for ax1, Fileobject in zip(AXS, Fileobjects):

    print('Nytt filobjekt')


    Filobjekt = N4.Dataset(Fileobject)


   

    TIMESTEPS = Filobjekt.variables['XTIME'][0:TIDSLANGD]
    print(N.shape(TIMESTEPS))
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
    T2M = Filobjekt.variables['T2'][0:TIDSLANGD, 1, 1]-273.15
    T2M[0]= T2M[1]
    QVAPOR_2D = Filobjekt.variables['QVAPOR'][0:TIDSLANGD, :, 1, 1]
    E_MATTNAD_2D = N.exp(N.log(611.2)+(17.62*TEMP_C_2D/(243.12+TEMP_C_2D)))
    QVAPOR_MATTNAD_2D = 0.622*E_MATTNAD_2D/(P_2D-E_MATTNAD_2D)
          


    RH_2D = QVAPOR_2D/QVAPOR_MATTNAD_2D
    RH_2D = N.transpose(RH_2D)

    TSK = Filobjekt.variables['TSK'][0:TIDSLANGD, 1, 1]-273.15
    EMISS = Filobjekt.variables['EMISS'][0:TIDSLANGD, 1, 1]

    CLDFRA = Filobjekt.variables['CLDFRA'][0:TIDSLANGD, 0, 1, 1]
    CLDFRA_2D = Filobjekt.variables['CLDFRA'][0:TIDSLANGD, :, 1, 1]


    QCLOUD = Filobjekt.variables['QCLOUD'][0:TIDSLANGD, 0, 1, 1]
    QCLOUD_2D = Filobjekt.variables['QCLOUD'][0:TIDSLANGD, :, 1, 1]


    QRAIN = Filobjekt.variables['QRAIN'][0:TIDSLANGD, 0, 1, 1]
    QRAIN_2D = Filobjekt.variables['QRAIN'][0:TIDSLANGD, :, 1, 1]


    QSNOW = Filobjekt.variables['QSNOW'][0:TIDSLANGD, 0, 1, 1]
    QSNOW_2D = Filobjekt.variables['QSNOW'][0:TIDSLANGD, :, 1, 1]


    QICE = Filobjekt.variables['QICE'][0:TIDSLANGD, 0, 1, 1]
    QICE_2D = Filobjekt.variables['QICE'][0:TIDSLANGD, :, 1, 1]


    QGRAUP = Filobjekt.variables['QGRAUP'][0:TIDSLANGD, 0, 1, 1]
    QGRAUP_2D = Filobjekt.variables['QGRAUP'][0:TIDSLANGD, :, 1, 1]

                                                                                                                 

    '''
    RTHRATEN = Filobjekt.variables['RTHRATEN'][0:TIDSLANGD, 0, 1, 1]
    RTHRATEN_2D = Filobjekt.variables['RTHRATEN'][0:TIDSLANGD, :, 1, 1]
    RTHRATEN_2D = N.transpose(RTHRATEN_2D)

    RTHBLTEN = Filobjekt.variables['RTHBLTEN'][0:TIDSLANGD, 0, 1, 1]
    RTHBLTEN_2D = Filobjekt.variables['RTHBLTEN'][0:TIDSLANGD, :, 1, 1]
    RTHBLTEN_2D = N.transpose(RTHBLTEN_2D)

    EXCH_H = Filobjekt.variables['EXCH_H'][0:TIDSLANGD, 0, 1, 1]
    EXCH_H_2D = Filobjekt.variables['EXCH_H'][0:TIDSLANGD, :, 1, 1]
    EXCH_H_2D = N.transpose(EXCH_H_2D)'''

    QTOTAL_2D = QVAPOR_2D + QCLOUD_2D + QICE_2D + QRAIN_2D + QSNOW_2D
    POISSON_C = 0.2854*(1-0.24*QVAPOR_2D)
    GAMMA = QVAPOR_2D*461/(1005.7+QTOTAL_2D*1.996e3)

    #THETA_L_2D = THETA_2D*((0.622*QVAPOR_2D)/(0.622+QTOTAL_2D))**POISSON_C * (QVAPOR_2D/QTOTAL_2D)**(-GAMMA) * N.exp(-2.257e6*(QCLOUD_2D+QRAIN_2D)/((1005.7+QTOTAL_2D*1.996e3)*TEMP_2D))  
    #THETA_L_2D = THETA_2D-(2.257e3/2.257e3)*(QCLOUD_2D+QRAIN_2D)

    THETA_L_2D = THETA_2D - ((2.26e6*THETA_2D)/(1004*TEMP_2D)* QCLOUD_2D)

        
    '''Väntar med att transponera dessa tills Theta_L har beräknats'''
    THETA_2D = N.transpose(THETA_2D)
    THETA_L_2D = N.transpose(THETA_L_2D)
    CLDFRA_2D = N.transpose(CLDFRA_2D)  #Transponerar för att det ska passa TIME och LEVELS nedan
    QCLOUD_2D = N.transpose(QCLOUD_2D)
    QRAIN_2D = N.transpose(QRAIN_2D)
    QSNOW_2D = N.transpose(QSNOW_2D)
    QICE_2D = N.transpose(QICE_2D)
    QGRAUP_2D = N.transpose(QGRAUP_2D)



                                                                                                                    

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
    GLW[0] = GLW[1]

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
    MODELLNIVAHOJD_TER = MODELLNIVAHOJD[0:TIDSLANGD, :] -MODELLNIVAHOJD[0, 0]


    MASSLEVELS = 0.5*(MODELLNIVAHOJD_TER[0:TIDSLANGD, :-1] + MODELLNIVAHOJD_TER[0:TIDSLANGD, 1:])  #Olika masslevels för varje tidssteg, om än mycket små skillnader!



    #################################### 
    #            PLOTTNING             # 
    ####################################################################
        
    
    


    VARIABLE_DICT = {'QCLOUD_2D': QCLOUD_2D}

    list(VARIABLE_DICT.keys())[0]

    '''Det vi nu vill göra är att få massnivåer som varierar i tiden och inte är statiskt lika som vid tidpunkt 1'''

    '''Använder massnivåerna vid tidpunkt 0 för att skapa matriser TIME och LEVELS med rätt dimension'''
    TIME, LEVELS = N.meshgrid(TIMESTEPS, MASSLEVELS[0, :])   #Behövs för contour Theta på massnivåerna(fulla nivåerna) samt contourf
    print(N.shape(TIMESTEPS))
    '''Nollställer LEVELS med med fortsatt samma dimension'''
    LEVELS = N.zeros_like(LEVELS)


    '''Skriver in massnivåer för de 7561 olika tidpunkterna. Antal rows är 90 stycken. 0-89 om man indexerar enligt Python'''
    for row in range(N.shape(LEVELS)[0]):
        LEVELS[row, :] = (N.transpose(MASSLEVELS))[row, :]
    
    

    '''Denna tas ggr 1000 för alla utom RH och CLDFRA för att göra om till g/kg i stället för kg/kg'''

    if list(VARIABLE_DICT.keys())[0] != 'CLDFRA_2D' and list(VARIABLE_DICT.keys())[0] != 'RH_2D':
        VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]] = VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]]*1000


    #Speciell plottning av RH

    if list(VARIABLE_DICT.keys())[0] == 'RH_2D':
    
        norm = mplc.Normalize(0, 1.1)
        A = 1.1/44
        myplot = ax1.contourf(TIME, LEVELS, VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]], N.arange(0, 1.1, A), norm = norm, extend = 'both', cmap ='Greens')
        contour_levels = N.arange(0, 1, 0.05)
        cp = plt.contour(TIME, LEVELS, VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]], contour_levels, colors='black', linestyles='dotted', linewidth='0.2' )
        #plt.clabel(cp, inline=True, fontsize=8)

        '''Plottar masslevels längs y-axeln'''   
        X = N.zeros_like(MASSLEVELS[0, :])+3     #Lägger till 3 för att plottarna inte ska hamna rakt på y-axeln.
        ax1.scatter(X, MASSLEVELS[0, :], marker='_', color='lime', s=8)
        plt.ylim(0,1000)
        plt.xlim(0, N.max(TIME))
    
    else:

        norm = mplc.Normalize(0, 5e-1)
        A = 5e-1/20
        myplot = ax1.contourf(TIME, LEVELS, VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]], N.arange(0, 5e-1, A), norm = norm, extend = 'both', cmap ='binary')
        #myplot = ax1.contourf(TIME, LEVELS, VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]], extend = 'both', cmap ='binary')

        '''Välj interpolering i fälten med contourf, eller plotta med scatter för att se exakt på vilka modellnivåer det finns moln'''
        #myplot = ax1.contourf(TIME, LEVELS, VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]], extend = 'both', cmap ='binary')
        #myplot = ax1.scatter(TIME, LEVELS, c=VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]], cmap='binary')
        #myplot = ax1.pcolormesh(TIME, LEVELS, VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]], vmin=-0.5, vmax=0.5, cmap = 'seismic')
    
        '''Plottar contours av Theta tillsammans med contourf'''
    
        contour_levels = N.arange(0, 330, 1.5)
        cp = ax1.contour(TIME, LEVELS, THETA_2D, contour_levels, linewidths= 0.5, colors='red', linestyles='solid')
        #cp = plt.contour(TIME, LEVELS, THETA_L_2D, contour_levels, colors='yellow', linestyles='dotted', linewidth='0.2' )
        plt.clabel(cp, contour_levels[::2], inline=True, fmt='%1.0f', fontsize=10)
    
        
        '''Plottar masslevels längs y-axeln''' 
        X = N.zeros_like(MASSLEVELS)[0, :]+3
        #ax1.scatter(X, MASSLEVELS[0, :], marker= '>', color='black', s=2)
        plt.xlim(0-100, N.max(TIME))
    
        ax1.text(700, 1600, Figure_letter_list[i], fontsize = 18, fontweight='bold')
        
        cbar_ax = fig1.add_axes([0.87, 0.30, 0.01, 0.55])       #Gemensam bar
        cbar = fig1.colorbar(myplot, cax=cbar_ax, format = "%8.1e")
    
        #cbar = fig1.colorbar(myplot, ax=ax1)           # För varje bar. Här måste  fig1.add_subplot(gs[4, :]) sättas till fig1.add_subplot(gs[4, :16])
        cbar.ax.tick_params(labelsize = 14)

          

        #ax1.set_title(str(i+1), fontsize = 15)
    '''

    VARIABLE = 'QCLOUD' #omdefinierar VARIABLE för titeln



    #Speciellt format på siffror(exponent) för colorbar för alla utom CLDFRA och RH

    if VARIABLE != 'CLDFRA' and VARIABLE != 'RH':
    
        cbar_ax = fig1.add_axes([0.92, 0.30, 0.01, 0.58])
        cbar = fig1.colorbar(myplot,cax=cbar_ax, format = "%8.1e")
    
    else:
    
    
        cbar_ax = fig1.add_axes([0.92, 0.30, 0.01, 0.58])
        cbar = fig1.colorbar(myplot, cax=cbar_ax)#, format = "%8.1e")


    #RH är en beräknad variabel och finns inget description-attribut för

    if VARIABLE == 'RH':

        ax1.set_title(VARIABLE + ' ' + TYPE_OF_RUN +' '+ DATE)

    elif VARIABLE == 'CLDFRA':

        ax1.set_title(Filobjekt.variables[VARIABLE].description + ' ' + TYPE_OF_RUN + ' ' + DATE, fontsize=18)
    

    elif VARIABLE == 'RTHRATEN' or VARIABLE == 'RTHBLTEN':
    
        ax1.set_title(Filobjekt.variables[VARIABLE].description + ' ' + TYPE_OF_RUN + ' ' + DATE + ' (K/s)', fontsize=18)
    

    elif VARIABLE == 'EXCH_H' or VARIABLE == 'EXCH_M':
    
        ax1.set_title(Filobjekt.variables[VARIABLE].description + ' ' + TYPE_OF_RUN + ' ' + DATE + ' (m2/s)', fontsize=18)
    
    
    else:

        #ax1.set_title(Filobjekt.variables[VARIABLE].description + ' ' + TYPE_OF_RUN + ' ' + DATE + ' (g/kg)', fontsize=18)
        ax1.set_title('2018-02-18 Reference')'''



    ax1.set_xticks(N.arange(0, len(ALLA_DATUM)/3, 180 ))  #X-axeln har enheten minuter, medan alla datum är var 20 sekund. Utgå därför från minuter,
                                                      # då ett "hopp" längs x-axeln är en minut, när den ska splittas upp i 3h-intervall för x-ticksen
    ax1.set_yticks([500, 1500, 2500])
    ax1.set_yticklabels([500, 1500, 2500], fontsize=16)
    

    ax1.tick_params(labelsize=10)

    ax1.set_xticklabels([])

        
    if i == 2:
        ax1.set_ylabel('              Height $(m)$', fontsize =18)
    ax1.tick_params(which='major', direction='in')

    #ax1.set_xlim(0, 360)

    #cbar.ax.set_title('                 ($gkg^-^1)')    #En cbar för varje
    cbar_ax.set_title(' ($gkg^{-1})$', fontsize = 14)            #Gemensam

    ax1.set_ylim(0, 2000)

    ax1.grid(linestyle="dotted")
    ax1.set_xlim(0-7, N.max(TIME))






####### 2nd AXIS #######

    '''Lagrar temperaturkurvor från olika körningar så att jag kan ha dem i samma legend'''

    ax2 = AXS[4]

    if i == 0:
        line = 'dotted'
        lab = 'T-Skin a'
    elif i == 1:
        line = 'dashed'
        lab = 'T-Skin b'
    elif i == 2:
        line = 'dashdot'
        lab = 'T-Skin c'
    else:
        line = 'solid'
        lab = 'T-Skin d'


    ax2.plot(TIMESTEPS[0:TIDSLANGD], TSK[0:TIDSLANGD], linewidth = 2,linestyle = line, label = lab)
    #ax2.plot(TIMESTEPS[0:TIDSLANGD], TEMP_C[0:TIDSLANGD], linewidth = 1, label = 'T 1st level')
    #ax2.plot(TIMESTEPS[0:TIDSLANGD], T2M[0:TIDSLANGD], linewidth = 1,linestyle = line, label = 'T2m')
    ax2.legend(fontsize = 'x-large', loc='center right', bbox_to_anchor=(1.15, 0.5))
    ax2.set_xticks(N.arange(0, len(ALLA_DATUM)/3, 180 ))
    ax2.set_xticklabels(ALLA_DATUM_3h, rotation=45, fontsize=16)   #, ha='right'
    ax2.set_xticklabels([])
    ax2.set_xlim(0-7, N.max(TIME))
    ax2.set_ylabel('T-skin $(^oC)$', fontsize=16) 
    #ax2.set_xlabel('Time (Hours)')
    #ax2.set_title('Longwave radiation down', fontsize=12)
    ax2.grid(linestyle="dotted")

    i+=1

ax2.text(700, -7.8, Figure_letter_list[4], fontsize = 18, fontweight='bold')

ax2.tick_params(which='major', direction='in')

ax2.set_xticks(N.arange(0, len(ALLA_DATUM)/3, 180 ))
ax2.set_xticklabels(ALLA_DATUM_3h, rotation=45, fontsize=15)   #, ha='right'
#ax2.set_xlabel('$Time$', fontsize=10)

'''
####### 3rd AXIS #######

ax3 = plt.subplot(414)


ax2.plot(TIMESTEPS[0:TIDSLANGD], GLW[0:TIDSLANGD], linewidth = 1, label = 'LW down')
'''



#################################### 
#            OBSERVATIONER         # 
####################################################################



###############
# OBS NC-FILE #
#################################################################################################################
print('Nc-obs')


Filobjekt_OBS = N4.Dataset('/home/sm_marha/sodankyla_yopp.nc', mode='r')

LONGWAVE_DOWN_SURFACE_OBS = Filobjekt_OBS.variables['rlds'][:]
LONGWAVE_UP_SURFACE_OBS = Filobjekt_OBS.variables['rlus'][:]

TSK_OBS_K = (LONGWAVE_UP_SURFACE_OBS/5.67e-8)**0.25

TSK_OBS_C = TSK_OBS_K - 273.15

TIMESTEPS_NC_OBS = Filobjekt_OBS.variables['time'][:]
time = Filobjekt_OBS.variables['time']     #OBS, ingen array
dates = N4.num2date(TIMESTEPS_NC_OBS, time.units, time.calendar)

ALL_DATES = []

'''Gör om tiden till ett visst datumformat.'''

for date in dates:
    datum = date.strftime('%y%m%d %Hz')
    ALL_DATES.append(datum)
    
for place, date in enumerate(ALL_DATES):
    if date == '180218 00z':                                #ÄNDRA
        print(place)
        datum_START = place
    if date == '180218 12z':                                #ÄNDRA
        print(place)
        datum_END = place




ax2.plot(TIMESTEPS[0:TIDSLANGD][::180],TSK_OBS_C[datum_START:datum_END+1],  color = 'k', linewidth = 2.5, label = 'OBS')
ax2.legend( loc='center right', bbox_to_anchor=(1.28, 0.5), fontsize = 'x-large')

fig1.savefig('/home/sm_marha/FIGURES/FINAL/20180218/20180218_SCM_SUBPLOTS_5_Corrected', dpi = 300)

plt.show()

