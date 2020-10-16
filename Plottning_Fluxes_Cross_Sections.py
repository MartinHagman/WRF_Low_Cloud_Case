import numpy as N
import scipy.io.netcdf as S
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import os
import datetime
import matplotlib.colors as mplc
import netCDF4 as N4



TYPE_OF_RUN = 'T_TO_Td'   #Kan vara 'REFERENCE', 'T_TO_Td', 'Wv_TO_Wvmax'
DATE = '18011800'
INTERPOLATION = 'SCM'   #Kan vara '_SCM', '_3D', ''   (Tomt!!)


'''Välj 46 eller 91 nivåer! Se till att ändra även var figurerna skall sparas utifrån det.'''


Filobjekt = N4.Dataset('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfout_d01_2018-01-18_00:00:00')
#Filobjekt = N4.Dataset('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY/test/em_scm_xy/wrfout_d01_2018-02-18_00:00:00')

TIDSLANGD = 7561 #Sätt hur mycket data som ska plockas ut. Kan även göras vid plottning genom att sätta variabeln t.

TIMESTEPS = Filobjekt.variables['XTIME'][0:TIDSLANGD]
time = Filobjekt.variables['XTIME']     #OBS, ingen array
dates = N4.num2date(TIMESTEPS, time.units)
ALLA_DATUM = N.empty(0)



'''Gör om tiden till ett visst datumformat.'''

for date in dates:
    datum = date.strftime('%d/%m %Hz')
    ALLA_DATUM = N.append(ALLA_DATUM, datum)

    


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
RH_2D = QVAPOR_2D/QVAPOR_MATTNAD_2D
RH_2D = N.transpose(RH_2D)

TSK = Filobjekt.variables['TSK'][0:TIDSLANGD, 1, 1]-273.15
EMISS = Filobjekt.variables['EMISS'][0:TIDSLANGD, 1, 1]

TEST = Filobjekt.variables['CLDFRA']
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


RAINNC = Filobjekt.variables['RAINNC'][0:TIDSLANGD, 1, 1]
SNOWNC = Filobjekt.variables['SNOWNC'][0:TIDSLANGD, 1, 1]


'''RTHRATEN = Filobjekt.variables['RTHRATEN'][0:TIDSLANGD, 0, 1, 1]
RTHRATEN_2D = Filobjekt.variables['RTHRATEN'][0:TIDSLANGD, :, 1, 1]
RTHRATEN_2D = N.transpose(RTHRATEN_2D)


RTHBLTEN = Filobjekt.variables['RTHBLTEN'][0:TIDSLANGD, 0, 1, 1]
RTHBLTEN_2D = Filobjekt.variables['RTHBLTEN'][0:TIDSLANGD, :, 1, 1]
RTHBLTEN_2D = N.transpose(RTHBLTEN_2D)


EXCH_H = Filobjekt.variables['EXCH_H'][0:TIDSLANGD, 0, 1, 1]
EXCH_H_2D = Filobjekt.variables['EXCH_H'][0:TIDSLANGD, :, 1, 1]
EXCH_H_2D = N.transpose(EXCH_H_2D)'''




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
MODELLNIVAHOJD_TER = MODELLNIVAHOJD[0:TIDSLANGD, :] -MODELLNIVAHOJD[0,0]


MASSLEVELS = 0.5*(MODELLNIVAHOJD_TER[0:TIDSLANGD,:-1] + MODELLNIVAHOJD_TER[0:TIDSLANGD, 1:])







#################################### 
#            PLOTTNING             # 
#################################### 


##########
#FIGURE 1#
##########
fig1 = plt.figure(1)

t=7561


ax1 = plt.subplot(211)

ax1.plot(TIMESTEPS[0:t], TSK[0:t], linewidth = 1, label = 'T-SKIN')
ax1.plot(TIMESTEPS[0:t], TEMP_C[0:t], linewidth = 1, label = 'T-1st LEVEL')

ax1.legend()
ax1.grid(linestyle="dotted")

plt.ylabel('Temperature (C)')
plt.title(TYPE_OF_RUN + ' '+ DATE)



ax2 = plt.subplot(212)

ax2.plot(TIMESTEPS, RAINNC, linewidth = 1, label = 'RAINNC')
ax2.plot(TIMESTEPS, SNOWNC, linewidth = 1, label = 'SNOWNC')

ax2.legend()
ax2.grid(linestyle="dotted")

plt.xlabel('Time (Minutes)')
plt.ylabel('Precipitation (mm)')



'''Spar på olika ställen beroende på om det är 46 eller 91 nivåer'''
#fig1.savefig('/home/sm_marha/FIGURES/SINGLE_COLUMN/'+DATE +'Z/' + TYPE_OF_RUN + '/' + DATE +'z_TSKIN_Precip' + INTERPOLATION)
#fig1.savefig('/home/sm_marha/FIGURES/SINGLE_COLUMN/'+DATE +'Z/91_LEVELS/' + TYPE_OF_RUN + '/' + DATE +'z_TSKIN_Precip' + INTERPOLATION+ '_NOT_Xu')


##########
#FIGURE 2#
##########

fig2 = plt.figure(2, figsize=(16,6)) 

ax3 = plt.subplot(111)

fig2.subplots_adjust(bottom=0.2)


VARIABLE_DICT = {'QCLOUD': QCLOUD_2D}

list(VARIABLE_DICT.keys())[0]


'''Det vi nu vill göra är att få massnivåer som varierar i tiden och inte är statiskt lika som vid tidpunkt 1'''

'''Använder massnivåerna vid tidpunkt 0 för att skapa matriser TIME och LEVELS med rätt dimension'''
TIME, LEVELS = N.meshgrid(TIMESTEPS, MASSLEVELS[0, :])   #Behövs för contour Theta på massnivåerna(fulla nivåerna) samt contourf

'''Nollställer LEVELS med med fortsatt samma dimension'''
LEVELS = N.zeros_like(LEVELS)


'''Skriver in massnivåer för de 7561!!!! olika tidpunkterna. Antal rows är 90 stycken. 0-89 om man indexerar enligt Python'''
for row in range(N.shape(LEVELS)[0]):
    LEVELS[row, :] = (N.transpose(MASSLEVELS))[row, :]




#Denna tas ggr 1000 för alla utom RH och CLDFRA för att göra om till g/kg i stället för kg/kg#

if list(VARIABLE_DICT.keys())[0] != 'CLDFRA' and list(VARIABLE_DICT.keys())[0] != 'RH':
    VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]] = VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]]*1000


#Speciell plottning av RH

if list(VARIABLE_DICT.keys())[0] == 'RH':
    
    norm = mplc.Normalize(0, 1.1)
    A = 1.1/44
    myplot = ax3.contourf(TIME, LEVELS, VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]], N.arange(0, 1.1, A), norm = norm, extend = 'both', cmap ='Greens')
    contour_levels = N.arange(0, 1, 0.05)
    cp = plt.contour(TIME, LEVELS, VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]], contour_levels, colors='black', linestyles='dotted', linewidth='0.2' )
    plt.clabel(cp, inline=True, fontsize=8)
    plt.ylim(0,2500)

    '''Plottar masslevels längs y-axeln'''   
    X = N.zeros_like(MASSLEVELS[0, :])+3     #Lägger till 3 för att plottarna inte ska hamna rakt på y-axeln.
    ax3.scatter(X, MASSLEVELS[0, :], marker='_', color='lime', s=8)
    plt.xlim(0, N.max(TIME))
    
else:

    '''norm = mplc.Normalize(0, 3e-1)
    A = 3e-1/20
    myplot = ax3.contourf(TIME, LEVELS, VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]], N.arange(0, 3e-1, A), norm = norm, extend = 'both', cmap ='binary')'''

    '''Välj interpolering i fälten med contourf, eller plotta med scatter för att se exakt på vilka modellnivåer det finns moln'''
    myplot = ax3.contourf(TIME, LEVELS, VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]], extend = 'both', cmap ='binary')
    #myplot = ax3.scatter(TIME, LEVELS, c=VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]], cmap='binary')
    #myplot = ax3.pcolormesh(TIME, LEVELS, VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]], vmin=-0.5, vmax=0.5, cmap = 'seismic')


    '''Plottar masslevels längs y-axeln''' 
    X = N.zeros_like(MASSLEVELS[0, :])+3
    ax3.scatter(X, MASSLEVELS[0, :], marker= '>', color='black', s=8)
    plt.xlim(0, N.max(TIME))
    

VARIABLE = 'QCLOUD' #omdefinierar VARIABLE för titeln



#Speciellt format på siffror(exponent) för colorbar för alla utom CLDFRA och RH

if VARIABLE != 'CLDFRA' and VARIABLE != 'RH':
    
    cbar = fig2.colorbar(myplot,format = "%8.1e")
    
else:
    
    cbar = fig2.colorbar(myplot)
    


#RH är en beräknad variabel och finns inget description-attribut för

if VARIABLE == 'RH':

    plt.title(VARIABLE + ' ' + TYPE_OF_RUN +' '+ DATE)
    

elif VARIABLE == 'CLDFRA':

    plt.title(Filobjekt.variables[VARIABLE].description + ' ' + TYPE_OF_RUN + ' ' + DATE)
    

elif VARIABLE == 'RTHRATEN' or VARIABLE == 'RTHBLTEN':
    
    plt.title(Filobjekt.variables[VARIABLE].description + ' ' + TYPE_OF_RUN + ' ' + DATE + ' (K/s)')
    

elif VARIABLE == 'EXCH_H' or VARIABLE == 'EXCH_M':
    
    plt.title(Filobjekt.variables[VARIABLE].description + ' ' + TYPE_OF_RUN + ' ' + DATE + ' (m2/s)')
    
    
else:

    plt.title(Filobjekt.variables[VARIABLE].description + ' ' + TYPE_OF_RUN + ' ' + DATE + ' (g/kg) REFERENCE')


plt.xlabel('Time (Minutes)')

plt.ylabel('Height (m)')

ax3.set_xticks(N.arange(0,(len(TIMESTEPS)/3), (len(TIMESTEPS)-1)/42))  #Först skapa array upp till 2520, då antal tidssteg är 7561 och var 20e sekund och då ENHETEN ÄR MINUTER(högsta värdet på x-axeln är 2520).
                                                                               #sedan plocka var 180e element, d v s var 3e timme

ax3.set_xticklabels(ALLA_DATUM[0::int((len(TIMESTEPS)-1)/14)], rotation=45, ha='right') #Välja ut från 0 till slutet av 'length timestep' vart 540e element = var 3e h


myplot.set_clim(vmin=0)

plt.ylim(0, 2500)
#plt.xlim(0,2500)   #Minuter....., inte antal tidssteg i vektorn

#plt.show()

'''Spar på olika ställen beroende på om det är 46 eller 91 nivåer'''
#fig2.savefig('/home/sm_marha/FIGURES/SINGLE_COLUMN/'+DATE +'Z/' + TYPE_OF_RUN + '/' + DATE +'z_CROSSSECTION' + VARIABLE + INTERPOLATION)
#fig2.savefig('/home/sm_marha/FIGURES/SINGLE_COLUMN/'+DATE +'Z/91_LEVELS/' + TYPE_OF_RUN + '/' + DATE +'z_CROSSSECTION' + VARIABLE + '_' + INTERPOLATION+ '_NOT_Xu_1e-4')



##########
#FIGURE 3#
##########


fig3 = plt.figure(3)

t=7561



ax1 = plt.subplot(211)

ax1.plot(TIMESTEPS[0:t], SWDOWN[0:t], linewidth = 1, label = 'SWDOWN')
ax1.plot(TIMESTEPS[0:t], GRDFLX[0:t], linewidth = 1, label = 'GRDFLX')
ax1.plot(TIMESTEPS[0:t], HFX[0:t], linewidth = 1, label = 'HFX')
ax1.plot(TIMESTEPS[0:t], LH[0:t], linewidth = 1, label = 'LH')
ax1.plot(TIMESTEPS[0:t], SWUP[0:t], linewidth = 1, label = 'SWUP')
ax1.plot(TIMESTEPS[0:t], LONGW_RESIDUAL[0:t], linewidth = 1, label = 'LW_RES')
ax1.plot(TIMESTEPS[0:t], RESIDUAL[0:t], linewidth = 1, label = 'RESIDUAL')
ax1.plot(TIMESTEPS[0:t], NOAHRES[0:t], linewidth = 1, label = 'NOAHRES')

ax1.grid(linestyle="dotted")
ax1.legend(loc=4, prop={'size': 8})
ax1.set_xticklabels([])


plt.title(TYPE_OF_RUN +' '+ DATE)
plt.ylabel('Flux (W/m2)')



ax2 = plt.subplot(212)

ax2.plot(TIMESTEPS[0:t], EMISS[0:t]*GLW[0:t], linewidth = 1, label = 'GLW')
ax2.plot(TIMESTEPS[0:t], GUPLW[0:t], linewidth = 1, label = 'GUPLW')


ax2.grid(linestyle="dotted")
ax2.legend(loc=4, prop={'size': 8})
plt.xlabel('Time (Minutes)')
plt.ylabel('Flux (W/m2)')

'''Spar på olika ställen beroende på om det är 46 eller 91 nivåer'''
#fig3.savefig('/home/sm_marha/FIGURES/SINGLE_COLUMN/' + DATE +'Z/' + TYPE_OF_RUN + '/' + DATE +'z_FLUXES' + INTERPOLATION)
#fig3.savefig('/home/sm_marha/FIGURES/SINGLE_COLUMN/' + DATE +'Z/91_LEVELS/' + TYPE_OF_RUN + '/NEW_FIGURES/' + DATE +'z_FLUXES_' + INTERPOLATION + '_NOT_Xu')

plt.show()



'''

##########
#FIGURE 4#
##########

fig4 = plt.figure(4)

t=1080



ax1 = plt.subplot(211)


ax1.plot(TIMESTEPS[0:t], FLX1[0:t], linewidth = 1, label = 'PRECIP-SNOW SFC')
ax1.plot(TIMESTEPS[0:t], FLX2[0:t], linewidth = 1, label = 'FREEZING RAIN')
ax1.plot(TIMESTEPS[0:t], FLX3[0:t], linewidth = 1, label = 'SNOWMELT')

ax1.grid(linestyle="dotted")
ax1.legend(loc=4, prop={'size': 8})
plt.title('REFERENCE')
plt.ylabel('Flux (W/m2)')

#fig4.savefig('/home/sm_marha/FIGURES/SINGLE_COLUMN/18011800Z/ADDING_4xQVAPOR/18011800z_FLUXES_FLX') '''



