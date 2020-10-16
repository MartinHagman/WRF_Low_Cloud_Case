import numpy as N
import scipy.io.netcdf as S
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import os
import datetime
import matplotlib.colors as mplc
import netCDF4 as N4
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
from pylab import *
import Picking_values_Sounding_Dat_file_version_2
import pandas as pd


######################
# CREATING AXES GRID #
#######################################################

rc('axes', linewidth=1.5)



fig1 = plt.figure(1,figsize=(12,6))
fig1.subplots_adjust(hspace=0.3)
fig1.subplots_adjust(wspace=0.5)
fig1.subplots_adjust(right=0.98)
fig1.subplots_adjust(left=0.08)
fig1.subplots_adjust(bottom=0.15)
fig1.subplots_adjust(top=0.90)

gs = gridspec.GridSpec(ncols=20, nrows=4)


AXS = [ fig1.add_subplot(gs[0:3, 0:11]), fig1.add_subplot(gs[3, 0:11]), fig1.add_subplot(gs[0:3,11:14]), fig1.add_subplot(gs[0:3, 14:17]), fig1.add_subplot(gs[0:3, 17:20])]

#cbar_ax = fig1.add_axes([0.92, 0.31, 0.01, 0.56])
cbar_ax = fig1.add_axes([0.64, 0.16, 0.34, 0.03])


AXS[0].set_xticklabels([])
AXS[2].set_yticklabels([])
AXS[3].set_yticklabels([])
AXS[4].set_yticklabels([])

AXS[0].tick_params(which='major', direction='in')
AXS[1].tick_params(which='major', direction='in')
AXS[2].tick_params(which='major', direction='in')
AXS[3].tick_params(which='major', direction='in')
AXS[4].tick_params(which='major', direction='in')
#cbar_ax.tick_params(which='major', direction='in')





#############
# SCM MODEL #
###########################################################


#Filobjekt = N4.Dataset('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfout_d01_2018-02-27_00:00:00_025T_LIQ_AND_VARIATION_YSU_TD0')
Filobjekt = N4.Dataset('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfout_d01_2018-02-27_00:00:00_T_LIQ_YSU_TOPD0_SOUNDING_SODANKYLA_QVAPOR_100')

TIDSLANGD = 2701   #7561 totalt....      540, 1080, 1620, 2160, 2701, 3241, 3781, 4321, 4861.......osv
TIDPUNKT = 1081

   

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


QCLOUD = Filobjekt.variables['QCLOUD'][0:TIDSLANGD, 0, 1, 1]*1000
QCLOUD_2D = Filobjekt.variables['QCLOUD'][0:TIDSLANGD, :, 1, 1]*1000


QRAIN = Filobjekt.variables['QRAIN'][0:TIDSLANGD, 0, 1, 1]*1000
QRAIN_2D = Filobjekt.variables['QRAIN'][0:TIDSLANGD, :, 1, 1]*1000


QSNOW = Filobjekt.variables['QSNOW'][0:TIDSLANGD, 0, 1, 1]*1000
QSNOW_2D = Filobjekt.variables['QSNOW'][0:TIDSLANGD, :, 1, 1]*1000


QICE = Filobjekt.variables['QICE'][0:TIDSLANGD, 0, 1, 1]*1000
QICE_2D = Filobjekt.variables['QICE'][0:TIDSLANGD, :, 1, 1]*1000


QGRAUP = Filobjekt.variables['QGRAUP'][0:TIDSLANGD, 0, 1, 1]
QGRAUP_2D = Filobjekt.variables['QGRAUP'][0:TIDSLANGD, :, 1, 1]

                                                                                                             


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
'''
TKE_PBL = Filobjekt.variables['TKE_PBL'][0:TIDSLANGD, 0, 1, 1]
TKE_PBL_2D = Filobjekt.variables['TKE_PBL'][0:TIDSLANGD, :, 1, 1]
TKE_PBL_2D = N.transpose(TKE_PBL_2D)

DTKE = Filobjekt.variables['DTKE'][0:TIDSLANGD, 0, 1, 1]                   #ÄNDRA TKE
QSHEAR = Filobjekt.variables['QSHEAR'][0:TIDSLANGD, 0, 1, 1]
QBUOY = Filobjekt.variables['QBUOY'][0:TIDSLANGD, 0, 1, 1]
QDISS = Filobjekt.variables['QDISS'][0:TIDSLANGD, 0, 1, 1]
QWT = Filobjekt.variables['QWT'][0:TIDSLANGD, 0, 1, 1]
'''



QTOTAL_2D = QVAPOR_2D + QCLOUD_2D + QICE_2D + QRAIN_2D + QSNOW_2D
POISSON_C = 0.2854*(1-0.24*QVAPOR_2D)
GAMMA = QVAPOR_2D*461/(1005.7+QTOTAL_2D*1.996e3)


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



U = Filobjekt.variables['U'][:]
V = Filobjekt.variables['V'][:]
U_10 = Filobjekt.variables['U10'][:]
V_10 = Filobjekt.variables['V10'][:]
XLAT = Filobjekt.variables['XLAT'][:]
XLONG = Filobjekt.variables['XLONG'][:]
X_LAT_U = Filobjekt.variables['XLAT_U'][:]
X_LAT_U = Filobjekt.variables['XLAT_U'][:]
X_LON_U = Filobjekt.variables['XLONG_U'][:]
X_LAT_V = Filobjekt.variables['XLAT_V'][:]
X_LONG_U = Filobjekt.variables['XLONG_V'][:]
SINALPHA = Filobjekt.variables['SINALPHA'][:]
COSALPHA = Filobjekt.variables['COSALPHA'][:]

SINALPHA = SINALPHA[TIDPUNKT, :, :]      #Detta är ett tillägg i SCM jämfört med 3D
COSALPHA = COSALPHA[TIDPUNKT, :, :]      #Detta är ett tillägg i SCM jämfört med 3D'''

U_3D = U[TIDPUNKT, :, :, :] 
V_3D = V[TIDPUNKT, :, :, :]

U_3D_UNSTAG = 0.5*(U_3D[:,:,:-1] + U_3D[:,:,1:])
V_3D_UNSTAG = 0.5*(V_3D[:,:-1,:] + V_3D[:,1:,:])

'''print(N.shape(U_3D_UNSTAG))
print(N.shape(V_3D_UNSTAG))'''

U_3D_UNSTAG_EARTH = U_3D_UNSTAG*COSALPHA - V_3D_UNSTAG*SINALPHA
V_3D_UNSTAG_EARTH = V_3D_UNSTAG*COSALPHA + U_3D_UNSTAG*SINALPHA

U_SOD = U_3D_UNSTAG_EARTH[:, 1, 1]
V_SOD = V_3D_UNSTAG_EARTH[:, 1, 1]

VELOCITY_SOD = N.sqrt(U_SOD**2+V_SOD**2)

RESIDUAL = EMISS*GLW + SWDOWN + SWUP - GUPLW + HFX - LH + GRDFLX - FLX1 - FLX2 - FLX3



ZNW = Filobjekt.variables['ZNW'][0:TIDSLANGD, :]

P = Filobjekt.variables['P'][0:TIDSLANGD, :, 1, 1]

PB = Filobjekt.variables['PB'][0:TIDSLANGD, :, 1, 1]



PH = Filobjekt.variables['PH'][:, :, 1, 1]
PHB = Filobjekt.variables['PHB'][:, :, 1, 1]
MODELLNIVAHOJD = (PH+PHB)/9.81
MODELLNIVAHOJD_TER = MODELLNIVAHOJD[0:TIDSLANGD, :] -MODELLNIVAHOJD[0, 0]


MASSLEVELS = 0.5*(MODELLNIVAHOJD_TER[0:TIDSLANGD, :-1] + MODELLNIVAHOJD_TER[0:TIDSLANGD, 1:])  #Olika masslevels för varje tidssteg, om än mycket små skillnader!


'''Det vi nu vill göra är att få massnivåer som varierar i tiden och inte är statiskt lika som vid tidpunkt 1'''

'''Använder massnivåerna vid tidpunkt 0 för att skapa matriser TIME och LEVELS med rätt dimension'''
TIME, LEVELS = N.meshgrid(TIMESTEPS, MASSLEVELS[0, :])   #Behövs för contour Theta på massnivåerna(fulla nivåerna) samt contourf

'''Nollställer LEVELS med med fortsatt samma dimension'''

LEVELS = N.zeros_like(LEVELS)


'''Skriver in massnivåer för de 7561 olika tidpunkterna. Antal rows är 90 stycken. 0-89 om man indexerar enligt Python'''
for row in range(N.shape(LEVELS)[0]):
    LEVELS[row, :] = (N.transpose(MASSLEVELS))[row, :]


'''Plottar linje vid nedersta molnvattnet'''

CLOUDBASE = N.zeros(N.shape(QCLOUD_2D[1])[0])

for i in range((N.shape(QCLOUD_2D))[1]):
    for place, k in enumerate(range((N.shape(QCLOUD_2D))[0])):
        if QCLOUD_2D[k, i] > 0 and k > 4:
            CLOUDBASE[i] = MASSLEVELS[i, k]
            break

        


################# PLOTTING SCM ########################
FIELD_2D = QCLOUD_2D
########## CLOUD WATER; THETA #######################

norm = mplc.Normalize(0, 4.2e-1)
A = 4.2e-1/40
myplot = AXS[0].contourf(TIME, LEVELS, FIELD_2D, N.arange(0, 4.2e-1, A), norm = norm, extend = 'both', cmap ='binary')
#myplot = AXS[0].contourf(TIME, LEVELS, FIELD_2D, extend = 'both', cmap ='binary')

'''Välj interpolering i fälten med contourf, eller plotta med scatter för att se exakt på vilka modellnivåer det finns moln'''
#myplot = AXS[0].contourf(TIME, LEVELS, VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]], extend = 'both', cmap ='binary')
#myplot = AXS[0].scatter(TIME, LEVELS, c=VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]], cmap='binary')
#myplot = AXS[0].pcolormesh(TIME, LEVELS, VARIABLE_DICT[list(VARIABLE_DICT.keys())[0]], vmin=-0.5, vmax=0.5, cmap = 'seismic')


'''Plottar contours av Theta tillsammans med contourf'''
contour_levels = N.arange(0, 330, 1)
cp = AXS[0].contour(TIME, LEVELS, THETA_2D, contour_levels, colors='red', linestyles='solid', linewidths=0.8)
#cp = plt.contour(TIME, LEVELS, THETA_L_2D, contour_levels, colors='yellow', linestyles='dotted', linewidth='0.2' )
plt.clabel(cp,contour_levels[::2],inline=True, fmt='%1.0f', fontsize=10)


'''Plottar masslevels längs y-axeln''' 
X = N.zeros_like(MASSLEVELS)[0, :]+3
#ax1.scatter(X, MASSLEVELS[0, :], marker= '>', color='black', s=2)
AXS[0].set_xlim(-4, N.max(TIME))
AXS[0].set_ylim(0, 2500)
AXS[0].grid(linestyle="dotted")
AXS[0].set_xticks(N.arange(0, len(ALLA_DATUM)/3, 180 ))
AXS[0].set_ylabel('Height (m)', fontsize = 16)
AXS[0].text(7, 2350, 'a', fontsize = 18, fontweight='bold')  

cbar = fig1.colorbar(myplot, cax=cbar_ax, orientation = 'horizontal', ticks=[0, 0.1, 0.2, 0.3, 0.4])#,  format = "%8.1e")
cbar.ax.tick_params(labelsize ='x-large')
cbar_ax.set_title(' ($gkg^{-1})$', fontsize = 16)

#cbar = fig1.colorbar(cp, ticks=[0.1, 0.2, 0.3, 0.4, 0.5])
#cbar.ax.set_yticklabels(['< -1', '0', '> 1'])  # vertically oriented colorbar




'''Plottar linje för var tvärsnittet läggs '''
Y = MASSLEVELS[TIDPUNKT, :]
X = N.ones(len(Y))*(TIDPUNKT-1)/3           # TIDPUNKT är antal steg i tiden, alltså antal 20 sekunders-steg. På x-axen är det dock minuter, därav delat på 3.
AXS[0].plot(X, Y, c='k', linestyle = 'dotted')


'''Plottar masslevels längs y-axeln''' 
XX = N.zeros_like(MASSLEVELS)[0, :]
AXS[0].scatter(XX, MASSLEVELS[0, :], marker= '>', color='black', s=3)


'''Plottar molnbas'''

AXS[0].plot(TIMESTEPS[0:TIDSLANGD], CLOUDBASE[0:TIDSLANGD],c='k', linestyle='dashed', linewidth = 1)
AXS[0].tick_params(axis="y", labelsize=12)









############# PLOTTING SKIN TEMPERATURE ######################


AXS[1].plot(TIMESTEPS[0:TIDSLANGD], TSK[0:TIDSLANGD], linewidth = 2, linestyle='dashed', c='C0')
AXS[1].plot(TIMESTEPS[0:TIDSLANGD], T2M[0:TIDSLANGD], linewidth = 2, c = 'C0')


AXS[1].set_xticklabels(ALLA_DATUM_3h, rotation=35, fontsize=12)
AXS[1].tick_params(axis="y", labelsize=12)


AXS[1].grid(linestyle="dotted")
AXS[1].set_xticks(N.arange(0, len(ALLA_DATUM)/3, 180 ))
AXS[1].set_xlim(-4, N.max(TIME))
AXS[1].set_ylabel('T $(^oC)$', fontsize = 15)
AXS[1].text(870, -7.2, 'e', fontsize = 18, fontweight='bold') 


############ PLOTTING CLOUD WATER, THETA ###########################


MASSLEVELS = N.transpose(MASSLEVELS)


lns1 = AXS[2].plot(QCLOUD_2D[:, TIDPUNKT], MASSLEVELS[:, TIDPUNKT], linewidth = 3, linestyle='dotted', c='C0', label = 'CW')
AXS[2].grid(linestyle="dotted")
#AXS[2].set_yscale('log')
AXS[2].set_ylim(0, 2500)
AXS[2].set_xlabel('CW $(gkg^{-1})$', fontsize=13, labelpad = -1)
#AXS[2].text(0, 2300, 'b', fontsize = 18, fontweight='bold')
AXS[2].tick_params(axis="x", labelsize=11)


AXS22 = AXS[2].twiny()
lns2 = AXS22.plot(THETA_2D[:, TIDPUNKT], MASSLEVELS[:, TIDPUNKT],linewidth = 2, c='C0', label = '$\Theta$' )
AXS22.set_ylim(0, 2500)
AXS22.set_xlim(260, 280)
AXS22.set_xlabel('$\Theta$ $(K)$', fontsize=12)
AXS22.text(260.1, 2330, 'b', fontsize = 18, fontweight='bold')
AXS22.tick_params(axis="x", labelsize=11)

############ PLOTTING WIND VELOCITY ###########################

AXS[3].plot(VELOCITY_SOD, MASSLEVELS[:, TIDPUNKT], linewidth = 2, label = 'U')
AXS[3].grid(linestyle="dotted")
#AXS[3].set_yscale('log')
AXS[3].set_ylim(0, 2500)
AXS[3].set_xlabel('Wind $(ms^{-1}$)', fontsize=13, labelpad = -1)
AXS[3].set_xlim(0, 20)
#AXS[3].text(1, 2350, 'c', fontsize = 18, fontweight='bold')
AXS[3].tick_params(axis="x", labelsize=11)
AXS[3].legend()


########### PLOTTING EXCHANGE COEFFICIENTS ###################

MODELLNIVAHOJD_TER = N.transpose(MODELLNIVAHOJD_TER)


lnsA = AXS[4].plot(EXCH_H_2D[:, TIDPUNKT], MODELLNIVAHOJD_TER[:, TIDPUNKT], linewidth = 2, label = 'K') #Ändra TKE
#lnsB = AXS[4].plot(EXCH_M_2D[:, TIDPUNKT], MODELLNIVAHOJD_TER[:, TIDPUNKT], linestyle='dashed', c='C0', linewidth = 1.4, label = 'Km')   
AXS[4].legend()
AXS[4].grid(linestyle="dotted")
AXS[4].set_ylim(0, 2500)
AXS[4].set_xlabel('$K$ $(m^2s^{-1})$', fontsize=13, labelpad = -1)
AXS[4].text(0.1, 2330, 'd', fontsize = 18, fontweight='bold')
AXS[4].tick_params(axis="x", labelsize=11)
#AXS[4].set_xlim(0, 15)

#AXS44 = AXS[4].twiny() #ÄNDRA TKE
#lnsC = AXS44.plot(TKE_PBL_2D[:, TIDPUNKT], MODELLNIVAHOJD_TER[:, TIDPUNKT], linestyle='dotted', c='C0', linewidth = 2, label = 'TKE')      #ÄNDRA TKE
#AXS44.set_ylim(0, 2500)
#AXS44.set_xlabel('$TKE$ $(m^2s^{-2})$', fontsize=12)        #ÄNDRA TKE
#AXS44.tick_params(axis="x", labelsize=11)

 
#AXS44.set_xticklabels(['0', '0', '2.5e-2', '  5.0e-2'])


lns = lnsA #+ lnsB + lnsC        #ÄNDRA TKE

labs = [l.get_label() for l in lns]

AXS[4].legend(lns, labs, loc='upper right', fontsize = 'large')



############## PLOTTING LOGWAVE RADIATION ######################

AXS11 = AXS[1].twinx()
AXS11.plot(TIMESTEPS[0:TIDSLANGD], GLW[0:TIDSLANGD], linewidth = 2.5, linestyle='dotted', c='C0')
AXS11.set_ylabel('$LW$ $(Wm^{-2})$     ', fontsize=11)
AXS11.tick_params(axis="y", labelsize=11)


#############
# SOUNDINGS #
###########################################################


j=58    #ÄNDRA (25 FEB)          #ÄNDRA

    
SOUNDING_TIME = '6'


'''Hämtar sonderingsdata med hjälp av funktionen 'pick dates' i Reading_Dat_Files'''


MATRIX_2, datum = Picking_values_Sounding_Dat_file_version_2.pick_dates(j, SOUNDING_TIME)


MATRIX_2 = (N.array(MATRIX_2)).astype(float)




'''Beräknar relativ fuktighet och specifik fuktighet'''



VV_SOND = 530  # 470   #2000



TEMP_C_SONDERING = MATRIX_2[0:VV_SOND, 4]

TEMP_K_SONDERING = TEMP_C_SONDERING + 273.15

DAGGPUNKT_SONDERING = MATRIX_2[0:VV_SOND, 7]

P_TEMPERATURE_SONDERING = MATRIX_2[0:VV_SOND, 6] * 100  #Gör om från hPa till Pa

P_DEWPOINT_SONDERING = MATRIX_2[0:VV_SOND, 9] * 100  #Gör om från hPa till Pa

E_MATTNAD_SONDERING = N.exp(N.log(611.2)+(17.62*TEMP_C_SONDERING/(243.12+TEMP_C_SONDERING)))

QVAPOR_MATTNAD_SONDERING = 0.622*E_MATTNAD_SONDERING/(P_TEMPERATURE_SONDERING-E_MATTNAD_SONDERING)

E_AKTUELL_SONDERING = N.exp(N.log(611.2)+(17.62*DAGGPUNKT_SONDERING/(243.12+DAGGPUNKT_SONDERING)))

RH_SONDERING = E_AKTUELL_SONDERING/E_MATTNAD_SONDERING

QVAPOR_SONDERING = 0.622*E_AKTUELL_SONDERING/(P_TEMPERATURE_SONDERING-E_AKTUELL_SONDERING)

THETA_SONDERING = TEMP_K_SONDERING *((1000/(P_TEMPERATURE_SONDERING/100))**0.286)




HEIGHT_TER_SONDERING = MATRIX_2[0:VV_SOND, 5]-180

WIND_VELOCITY_SONDERING = MATRIX_2[0:VV_SOND, 0]

WIND_VELOCITY_SONDERING = MATRIX_2[0:VV_SOND, 0]

WIND_DIRECTION_SONDERING = MATRIX_2[0:VV_SOND, 2]

U_WIND_SONDERING = -WIND_VELOCITY_SONDERING*N.sin(pi/180*WIND_DIRECTION_SONDERING)

V_WIND_SONDERING = -WIND_VELOCITY_SONDERING*N.cos(pi/180*WIND_DIRECTION_SONDERING)

ALTITUDE_WIND_SONDERING = MATRIX_2[0:VV_SOND, 1]


#####################  PLOTTING SOUNDINGS  ###############################



lns3 = AXS22.plot(THETA_SONDERING[0:VV_SOND], HEIGHT_TER_SONDERING[0:VV_SOND],linewidth = 3, c = 'k', label = '$\Theta$ OBS' )

lns = lns1 + lns2 + lns3

labs = [l.get_label() for l in lns]

AXS[2].legend(lns, labs, loc='upper left')#, fontsize = 'large')



AXS[3].plot(WIND_VELOCITY_SONDERING[0:VV_SOND], HEIGHT_TER_SONDERING[0:VV_SOND],linewidth = 3, c = 'k', label='U OBS')

AXS[3].legend(loc=1)

AXS[3].text(1, 2350, 'c', fontsize = 18, fontweight='bold')



#####################
# AUTOMATIC STATION #
#################################################################################################################
print('Automatic station')


df = pd.read_csv(r'/home/sm_marha/TEXTFILER/Automatic_Station_jan_feb_all_parameters.txt', encoding='latin-1', delimiter=',')

df = df.loc[:,['DATE TIME','T']]

df_numpy_array = df.values

TIDSVEKTOR_OBS = df_numpy_array[2:, 0]

TIDSVEKTOR_MINUTER_OBS = []

'''Letar upp platserna i arrayen för ändpunkterna på de datum jag är intresserad av.'''


for place, date in enumerate(TIDSVEKTOR_OBS):
    if date == '2018-02-27 00:00:00+00':                #ÄNDRA
        print(place)
        date_START = place
    if date == '2018-02-27 15:00:00+00':                #ÄNDRA
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

for place, temp in enumerate(T2_TIDSSERIE_OBS):
    T2_TIDSSERIE_OBS[place] = float(temp)


for place, temp in enumerate(T2_GLES_OBS):
    T2_GLES_OBS[place] = float(temp)


TIDSVEKTOR_NUMBER_OBS = N.arange(0,(date_END-date_START)*10+1, 10)



#####################  Plotting Automatic station T2m  ###############################

AXS[1].plot(TIDSVEKTOR_NUMBER_OBS, T2_TIDSSERIE_OBS, linewidth=3, c = 'k')     

#AXS[0].set_title('27, 005T-LIQ(And variation), YSU-0 ')







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
    if date == '180227 00z':                                #ÄNDRA
        print(place)
        datum_START = place
    if date == '180227 15z':                                #ÄNDRA
        print(place)
        datum_END = place


''' Anpassa strålningsobservationerna till en xaxel i minuter'''

TID = datum_END - datum_START

TIMESTEPS_NC = N.arange(0, 901, 60)

################### PLOTTING NC-FILE Eirikur ###############################################


AXS11.plot(TIMESTEPS_NC, LONGWAVE_DOWN_SURFACE_OBS[datum_START:datum_END+1], linestyle='dotted', linewidth=3, c = 'k')

AXS11.set_xlim(-4, N.max(TIMESTEPS_NC))


#fig1.savefig('/home/sm_marha/FIGURES/FINAL/20180227/27_YSU_VARIATION_of_Theta_liq')
#fig1.savefig('/home/sm_marha/FIGURES/FINAL/20180227/27_YSU_VARIATION_of_Theta_liq_plus_Sodankyla_temp_from_sounding')

plt.show()

