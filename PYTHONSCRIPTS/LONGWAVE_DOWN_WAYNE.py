import numpy as N
import netCDF4 as N4
import matplotlib.pyplot as plt
import datetime
from matplotlib.colors import LogNorm

Filobjekt = N4.Dataset('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfout_d01_2018-02-26_00:00:00_QC_QI_REFERENCE')
Filobjekt2 = N4.Dataset('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfout_d01_2018-02-26_00:00:00_QC_QI_REFERENCE_Icloud2')

TIDSLANGD = 7561

V_LEVELS = 35

LONGWAVE_DOWN_SURFACE_WRF_Icloud1 = Filobjekt.variables['GLW'][:, 1, 1]
LONGWAVE_DOWN_SURFACE_WRF_Icloud2 = Filobjekt2.variables['GLW'][:, 1, 1]

CLOUDWATER_Icloud1 = Filobjekt.variables['QCLOUD'][:, :V_LEVELS, 1, 1]
CLOUDWATER_Icloud2 = Filobjekt.variables['QCLOUD'][:, :V_LEVELS, 1, 1]

CLOUDWATER_Icloud1 = N.transpose(CLOUDWATER_Icloud1)
CLOUDWATER_Icloud2 = N.transpose(CLOUDWATER_Icloud2)

PH = Filobjekt.variables['PH'][:, :, 1, 1]
PHB = Filobjekt.variables['PHB'][:, :, 1, 1]
MODELLNIVAHOJD = (PH+PHB)/9.81
MODELLNIVAHOJD_TER = MODELLNIVAHOJD[0:TIDSLANGD, :] -MODELLNIVAHOJD[0, 0]






TIMESTEPS = Filobjekt.variables['XTIME'][0:TIDSLANGD]
print(N.shape(TIMESTEPS))
time = Filobjekt.variables['XTIME']
dates = N4.num2date(TIMESTEPS, time.units)

ALLA_DATUM=[]

for date in dates:
    datum = date.strftime('%d/%-m %Hz')
    ALLA_DATUM.append(datum)

ALLA_DATUM_3h = ALLA_DATUM[0:len(ALLA_DATUM):540]



MASSLEVELS = 0.5*(MODELLNIVAHOJD_TER[0:TIDSLANGD, :-1] + MODELLNIVAHOJD_TER[0:TIDSLANGD, 1:])       #Olika masslevels för varje tidssteg, om än mycket små skillnader!

'''Det vi nu vill göra är att få massnivåer som varierar i tiden och inte är statiskt lika som vid tidpunkt 1'''

'''Använder massnivåerna vid tidpunkt 0 för att skapa matriser TIME och LEVELS med rätt dimension'''
TIME, LEVELS = N.meshgrid(TIMESTEPS, MASSLEVELS[0, 0:V_LEVELS])   #Behövs för contour Theta på massnivåerna(fulla nivåerna) samt contourf
print(N.shape(TIMESTEPS))

'''Nollställer LEVELS med med fortsatt samma dimension'''
LEVELS = N.zeros_like(LEVELS)


'''Skriver in massnivåer för de 7561 olika tidpunkterna. Antal rows är 90 stycken. 0-89 om man indexerar enligt Python'''
for row in range(N.shape(LEVELS)[0]):
    LEVELS[row, :] = (N.transpose(MASSLEVELS))[row, :]





fig1 = plt.figure(1, figsize=(12,5))
fig1.subplots_adjust(bottom=0.3)

t=7561


ax1 = plt.subplot(111)



ax1.plot(TIMESTEPS[0:t], LONGWAVE_DOWN_SURFACE_WRF_Icloud1[0:t], linewidth = 1, label = 'Icloud 1')
ax1.plot(TIMESTEPS[0:t], LONGWAVE_DOWN_SURFACE_WRF_Icloud2[0:t], linewidth = 1, label = 'Icloud 2')

ax1.legend(fontsize = 'xx-large')
ax1.grid(linestyle="dotted")
ax1.set_xticks(N.arange(0, len(ALLA_DATUM)/3, 180 ))
ax1.set_xticklabels(ALLA_DATUM_3h, rotation=45, ha='right', fontsize=15)
ax1.set_ylabel('Downwelling Longwave Radiation')
ax1.set_ylim(190, 270)
ax1.set_xlabel('Time', fontsize=18)
ax1.set_title('Comparing Icloud1 and Icloud2', fontsize=18)






fig2 = plt.figure(2, figsize=(12,5))
fig2.subplots_adjust(bottom=0.3)

t=7561


ax1 = plt.subplot(111)


myplot = ax1.contourf(TIME, LEVELS,CLOUDWATER_Icloud1 , extend = 'both', cmap ='binary')

cbar_ax = fig2.add_axes([0.92, 0.30, 0.01, 0.58])
cbar = fig2.colorbar(myplot, cax=cbar_ax)



ax1.grid(linestyle="dotted")
ax1.set_xticks(N.arange(0, len(ALLA_DATUM)/3, 180 ))
ax1.set_xticklabels(ALLA_DATUM_3h, rotation=45, ha='right', fontsize=15)
ax1.set_ylabel('Cloud water', fontsize=18)
#ax1.set_ylim(190, 270)
ax1.set_xlabel('Time', fontsize=18)
ax1.set_title('Icloud1', fontsize=18)





fig3 = plt.figure(3, figsize=(12,5))
fig3.subplots_adjust(bottom=0.3)

t=7561


ax1 = plt.subplot(111)


myplot = ax1.contourf(TIME, LEVELS, CLOUDWATER_Icloud2 , extend = 'both', cmap ='binary')

cbar_ax = fig3.add_axes([0.92, 0.30, 0.01, 0.58])
cbar = fig3.colorbar(myplot, cax=cbar_ax)



ax1.grid(linestyle="dotted")
ax1.set_xticks(N.arange(0, len(ALLA_DATUM)/3, 180 ))
ax1.set_xticklabels(ALLA_DATUM_3h, rotation=45, ha='right', fontsize=15)
ax1.set_ylabel('Cloud water', fontsize=18)
#ax1.set_ylim(190, 270)
ax1.set_xlabel('Time', fontsize=18)
ax1.set_title('Icloud2', fontsize=18)


plt.show()
