import numpy as N
import scipy.io.netcdf as S
import matplotlib.pyplot as plt





File_object = S.netcdf_file('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY/test/em_scm_xy/wrfout_d01_2018-02-18_00:00:00_REFERENCE')  #Reading NetCDF file

TIME_RANGE = 2160      #Setting time range to plot.

TIMESTEPS = File_object.variables['XTIME'][0:TIME_RANGE]



TSK = TSK = File_object.variables['TSK'][0:TIME_RANGE, 1, 1]-273.15    #Change from Kelvin to Celsius
EMISS = File_object.variables['EMISS'][0:TIME_RANGE, 1, 1]


'''Extracting variables from netcdf file'''

GLW = File_object.variables['GLW'][0:TIME_RANGE, 1, 1]
SWDOWN = File_object.variables['SWDOWN'][0:TIME_RANGE, 1, 1]
GRDFLX = File_object.variables['GRDFLX'][0:TIME_RANGE, 1, 1]
HFX = File_object.variables['HFX'][0:TIME_RANGE, 1, 1]
LH = File_object.variables['LH'][0:TIME_RANGE, 1, 1]
EMISS = File_object.variables['EMISS'][0:TIME_RANGE, 1, 1]
ALBEDO = File_object.variables['ALBEDO'][0:TIME_RANGE, 1, 1]
NOAHRES = File_object.variables['NOAHRES'][0:TIME_RANGE, 1, 1]
FLX1 = File_object.variables['FLX1'][0:TIME_RANGE, 1, 1]
FLX2 = File_object.variables['FLX2'][0:TIME_RANGE, 1, 1]
FLX3 = File_object.variables['FLX3'][0:TIME_RANGE, 1, 1]
        
SWUP = ALBEDO*SWDOWN
GUPLW = EMISS*5.67e-8*(TSK+273.15)**4 #Calculating upwelling longwave with the help of TSKIN
GUPLW[0] = 0    #Setting my own calculated round upwelling long wave equal to zero at time=0, because the other fluxes is zero at time_step zero          

LONGW_RESIDUAL = GLW - GUPLW


RESIDUAL = EMISS*GLW + SWDOWN - SWUP - GUPLW - HFX - LH + GRDFLX - FLX1 - FLX2 - FLX3      #Calculating residual in the same way that it is done in Noah driver


'''This is how NOAHRES is calculated in module_sf_noahdrv.F

noahres(i,j) = ( solnet + lwdn ) - sheat + ssoil - eta &
              - ( emissi * STBOLT * (t1**4) ) - flx1 - flx2 - flx3 - flx4'''


#################################### 
#            PLOTTING              # 
#################################### 


fig1 = plt.figure(1)

t=2160





#Plotting different fluxes and radiation terms.

fig1 = plt.figure(1)

ax1 = plt.subplot(111)

ax1.plot(TIMESTEPS[0:t], SWDOWN[0:t], linewidth = 1, label = 'SWDOWN')
ax1.plot(TIMESTEPS[0:t], GRDFLX[0:t], linewidth = 1, label = 'GRDFLX')
ax1.plot(TIMESTEPS[0:t], HFX[0:t], linewidth = 1, label = 'HFX')
ax1.plot(TIMESTEPS[0:t], LH[0:t], linewidth = 1, label = 'LH')
ax1.plot(TIMESTEPS[0:t], SWUP[0:t], linewidth = 1, label = 'SWUP')
ax1.plot(TIMESTEPS[0:t], LONGW_RESIDUAL[0:t], linewidth = 1, label = 'LW_RES')
ax1.grid(linestyle="dotted")
ax1.legend(loc=4, prop={'size': 8})
plt.title('Different terms in energy buget')



#Plotting my own calculated NOAHRES=RESIDUAL and models calculation NOAHRES

fig2 = plt.figure(2)

ax2 = plt.subplot(211)

ax2.plot(TIMESTEPS[0:t], NOAHRES[0:t], linewidth = 1, label = 'NOAHRES')
ax2.plot(TIMESTEPS[0:t], RESIDUAL[0:t], linewidth = 1, label = 'RESIDUAL')
ax2.grid(linestyle="dotted")
ax2.legend(loc=4, prop={'size': 8})
ax2.set_xticklabels([])
plt.title('My residual vs model residual')
plt.ylabel('Flux (W/m2)')



#Plotting upwelling and downwelling longwave

ax3 = plt.subplot(212)

ax3.plot(TIMESTEPS[0:t], GLW[0:t], linewidth = 1, label = 'GLW')
ax3.plot(TIMESTEPS[0:t], GUPLW[0:t], linewidth = 1, label = 'GUPLW')
ax3.grid(linestyle="dotted")
ax3.legend(loc=4, prop={'size': 8})
plt.xlabel('Time (Minutes)')
plt.ylabel('Flux (W/m2)')
plt.title('Long wave upwelling vs long wave downwelling')


plt.show()






