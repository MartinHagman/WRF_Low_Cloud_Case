import numpy as N
import scipy.io.netcdf as S
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import os
import datetime
import netCDF4 as N4



Filobjekt = N4.Dataset('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfinput_d01_20180227_QC_025THETA_LIQUID_AND_VARIATION_SOUNDING_YSU')


CLOUD_WATER_1D = Filobjekt.variables['QCLOUD'][0, :, 1, 1]*1000

THETA_PERT_1D = Filobjekt.variables['T'][0, :, 1, 1]
P_PERT_1D = Filobjekt.variables['P'][0, :, 1, 1]
P_BASE_1D = Filobjekt.variables['PB'][0, :, 1, 1]   

P_TOT_1D = P_BASE_1D + P_PERT_1D
THETA_1D = THETA_PERT_1D+300
TEMP_K_1D = THETA_1D/((1000/(P_TOT_1D/100))**0.286)
T_C_1D = TEMP_K_1D-273.15

QVAPOR_1D = Filobjekt.variables['QVAPOR'][0, :, 1, 1]
E_MATTNAD_1D = N.exp(N.log(611.2)+(17.62*T_C_1D/(243.12+T_C_1D)))
QVAPOR_MATTNAD_1D = 0.622*E_MATTNAD_1D/(P_TOT_1D-E_MATTNAD_1D)
RH_1D = QVAPOR_1D/QVAPOR_MATTNAD_1D
T_DAGGPUNKT_1D = (243.12*N.log(611.2)-243.12*N.log(RH_1D*E_MATTNAD_1D)) / (N.log(RH_1D*E_MATTNAD_1D)-17.62-N.log(611.2))


PH_WRF_3D_WRFOUT = Filobjekt.variables['PH'][0, :, 1, 1]
PHB_WRF_3D_WRFOUT = Filobjekt.variables['PHB'][0, :, 1, 1]
MODELLNIVAHOJD_WRF_3D_WRFOUT = (PH_WRF_3D_WRFOUT+PHB_WRF_3D_WRFOUT)/9.81
MODELLNIVAHOJD_TER_WRF_3D_WRFOUT = MODELLNIVAHOJD_WRF_3D_WRFOUT[:] -MODELLNIVAHOJD_WRF_3D_WRFOUT[0]




MASSLEVELS_WRF_3D_WRFOUT = 0.5*(MODELLNIVAHOJD_TER_WRF_3D_WRFOUT[:-1] + MODELLNIVAHOJD_TER_WRF_3D_WRFOUT[1:])

plt.figure(1)
plt.plot(CLOUD_WATER_1D[0: 30], MASSLEVELS_WRF_3D_WRFOUT[0: 30])
plt.scatter(CLOUD_WATER_1D[0: 30], MASSLEVELS_WRF_3D_WRFOUT[0: 30])
plt.title('CLW TLIQ 100% from sounding and T sounding')

plt.figure(2)
plt.plot(T_C_1D[0: 30], MASSLEVELS_WRF_3D_WRFOUT[0: 30])
plt.scatter(T_C_1D[0: 30], MASSLEVELS_WRF_3D_WRFOUT[0: 30])
plt.plot(T_DAGGPUNKT_1D[0: 30], MASSLEVELS_WRF_3D_WRFOUT[0: 30], linestyle = 'dotted')
plt.scatter(T_DAGGPUNKT_1D[0: 30], MASSLEVELS_WRF_3D_WRFOUT[0: 30])
plt.title('T and DP TLIQ 100% from sounding and T sounding')

plt.show()
