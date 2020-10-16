import numpy as N
import scipy.io.netcdf as S
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import os
import datetime
import netCDF4 as N4


Filobjekt = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-27_00z_T_to_Td_0_Icloud2/run/wrfinput_d01')
#Filobjekt = N4.Dataset('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfinput_d01_0218_QC_QI_REFERENCE')

LAT = 500

LON = 250


CLOUD_WATER_1D = Filobjekt.variables['QCLOUD'][0, :, LAT, LON]*1000

THETA_PERT_1D = Filobjekt.variables['T'][0, :, LAT, LON]
P_PERT_1D = Filobjekt.variables['P'][0, :, LAT, LON]
P_BASE_1D = Filobjekt.variables['PB'][0, :, LAT, LON]   

P_TOT_1D = P_BASE_1D + P_PERT_1D
THETA_1D = THETA_PERT_1D+300
TEMP_K_1D = THETA_1D/((1000/(P_TOT_1D/100))**0.286)
T_C_1D = TEMP_K_1D-273.15

QVAPOR_1D = Filobjekt.variables['QVAPOR'][0, :, LAT, LON]
E_MATTNAD_1D = N.exp(N.log(611.2)+(17.62*T_C_1D/(243.12+T_C_1D)))
QVAPOR_MATTNAD_1D = 0.622*E_MATTNAD_1D/(P_TOT_1D-E_MATTNAD_1D)
RH_1D = QVAPOR_1D/QVAPOR_MATTNAD_1D
T_DAGGPUNKT_1D = (243.12*N.log(611.2)-243.12*N.log(RH_1D*E_MATTNAD_1D)) / (N.log(RH_1D*E_MATTNAD_1D)-17.62-N.log(611.2))


PH_WRF_3D_WRFOUT = Filobjekt.variables['PH'][0, :, 563, 352]
PHB_WRF_3D_WRFOUT = Filobjekt.variables['PHB'][0, :, 563, 352]
MODELLNIVAHOJD_WRF_3D_WRFOUT = (PH_WRF_3D_WRFOUT+PHB_WRF_3D_WRFOUT)/9.81
MODELLNIVAHOJD_TER_WRF_3D_WRFOUT = MODELLNIVAHOJD_WRF_3D_WRFOUT[:] -MODELLNIVAHOJD_WRF_3D_WRFOUT[0]




MASSLEVELS_WRF_3D_WRFOUT = 0.5*(MODELLNIVAHOJD_TER_WRF_3D_WRFOUT[:-1] + MODELLNIVAHOJD_TER_WRF_3D_WRFOUT[1:])

plt.figure(1)

ax1 = plt.subplot(111)
ax1.plot(CLOUD_WATER_1D[0: 30], MASSLEVELS_WRF_3D_WRFOUT[0: 30])
ax1.grid(linestyle="dotted")

plt.figure(2)
ax1 = plt.subplot(111)
ax1.plot(T_C_1D[0: 30], MASSLEVELS_WRF_3D_WRFOUT[0: 30])
#ax1.scatter(T_C_1D[0: 30], MASSLEVELS_WRF_3D_WRFOUT[0: 30])

ax1.plot(T_DAGGPUNKT_1D[0: 30], MASSLEVELS_WRF_3D_WRFOUT[0: 30], linestyle = 'dotted')
ax1.grid(linestyle="dotted")

plt.show()
