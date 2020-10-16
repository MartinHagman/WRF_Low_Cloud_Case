import numpy as N
import scipy.io.netcdf as S
import matplotlib as mpl
import matplotlib.pyplot as plt





Filobjekt = S.netcdf_file('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM/test/em_scm_xy/wrfinput_d01')

P_PERT = Filobjekt.variables['P'][0, :, 1, 1]
P_BASE = Filobjekt.variables['PB'][0, :, 1, 1]
P = P_BASE+P_PERT
PSFC = Filobjekt.variables['PSFC'][0, 1, 1]

PH = Filobjekt.variables['PH'][0, :, 1, 1]
PHB = Filobjekt.variables['PHB'][0, :, 1, 1]
MODELLNIVAHOJD = (PH+PHB)/9.81
#MODELLNIVAHOJD_TER = MODELLNIVAHOJD -MODELLNIVAHOJD[0]


MASSLEVELS = 0.5*(MODELLNIVAHOJD[:-1] + MODELLNIVAHOJD[1:])
