import numpy as N
import scipy.io.netcdf as S
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import os
import datetime
import pandas as pd
import netCDF4 as N4

#KORNING = '2018-01-18_00z'



###########
#VARIABLES#
###########


Filobjekt = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-20_00z_91_LEVELS_42h/run/wrfinput_d01', mode='r')



TSK = Filobjekt.variables['TSK'][0, 563, 352]

P = Filobjekt.variables['P'][0,:, 563, 352]
PB = Filobjekt.variables['PB'][0, :, 563, 352]
P_TOT = P+PB

PH = Filobjekt.variables['PH'][0, :, 563, 352]
PHB =Filobjekt.variables['PHB'][0, :, 563, 352]

PERT_THETA = Filobjekt.variables['T'][0, :, 563, 352]
THETA = PERT_THETA + 300
T_K = THETA/((1000/(P_TOT/100))**0.286)
T_C = T_K - 273.15

QVAPOR = Filobjekt.variables['QVAPOR'][0, :, 563, 352]
QCLOUD = Filobjekt.variables['QCLOUD'][0, :, 563, 352]
QICE = Filobjekt.variables['QICE'][0, :, 563, 352]
#CLDFRA = Filobjekt.variables['CLDFRA'][0, :, 563, 352]

U = Filobjekt.variables['U'][0, :, 563, 352]
V = Filobjekt.variables['V'][0, :, 563, 352]
WIND_SPEED = N.sqrt(U**2+V**2)

XLONG = Filobjekt.variables['XLONG'][0, 563, 352]
XLAT = Filobjekt.variables['XLAT'][0, 563, 352]

PH = Filobjekt.variables['PH'][0, :, 563, 352]
PHB = Filobjekt.variables['PHB'][0, :, 563, 352]
Z_FULL = (PH + PHB)/9.81
Z_MASS = 0.5*(Z_FULL[:-1] + Z_FULL[1:])

E_MATTNAD = N.exp(N.log(611.2)+(17.62*T_C/(243.12+T_C)))
QVAPOR_MATTNAD = 0.622*E_MATTNAD/(P_TOT-E_MATTNAD)
RH = QVAPOR/QVAPOR_MATTNAD
SNOWH = Filobjekt.variables['SNOWH'][0, 563, 352]
SNOW = Filobjekt.variables['SNOW'][0, 563, 352]
SNOALB = Filobjekt.variables['SNOALB'][0, 563, 352]


HGT = Filobjekt.variables['HGT'][0, 563, 352]
U_TEN = Filobjekt.variables['U10'][0, 563, 352]
V_TEN = Filobjekt.variables['V10'][0, 563, 352]
T_TVA = Filobjekt.variables['T2'][0, 563, 352]
Q_TVA = Filobjekt.variables['Q2'][0, 563, 352]
P_SURFACE = Filobjekt.variables['PSFC'][0, 563, 352]

Filobjekt.close()




###############################
#MATRIX OF UPPER AIR VARIABLES#
###############################

MATRIX = N.column_stack((Z_MASS, U, V, THETA, QVAPOR))

#FIRST_ROW_ARRAY = N.array([HGT, U_TEN, V_TEN, T_TVA, Q_TVA, P_SURFACE])




##################################
#ROUNDING VARIABLES AT 2M AND 10M#
##################################

float_formatter = lambda x: "%8.2f" % x
HGT = float_formatter(HGT)

float_formatter = lambda x: "%4.1f" % x
U_TEN = float_formatter(U_TEN)

float_formatter = lambda x: "%4.1f" % x
V_TEN = float_formatter(V_TEN)

float_formatter = lambda x: "%6.2f" % x
T_TVA = float_formatter(T_TVA)

float_formatter = lambda x: "%8.5e" % x
Q_TVA = float_formatter(Q_TVA)

float_formatter = lambda x: "%6.1f" % x
P_SURFACE = float_formatter(P_SURFACE)




##################################
#DATAFRAME OF UPPER AIR VARIABLES#
##################################

df = pd.DataFrame(MATRIX)
#df.loc[-1] = FIRST_ROW_ARRAY
#df.index = df.index + 1
#df = df.sort_index()




###########################################
#ROUNDING OF DATAFRAME UPPER AIR VARIABLES#
###########################################

df.iloc[:, 0] = df.iloc[:, 0].map('{:8.2f}'.format)
df.iloc[:, 1] = df.iloc[:, 1].map('{:4.1f}'.format)
df.iloc[:, 2] = df.iloc[:, 2].map('{:4.1f}'.format)
df.iloc[:, 3] = df.iloc[:, 3].map('{:6.2f}'.format)
df.iloc[:, 4] = df.iloc[:, 4].map('{:8.5e}'.format)





##################
#WRITING TEXTFILE#
##################

df.to_csv('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY/test/em_scm_xy/N_input_sounding.txt', sep='\t', index=False, header=False)






#####################################################
#OPENING AND READING TEXTFILE AND STORE IT IN "data"#
#####################################################

textfile = open('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY/test/em_scm_xy/N_input_sounding.txt', 'r')
data = textfile.read()
textfile.close()




#####################################################
#OPENING NEW TEXTFILE AND WRITE 2M and 10M VARIABLES#
#####################################################

textfile = open('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY/test/em_scm_xy/2018_02_20_00z_input_sounding_NEW.txt', 'w')
textfile.write(str(HGT)+'\t'+str(U_TEN)+'\t'+str(V_TEN)+'\t'+str(T_TVA)+'\t'+str(Q_TVA)+'\t'+str(P_SURFACE)+'\n')



#############################################
#THEN WRITNG BACK VARIABLES STORED IN "data"#
#############################################

textfile.write(data)

textfile.close()



############
#INPUT SOIL#
############

Filobjekt_TVA = N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-20_00z_91_LEVELS_42h/run/wrfinput_d01', mode='r')

ZS = Filobjekt_TVA.variables['ZS'][0, :]
TMN = Filobjekt_TVA.variables['TMN'][:, 563, 352]           #Soil temperature at lower boundary
TSK = Filobjekt_TVA.variables['TSK'][0, 563, 352]

TSLB = Filobjekt_TVA.variables['TSLB'][0, :, 563, 352]      #Soil temperature
SMOIS = Filobjekt_TVA.variables['SMOIS'][0, :, 563, 352]    #Soil moisture

ARRAY = N.array([0.0, TMN, TSK])

MATRIX_2 = N.column_stack(( ZS, TSLB, SMOIS))

MATRIX_2 = N.vstack([ARRAY, MATRIX_2])

df_2 = pd.DataFrame(MATRIX_2)



df_2.iloc[:, 0] = df_2.iloc[:, 0].map('{:4.2f}'.format)
df_2.iloc[:, 1] = df_2.iloc[:, 1].map('{:6.2f}'.format)
df_2.iloc[:, 2] = df_2.iloc[:, 2].map('{:6.2f}'.format)



df_2.to_csv('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY/test/em_scm_xy/2018_02_20_00z_input_soil_NEW.txt', sep='\t', index=False, header=False)


print('SNOW ALBEDO= ' ,SNOALB)
print('SNOW WATER EQUIVALENT= ' ,SNOW)
print('SNOW DEPTH= ' ,SNOWH)





