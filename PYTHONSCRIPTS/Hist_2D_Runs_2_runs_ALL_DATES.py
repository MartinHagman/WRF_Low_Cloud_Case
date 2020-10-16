import numpy as N
import scipy.io.netcdf as S
import netCDF4 as N4
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import os
import pandas as pd
from math import *
import sys
from matplotlib.colors import LogNorm
from pylab import *



rc('axes', linewidth=1.5) 

fig1, axs = plt.subplots(1, 4, figsize=(11,3))

fig1.subplots_adjust(bottom=0.2)
fig1.subplots_adjust(wspace=0.1)
fig1.subplots_adjust(right=0.92)
fig1.subplots_adjust(left=0.07)


Figure_letter_list = ['a', 'b', 'c', 'd']


#18 Feb

Dirs_1 = [N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-18_00z_91_LEVELS_42h/run/wrfout_d01_2018-02-18_06:00:00', mode='r'),
        N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-18_00z_91_LEVELS_QC_QI_T_to_Td_1e-4_Icloud_2/run/wrfout_d01_2018-02-18_06:00:00', mode='r')]


#20 feb

Dirs_2 = [N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-20_00z_91_LEVELS_15h/run/wrfout_d01_2018-02-20_06:00:00', mode='r'),
        N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-20_00z_91_LEVELS_QC_QI_18h_T-to_Td_0_Icloud2/run/wrfout_d01_2018-02-20_06:00:00', mode='r')]


#26 feb

Dirs_3 = [N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-26_00z_91_LEVELS_REFERENCE/run/wrfout_d01_2018-02-26_06:00:00', mode='r'),
        N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-26_00z_91_LEVELS_QC_QI_24h_T_to_Td_Whole_domain_0_icloud2/run/wrfout_d01_2018-02-26_06:00:00', mode='r')]


#27 feb

Dirs_4 = [N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-27_00z_91_LEVELS/run/wrfout_d01_2018-02-27_06:00:00', mode='r'),
        N4.Dataset('/nobackup/smhid12/sm_marha/2018-02-27_00z_T_to_Td_1e-8_lower_0_Icloud2/run/wrfout_d01_2018-02-27_06:00:00', mode='r')]



DIRECTORIES = [Dirs_1, Dirs_2, Dirs_3, Dirs_4]



b = 0
for ax, dirs in zip(axs, DIRECTORIES):
    Filobjekt_Ref = dirs[0]
    Filobjekt_Ass = dirs[1]

    ################
    # REFERENCE WRF#
    ############################################################################


    VARIABEL = 'CLDFRA'

    VARIABELVARDE_REF = Filobjekt_Ref.variables[VARIABEL][:]
    VARIABELVARDE_3D_REF = VARIABELVARDE_REF[0, :, :, :]
    P_PERT_REF = Filobjekt_Ref.variables['P'][:]
    P_BASE_REF = Filobjekt_Ref.variables['PB'][:]
    P_REF = P_BASE_REF+P_PERT_REF
    TEMP_REF = VARIABELVARDE_REF/((1000/(P_REF/100))**0.286)
    TEMP_REF = TEMP_REF-273.15
    PH_REF = Filobjekt_Ref.variables['PH'][:]
    PHB_REF = Filobjekt_Ref.variables['PHB'][:]
    MODELLNIVAHOJD_REF = (PH_REF+PHB_REF)/9.81
    MASSLEVELS_REF = 0.5*(MODELLNIVAHOJD_REF[0,:-1, :, :] + MODELLNIVAHOJD_REF[0, 1:, :, :])
    LON_REF = Filobjekt_Ref.variables['XLONG'][0, :, :]
    LAT_REF = Filobjekt_Ref.variables['XLAT'][0, :, :]
    lon_units = Filobjekt_Ref.variables['XLONG'].units
    lat_units = Filobjekt_Ref.variables['XLAT'].units
    LANDMASK_REF = Filobjekt_Ref.variables['LANDMASK'][0, :, :]




    CLOUDBASE_REF = N.zeros(N.shape(LAT_REF))
    #CLOUDBASE_REF = CLOUDBASE_REF.astype(float)

    for j in range(500, N.shape(VARIABELVARDE_3D_REF)[1]):
        for i in range(250, N.shape(VARIABELVARDE_3D_REF)[2]):
            if LANDMASK_REF[j, i] == 1:
                for k in range(0, 30):
                #for k in range(N.shape(VARIABELVARDE_3D_REF)[0]):
                    if VARIABELVARDE_3D_REF[k,j,i] >=0.8:
                        CLOUDBASE_REF[j, i] = MASSLEVELS_REF[k,j,i]-MODELLNIVAHOJD_REF[0,0,j,i]
                        break
                    else:
                        CLOUDBASE_REF[j, i] = -200


            else:
                 CLOUDBASE_REF[j, i] = N.nan


    CLOUDBASE_REF = CLOUDBASE_REF[500:, 250:]
    CLOUDBASE_REF = CLOUDBASE_REF.reshape(1, -1)
    CLOUDBASE_REF = CLOUDBASE_REF[0, :]
    CLOUDBASE_REF = CLOUDBASE_REF[N.logical_not(N.isnan(CLOUDBASE_REF))]
    print(len(CLOUDBASE_REF))


    ###################
    # ASSIMILATION WRF#
    ######################################################################


    VARIABEL = 'CLDFRA'

    VARIABELVARDE_ASS = Filobjekt_Ass.variables[VARIABEL][:]
    VARIABELVARDE_3D_ASS = VARIABELVARDE_ASS[0, :, :, :]
    P_PERT_ASS = Filobjekt_Ass.variables['P'][:]
    P_BASE_ASS = Filobjekt_Ass.variables['PB'][:]
    P_ASS = P_BASE_ASS+P_PERT_ASS
    TEMP_ASS = VARIABELVARDE_ASS/((1000/(P_ASS/100))**0.286)
    TEMP_ASS = TEMP_ASS-273.15
    PH_ASS = Filobjekt_Ass.variables['PH'][:]
    PHB_ASS = Filobjekt_Ass.variables['PHB'][:]
    MODELLNIVAHOJD_ASS = (PH_ASS+PHB_ASS)/9.81
    MASSLEVELS_ASS = 0.5*(MODELLNIVAHOJD_ASS[0,:-1, :, :] + MODELLNIVAHOJD_ASS[0, 1:, :, :])
    LON_ASS = Filobjekt_Ass.variables['XLONG'][0, :, :]
    LAT_ASS = Filobjekt_Ass.variables['XLAT'][0, :, :]
    lon_units = Filobjekt_Ass.variables['XLONG'].units
    lat_units = Filobjekt_Ass.variables['XLAT'].units
    LANDMASK_ASS = Filobjekt_Ass.variables['LANDMASK'][0, :, :]




    CLOUDBASE_ASS = N.zeros(N.shape(LAT_ASS))
    #CLOUDBASE_ASS = CLOUDBASE_ASS.astype(float)


    for j in range(500, N.shape(VARIABELVARDE_3D_ASS)[1]):
        for i in range(250, N.shape(VARIABELVARDE_3D_ASS)[2]):
            if LANDMASK_ASS[j, i] == 1:
                for k in range(0,30):
                #for k in range(N.shape(VARIABELVARDE_3D_ASS)[0]):
                    if VARIABELVARDE_3D_ASS[k,j,i] >=0.8:
                        CLOUDBASE_ASS[j, i] = MASSLEVELS_ASS[k,j,i]-MODELLNIVAHOJD_ASS[0,0,j,i]
                        break
                    else:
                        CLOUDBASE_ASS[j, i] = -200
            else:
                 CLOUDBASE_ASS[j, i] = N.nan
                

    CLOUDBASE_ASS = CLOUDBASE_ASS[500:, 250:]
    CLOUDBASE_ASS = CLOUDBASE_ASS.reshape(1, -1)
    CLOUDBASE_ASS = CLOUDBASE_ASS[0, :]
    print(CLOUDBASE_ASS)
    print(len(CLOUDBASE_ASS))
    CLOUDBASE_ASS = CLOUDBASE_ASS[N.logical_not(N.isnan(CLOUDBASE_ASS))]
    print(len(CLOUDBASE_ASS))
    print(CLOUDBASE_ASS)


    ###########
    # PLOTTING#
    ############################################################################
    bins = (20, 20)              

    if b == 3:
        bins = (17,17)
        

    hsist = ax.hist2d(CLOUDBASE_REF, CLOUDBASE_ASS, bins=bins ,norm = LogNorm(), cmap = 'PuRd')

    #cbar = plt.colorbar()

    '''if b == 3:
        cbar = fig1.colorbar(hsist[3], ax = ax)'''

    X = N.arange(-15,10000)
    Y = N.arange(-15,10000)
    X_0 = N.zeros(10015)
    Y_0 = N.zeros(10015)



    ax.plot(X, Y , 'r', linewidth = 2)
    ax.plot(X_0-15, Y , 'k', linewidth = 2)
    ax.plot(X, Y_0-15 , 'k', linewidth = 2)


    #ax.set_title('Cloudbase Ref vs Ass Sodankylä Jan 18 2018 1e-4 Icloud2', fontsize = 16)
    ax.tick_params(which='major', direction='in', width =2, labelsize = '12')
    ax.set_xlabel('Reference (m)', fontsize = 15)

    if b == 0:
        ax.set_ylabel('Assimilation (m)', fontsize = 15)

    if b != 0:
        ax.set_yticklabels([])
        
    ax.set_xlim(-300, 1700)
    ax.set_ylim(-300, 1700)
    ax.text(-250, 1550, Figure_letter_list[b], fontsize = 12, fontweight='bold')
    #ax.set_xticks([0, 250, 500, 750, 1000, 1250, 1500])
    #ax.set_yticks([0, 250, 500, 750, 1000, 1250, 1500])
    #ax.set_xticklabels([0, 250, 500, 750, 1000, 1250, 1500])
    #ax.set_yticklabels([0, 250, 500, 750, 1000, 1250, 1500])
    #ax.set_tick_params(labelsize=10)

    '''if b ==3:
        cbar.set_label('Number of hits', fontsize = 10)
        cbar.ax.tick_params(labelsize=10)'''

    b+=1

cbar_ax = fig1.add_axes([0.94, 0.19, 0.01, 0.69])       #Gemensam bar
cbar = fig1.colorbar(hsist[3], cax=cbar_ax)
cbar.ax.tick_params(labelsize=15)
    
plt.savefig('/home/sm_marha/FIGURES/FINAL/20180218/Run_2_Run_ALL_DATES_30_vertical_levels_northern_parts')

plt.show()
