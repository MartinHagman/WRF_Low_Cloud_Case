import numpy as N
import scipy.io.netcdf as S
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import netCDF4 as N4


#18 FEB

Korningar = ['/nobackup/smhid12/sm_marha/2018-02-18_00z_91_LEVELS_42h/run/wrfout_d01_2018-02-18_03:00:00',
             '/nobackup/smhid12/sm_marha/2018-02-18_00z_91_LEVELS_Icloud_2/run/wrfout_d01_2018-02-18_03:00:00',
             '/nobackup/smhid12/sm_marha/2018-02-18_00z_91_LEVELS_QC_QI_T_to_Td_0_Whole_Dom/run/wrfout_d01_2018-02-18_03:00:00',
             '/nobackup/smhid12/sm_marha/2018-02-18_00z_91_LEVELS_QC_QI_T_to_Td_0_Whole_Dom_Icloud2/run/wrfout_d01_2018-02-18_03:00:00']


#26 FEB
'''
Korningar = ['/nobackup/smhid12/sm_marha/2018-02-26_00z_91_LEVELS_REFERENCE/run/wrfout_d01_2018-02-26_06:00:00',
             '/nobackup/smhid12/sm_marha/2018-02-26_00z_91_LEVELS_Icloud2/run/wrfout_d01_2018-02-26_06:00:00',
             '/nobackup/smhid12/sm_marha/2018-02-26_00z_91_LEVELS_QC_QI_24h_T_to_Td_Whole_domain_0_icloud1/run/wrfout_d01_2018-02-26_06:00:00',
             '/nobackup/smhid12/sm_marha/2018-02-26_00z_91_LEVELS_QC_QI_24h_T_to_Td_Whole_domain_0_icloud2/run/wrfout_d01_2018-02-26_06:00:00']'''


# 20 FEB
'''
Korningar = ['/nobackup/smhid12/sm_marha/2018-02-20_00z_91_LEVELS_15h/run/wrfout_d01_2018-02-20_06:00:00',
             '/nobackup/smhid12/sm_marha/2018-02-20_00z_91_LEVELS_QC_QI_42h/run/wrfout_d01_2018-02-20_06:00:00',
             '/nobackup/smhid12/sm_marha/2018-02-20_00z_91_LEVELS_QC_QI_Icloud_2/run/wrfout_d01_2018-02-20_06:00:00',
             '/nobackup/smhid12/sm_marha/2018-02-20_00z_91_LEVELS_QC_QI_18h_T-to_Td_0_Icloud2/run/wrfout_d01_2018-02-20_06:00:00']'''


# 27 FEB

'''Korningar = ['/nobackup/smhid12/sm_marha/2018-02-27_00z_91_LEVELS/run/wrfout_d01_2018-02-27_21:00:00',
             '/nobackup/smhid12/sm_marha/2018-02-27_00z_91_LEVELS_Icloud2/run/wrfout_d01_2018-02-27_21:00:00',
             '/nobackup/smhid12/sm_marha/2018-02-27_00z_91_LEVELS_Clouds_Where_Clouds_In_Sounding_YSU_1_100/run/wrfout_d01_2018-02-27_21:00:00',
             '/nobackup/smhid12/sm_marha/2018-02-27_00z_91_LEVELS_Clouds_Where_Clouds_In_Sounding_YSU_TOPD1_CORRECT_CLW_100_dryer_below250/run/wrfout_d01_2018-02-27_21:00:00']'''



fig = plt.figure(1, figsize = (11, 4))

fig.subplots_adjust(right=0.90)
fig.subplots_adjust(left=0.01)

#cbar_ax = fig.add_axes([0.92, 0.16, 0.01, 0.67])

cbar_ax = fig.add_axes([0.92, 0.11, 0.01, 0.75])


Panel_letter_list = ['a', 'b', 'c', 'd']


Fig_No = 1

b = 0

for korning in Korningar:


    Filobjekt = N4.Dataset(korning, mode='r')

    VARIABEL = 'CLDFRA'

    VARIABELVARDE = Filobjekt.variables[VARIABEL][:]
    VARIABELVARDE_3D = VARIABELVARDE[0, :, :, :]
    P_PERT = Filobjekt.variables['P'][:]
    P_BASE = Filobjekt.variables['PB'][:]
    P = P_BASE+P_PERT
    TEMP = VARIABELVARDE/((1000/(P/100))**0.286)
    TEMP = TEMP-273.15
    PH = Filobjekt.variables['PH'][:]
    PHB = Filobjekt.variables['PHB'][:]
    MODELLNIVAHOJD = (PH+PHB)/9.81
    MASSLEVELS = 0.5*(MODELLNIVAHOJD[0,:-1, :, :] + MODELLNIVAHOJD[0, 1:, :, :])
    LON = Filobjekt.variables['XLONG'][0, :, :]
    LAT = Filobjekt.variables['XLAT'][0, :, :]
    lon_units = Filobjekt.variables['XLONG'].units
    lat_units = Filobjekt.variables['XLAT'].units




    CLOUDBASE = N.zeros(N.shape(LAT))

    for j in range(N.shape(VARIABELVARDE_3D)[1]):
        for i in range(N.shape(VARIABELVARDE_3D)[2]):
            for k in range(N.shape(VARIABELVARDE_3D)[0]):
                if VARIABELVARDE_3D[k,j,i] >=0.8:
                    CLOUDBASE[j, i] = MASSLEVELS[k,j,i]-MODELLNIVAHOJD[0,0,j,i]
                    break
                else:
                    CLOUDBASE[j, i] = -999    
                





    #plt.figure(1, figsize = (5,10))

    plt.subplot(1, 4, Fig_No)

    mapproj = bm.Basemap(width=1284000, height=2103000, resolution='h', rsphere=(6378137.00,6356752.3142),projection='lcc',lat_1=55.0,
                     lat_2=70.0, lat_0=62.5, lon_0=15.0)

    mapproj.drawcoastlines(linewidth=0.2)
    lonproj, latproj = mapproj(LON, LAT)
    levels = [0, 90, 180, 300, 450, 750, 1000, 1500, 1800, 5100, 9000, 15000, 18000]


    #mymap = plt.contour(lonproj, latproj, CLOUDBASE,linewidths=0.5, colors='k')
    #mymap = plt.contour(lonproj, latproj, LON,linewidths=0.5, colors='k')
    #mymap = plt.contour(lonproj, latproj, LAT,linewidths=0.5, colors='k')
    mymapf = plt.contourf(lonproj, latproj, CLOUDBASE, levels, colors=('r', '#A52A2A', '#D2691E', '#FFFF00','#98FB98','#00FF00', '#228B22', 'g', '#006400', '#0000CD', 'b', '#4169E1', '#B0C4DE' ))#cmap=plt.cm.RdYlBu)
    #plot.contour(lonproj, latproj, temperature_2D_C, 10, colors='k')

    plt.plot(lonproj[563, 352], latproj[563, 352], 'k^')
    plt.text(lonproj[630, 0], latproj[630, 0], Panel_letter_list[b], fontsize = 18, fontweight='bold')#, backgroundcolor = 'w')

    #plt.clabel(mymap, fontsize=3)
    #plt.axis([0, 360, -90, 90])
    #plt.suptitle('Cloud base at 0600 UTC 20 February 2018')
    '''if Fig_No == 3:
        plt.title('TEMF', fontweight="bold", fontsize=16)
    else:
        plt.title(korning[42:45], fontweight="bold", fontsize=16)'''

    
    #plt.xlabel('Longitude [' + str(lon_units) + ']')
    #plt.ylabel('Latitude [' + str(lat_units) + ']')
    
    Fig_No +=1
    b+=1
    
   

    #plt.colorbar(mymapf, orientation='vertical')
cbar = fig.colorbar(mymapf, cax=cbar_ax)
cbar.set_ticks([0, 90, 180, 300, 450, 750, 1000, 1500, 1800, 5100, 9000, 15000, 18000])
cbar.set_ticklabels([0, 90, 180, 300, 450, 750, 1000, 1500, 1800, 5100, 9000, 15000, 18000])
cbar.ax.tick_params(labelsize=14)
cbar.ax.set_title('(m)', fontsize=14)

plt.savefig('/home/sm_marha/FIGURES/FINAL/20180218/CLOUDBASE_BEFORE_AND_AFTER_20180218_00_+_03')
plt.show()

                    






