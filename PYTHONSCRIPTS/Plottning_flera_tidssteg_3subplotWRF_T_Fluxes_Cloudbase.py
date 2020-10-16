import numpy as N
import scipy.io.netcdf as S
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import os
import datetime


file_list = os.listdir('/nobackup/smhid13/fm/archive/output/wrf/2018/02/2018021800ECth1op')
#file_list = os.listdir('/nobackup/smhid12/sm_marha/2018-01-18_00z/run')
#file_list = os.listdir('/nobackup/smhid12/sm_marha/2018-02-17_18z_radt1/run/')

wrfout_list = []

for file in file_list:
    if file.startswith('wrfout_d01_2018-02-18'):
        wrfout_list.append(file)

wrfout_list.sort()




'''Skriver ut alla wrf-filer'''


for file in wrfout_list:
    print (file)

    


'''Listar diverse variabler'''        


VARIABEL_9 = 'SMOIS'

VARIABEL_10 = 'PBLH'

VARIABEL_11 = 'SNOWH'

VARIABEL_12 = 'SNOW'




wrfout_list_datum = []

TSK_TIDSSERIE = []

T2_TIDSSERIE = []

TEMP_C_TIDSSERIE = []

SNOWH_TIDSSERIE = []

SNOW_TIDSSERIE =[]

SOILT_1_TIDSSERIE = []

SOILT_2_TIDSSERIE = []

SOILT_3_TIDSSERIE = []

SOILT_4_TIDSSERIE = []

PBLH_TIDSSERIE = []
        

GLW_TIDSSERIE = []

SWDOWN_TIDSSERIE = []

GRDFLX_TIDSSERIE = []

HFX_TIDSSERIE = []

LH_TIDSSERIE = []

SWUP_TIDSSERIE = []

GUPLW_TIDSSERIE = []

CLOUDBASE_TIDSSERIE = []

NOAHRES_TIDSSERIE = []



TIDSVEKTOR_MINUTER = range(len(wrfout_list))
TIDSVEKTOR_DATUM = []


for file in wrfout_list:
        Filobjekt = S.netcdf_file('/home/sm_marha/WRF_JAN_FEB_KORNINGAR/ALLA_00_KORNINGAR/'+file, mode='r')
        print(file)
        #Filobjekt = S.netcdf_file('/nobackup/smhid12/sm_marha/2018-01-18_00z/run/'+file, mode='r')

        TSK = Filobjekt.variables['TSK'][0, 563, 352]-273.15
        T2 = Filobjekt.variables['T2'][0, 563, 352]-273.15

        
        
        SOILT_1 = Filobjekt.variables['TSLB'][0, 0, 563, 352]-273.15
        SOILT_2 = Filobjekt.variables['TSLB'][0, 1, 563, 352]-273.15
        SOILT_3 = Filobjekt.variables['TSLB'][0, 2, 563, 352]-273.15
        SOILT_4 = Filobjekt.variables['TSLB'][0, 3, 563, 352]-273.15

        
        
        THETA_PERT = Filobjekt.variables['T'][0, 0, 563, 352]
        P_PERT = Filobjekt.variables['P'][0, 0, 563, 352]
        P_BASE = Filobjekt.variables['PB'][0, 0, 563, 352]   
        
        P_TOT = P_BASE + P_PERT
        THETA = THETA_PERT+300
        TEMP_K = THETA/((1000/(P_TOT/100))**0.286)
        TEMP_C = TEMP_K-273.15

        

        SNOWH = Filobjekt.variables['SNOWH'][0, 563, 352]
        SNOW = Filobjekt.variables['SNOW'][0, 563, 352]
        PBLH = Filobjekt.variables['PBLH'][0, 563, 352]

        

        GLW = Filobjekt.variables['GLW'][0, 563, 352]
        SWDOWN = Filobjekt.variables['SWDOWN'][0, 563, 352]
        GRDFLX = Filobjekt.variables['GRDFLX'][0, 563, 352]
        HFX = Filobjekt.variables['HFX'][0, 563, 352]
        HFX = -HFX
        LH = Filobjekt.variables['LH'][0, 563, 352]
        EMISS = Filobjekt.variables['EMISS'][0, 563, 352]
        ALBEDO = Filobjekt.variables['ALBEDO'][0, 563, 352]
        GLW = EMISS*GLW                                      #How much from the downwelling long wave flux does the surface take up? 
        SWUP = ALBEDO*SWDOWN
        GUPLW = EMISS*5.67e-8*(TSK+273.15)**4
        NOAHRES = Filobjekt.variables['NOAHRES'][0, 563, 352]
        CLDFRA_1D = Filobjekt.variables['CLDFRA'][0, :, 563, 352]

        PH = Filobjekt.variables['PH'][:]
        PHB = Filobjekt.variables['PHB'][:]
        MODELLNIVAHOJD = (PH+PHB)/9.81
        MASSLEVELS = 0.5*(MODELLNIVAHOJD[0,:-1, 563, 352] + MODELLNIVAHOJD[0, 1:, 563, 352])


        for k in range(len(CLDFRA_1D)):
            if CLDFRA_1D[k] >=0.8:
                CLOUDBASE = MASSLEVELS[k]-MODELLNIVAHOJD[0,0,563,352]
                break
            else:
                CLOUDBASE = -999


        #CLOUD_ICE= Filobjekt.variables['QICE'][0, 563, 352]
        #CLOUD_WATER = Filobjekt.variables['QWATER'][0, 563, 352]

        
        
        TSK_TIDSSERIE.append(TSK)
        T2_TIDSSERIE.append(T2)
        SOILT_1_TIDSSERIE.append(SOILT_1)
        SOILT_2_TIDSSERIE.append(SOILT_2)
        SOILT_3_TIDSSERIE.append(SOILT_3)
        SOILT_4_TIDSSERIE.append(SOILT_4)
        TEMP_C_TIDSSERIE.append(TEMP_C)
        SNOWH_TIDSSERIE.append(SNOWH)
        SNOW_TIDSSERIE.append(SNOW)
        PBLH_TIDSSERIE.append(PBLH)
        GLW_TIDSSERIE.append(GLW)
        SWDOWN_TIDSSERIE.append(SWDOWN) 
        GRDFLX_TIDSSERIE.append(GRDFLX)
        HFX_TIDSSERIE.append(HFX)
        LH_TIDSSERIE.append(LH)
        SWUP_TIDSSERIE.append(SWUP)
        GUPLW_TIDSSERIE.append(GUPLW)
        CLOUDBASE_TIDSSERIE.append(CLOUDBASE)
        NOAHRES_TIDSSERIE.append(NOAHRES)

        RESIDUAL_TIDSSERIE = N.array(GLW_TIDSSERIE) + N.array(SWDOWN_TIDSSERIE) - N.array(SWUP_TIDSSERIE) - N.array(GUPLW_TIDSSERIE) + N.array(HFX_TIDSSERIE) - N.array(LH_TIDSSERIE) + N.array(GRDFLX_TIDSSERIE)
        LONGWAVE_RESIDUAL_TIDSSERIE = N.array(GLW_TIDSSERIE)-N.array(GUPLW_TIDSSERIE)

        
        wrfout_list_datum.append(file[11:])
        datum = datetime.datetime.strptime(str(file[11:]), '%Y-%m-%d_%H:%M:%S').strftime('%d/%m %H:%M')
        TIDSVEKTOR_DATUM.append(datum)

GUPLW_TIDSSERIE[0]=0
RESIDUAL_TIDSSERIE[0] = 0
LONGWAVE_RESIDUAL_TIDSSERIE[0] =0



'''Plottar data'''


fig1 = plt.figure(1, figsize=(10,7))

plt.suptitle(str(wrfout_list[1][11:21]) + ' 00z-run')


ax1 = plt.subplot(311)
ax1.plot(TIDSVEKTOR_MINUTER, TSK_TIDSSERIE, linewidth = 1, label = 'T-SKIN')
ax1.plot(TIDSVEKTOR_MINUTER, T2_TIDSSERIE, linewidth = 1, label = 'T2')
ax1.plot(TIDSVEKTOR_MINUTER, TEMP_C_TIDSSERIE, linewidth = 1, label = 'T 1st Level')

ax1.set_xticklabels([])
ax1.grid(linestyle="dotted")
ax1.legend(loc=1,prop={'size': 12})
plt.ylabel('T (C)')

'''
ax2 = plt.subplot(612)
ax2.plot(TIDSVEKTOR_MINUTER, SOILT_1_TIDSSERIE, linewidth = 1, label = 'T SOILLAYER 1')
ax2.plot(TIDSVEKTOR_MINUTER, SOILT_2_TIDSSERIE, linewidth = 1, label = 'T SOILLAYER 2')
ax2.plot(TIDSVEKTOR_MINUTER, SOILT_3_TIDSSERIE, linewidth = 1, label = 'T SOILLAYER 3')
ax2.plot(TIDSVEKTOR_MINUTER, SOILT_4_TIDSSERIE, linewidth = 1, label = 'T SOILLAYER 4')


ax2.set_xticklabels([])
ax2.grid(linestyle="dotted")
ax2.legend(loc=4, prop={'size': 8})
plt.ylabel('T (C)')


ax3 = plt.subplot(313)
ax3.plot(TIDSVEKTOR_MINUTER, PBLH_TIDSSERIE, linewidth = 1, label = 'PBLH')
#ax1.plot(TIDSVEKTOR_MINUTER, T2_TIDSSERIE, linewidth = 1, label = 'T2')
#ax1.plot(TIDSVEKTOR_MINUTER, TEMP_C_TIDSSERIE, linewidth = 1, label = 'T_1st Level')


ax3 = plt.subplot(312)

ax3.plot(TIDSVEKTOR_MINUTER, LONGWAVE_RESIDUAL_TIDSSERIE, linewidth = 1, label = 'LONGWAVE_RES')
ax3.plot(TIDSVEKTOR_MINUTER, RESIDUAL_TIDSSERIE, linewidth = 1, label = 'RESIDUAL')
ax3.plot(TIDSVEKTOR_MINUTER, GRDFLX_TIDSSERIE, linewidth = 1, label = 'GRDFLX')
ax3.plot(TIDSVEKTOR_MINUTER, HFX_TIDSSERIE, linewidth = 1, label = 'HFX')
ax3.plot(TIDSVEKTOR_MINUTER, LH_TIDSSERIE, linewidth = 1, label = 'LH')
ax3.plot(TIDSVEKTOR_MINUTER, NOAHRES_TIDSSERIE, linewidth = 1, label = 'NOAHRES')


ax3.grid(linestyle="dotted")
ax3.legend(loc=1, prop={'size': 10})
ax3.set_xticklabels([])
plt.ylabel('Flux (W/m2)')



ax4 = plt.subplot(614)

ax4.plot(TIDSVEKTOR_MINUTER, GLW_TIDSSERIE, linewidth = 1, label = 'GLW')
ax4.plot(TIDSVEKTOR_MINUTER, GUPLW_TIDSSERIE, linewidth = 1, label = 'GUPLW')




ax4.grid(linestyle="dotted")
ax4.legend(loc=4, prop={'size': 8})
ax4.set_xticklabels([])
plt.ylabel('Flux (W/m2)')


ax5 = plt.subplot(615)

ax5.plot(N.array(TIDSVEKTOR_MINUTER), RESIDUAL_TIDSSERIE, linewidth = 1, label = 'RESIDUAL')





ax5.grid(linestyle="dotted")
ax5.legend(loc=4, prop={'size': 8})
ax5.set_xticklabels([])
plt.ylabel('Flux (W/m2)')'''



ax6 = plt.subplot(313)

ax6.scatter(TIDSVEKTOR_MINUTER, CLOUDBASE_TIDSSERIE, label = 'CLOUDBASE')

ax6.grid(linestyle="dotted")
ax6.legend(loc=1, prop={'size': 12})
plt.xticks(TIDSVEKTOR_MINUTER, TIDSVEKTOR_DATUM  , rotation= 45, size=6)
plt.ylabel('Cloudbase (m)')
plt.ylim(0,max(CLOUDBASE_TIDSSERIE)+ 10)

#plt.plot(TIDSVEKTOR_MINUTER, VARIABEL_TIDSSERIE, linewidth=2)



#plt.xticks(TIDSVEKTOR_TIMMAR, TIDSVEKTOR_DATUM  , rotation= 45, size=6)

#plt.yticks(size=8)



#plt.tick_params(which = 'both', direction = 'in', pad = 10)

#plt.legend(VARIABEL)


#plt.title(VARIABEL + '  Sodankylä')

#fig1.savefig('/home/sm_marha/FIGURES/FOKUS_17_27_FEB/FLUXES_TEMPERATURE_RESIDUAL/' + wrfout_list[1][11:21] + '00z')

plt.show()
        
