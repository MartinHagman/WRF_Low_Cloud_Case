import numpy as N
import scipy.io.netcdf as S
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm

'''Öppnar filobjektet'''

Filobjekt = S.netcdf_file('/nobackup/smhid12/sm_marha/2018-01-18_00z/run/wrfout_d01_2018-01-18_00:00:00', mode='r')



'''Väljer variabel'''

VARIABEL = 'T'



'''Skriver ut alla variablerna samt dimensionerna'''

#print(Filobjekt.variables)
#print(Filobjekt.dimensions)



'''Beräknar diverse parametrar'''

VARIABELVARDE = Filobjekt.variables[VARIABEL][:]
ETANIVAER_MASS = Filobjekt.variables['ZNU'][:]
ETANIVAER_FULL = Filobjekt.variables['ZNW'][:]
ZETATOP = Filobjekt.variables['ZETATOP'][:]
HGT = Filobjekt.variables['HGT'][:]
PH = Filobjekt.variables['PH'][:]
PHB = Filobjekt.variables['PHB'][:]
GEOPOTENTIAL = PH + PHB
MODELLNIVAHOJD = (PH+PHB)/9.81
W = Filobjekt.variables['W'][:]
'''PERT_THETA är den variabeln som i detta fall är vald längst upp. I många andra fall behövs inga omräkningar av variabeln utan
   VARIABELVARDE är den viktiga variabeln genom hela skriptet.'''
PERT_THETA = Filobjekt.variables['T'][:]
THETA = PERT_THETA+300
P_PERT = Filobjekt.variables['P'][:]
P_BASE = Filobjekt.variables['PB'][:]
P = P_BASE+P_PERT
TEMP = THETA/((1000/(P/100))**0.286)
TEMP_C = TEMP-273.15



'''Skriver ut matriserna för diverse parametrar'''

BF =  Filobjekt.variables['BF'][:]
C1H = Filobjekt.variables['C1H'][:]
C2H = Filobjekt.variables['C2H'][:]
BH = Filobjekt.variables['BH'][:]
C1F = Filobjekt.variables['C1F'][:]
C2F = Filobjekt.variables['C2F'][:]
PSFC = Filobjekt.variables['PSFC'][:]



'''Skriver ut shapen för olika variabler'''


print(VARIABEL,'= ', N.shape(VARIABELVARDE))
print('ZNU(MASS) = ', N.shape(ETANIVAER_MASS))
print('ZNU(FULL) = ',N.shape(ETANIVAER_FULL))
print('ZETATOP = ',N.shape(ZETATOP))
print('HGT = ', N.shape(HGT))
print('PH = ', N.shape(PH))
print('PHB = ', N.shape(PHB))
print('GEOPOTENTIAL = ', N.shape(GEOPOTENTIAL))
print('MODELLNIVAHOJD = ', N.shape(MODELLNIVAHOJD))
print('PSFC = ', N.shape(PSFC))
print('W = ', N.shape(W))
print('TEMP = ', N.shape(TEMP))
print('TEMP_C = ', N.shape(TEMP_C))
print('THETA = ', N.shape(THETA))




'''Läser ut longitud, latitud och temperatur ur NetCDF-filen. Deklarerar även variabler
för enheter.'''

lon = Filobjekt.variables['XLONG'][0, :, :]
lon_units = Filobjekt.variables['XLONG'].units
lat = Filobjekt.variables['XLAT'][0, :, :]
lat_units = Filobjekt.variables['XLAT'].units
VARIABELVARDE_3D = VARIABELVARDE[0,:, :, :]
'''Omräkning av Kelvin till Celsius'''
#VARIABELVARDE_2D_C = VARIABELVARDE_2D-273.15

TEMP_C_3D = TEMP_C[0, :, :, :]



'''Dubbelkollar shapen på VARIABELVARDEt'''

print('VARIABEL ', N.shape(VARIABELVARDE_3D))




'''Färdiga rader för att skriva ut olika matriser'''

#print('VARIABEL ', '= ',  VARIABELVARDE_3D)
#print('ZNU(MASS) = ', ETANIVAER_MASS)
#print('ZNU(FULL) = ', ETANIVAER_FULL)
#print('ZETATOP = ', ZETATOP)
#print('HGT = ', HGT)
#print('PH = ', PH)
#print('PHB = ', PHB)
#print('GEOPOTENTIAL = ', GEOPOTENTIAL)
#print('MODELLNIVAHOJD = ', MODELLNIVAHOJD)
#print('PSFC = ', PSFC)
#print('W = ', W)
#print('TEMP = ', TEMP)
#print('THETA = ', THETA)
#print('P = ', P)
#print('TEMP_C = ', TEMP_C)

#print('BF = ', BF)
#print('C1H = ', C1H)
#print('C2H = ', C2H)
#print('BH = ', BH)
#print('C1F = ', C1F)
#print('C2F = ', C2F)




'''I dett fall räknar vi på TEMP_C_3D, men normalt mär inga omräkningar görs kommer VARIABELVARDEt att stå i stället för TEMP_C_3D
Vi skapar endimensionella arrayer i de omgivande punkterna till Sodankylä'''

VARIABELPROFILER = N.array([[TEMP_C_3D[:, 564, 351]],
                             [TEMP_C_3D[:, 564, 352]],
                             [TEMP_C_3D[:, 564, 353]],
                             [TEMP_C_3D[:, 563, 351]],
                             [TEMP_C_3D[:, 563, 353]],
                             [TEMP_C_3D[:, 562, 351]],
                             [TEMP_C_3D[:, 562, 352]],
                             [TEMP_C_3D[:, 562, 353]],
                             [TEMP_C_3D[:, 563, 352]]])


'''Skapar ny array som är tom. samt sätter iterationsvariabeln till 0 för den for-slinga som följer. Därefter skapas två listor
   med x- och y-koordinaterna för de olika punkterna runt Sodankylä. Modellnivåhöjden som tidigare räknats ut från PH och PHB,
   och där nedersta värdet är höjden på terrängen görs om till höjden till massnivåerna. På de två första raderna i loopen
   skapas 1-dimensionella höjdvektorer för massnivåerna och terränghöjden dras bort. Då blir den lägsta nivån 0 m. Efter detta
   skapas matrisen MASSLEVELS, där vi kommer använda i while-satsen. Matrisen ska vara ett element kortare än den ursprungliga
   modellnivåhöjdvektorn, då ju massnivåerna är en färre. I whilesatsen skapas sedan MASSLEVELS för en punkt tills den är full.
   Sedan lagras den i ytterligare en matris, MASSLEVELSARRAY. OBS att len() alltid är en större än sista indexet i matrisen.
   Därav len(MODELLNIVAHOJD_1D)-1. När while-satsen är klar startar nästa steg i for-loopen.När alla steg är klara så splitrar
   jag upp MASSREVELSARRAY I MASSLEVELSARRAY_RESHAPED så att jag får alla punkterna i olika vektorer'''

MASSLEVELSARRAY = N.array([])

i = 0

for i in range(9):

    Y = [564, 564, 564, 563, 563, 562, 562, 562, 563]
    X = [351, 352, 353, 351, 353, 351, 352, 353, 352]

#    VARIABELVARDE_1D = VARIABELVARDE_3D[:, Y[i], X[i]] Denna behövs ej just nu då denna bara innehåller potentiell temperatur.
    MODELLNIVAHOJD_1D = MODELLNIVAHOJD[0, :, Y[i], X[i]]
    MODELLNIVAHOJD_1D_MINUS_TER = MODELLNIVAHOJD_1D-MODELLNIVAHOJD_1D[0]


    MASSLEVELS = N.zeros(len(MODELLNIVAHOJD_1D)-1)
    j = 0
    
    while j<(len(MODELLNIVAHOJD_1D_MINUS_TER)-1):
        MASSLEVELS[j] = (MODELLNIVAHOJD_1D_MINUS_TER[j]+ MODELLNIVAHOJD_1D_MINUS_TER[j+1])/2        
        j+=1
        

    MASSLEVELSARRAY = N.concatenate((MASSLEVELSARRAY, MASSLEVELS))

    
MASSLEVELSARRAY_RESHAPED = N.reshape(MASSLEVELSARRAY,(9, 45))





''' Kontrollerar att löängden på massnivåerna och Variabelprofilerna är lika och har samma dimension. Annars går de ej att plotta'''


print(N.shape(VARIABELPROFILER[0]))
print(VARIABELPROFILER[0])
print(N.shape(MASSLEVELSARRAY_RESHAPED[0]))
print(MASSLEVELSARRAY_RESHAPED[0])



'''Använder for-slinga för att skriva ut alla "sonderingarna" i samma figur. Sodankylä i rött, omgivande i grönt.'''


for i in range(8):
    
    plt.plot(VARIABELPROFILER[i][0,0:17], MASSLEVELSARRAY_RESHAPED[i][0:17] ,'g--')
#plt.axis([-85, -15], [0, 6000])


plt.plot(VARIABELPROFILER[8][0,0:17], MASSLEVELSARRAY_RESHAPED[8][0:17] ,'r-')
plt.xlabel('Temperature (C)')
plt.ylabel('Height (m)')
plt.title('Nearest gridpoint Sodankylä(red), Adjacent gridpoints (green)') 


'''Definierar variabeln, höjden till lägsta nivån'''

HGT_lagsta_nivan = HGT[0, 563, 352]






'''Plottar vektorn modellnivåhöjd. Första värdet i denna ska vara lika med första värdet i HGT'''

print('MODELLNIVAHOJD_1D', MODELLNIVAHOJD_1D)
#print('lon = ', lon[563,352])
#print('lat = ', lat[563,352])
print(HGT_lagsta_nivan)

#plt.plot(, A1,'r-')


'''Sparar figuren och visar plotten'''


plt.savefig('Temperature_Celsius')

plt.show()

Filobjekt.close()
