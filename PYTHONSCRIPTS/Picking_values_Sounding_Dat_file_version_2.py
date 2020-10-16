import numpy as N
import scipy.io.netcdf as S
import netCDF4 as N4
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys 
import pandas as pd
import datetime

print('Python version ' + sys.version)
print('Pandas version ' + pd.__version__)
print('Matplotlib version ' + matplotlib.__version__)


#df = pd.read_csv('/home/sm_marha/TEXTFILER/luotaus_martin4.dat', encoding='latin-1', delimiter='\s+')
df = pd.read_csv('/home/sm_marha/TEXTFILER/YOPP2018FM.dat', encoding='latin-1', delimiter='\s+')


df_VARIABLES = df.loc[:, ['DATE', 'HOUR', 'ALTITUDE', 'PRESSURE', 'PARAMETER', 'DATA_VALUE'] ]

DATA_MATRIX = df_VARIABLES.values  #Görs om till Numpy array

print(DATA_MATRIX)
    

DATA_VALUES = DATA_MATRIX[1:, 5]





DATE = DATA_MATRIX[1:, 0]  

HOUR = DATA_MATRIX[1:, 1]

PRESSURE = DATA_MATRIX[1:, 3]




def pick_dates(d, sounding_time):



    DATUM = []
    datum = 20180100

    if d > 31:
        d = d-31
        datum = 20180200
    

    for i in range(len(DATE)):
        if DATA_MATRIX[i,0] == str(datum+d) and  DATA_MATRIX[i,1] == sounding_time:
            DATUM.append(DATA_MATRIX[i, :])
        
        

    print(str(datum+d) + sounding_time.zfill(2))
    print(N.shape(DATUM))
    #print(DATUM[0][0])
    


    WIND_VELOCITY = []

    ALTITUDE_WIND_VELOCITY = []

    for i in range(len(DATUM)):

        if DATUM[i][4] == str(332):
            WIND_VELOCITY.append(DATUM[i][5])
            ALTITUDE_WIND_VELOCITY.append(DATUM[i][2])



    WIND_DIRECTION = []

    ALTITUDE_WIND_DIRECTION = []


    for i in range(len(DATUM)):

        if DATUM[i][4] == str(333):
            WIND_DIRECTION.append(DATUM[i][5])
            ALTITUDE_WIND_DIRECTION.append(DATUM[i][2])




    TEMPERATURE = []

    ALTITUDE_TEMPERATURE = []

    PRESSURE_TEMPERATURE = []


    for i in range(len(DATUM)):

        if DATUM[i][4] == str(336):
            TEMPERATURE.append(DATUM[i][5])
            ALTITUDE_TEMPERATURE.append(DATUM[i][2])
            PRESSURE_TEMPERATURE.append(DATUM[i][3])
    

    DEWPOINT = []

    ALTITUDE_DEWPOINT = []

    PRESSURE_DEWPOINT = []


    for i in range(len(DATUM)):

        if DATUM[i][4] == str(337):
            DEWPOINT.append(DATUM[i][5])
            ALTITUDE_DEWPOINT.append(DATUM[i][2])
            PRESSURE_DEWPOINT.append(DATUM[i][3])

    MATRIX_2 = [WIND_VELOCITY, ALTITUDE_WIND_VELOCITY, WIND_DIRECTION, ALTITUDE_WIND_DIRECTION, TEMPERATURE, ALTITUDE_TEMPERATURE, PRESSURE_TEMPERATURE, DEWPOINT, ALTITUDE_DEWPOINT, PRESSURE_DEWPOINT]
            
    #print(MATRIX_2)

    print(N.shape(MATRIX_2))

    MATRIX_2 = N.transpose(MATRIX_2)
    
    print(N.shape(MATRIX_2))

    print('\n')

    datum = str(datum+d)[6:8] + '/' + str(datum+d)[5:6] + ' ' + sounding_time.zfill(2) + 'z'
        
    return (MATRIX_2, datum)
   
