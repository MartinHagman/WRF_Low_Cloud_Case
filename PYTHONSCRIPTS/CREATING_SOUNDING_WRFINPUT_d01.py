import numpy as N
import scipy.io.netcdf as S
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import os
import datetime
import matplotlib.colors as mplc
import netCDF4 as N4
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
from pylab import *
import Picking_values_Sounding_Dat_file_version_2
import pandas as pd






Filobjekt = N4.Dataset('/home/sm_marha/WRF_SCM/WRFV3.9.1.1_FM_FLX_KELLY_TEND_91_LEVELS/test/em_scm_xy/wrfout_d01_2018-02-27_00:00:00')


PERT_THETA = Filobjekt.variables['T'][0, :, 1, 1]


THETA = PERT_THETA+300
