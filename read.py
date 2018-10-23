import numpy as np
from netCDF4 import Dataset
from scipy.io import netcdf

def read_nc(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    print(rootgrp.dimensions)
    for kk in rootgrp.variables.keys():
        print('key=',kk)
        if rootgrp.variables[kk].typecode()=='c':
            if len(rootgrp.variables[kk].shape)==1:
                print(''.join([ii.decode('UTF-8') for ii in list(rootgrp.variables[kk])]))
            else:
                print([''.join([ii.decode('UTF-8') for ii in k]) for k in list(rootgrp.variables[kk])])
    
        else:
            print(rootgrp.variables[kk].data)
    rootgrp.close()

def read_nc_head(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    vers = [''.join([ii.decode('UTF-8') for ii in list(rootgrp.variables['vgosDB_Version'])])]
    stat = [''.join([ii.decode('UTF-8') for ii in k]) for k in list(rootgrp.variables['StationList'])]
    sour = [''.join([ii.decode('UTF-8') for ii in k]) for k in list(rootgrp.variables['SourceList'])]
    rootgrp.close()
    return vers, stat, sour

def read_nc_stat(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    station = [''.join([ii.decode('UTF-8') for ii in k]) for k in list(rootgrp.variables['AprioriStationList'])]
    coord = rootgrp.variables['AprioriStationXYZ'].data
    rootgrp.close()
    return coord
    
def read_nc_sour(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    source = [''.join([ii.decode('UTF-8') for ii in k]) for k in list(rootgrp.variables['AprioriSourceList'])]
    coord = rootgrp.variables['AprioriSource2000RaDec'].data # in radians
    rootgrp.close()
    return coord
    
def read_nc_delay(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    group_delay = rootgrp.variables['GroupDelay'].data # in radians
    group_delay_sig = rootgrp.variables['GroupDelaySig'].data # in radians
    rootgrp.close()
    return group_delay,group_delay_sig
    
def read_nc_delay_rate(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    group_delay = rootgrp.variables['GroupRate'].data # in radians
    group_delay_sig = rootgrp.variables['GroupRateSig'].data # in radians
    rootgrp.close()
    return group_delay,group_delay_sig
    
def read_nc_CR(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    obs2stat1_stat2 = rootgrp.variables['Obs2Baseline'].data # in radians
    obs2scan = rootgrp.variables['Obs2Scan'].data # in radians
    rootgrp.close()
    return obs2scan, obs2stat1_stat2

def create_NGS(file, version,stations,sources, delay,delay_sigma, delay_rate,delay_rate_sigma, obs2scan, obs2stat1_stat2):
    out=open(file,'w')
    for i in range(len(delay)):
        n_stat1=obs2stat1_stat2[obs2scan[i]-1][0]
        n_stat2=obs2stat1_stat2[obs2scan[i]-1][1]
        #card 1
        out.write('{:s}{:s}01\n'.format(stations[n_stat1-1],stations[n_stat2-1]))
    out.close()
    
path='in/18APR01XK'
out='out/18APR01XK.ngs'
    
version,stations,sources = read_nc_head(path+"/Head.nc")

# Apriori
coord_stations=read_nc_stat(path+'/Apriori/Station.nc')
coord_sources=read_nc_sour(path+'/Apriori/Source.nc')

# Observables
delay,delay_sigma=read_nc_delay(path+'/Observables/GroupDelay_bX.nc')
delay_rate,delay_rate_sigma=read_nc_delay_rate(path+'/Observables/GroupRate_bX.nc')

# CrossReference
obs2scan, obs2stat1_stat2 = read_nc_CR(path+"/CrossReference/ObsCrossRef.nc")
    
# write NGS
create_NGS(out, version,stations,sources, delay,delay_sigma, delay_rate,delay_rate_sigma, obs2scan, obs2stat1_stat2)