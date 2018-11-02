import numpy as np
from scipy.io import netcdf
import os
import sys

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
    #DataFlag=rootgrp.variables['DataFlag'].data # in radians
    rootgrp.close()
    return group_delay,group_delay_sig
    
def read_nc_CR(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    obs2stat1_stat2 = rootgrp.variables['Obs2Baseline'].data
    obs2scan = rootgrp.variables['Obs2Scan'].data
    rootgrp.close()
    return obs2scan, obs2stat1_stat2
    
def read_nc_CR_sour(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    obs2sour = rootgrp.variables['Scan2Source'].data
    rootgrp.close()
    return obs2sour
    
def read_nc_O_time(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    YMDHM = rootgrp.variables['YMDHM'].data
    second = rootgrp.variables['Second'].data
    rootgrp.close()
    return YMDHM, second
    
def read_nc_A_ant(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    axis_type = rootgrp.variables['AntennaAxisType'].data
    axis_offset = rootgrp.variables['AntennaAxisOffset'].data
    rootgrp.close()
    return axis_type, axis_offset
    
def read_nc_O_RF(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    RefFreq = rootgrp.variables['RefFreq'].data
    rootgrp.close()
    return RefFreq

def read_nc_St(stations,path):
    CableCal={}; TempC={}; AtmPres={}; RelHum={}
    for k in stations:
        rootgrp = netcdf.NetCDFFile(path+'/'+k.replace(' ','')+'/Cal-Cable.nc', "r")
        CableCal[k] = rootgrp.variables['Cal-Cable'].data
        rootgrp = netcdf.NetCDFFile(path+'/'+k.replace(' ','')+'/Met.nc', "r")
        TempC[k] = rootgrp.variables['TempC'].data
        AtmPres[k] = rootgrp.variables['AtmPres'].data
        RelHum[k] = rootgrp.variables['RelHum'].data
    return CableCal, TempC, AtmPres, RelHum
    
def create_NGS(file, version,stations,sources, delay,delay_sigma, delay_rate,delay_rate_sigma, tau_ion,tau_r_ion, YMDHM, second,RefFreq,\
           obs2scan, obs2stat1_stat2,scan2sour,coord_stations,coord_sources,axis_type,axis_offset,\
           CableCal, TempC, AtmPres, RelHum):
    out=open(file,'w')
    out.write('DATA IN NGS FORMAT FROM DATABASE '+file[-16:-7]+'_V'+version[0][2:5]+'\n')
    out.write('Observed delays and rates in card #2, modified errors in card #9\n')
    #Site cards
    dic_axis_type={1:'EQ  ', 2:'XY  ', 3:'AZEL'}
    for i in range(len(stations)):
        out.write('{:8s}  {:15.5f}{:15.5f}{:15.5f} {:4s}{:10.5f}\n'.format(stations[i],
                  coord_stations[i][0],coord_stations[i][1],coord_stations[i][2],
                  dic_axis_type[axis_type[i]],axis_offset[i]))
    out.write('$END\n')
    # Radio source position cards
    for i in range(len(sources)):
        h=int((coord_sources[i][0]*12/np.pi)//1)
        mh=int(((coord_sources[i][0]*12/np.pi)%1)*60//1)
        sh=((coord_sources[i][0]*12/np.pi)%1)*60%1*60
        d=int((coord_sources[i][1]*180/np.pi)//1)
        md=int(((coord_sources[i][1]*180/np.pi)%1)*60//1)
        sd=((coord_sources[i][1]*180/np.pi)%1)*60%1*60
        out.write('{:8s}  {:2d} {:2d} {:10f} {:2d} {:2d} {:10f}\n'.format(sources[i],h,mh,sh,d,md,sd))
    out.write('$END\n')
    # Auxiliary parameters
    out.write('{:15.10e}            GR PH\n'.format(RefFreq[0]))
    out.write('$END\n')
    # Data cards
    for i in range(len(delay)):
        if len(stations)>2:
            n_stat1=obs2stat1_stat2[i][0]
            n_stat2=obs2stat1_stat2[i][1]
        else:
            n_stat1=obs2stat1_stat2[0]
            n_stat2=obs2stat1_stat2[1]
        n_sour=scan2sour[obs2scan[i]-1]-1
        #card 1
        out.write('{:10s}{:10s}{:8s}{:5d} {:02d} {:02d} {:02d} {:02d}{:15.10f}          {:8d}01\n'.format(stations[n_stat1-1],
                  stations[n_stat2-1],sources[n_sour],
                  YMDHM[i][0],YMDHM[i][1],YMDHM[i][2],YMDHM[i][3],YMDHM[i][4],float(second[i]),i+1))
        #card 2
        out.write('{:20.7f}{:10.5f}{:20.10f}{:10.5f} 0      I {:8d}02\n'.format(delay[i]*10**9,delay_sigma[i]*10**9,
                  delay_rate[i]*10**9,delay_rate_sigma[i]*10**9,i+1))
        #card 3
        out.write('    .00205    .00000    .00000    .00000    .000000000000000        0.'+'{:8d}'.format(i+1)+'03\n')
        #out.write('{:s}{:8d}03\n'.format(70*' ',i+1))
        #card 4
        out.write('       .00   .0       .00   .0       .00   .0       .00   .0          '+'{:8d}'.format(i+1)+'04\n')
        #out.write('{:s}{:8d}04\n'.format(70*' ',i+1))
        #card 5
        out.write('{:10.5f}{:10.5f}    .00000    .00000    .00000    .00000          {:8d}05\n'.format(10**9*CableCal[stations[n_stat1-1]][obs2scan[i]-1],\
                  10**9*CableCal[stations[n_stat2-1]][obs2scan[i]-1],i+1))
        #card 6
        out.write('{:10.3f}{:10.3f}{:10.3f}{:10.3f}{:10.3f}{:10.3f} 0 0      {:8d}06\n'.format(TempC[stations[n_stat1-1]][obs2scan[i]-1],\
                  TempC[stations[n_stat2-1]][obs2scan[i]-1],\
                  AtmPres[stations[n_stat1-1]][obs2scan[i]-1],AtmPres[stations[n_stat2-1]][obs2scan[i]-1],\
                  RelHum[stations[n_stat1-1]][obs2scan[i]-1]*100,RelHum[stations[n_stat2-1]][obs2scan[i]-1]*100,i+1))
        #card 8
        out.write('{:20.10f}{:10f}{:20.10f}{:10f}  0       {:8d}08\n'.format(-tau_ion[i]*10**9,0,-tau_r_ion[i]*10**9,0,i+1))
        #card 9
        out.write('{:s}{:8d}09\n'.format(70*' ',i+1))
    out.close()
    
path=sys.argv[1]
out=sys.argv[2]
    
version,stations,sources = read_nc_head(path+"/Head.nc")

# Apriori
coord_stations=read_nc_stat(path+'/Apriori/Station.nc')
coord_sources=read_nc_sour(path+'/Apriori/Source.nc')
axis_type, axis_offset=read_nc_A_ant(path+"/Apriori/Antenna.nc")

# Observables X
delay_x,delay_sigma_x=read_nc_delay(path+'/Observables/GroupDelay_bX.nc')
delay_rate_x,delay_rate_sigma_x=read_nc_delay_rate(path+'/Observables/GroupRate_bX.nc')
RefFreq_x = read_nc_O_RF(path+'/Observables/RefFreq_bX.nc')

# Observables S
delay_s,delay_sigma_s=read_nc_delay(path+'/Observables/GroupDelay_bS.nc')
delay_rate_s,delay_rate_sigma_s=read_nc_delay_rate(path+'/Observables/GroupRate_bS.nc')
RefFreq_s = read_nc_O_RF(path+'/Observables/RefFreq_bS.nc')

# Observables
YMDHM, second = read_nc_O_time(path+"/Observables/TimeUTC.nc")
tau_ion=[]; tau_r_ion=[]
for i in range(len(delay_x)):
    tau_ion.append((delay_x[i]-delay_s[i])*RefFreq_s[0]*RefFreq_s[0]/(RefFreq_x[0]*RefFreq_x[0]-RefFreq_s[0]*RefFreq_s[0]))
    tau_r_ion.append((delay_rate_x[i]-delay_rate_s[i])*RefFreq_s[0]*RefFreq_s[0]/(RefFreq_x[0]*RefFreq_x[0]-RefFreq_s[0]*RefFreq_s[0]))

# S or X
if path[-2]=='X':
    delay=delay_x; delay_sigma=delay_sigma_x
    delay_rate=delay_rate_x; delay_rate_sigma=delay_rate_sigma_x
    RefFreq=RefFreq_x
if path[-2]=='S':
    delay=delay_s; delay_sigma=delay_sigma_s
    delay_rate=delay_rate_s; delay_rate_sigma=delay_rate_sigma_s
    RefFreq=RefFreq_s

# CrossReference
obs2scan, obs2stat1_stat2 = read_nc_CR(path+"/CrossReference/ObsCrossRef.nc")
scan2sour = read_nc_CR_sour(path+"/CrossReference/SourceCrossRef.nc")

# station
CableCal, TempC, AtmPres, RelHum=read_nc_St(stations,path)
    
# write NGS
create_NGS(out, version,stations,sources, delay,delay_sigma, delay_rate,delay_rate_sigma, tau_ion,tau_r_ion, YMDHM, second,RefFreq,\
           obs2scan, obs2stat1_stat2,scan2sour,coord_stations,coord_sources,axis_type,axis_offset,\
           CableCal, TempC, AtmPres, RelHum)

#read_nc(path+"/ISHIOKA/Met.nc")