import numpy as np
from scipy.io import netcdf
import os
import sys
import tarfile

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
    #rootgrp=file
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

def read_nc_St(stations,path,members):
    CableCal={}; TempC={}; AtmPres={}; RelHum={}
    
    f_C=False; f_M=False
    for m in members:
        if 'Cal-Cable.nc' in m.name:
            f_C=True
        if 'Met.nc' in m.name:
            f_M=True

    for k in stations:
        if f_C:
            rr=tar.extractfile(path[-16:-7]+'/'+k.replace(' ','')+'/Cal-Cable.nc')
            rootgrp = netcdf.NetCDFFile(rr, "r")
            CableCal[k] = rootgrp.variables['Cal-Cable'].data
        else:
            CableCal[k] = [0,]
        if f_M:
            rr=tar.extractfile(path[-16:-7]+'/'+k.replace(' ','')+'/Met.nc')
            rootgrp = netcdf.NetCDFFile(rr, "r")
            TempC[k] = rootgrp.variables['TempC'].data
            AtmPres[k] = rootgrp.variables['AtmPres'].data
            RelHum[k] = rootgrp.variables['RelHum'].data
        else:
            TempC[k] = [0,]
            AtmPres[k] = [0,]
            RelHum[k] = [0,]
    return CableCal, TempC, AtmPres, RelHum
    
def create_NGS(name,file, version,stations,sources, delay,delay_sigma, delay_rate,delay_rate_sigma, \
           tau_ion, d_tau_ion, tau_r_ion, d_tau_r_ion, YMDHM, second,RefFreq,\
           obs2scan, obs2stat1_stat2,scan2sour,coord_stations,coord_sources,axis_type,axis_offset,\
           CableCal, TempC, AtmPres, RelHum):
    # correct data
    #RelHum_i={}; CableCal_i={}
    if len(second)==1:
        second=list(second)
        for i in range(len(YMDHM)):
            second.append(second[0])

    for i in stations:
        #RelHum_i[i]=[]; CableCal_i[i]=[]
        if len(RelHum[i])==1:
            RelHum[i]=[RelHum[i][0],]
            for j in range(len(delay)):
                RelHum[i].append(RelHum[i][0])
        if len(CableCal[i])==1:
            CableCal[i]=[CableCal[i][0],]
            for j in range(len(delay)):
                CableCal[i].append(CableCal[i][0])
        if len(TempC[i])==1:
            TempC[i]=[TempC[i][0],]
            for j in range(len(delay)):
                TempC[i].append(TempC[i][0])
        if len(AtmPres[i])==1:
            AtmPres[i]=[AtmPres[i][0],]
            for j in range(len(delay)):
                AtmPres[i].append(AtmPres[i][0])
     
    scan2stat1_stat2=[]
    for i in range(obs2scan[-1]):
        scan2stat1_stat2.append([])
    
    if len(obs2stat1_stat2)>2:
        for ni,i in enumerate(obs2stat1_stat2):
            #print(obs2scan[ni]-1)
            if i[0]-1 not in scan2stat1_stat2[obs2scan[ni]-1]:
                scan2stat1_stat2[obs2scan[ni]-1].append(i[0]-1)
            if i[1]-1 not in scan2stat1_stat2[obs2scan[ni]-1]:
                scan2stat1_stat2[obs2scan[ni]-1].append(i[1]-1)
    #print('after exit',len(scan2stat1_stat2),scan2stat1_stat2[:10])
            
    n_data={}
    if len(obs2stat1_stat2)>2:
        for i in stations:
            n_data[i]=[]
        for ni,i in enumerate(scan2stat1_stat2):
            for nj, j in enumerate(stations):
                if nj in i:
                    #n_data[j].append(obs2scan[ni]-1)
                    n_data[j].append(ni)
    else:
        for i in stations:
            n_data[i]=[]
        for i in obs2scan:
            n_data[stations[0]].append(i-1)
            n_data[stations[1]].append(i-1)
    
    out=open(file,'w')
    out.write('DATA IN NGS FORMAT FROM DATABASE '+name+'_V'+version[0][2:5]+'\n')
    out.write('Observed delays and rates in card #2, modified errors in card #9\n')
    #Site cards
    dic_axis_type={1:'EQUA', 2:'XY  ', 3:'AZEL', 4:'RICH', 5:'X-YE', 6:'X-YN'}
    axis_type_i=[]
    if len(axis_type)==1:
        for i in range(len(stations)):
            axis_type_i.append(axis_type[0])
    else:
        axis_type_i=[i for i in axis_type]
    #print(len(axis_type))
    for i in range(len(stations)):
        out.write('{:8s}  {:15.5f}{:15.5f}{:15.5f} {:4s}{:10.5f}\n'.format(stations[i],
                  coord_stations[i][0],coord_stations[i][1],coord_stations[i][2],
                  dic_axis_type[axis_type_i[i]],axis_offset[i]))
    out.write('$END\n')
    # Radio source position cards
    for i in range(len(sources)):
        h=int((coord_sources[i][0]*12/np.pi)//1)
        mh=int(((coord_sources[i][0]*12/np.pi)%1)*60//1)
        sh=((coord_sources[i][0]*12/np.pi)%1)*60%1*60
        sign=' '
        d=int((coord_sources[i][1]*180/np.pi)//1)
        md=int(((coord_sources[i][1]*180/np.pi)%1)*60//1)
        sd=((coord_sources[i][1]*180/np.pi)%1)*60%1*60
        if (d<0):sign='-'
        out.write('{:8s}  {:2d} {:2d} {:10f} {:1s}{:2d} {:2d} {:10f}\n'.format(sources[i],h,mh,sh,sign,abs(d),abs(md),abs(sd)))
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

        #print(stations[n_stat1-1],stations[n_stat2-1],i,obs2scan[i]-1)
        #print(n_data[stations[n_stat1-1]][0],n_data[stations[n_stat2-1]][0])
        n_data1=n_data[stations[n_stat1-1]].index(obs2scan[i]-1)
        n_data2=n_data[stations[n_stat2-1]].index(obs2scan[i]-1)
        #print('n_data',n_data1,n_data2,stations[n_stat1-1],stations[n_stat2-1],obs2scan[i]-1)

        #card 1
        out.write('{:10s}{:10s}{:8s}{:5d} {:02d} {:02d} {:02d} {:02d}{:15.10f}          {:8d}01\n'.format(stations[n_stat1-1],
                  stations[n_stat2-1],sources[n_sour],
                  YMDHM[i][0],YMDHM[i][1],YMDHM[i][2],YMDHM[i][3],YMDHM[i][4],float(second[i]),i+1))
        #card 2
        out.write('{:20.7f}{:10.5f}{:20.10f}{:10.5f} 0      I {:8d}02\n'.format(delay[i]*10**9,delay_sigma[i]*10**9,
                  delay_rate[i]*10**12,delay_rate_sigma[i]*10**12,i+1))
        #card 3
        out.write('    .00205    .00000    .00000    .00000    .000000000000000        0.'+'{:8d}'.format(i+1)+'03\n')
        #out.write('{:s}{:8d}03\n'.format(70*' ',i+1))
        #card 4
        out.write('       .00   .0       .00   .0       .00   .0       .00   .0          '+'{:8d}'.format(i+1)+'04\n')
        #out.write('{:s}{:8d}04\n'.format(70*' ',i+1))
        #card 5
        #print('card5',stations[n_stat1-1],len(CableCal[stations[n_stat1-1]]),n_data1,i)
        #print('card5',stations[n_stat2-1],len(CableCal[stations[n_stat2-1]]),n_data2,i)
        out.write('{:10.5f}{:10.5f}    .00000    .00000    .00000    .00000          {:8d}05\n'.format(10**9*CableCal[stations[n_stat1-1]][n_data1],\
                  10**9*CableCal[stations[n_stat2-1]][n_data2],i+1))
        #card 6
        #print('card6 ',i,n_stat1-1,n_data1,obs2scan[n_data1]-1,len(TempC[stations[n_stat1-1]]),TempC[stations[n_stat1-1]][obs2scan[n_data1]-1])
        #print('card6',i,n_stat2-1,n_data2,obs2scan[n_data2]-1,len(TempC[stations[n_stat2-1]]),TempC[stations[n_stat2-1]][obs2scan[n_data2]-1])
        out.write('{:10.3f}{:10.3f}{:10.3f}{:10.3f}{:10.3f}{:10.3f} 0 0      {:8d}06\n'.format(TempC[stations[n_stat1-1]][n_data1],\
                  TempC[stations[n_stat2-1]][n_data2],\
                  AtmPres[stations[n_stat1-1]][n_data1],AtmPres[stations[n_stat2-1]][n_data2],\
                  RelHum[stations[n_stat1-1]][n_data1]*100,RelHum[stations[n_stat2-1]][n_data2]*100,i+1))
        #card 8
        out.write('{:20.10f}{:10f}{:20.10f}{:10f}  0       {:8d}08\n'.format(-tau_ion[i]*10**9,d_tau_ion[i]*10**12,-tau_r_ion[i]*10**9,d_tau_r_ion[i]*10**12,i+1))
        #card 9
        out.write('{:s}{:8d}09\n'.format(70*' ',i+1))
    out.close()
    
path=sys.argv[1]
out= sys.argv[2]

tar = tarfile.open(path, "r")
members = tar.getmembers()

rr=tar.extractfile(path[-16:-7]+"/Head.nc")
version,stations,sources = read_nc_head(rr)

# Apriori
rr=tar.extractfile(path[-16:-7]+'/Apriori/Station.nc')
coord_stations=read_nc_stat(rr)
rr=tar.extractfile(path[-16:-7]+'/Apriori/Source.nc')
coord_sources=read_nc_sour(rr)
rr=tar.extractfile(path[-16:-7]+'/Apriori/Antenna.nc')
axis_type, axis_offset=read_nc_A_ant(rr)

# Observables X
rr=tar.extractfile(path[-16:-7]+'/Observables/GroupDelay_bX.nc')
delay_x,delay_sigma_x=read_nc_delay(rr)
rr=tar.extractfile(path[-16:-7]+'/Observables/GroupRate_bX.nc')
delay_rate_x,delay_rate_sigma_x=read_nc_delay_rate(rr)
rr=tar.extractfile(path[-16:-7]+'/Observables/RefFreq_bX.nc')
RefFreq_x = read_nc_O_RF(rr)

# Observables S
f_GDS=False; f_GRS=False; f_RFS=False
for m in members:
    if 'GroupDelay_bS.nc' in m.name:
        rr=tar.extractfile(path[-16:-7]+'/Observables/GroupDelay_bS.nc')
        delay_s,delay_sigma_s=read_nc_delay(rr)
        f_GDS=True
    if 'GroupRate_bS.nc' in m.name:
        rr=tar.extractfile(path[-16:-7]+'/Observables/GroupRate_bS.nc')
        delay_rate_s,delay_rate_sigma_s=read_nc_delay_rate(rr)
        f_GRS=True
    if 'RefFreq_bS.nc' in m.name:
        rr=tar.extractfile(path[-16:-7]+'/Observables/RefFreq_bS.nc')
        RefFreq_s = read_nc_O_RF(rr)
        f_RFS=True
if not f_GDS:
    delay_s=[0,]; delay_sigma_s=[0,]
if not f_GRS:
    delay_rate_s=[0,]; delay_rate_sigma_s=[0,]
if not f_RFS:
    RefFreq_s=[0,]

# Observables
rr=tar.extractfile(path[-16:-7]+'/Observables/TimeUTC.nc')
YMDHM, second = read_nc_O_time(rr)
tau_ion=[]; tau_r_ion=[]
d_tau_ion=[]; d_tau_r_ion=[]
# what is d_RefFreq_s and d_RefFreq_x
d_RefFreq_s=0; d_RefFreq_x=0
for i in range(len(delay_x)):
    if len(delay_s)>1:
        tau_ion_i=(delay_x[i]-delay_s[i])*RefFreq_s[0]*RefFreq_s[0]/(RefFreq_x[0]*RefFreq_x[0]-RefFreq_s[0]*RefFreq_s[0])
        d_tau_ion_i= tau_ion_i * ( (delay_sigma_s[i]+delay_sigma_x[i])/(delay_s[i]+delay_x[i]) + \
                                    2*d_RefFreq_s/RefFreq_s + \
                                    2*(RefFreq_s*d_RefFreq_s+RefFreq_x*d_RefFreq_x)/(RefFreq_s*RefFreq_s+RefFreq_x*RefFreq_x) )
        tau_ion.append(tau_ion_i)
        d_tau_ion.append(d_tau_ion_i[0])
    else:
        tau_ion.append(0)
        d_tau_ion.append(0)
    if len(delay_rate_s)>1:
        tau_r_ion_i=(delay_rate_x[i]-delay_rate_s[i])*RefFreq_s[0]*RefFreq_s[0]/(RefFreq_x[0]*RefFreq_x[0]-RefFreq_s[0]*RefFreq_s[0])
        d_tau_r_ion_i=tau_r_ion_i * ( (delay_rate_sigma_s[i]+delay_rate_sigma_x[i])/(delay_s[i]+delay_x[i]) + \
                                    2*d_RefFreq_s/RefFreq_s + \
                                    2*(RefFreq_s*d_RefFreq_s+RefFreq_x*d_RefFreq_x)/(RefFreq_s*RefFreq_s+RefFreq_x*RefFreq_x) )
        tau_r_ion.append(tau_r_ion_i)
        d_tau_r_ion.append(d_tau_r_ion_i[0])
    else:
        tau_r_ion.append(0)
        d_tau_r_ion.append(0)

# S or X
if path[-9]=='X':
    delay=delay_x; delay_sigma=delay_sigma_x
    delay_rate=delay_rate_x; delay_rate_sigma=delay_rate_sigma_x
    RefFreq=RefFreq_x
if path[-9]=='S':
    delay=delay_s; delay_sigma=delay_sigma_s
    delay_rate=delay_rate_s; delay_rate_sigma=delay_rate_sigma_s
    RefFreq=RefFreq_s

# CrossReference
rr=tar.extractfile(path[-16:-7]+'/CrossReference/ObsCrossRef.nc')
obs2scan, obs2stat1_stat2 = read_nc_CR(rr)
rr=tar.extractfile(path[-16:-7]+'/CrossReference/SourceCrossRef.nc')
scan2sour = read_nc_CR_sour(rr)

# station
print(path)
CableCal, TempC, AtmPres, RelHum=read_nc_St(stations,path,members)
    
# write NGS
create_NGS(path[-16:-7],out, version,stations,sources, delay,delay_sigma, delay_rate,delay_rate_sigma, \
           tau_ion, d_tau_ion, tau_r_ion, d_tau_r_ion, YMDHM, second,RefFreq,\
           obs2scan, obs2stat1_stat2,scan2sour,coord_stations,coord_sources,axis_type,axis_offset,\
           CableCal, TempC, AtmPres, RelHum)

#read_nc(path+"/ISHIOKA/Met.nc")