import traceback
import argparse
import logging
import os
import sys
import tarfile
import fnmatch
#from pathlib import Path
import numpy as np
from scipy.io import netcdf
np.seterr(all='raise')

logging.basicConfig(
    format=u'%(levelname)-8s [%(asctime)s] %(message)s',
    level=logging.DEBUG,
    filename='VGOS2ngs.log')


def main():
    '''
    Read VGOS_DB format files and write ngs-files
    '''
    args = parse_args()
    if args.archive:
        translate(args.archive)
    elif args.archives_dir:
        for archive_name in os.listdir(args.archives_dir):
            if archive_name.endswith('.tar.gz') or archive_name.endswith(
                    '.tgz'):
                try:
                    translate(os.path.join(args.archives_dir, archive_name))
                except Exception as err:
                    print('Error processing {}'.format(archive_name))
                    traceback.print_exc()
                    logging.error(u'Error processing {}'.format(archive_name))

def parse_args():
    args = argparse.ArgumentParser()
    args.add_argument('--archive')
    args.add_argument('--archives-dir')
    return args.parse_args()

def translate(vigos_archive):
    name = os.path.basename(vigos_archive).split('.')[0].upper()
    #path='in/18JAN25XE.tar.gz'#sys.argv[1]
    out = os.path.join('out', name)  #'out/18JAN25XE'#sys.argv[2]

    tar = tarfile.open(vigos_archive, "r")
    members = tar.getmembers()

    try:
        rr = tar.extractfile(name + '/Head.nc')
        version, stations, sources = read_nc_head(rr)
    except KeyError:
        print('Head file for {} is missing, skipping'.format(name))
        logging.error(u'Head file for {} is missing, skipping'.format(name))
        return

    try:
        rr = tar.extractfile(name + '/Apriori/Station.nc')
        coord_stations = read_nc_stat(rr)
        rr = tar.extractfile(name + '/Apriori/Source.nc')
        coord_sources = read_nc_sour(rr)
        rr = tar.extractfile(name + '/Apriori/Antenna.nc')
        axis_type, axis_offset = read_nc_A_ant(rr)
    except KeyError:
        print('Apriori file for {} is missing, skipping'.format(name))
        logging.error(u'Apriori file for {} is missing, skipping'.format(name))
#    else:  
        coord_stations = [ [0] * 3  for i in range(len(stations))]
        coord_sources = [ [0] * 2  for i in range(len(sources))]
        axis_type = [7]*len(stations)
        axis_offset = [0]*len(stations) 

    # Observables X

    rr = tar.extractfile(name + '/Observables/GroupDelay_bX.nc')
    delay_x, delay_sigma_x = read_nc_delay(rr)
    rr = tar.extractfile(name + '/Observables/GroupRate_bX.nc')
    delay_rate_x, delay_rate_sigma_x = read_nc_delay_rate(rr)
    rr = tar.extractfile(name + '/Observables/RefFreq_bX.nc')
    RefFreq_x = read_nc_O_RF(rr)
    
    rr = tar.extractfile(name + '/Observables/GroupDelay_bS.nc')
    delay_s, delay_sigma_s = read_nc_delay(rr)
    rr = tar.extractfile(name + '/Observables/GroupRate_bS.nc')
    delay_rate_s, delay_rate_sigma_s = read_nc_delay_rate(rr)
    rr = tar.extractfile(name + '/Observables/RefFreq_bS.nc')
    RefFreq_s = read_nc_O_RF(rr)

    try:
        file_names = [m.name for m in members]
        matched = fnmatch.filter( file_names, name + '/ObsEdit/GroupDelayFull_*bX*.nc')
        rr = tar.extractfile(matched[-1])
#        rr = tar.extractfile(name + '/ObsEdit/GroupDelayFull_iIVS_bX.nc') 
        delay_x_full = read_nc_delayfull(rr)
#        print('x_full', 'x')
#        for i in range (len(delay_x_full) ):
#            print (delay_x_full[i], delay_x[i])
    except IndexError:
        print('/ObsEdit file for {} is missing, skipping'.format(name))
        logging.error(u'/ObsEdit for {} is missing, skipping'.format(name))
        ambX = [0]*len(delay_x)
        rr = tar.extractfile(name + '/Observables/AmbigSize_bX.nc')
        ambX = read_ambiguity(rr)
        if len(ambX) == 1:
            ambX = ambX[0]*len(delay_x)
   
        delay_x_full = delay_x + ambX
        print(ambX)
  
    try:
        file_names = [m.name for m in members]
        matched = fnmatch.filter( file_names, name + '/ObsEdit/GroupDelayFull_*bS*.nc')
        rr = tar.extractfile(matched[-1])
#        rr = tar.extractfile(name + '/ObsEdit/GroupDelayFull_bS.nc')
#        rr = tar.extractfile(name + '/ObsEdit/GroupDelayFull_iIVS_bS.nc')
        delay_s_full = read_nc_delayfull(rr)
    except IndexError:
        print('/ObsEdit file for {} is missing, skipping'.format(name))
        logging.error(u'/ObsEdit for {} is missing, skipping'.format(name))  
        rr = tar.extractfile(name + '/Observables/AmbigSize_bS.nc')
        ambS = [0]*len(delay_s)
        ambS = read_ambiguity(rr)
        if len(ambS) == 1:
            ambS = ambS[0]*len(delay_s)

        delay_s_full = delay_s+ ambS
        print(ambS)
#        print(len(ambS))
  
        
   # S or X
    if name[-2] == 'X':
        delay = delay_x_full
        delay_sigma = delay_sigma_x
        delay_rate = delay_rate_x
        delay_rate_sigma = delay_rate_sigma_x
        RefFreq = RefFreq_x
    if name[-2] == 'S':
        delay = delay_s_full
        delay_sigma = delay_sigma_s
        delay_rate = delay_rate_s
        delay_rate_sigma = delay_rate_sigma_s
        RefFreq = RefFreq_s
      
    Nobs = len(delay)    
                
#    data_quality_flag = [0] * len(delay)
#  Correlation for ngs card 3
    try:
        rr = tar.extractfile(name + '/Observables/Correlation_bX.nc')
        Correlation = read_cor(rr)
    except KeyError:
        print('Correlation_bX file for {} is missing, skipping'.format(name))
        logging.error(u'Correlation_bX file for {} is missing, skipping'.format(name))

#  Data Quality Flag  for ngs card 2
    dataQualityFlag = [0]*Nobs
    try:
        file_names = [m.name for m in members]
        matched = fnmatch.filter( file_names, name + '/ObsEdit/Edit*.nc')
        rr = tar.extractfile(matched[-1])
#        rr = tar.extractfile(name + '/ObsEdit/Edit.nc')
#        rr = tar.extractfile(name + '/ObsEdit/Edit_iIVS.nc')
        delayflag = read_DelayFlag(rr)
        dataQualityFlag = delayflag[0]

#        print(dataQualityFlag)  
#        print(len(dataQualityFlag))
        if len(dataQualityFlag) == 1:
            dataQualityFlag = [0]*Nobs
#            for i in range(Nobs):
#               dataQualityFlag[i] = dataQualityFlag[0]
       #        print('222:',len(dataQualityFlag))
    except IndexError:
        print('ObsEdit/Edit_iIVS.nc for {} is missing, skipping'.format(name))
        logging.error(u'ObsEdit/Edit_iIVS.nc for {} is missing, skipping'.format(name))
        
        rr = tar.extractfile(name + '/Observables/QualityCode_bX.nc')
        QCodeX = read_qcode(rr)
    
        rr = tar.extractfile(name + '/Observables/QualityCode_bS.nc')
        QCodeS = read_qcode(rr)

        dataQualityFlag = read_DataQualityFlag_fromDataQualityCode(Nobs,QCodeX,QCodeS)

# Observables
    rr = tar.extractfile(name + '/Observables/TimeUTC.nc')
    YMDHM, second = read_nc_O_time(rr)
    rr = tar.extractfile(name + '/Observables/Phase_bX.nc')
    Phase, PhaseSig = read_Phase(rr)
    
# Ionospheric delay    
    delay_ion = [0]*Nobs
    delay_ion_r = [0]*Nobs
    sigma_delay_ion = [0]*Nobs
    sigma_delay_ion_r = [0]*Nobs

    try:
        rr = tar.extractfile(name + '/ObsDerived/Cal-SlantPathIonoGroup_bX.nc')
        tau_ion, d_tau_ion = read_del_IONO(rr)
     
        for i in range(Nobs):
            delay_ion[i] = tau_ion[i][0] * 10**9
            delay_ion_r[i] = tau_ion[i][1] * 10**12
            sigma_delay_ion[i] = d_tau_ion[i][0] * 10**9
            sigma_delay_ion_r[i] = d_tau_ion[i][1] * 10**12
    except KeyError:

        logging.error(u'Cal-SlantPathIonoGroup_bX.nc for {} is missing, skipping'.format(name))
        
        rr = tar.extractfile(name + '/Observables/ChannelInfo_bX.nc')
        numChannelsX,bandX,ChannelFreqX = read_channalinfo(rr)
        rr = tar.extractfile(name + '/Observables/ChannelInfo_bS.nc')
        numChannelsS,bandS,ChannelFreqS = read_channalinfo(rr)
          
        delay_ion, delay_ion_r,sigma_delay_ion, sigma_delay_ion_r, dataQualityFlag =\
        calc_del_IONO(name,numChannelsX,ChannelFreqX,numChannelsS,ChannelFreqS,\
        delay_x,delay_s,delay_sigma_x,delay_sigma_s,delay_rate_x,delay_rate_s,\
        delay_rate_sigma_x, delay_rate_sigma_s,Nobs,dataQualityFlag)

    # CrossReference
    rr = tar.extractfile(name + '/CrossReference/ObsCrossRef.nc')
    obs2scan, obs2stat1_stat2 = read_nc_CR(rr)
#    print('obs2scan: ', obs2scan)
#    print(len(obs2scan))
    
#    print('obs2stat1_stat2:', obs2stat1_stat2)
#    print(len(obs2stat1_stat2))
    
    rr = tar.extractfile(name + '/CrossReference/SourceCrossRef.nc')
    scan2sour = read_nc_CR_sour(rr)
#    print('scan2sour:', scan2sour)
#    print(len(scan2sour))
    
    rr = tar.extractfile(name + '/CrossReference/StationCrossRef.nc')
    scan2sta, sta2scan,numscanpersta = read_nc_CR_sta(rr)
#    print('sta2scan:', sta2scan)
#    print(len(sta2scan))
#    print( 'scan2sta', scan2sta)
#    print(len(scan2sta))
#    print(len(numscanpersta))

    # station
    CableCal, TempC, AtmPres, RelHum = read_nc_St(name, stations, members, tar)
    # write NGS
    
    create_NGS(name,out,version,stations,sources,delay,delay_sigma, delay_rate,delay_rate_sigma,\
               delay_ion, delay_ion_r, sigma_delay_ion, sigma_delay_ion_r,\
               YMDHM, second,RefFreq,obs2scan, obs2stat1_stat2,scan2sour,\
               coord_stations,coord_sources,axis_type,axis_offset,\
               CableCal, TempC, AtmPres, RelHum,Correlation,Phase,PhaseSig,dataQualityFlag)
    print('Translation of {} is complete'.format(name))
    logging.info(u'Translation of {} is complete'.format(name))


def read_nc(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    print('dimensions',rootgrp.dimensions)
    for kk in rootgrp.variables.keys():
        print('key=',kk)
        if rootgrp.variables[kk].typecode() == 'c':
            if len(rootgrp.variables[kk].shape) == 1:
                pass
            else:
                pass
        else:
            pass
    rootgrp.close()


def read_nc_head(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    vers = [
        ''.join([
            ii.decode('UTF-8')
            for ii in list(rootgrp.variables['vgosDB_Version'])
        ])
    ]
    stations = [
        ''.join([ii.decode('UTF-8') for ii in k])
        for k in list(rootgrp.variables['StationList'])
    ]
    stations = [station.strip() for station in stations]
    sour = [
        ''.join([ii.decode('UTF-8') for ii in k])
        for k in list(rootgrp.variables['SourceList'])
    ]
    rootgrp.close()
    return vers, stations, sour


def read_nc_stat(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    station = [
        ''.join([ii.decode('UTF-8') for ii in k])
        for k in list(rootgrp.variables['AprioriStationList'])
    ]
    coord = rootgrp.variables['AprioriStationXYZ'].data
    rootgrp.close()
    return coord


def read_nc_sour(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    source = [
        ''.join([ii.decode('UTF-8') for ii in k])
        for k in list(rootgrp.variables['AprioriSourceList'])
    ]
    coord = rootgrp.variables['AprioriSource2000RaDec'].data  # in radians
    rootgrp.close()
    return coord


def read_nc_delay(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    group_delay = rootgrp.variables['GroupDelay'].data  # in sec
    group_delay_sig = rootgrp.variables['GroupDelaySig'].data  # in sec
    rootgrp.close()
    return group_delay, group_delay_sig
    
def read_nc_delayfull(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    group_delay_full = rootgrp.variables['GroupDelayFull'].data  # in sec
    rootgrp.close()
    return group_delay_full


def read_nc_delay_rate(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    group_delay = rootgrp.variables['GroupRate'].data  # in 
    group_delay_sig = rootgrp.variables['GroupRateSig'].data  # in 
    #DataFlag=rootgrp.variables['DataFlag'].data # in radians
    rootgrp.close()
    return group_delay, group_delay_sig


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
    
   
def read_nc_CR_sta(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
#    numsta = rootgrp.variables['NumStation'].data
#    numscan = rootgrp.variables['NumScan'].data
    numscanpersta = rootgrp.variables['NumScansPerStation'].data
    scan2sta = rootgrp.variables['Scan2Station'].data
    sta2scan = rootgrp.variables['Station2Scan'].data
    rootgrp.close()
    return scan2sta, sta2scan,numscanpersta


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
    
def read_cor(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    print(rootgrp.variables.keys())
    cor = rootgrp.variables['Correlation'].data
    rootgrp.close()
    return cor
    return rootgrp.variables.keys()
    
def read_dataflag(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    print(rootgrp.variables.keys())
    dataflag = rootgrp.variables['DataFlag'].data
    rootgrp.close()
    return dataflag
    return rootgrp.variables.keys()       
  
def read_channalinfo(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    print(rootgrp.variables.keys())
    numChannels = rootgrp.variables['NumChannels'].data
    band = rootgrp.variables['Band'].data
#    sampleRate = rootgrp.variables['SampleRate'].data
#    BITSAMPL = rootgrp.variables['BITSAMPL'].data
#    ChannelID = rootgrp.variables['ChannelID'].data
#    Polarization = rootgrp.variables['Polarization'].data
    ChannelFreq = rootgrp.variables['ChannelFreq'].data 
#    NumAp = rootgrp.variables['NumAp'].data 
#    ChanAmpPhase = rootgrp.variables['ChanAmpPhase'].data 
#    NumSamples = rootgrp.variables['NumSamples'].data 
    rootgrp.close()
#    return numChannels,band,sampleRate,BITSAMPL,ChannelID,Polarization,ChannelFreq,NumAp,ChanAmpPhase,NumSamples
    return numChannels,band,ChannelFreq
    return rootgrp.variables.keys()  

def read_qcode(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    qcode = rootgrp.variables['QualityCode'].data
    rootgrp.close()
    return qcode
    return rootgrp.variables.keys()

def read_Phase(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    Phase = rootgrp.variables['Phase'].data
    PhaseSig = rootgrp.variables['PhaseSig'].data
    rootgrp.close()
    return Phase, PhaseSig
    return rootgrp.variables.keys()

def read_DelayFlag(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    print(rootgrp.variables.keys())
    delayflag = rootgrp.variables['DelayFlag'].data
    rateflag = rootgrp.variables['RateFlag'].data
    rootgrp.close()
    return delayflag, rateflag
    return rootgrp.variables.keys()    

def read_nc_Eff_Fr(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    Eff_Fr = rootgrp.variables['FreqGroupIono'].data
    rootgrp.close()
    return Eff_Fr
    return rootgrp.variables.keys()

def read_del_IONO(file):
    rootgrp = netcdf.NetCDFFile(file, "r")
    tau_ion = rootgrp.variables['Cal-SlantPathIonoGroup'].data
    d_tau_ion = rootgrp.variables['Cal-SlantPathIonoGroupSigma'].data
    rootgrp.close()
    return tau_ion, d_tau_ion
    return rootgrp.variables.keys()

def read_ambiguity(file):  
    rootgrp = netcdf.NetCDFFile(file, "r")
    print(rootgrp.variables.keys())
#    ambig = rootgrp.variables['GPDLAMBG'].data
    ambig = rootgrp.variables['AmbigSize'].data
    rootgrp.close()
    return ambig
    return rootgrp.variables.keys()
    

def read_nc_St(name, stations, members, tar):
    CableCal = {}
    TempC = {}
    AtmPres = {}
    RelHum = {}

    for k in stations:
        try:
            rr = tar.extractfile(name + '/' + k.replace(' ', '') +
                                 '/Cal-Cable.nc')
            rootgrp = netcdf.NetCDFFile(rr, "r")
            # FIXME: CAl-Cable file might be missing
            CableCal[k] = rootgrp.variables['Cal-Cable'].data
        except KeyError:
            CableCal[k] = [
                0,
            ]
        try:
            rr = tar.extractfile(name + '/' + k.replace(' ', '') + '/Met.nc')
            rootgrp = netcdf.NetCDFFile(rr, "r")
            TempC[k] = rootgrp.variables['TempC'].data
            AtmPres[k] = rootgrp.variables['AtmPres'].data
            RelHum[k] = rootgrp.variables['RelHum'].data
        except KeyError:
            print('/No meteo data  for {}, skipping'.format(name))
            logging.error(u'/No meteo data  for {}, skipping'.format(name))
            TempC[k] = [
                10,
            ]
            AtmPres[k] = [
                1000,
            ]
            RelHum[k] = [
                0.50,
            ]
    return CableCal, TempC, AtmPres, RelHum
    
def read_DataQualityFlag_fromDataQualityCode(N,QCodeX,QCodeS):
#    GoodQualityCode = {'5', '6', '7', '8','9','H'} 
    GoodQualityCode = {'5', '6', '7', '8','9'}
    DataQualityFlagX = [0]*N
    DataQualityFlagS = [0]*N
    DataQualityFlag = [0]*N

    if len(QCodeS) == 1:
        if QCodeS[0] in GoodQualityCode: DataQualityFlagS[0] = 0
        for i in range(N):
            DataQualityFlagS.append(DataQualityFlagS[0])
    else:
        for i in range(N):
            if QCodeS[i] in GoodQualityCode:
                DataQualityFlagS[i] = 0 
            else:
                DataQualityFlagS[i] = 1
    if len(QCodeX) == 1:
        if QCodeX[0] in GoodQualityCode: DataQualityFlagX[0] = 0
        for i in range(N):
            DataQualityFlagX.append(DataQualityFlagX[0])
    else:
        for i in range(N):
            if QCodeX[i] in GoodQualityCode:
                DataQualityFlagX[i] = 0 
            else:
                DataQualityFlagX[i] = 1

        DataQualityFlag[i] = DataQualityFlagS[i] + DataQualityFlagX[i]  
    return DataQualityFlag

    
def calc_del_IONO(name,numChannelsX,ChannelFreqX,numChannelsS,ChannelFreqS,\
                          tauX,tauS,sigma_tauX,sigma_tauS,\
                  tauX_r,tauS_r,sigma_tauX_r,sigma_tauS_r,N, data_quality_flag):
#  N - number of observations
    
    Car_Freq_S = [0]*N
    Car_Freq_X = [0]*N
    C_ion = [0]*N
    tau_ion =  [0]*N 
    tau_ion_r =  [0]*N  
    sigma_tau_ion = [0]*N
    sigma_tau_ion_r = [0]*N
    
    if len(numChannelsX) == 1:
        S = sum(ChannelFreqX)
        NC = numChannelsX[0]
        CarFreqX = S/NC
        for i in range(N):
#            Car_Freq_X.insert(i,CarFreqX)
            Car_Freq_X[i] = CarFreqX
    else:
        for i in range(N):
            try:
                S = sum(ChannelFreqX[i])
                Car_Freq_X[i] = S/numChannelsX[i]
            except FloatingPointError:
                if numChannelsX[i] == 0:
                   Car_Freq_X[i] =  8562.99
                   data_quality_flag[i] =  2
    if len(numChannelsS) == 1:
        try:
            S = sum(ChannelFreqS)
            CarFreqS = S/numChannelsS[0]
        except FloatingPointError:
            Car_Freq_S[0] = 2265.99
            data_quality_flag[i] =  2
        for i in range(N):
            Car_Freq_S.insert(i,CarFreqS)
    else:
        for i in range(N):
            try:
                S = sum(ChannelFreqS[i])
                S = S/ numChannelsS[i]
                Car_Freq_S[i] = S
            except FloatingPointError:
                if numChannelsS[i]==0:
                    Car_Freq_S[i] = 2265.99
                    data_quality_flag[i] = 8

    for i in range(N):            
        try:
            C_ion[i] = Car_Freq_S[i]**2 / (Car_Freq_X[i]**2 - Car_Freq_S[i]**2)
                       
            tau_ion[i] = C_ion[i] * (tauS[i] - tauX[i])*10**9
            sigma_tau_ion[i]  =  10**9 * C_ion[i]*(sigma_tauS[i]**2 + sigma_tauX[i]**2)**0.5
            
            tau_ion_r[i] =  10**12 *  C_ion[i] * (tauS_r[i] - tauX_r[i])
            sigma_tau_ion_r[i]  = 10**12 *  C_ion[i] * (sigma_tauS_r[i]**2 + sigma_tauX_r[i]**2)**0.5
        except FloatingPointError:
            C_ion[i] = 0.076
            tau_ion[i] = C_ion[i] * (tauS[i] - tauX[i])*10**9
            sigma_tau_ion[i]  = 10**9 * C_ion[i]*(sigma_tauS[i]**2 + sigma_tauX[i]**2)**0.5
            
            tau_ion_r[i] =  10**12 * C_ion[i] * (tauS_r[i] - tauX_r[i])
            sigma_tau_ion_r[i]  =  10**12 * C_ion[i] * (sigma_tauS_r[i]**2 + sigma_tauX_r[i]**2)**0.5
            data_quality_flag[i] =  2
    return (tau_ion,tau_ion_r,sigma_tau_ion,sigma_tau_ion_r,data_quality_flag)

def create_NGS(name,file,version,stations,sources,delay,delay_sigma,delay_rate,delay_rate_sigma,\
           delay_ion, delay_ion_r, sigma_delay_ion, sigma_delay_ion_r,\
           YMDHM, second,RefFreq,obs2scan, obs2stat1_stat2,scan2sour,\
           coord_stations,coord_sources,axis_type,axis_offset,\
           CableCal, TempC, AtmPres, RelHum,Correlation,Phase,PhaseSig,data_quality_flag):
    
        # correct data
    #RelHum_i={}; CableCal_i={}
    if len(second) == 1:
        second = list(second)
        for i in range(len(YMDHM)):
            second.append(second[0])     
    for station in stations:
        if len(RelHum[station]) == 1:
            RelHum[station] = [
                RelHum[station][0],
            ]
            for _ in range(len(delay)):
                RelHum[station].append(RelHum[station][0])
        if len(CableCal[station]) == 1:
            CableCal[station] = [
                CableCal[station][0],
            ]
            for _ in range(len(delay)):
                CableCal[station].append(CableCal[station][0])
        if len(TempC[station]) == 1:
            TempC[station] = [
                TempC[station][0],
            ]
            for _ in range(len(delay)):
                TempC[station].append(TempC[station][0])
        if len(AtmPres[station]) == 1:
            AtmPres[station] = [
                AtmPres[station][0],
            ]
            for _ in range(len(delay)):
                AtmPres[station].append(AtmPres[station][0])

    scan2stat1_stat2 = []
   
    for i in range(obs2scan[-1]):
        scan2stat1_stat2.append([])
                   
    if len(obs2stat1_stat2)>2:
        for ni,i in enumerate(obs2stat1_stat2):
            if i[0]-1 not in scan2stat1_stat2[obs2scan[ni]-1]:
                scan2stat1_stat2[obs2scan[ni]-1].append(i[0]-1)
            if i[1]-1 not in scan2stat1_stat2[obs2scan[ni]-1]:
                scan2stat1_stat2[obs2scan[ni]-1].append(i[1]-1)
                             
    out=open('out.txt','w')            
    n_data={}
    if len(obs2stat1_stat2)>2:
        for i in stations:
            n_data[i]=[]
        for ni,i in enumerate(scan2stat1_stat2):
            out.write('scan2stat1_stat2 ni='+str(ni)+' i='+str(i)+'\n')
            for nj, j in enumerate(stations):
                out.write('stations nj='+str(nj)+' j='+str(j)+'\n')
                if nj in i:
                    out.write('n_data['+str(j)+'.append('+str(ni)+')\n')
                    n_data[j].append(ni)
    else:
        for i in stations:
            n_data[i]=[]
        for i in obs2scan:
            n_data[stations[0]].append(i-1)
            n_data[stations[1]].append(i-1)
    out.close()  
 
    out = open(file+'_N'+version[0][2:5],'w')
    out.write('DATA IN NGS FORMAT FROM DATABASE {name}_V{version}\n'.format(
        name=name,version=version[0][2:5]))
    out.write(
        'Observed delays and rates in card #2\n')
#      'Observed delays and rates in card #2, modified errors in card #9\n')
    #Site cards
    dic_axis_type = {
        1: 'EQUA',
        2: 'XY  ',
        3: 'AZEL',
        4: 'RICH',
        5: 'X-YE',
        6: 'X-YN',
        7: '    ' 
    }
    axis_type_i = [0]*len(stations)
    
    if len(axis_type) == 1:
        for i in range(len(stations)):
            axis_type_i[i] = axis_type[0]
    else:
        for station_idx, station_name in enumerate(stations):
            try:
                axis_type_i = [i for i in axis_type]
            except IndexError:
                axis_type_i = 7
    for station_idx, station_name in enumerate(stations):
        try:
            offset = axis_offset[station_idx]
        except IndexError:
            offset = 0.0
            
        coords_for_station = coord_stations[station_idx]

        axis_type_name = dic_axis_type[axis_type_i[station_idx]]

        out.write('{:8s}  {:15.5f}{:15.5f}{:15.5f} {:4s}{:10.5f}\n'.format(
            station_name, coords_for_station[0], coords_for_station[1],
            coords_for_station[2], axis_type_name, offset))
    out.write('$END\n')

    # Radio source position cards
    for station_idx, source_name in enumerate(sources):
        coord_ra, coord_dec = coord_sources[station_idx]
        h = int((coord_ra * 12 / np.pi) // 1)
        mh = int(((coord_ra * 12 / np.pi) % 1) * 60 // 1)
        sh = ((coord_ra * 12 / np.pi) % 1) * 60 % 1 * 60
        sign = ' '
        d = int((coord_dec * 180 / np.pi) // 1)
        md = int(((coord_dec * 180 / np.pi) % 1) * 60 // 1)
        sd = ((coord_dec * 180 / np.pi) % 1) * 60 % 1 * 60
        if (d < 0): sign = '-'
        out.write('{:8s}  {:2d} {:2d} {:10f} {:1s}{:2d} {:2d} {:10f}\n'.format(
            source_name, h, mh, sh, sign, abs(d), md, sd))
    out.write('$END\n')
    # Auxiliary parameters
    out.write('{:15.10e}            GR PH\n'.format(RefFreq[0]))
    out.write('$END\n')
    # Data cards
    for i in range(len(delay)):
        
        if len(stations) > 2:
            n_stat1 = obs2stat1_stat2[i][0]
            n_stat2 = obs2stat1_stat2[i][1]
        else:
            n_stat1 = obs2stat1_stat2[0]
            n_stat2 = obs2stat1_stat2[1]

        n_sour = scan2sour[obs2scan[i] - 1] - 1

        n_data1 = n_data[stations[n_stat1 - 1]].index(obs2scan[i] - 1)
        n_data2 = n_data[stations[n_stat2 - 1]].index(obs2scan[i] - 1)

        if YMDHM[i][0] < 50:  YMDHM[i][0] =  YMDHM[i][0] + 2000
        #card 1
        out.write(
            '{:10s}{:10s}{:8s}{:5d} {:02d} {:02d} {:02d} {:02d}{:15.10f}          {:8d}01\n'\
            .format(stations[n_stat1 - 1], stations[n_stat2 - 1],\
                    sources[n_sour], YMDHM[i][0], YMDHM[i][1], YMDHM[i][2],\
                    YMDHM[i][3], YMDHM[i][4], float(second[i]), i + 1))
        #card 2

        if data_quality_flag[i] > 10: data_quality_flag[i] = 8 
        out.write(
            '{:20.8f}{:10.5f}{:20.10f}{:10.5f} {:1d}      I {:8d}02\n'.format(\
                delay[i] * 10**9, delay_sigma[i] * 10**9,
                delay_rate[i] * 10**12, delay_rate_sigma[i] * 10**12, data_quality_flag[i], i + 1))
        #card 3
        cor = Correlation[i]
        Phas = Phase[i]
        SigPhas = PhaseSig[i]

        out.write(
            '{:10.5f}    .00000    .00000    .00000{:20.15f}{:10.1f}{:8d}03\n'.\
            format(cor, Phas, SigPhas, i + 1))
        #card 4
        out.write(
            '       .00   .0       .00   .0       .00   .0       .00   .0          '
            + '{:8d}'.format(i + 1) + '04\n')
        #card 5
        out.write('{:10.5f}{:10.5f}    .00000    .00000    .00000    .00000          {:8d}05\n'.format(10**9*CableCal[stations[n_stat1-1]][n_data1],\
                  10**9*CableCal[stations[n_stat2-1]][n_data2],i+1))
        #card 6

        out.write('{:10.3f}{:10.3f}{:10.3f}{:10.3f}{:10.3f}{:10.3f} 0 0      {:8d}06\n'.format(TempC[stations[n_stat1-1]][n_data1],\
                  TempC[stations[n_stat2-1]][n_data2],\
                  AtmPres[stations[n_stat1-1]][n_data1],AtmPres[stations[n_stat2-1]][n_data2],\
                  RelHum[stations[n_stat1-1]][n_data1]*100,RelHum[stations[n_stat2-1]][n_data2]*100,i+1))
        #card 8

        out.write('{:20.10f}{:10.5f}{:20.5f}{:10.5f}  0       {:8d}08\n'.format\
                  (delay_ion[i],  sigma_delay_ion[i],\
                   delay_ion_r[i],  sigma_delay_ion_r[i],i + 1)) 
#        #card 9
#        #        out.write('{:s}{:8d}09\n'.format(70*' ',i+1))
#        out.write(
#            '{:20.8f}{:10.5f}{:20.10f}{:10.5f} 0      I {:8d}09\n'.format(
#                delay[i] * 10**9, delay_sigma[i] * 10**9,
#                delay_rate[i] * 10**12, delay_rate_sigma[i] * 10**12, i + 1))

    out.close()

if __name__ == '__main__':
    main()
