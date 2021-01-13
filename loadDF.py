#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 19:12:20 2021

@author: tadams@met.co.nz
"""

import shapefile
import numpy as np
import matplotlib.path as mpltPath
import datetime
import pygrib
import sys


#import matplotlib as mpl
#import matplotlib.pyplot as plt
#import pdb



def loadNWP(loadGrib,NWPfn,NWPfn2,mdlrun,progDay,data_timestep,nwpfn,maplims):
    # load NWP
    if loadGrib == False:
        TP = np.load(NWPfn)
        lats, lons = np.load(NWPfn2)
        daytxt =  'Wed 13 Jan'#progDT.strftime( '%a %d %b')
    else:
            
        runyr = int(mdlrun[:4])
        runmth = int(mdlrun[4:6])
        runday = int(mdlrun[6:8])
        runhr = int(mdlrun[8:10])
        
        runDT = datetime.datetime(runyr,runmth,runday,runhr)
        
        
        fields = ['Total precipitation','Convective precipitation','10 metre U wind component',
                  '10 metre V wind component','Low cloud cover','Medium cloud cover',
                  'High cloud cover'] 
        
        # Load a day's worth of NWP, in 3hr timesteps
        
        for i,p in enumerate(np.arange(0,27,data_timestep)):
            ntsteps = (24.0/data_timestep)+1
            print('Loading data '+str(int(100*(i/ntsteps)))+'%')
            prog = (progDay-1)*24 + p
            progDT = runDT + datetime.timedelta(hours = prog)
            if prog == 0: # Note there is a special analysis grib with land-sea mask etc.
                progfn = nwpfn+mdlrun+'/N1D'+runDT.strftime('%y%m%d%H0001')+progDT.strftime('%d%H011')        
            else:
                progfn = nwpfn+mdlrun+'/N1D'+runDT.strftime('%y%m%d%H0001')+progDT.strftime('%d%H001')
                       
            grbs = pygrib.open(progfn) 
            grbs.seek(0) 
            TPg = grbs.select(name=fields[0])[0]
            #CPg = grbs.select(name=fields[1])[0]
            #Ug = grbs.select(name=fields[2])[0]
            #Vg = grbs.select(name=fields[3])[0]
            #LCCg = grbs.select(name=fields[4])[0]
            #MCCg = grbs.select(name=fields[5])[0]
            #HCCg = grbs.select(name=fields[6])[0]
            #
            #for g in grbs:
            #    print(g.name)
                  
            TP_acc, lats, lons = TPg.data(lat1 = maplims[3], lat2 = maplims[2], lon1 = maplims[0], lon2 = maplims[1])
            #CP, lats, lons = CPg.data(lat1 = maplims[3], lat2 = maplims[2], lon1 = maplims[0], lon2 = maplims[1])
            #U, lats, lons = Ug.data(lat1 = maplims[3], lat2 = maplims[2], lon1 = maplims[0], lon2 = maplims[1])
            #V, lats, lons = Vg.data(lat1 = maplims[3], lat2 = maplims[2], lon1 = maplims[0], lon2 = maplims[1])
            #LCC, lats, lons = LCCg.data(lat1 = maplims[3], lat2 = maplims[2], lon1 = maplims[0], lon2 = maplims[1])
            #MCC, lats, lons = MCCg.data(lat1 = maplims[3], lat2 = maplims[2], lon1 = maplims[0], lon2 = maplims[1])
            #HCC, lats, lons = HCCg.data(lat1 = maplims[3], lat2 = maplims[2], lon1 = maplims[0], lon2 = maplims[1])
            
            grbs.close()
            if i == 0:
                TP = np.zeros([TP_acc.shape[0],TP_acc.shape[1],8])
                prevTP = TP_acc
            else:
                TP[:,:,i-1] = (TP_acc - prevTP) * (1000/data_timestep) # Convert m to mm per hour
        daytxt =  progDT.strftime( '%a %d %b')    

    np.save(NWPfn,TP)
    #np.save(NWPfn2,[lats,lons])
    return lats,lons,TP,daytxt
    # Define temporal density functions

def loadTDF(data_timestep):   # 1 is X, 0 is allowing Y 
    
    nTDF = 6 # Number of functions
    TDF = np.zeros([nTDF,24/data_timestep]).astype('bool')
    TDF_text = np.zeros(nTDF).astype('U64')
    # UTC hours for TDF = 15,18,21,0,3,6,9,12
    TDF_text[0] = 'XXX' # constant
    TDF_text[1] = 'XXX but YYY' # constant with spatial split
    TDF_text[2] = 'XXX becoming YYY for a time in the afternoon and evening' # Replace XXX and YYY with phenomena
    TDF_text[3] = 'XXX, but becoming YYY in the morning'
    TDF_text[4] = 'XXX in the morning, then YYY from afternoon'
    TDF_text[5] = 'XXX, then YYY from evening' 
    # Values are 1 is XXX, 0 is YYY # Must be boolean if using spatial subregions
    #                15 ,18 ,21 ,  0, 3,  6,  9, 12
    if data_timestep==3:
        TDF[0,:] = [ 1,  1,  1,   1, 1,  1,  1, 1]
        TDF[1,:] = [ 0,  0,  0,   0, 0,  0,  0, 0]
        TDF[2,:] = [ 1,  1,  1,   1, 0,  0,  1, 1]
        TDF[3,:] = [ 1,  1,  0,   0, 0,  0,  0, 0]
        TDF[4,:] = [ 1,  1,  1,   1, 1,  0,  0, 0]
        TDF[5,:] = [ 1,  1,  1,   1, 1,  1,  0, 0]
    if data_timestep==6:
    #              18 ,  0,  6,  12
        TDF[0,:] = [ 1,  1,  1,   1]
        TDF[1,:] = [ 0,  0,  0,   0]
        TDF[2,:] = [ 1,  1,  0,   1]
        TDF[3,:] = [ 1,  0,  0,   0]
        TDF[4,:] = [ 1,  1,  0,   0]
        TDF[5,:] = [ 1,  1,  1,   0]
    
    #fig = plt.figure(figsize=(20, 20), facecolor='w', edgecolor='w')
    #gs = mpl.gridspec.GridSpec(3, 2)
    #x = np.arange(8)
    #utc_labels = np.concatenate([(np.arange(15,24,3)),(np.arange(0,15,3))]).astype(str)
    #for i in range(5):
    #    ax=plt.subplot(gs[i])    
    #    ax.plot(x,TDF[i,:],'b-')
    #    ax.set_title(TDF_text[i],fontsize=14)
    #    ax.set_xticklabels(utc_labels)
    #    ax.set_xlabel('UTC time')
    #    ax.set_ylim([0,1])
    #    ax.set_yticks([0,1])
    #    ax.set_yticklabels(['XXX','YYY'])
    #
    #fig.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.2, hspace=0.5)
    return nTDF,TDF,TDF_text
    
    # Define phenomena density functions
def loadDFs():    
    nDF = 10 # Number of phenomena 
    nParams = 1# Number of parameters considered
    TP_bins = [0  ,0.1, 0.5,  1,  5, 10, 100] # precip bins in timestep, average mm /hr 
    DF_unnorm = np.zeros([nDF,len(TP_bins)-1,nParams])
    DF = np.zeros([nDF,len(TP_bins)-1,nParams])
    DF_text = np.zeros(nDF).astype('U64')
    # UTC hours for TDF = 15,18,21,0,3,6,9,12
    DF_text[0] = 'fine'
    DF_text[1] = 'patchy drizzle'
    DF_text[2] = 'drizzle'
    DF_text[3] = 'isolated showers'    
    DF_text[4] = 'a few showers'
    DF_text[5] = 'scattered rain'
    DF_text[6] = 'showers'
    DF_text[7] = 'showers, some heavy'
    DF_text[8] = 'rain'
    DF_text[9] = 'rain with heavy falls'
    
    
    #TP_bins =         [0.1 , 0.5,   1,   5,  10, 100] # upper limits mm / hr (ignore 0)
    DF_unnorm[0,:,0] = [ 0.9, 0.1,1e-4,1e-5,1e-5,1e-5] # Fine
    DF_unnorm[1,:,0] = [ 0.6, 0.2,1e-3,1e-3,1e-5,1e-5] # Patchy drizzle
    DF_unnorm[2,:,0] = [ 0.2, 0.6, 0.1,1e-3,1e-5,1e-5] # Drizzle
    DF_unnorm[3,:,0] = [ 0.6, 0.1, 0.2, 0.1,1e-3,1e-5] # Isolated showers
    DF_unnorm[4,:,0] = [0.25,0.25, 0.3, 0.2,1e-3,1e-5] # A few showers
    DF_unnorm[5,:,0] = [   0, 0.4, 0.4, 0.2,0.01,1e-3] # scattered rain
    DF_unnorm[6,:,0] = [ 0.2, 0.2, 0.2, 0.3, 0.1,1e-3] # showers
    DF_unnorm[7,:,0] = [ 0.1, 0.1, 0.1, 0.4, 0.2, 0.1] # showers, some heavy
    DF_unnorm[8,:,0] = [1e-3, 0.1, 0.3, 0.5,1e-3,1e-3] # rain
    DF_unnorm[9,:,0] = [1e-3,1e-3, 0.2, 0.5, 0.3, 0.1] # rain with heavy falls
    
    # Normalise
    for i in range(nDF):
        DF[i,:,0] = DF_unnorm[i,:,0] / np.sum(DF_unnorm[i,:,0])
      
    banned_list_precip = {0:[0],1:[1,2,3],2:[1,2],3:[3,4,5],4:[3,4,5],5:[3,4,5],6:[6,8],7:[7,9],8:[6,8],9:[7,9]}
    #e.g. 4 (a few showers) is banned with a few showers, and also scattered rain
    #
#    fig = plt.figure(figsize=(20, 20), facecolor='w', edgecolor='w')
#    gs = mpl.gridspec.GridSpec(3, 3)
#    nlabels = len(TP_bins)-1  
#    precip_labels = np.zeros(nlabels+2).astype(str)
#    precip_labels[0]=''
#    precip_labels[-1]=''
#    for i in range(nlabels):
#        precip_labels[i+1] = str(TP_bins[i])+'-'+str(TP_bins[i+1])
#        
#    x = np.arange(nlabels)
#    for i in range(9):
#        ax=plt.subplot(gs[i])    
#        ax.plot(x,DF[i,:,0],'b-')
#        ax.set_title(DF_text[i],fontsize=14)
#        ax.set_yticks(x)
#        ax.set_xlim([-1,6])
#        ax.set_xticklabels(precip_labels)
#        ax.set_xlabel('Precip amount / mm')
#        ax.set_ylim([0,1])
#        ax.set_yticks([0,1])
#    fig.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.2, hspace=0.5)
#    pdb.set_trace()
    return nDF,DF,DF_text,banned_list_precip,TP_bins
    # Load whole region
def loadSDF(lats,lons,shpfn):
    pts = np.array(zip(list(lons.flatten()), list(lats.flatten()))).astype('float')
    sf = shapefile.Reader(shpfn)
    
    foundAll = False
    c = -1
    for shape in sf.shapeRecords():
        c +=1
        if shape.record[0] == 'All':
            foundAll = True
            if len(shape.shape.parts)==1:                           # Single part polygon
                i_start = shape.shape.parts[0]
                i_end = len(shape.shape.points)
                polygon = shape.shape.points[i_start:i_end]
                path = mpltPath.Path(polygon)
                inside = path.contains_points(pts)
                inArea = np.reshape(inside,lats.shape)
            else:                                                   # Deal with multiple parts to polygon one by one
                inSubTemp = np.zeros(lats.shape)
                for p in range(len(shape.shape.parts)):
                    i_start = shape.shape.parts[p]
                    if p<(len(shape.shape.parts)-1):
                        i_end = shape.shape.parts[p+1]
                    else:
                        i_end=len(shape.shape.points)
                    polygon = shape.shape.points[i_start:i_end]
                    path = mpltPath.Path(polygon)
                    inside = path.contains_points(pts)
                    inSubTemp += np.reshape(inside,lats.shape)
                inArea = inSubTemp>0
    
    if foundAll == False:
        print('No All region - check shapefile')
        sys.exit()
        
    # Load subregions
    nSDF = int(0) # number of subregions
    for shape in sf.shapeRecords():
        nSDF+=1
        if shape.record[1]!='None': # Has an inverse
            nSDF+=1
            
    SDF = np.zeros([lats.shape[0],lats.shape[1],nSDF]).astype(bool) # Dimensions lat, lon, subregion
    SDF_text = np.zeros(nSDF).astype('str')
    
    c = -1    
    for shape in sf.shapeRecords():
        c+=1
        SDF_text[c] = shape.record[0]
        if SDF_text[c] == 'All':
            wholeArea = c
        inverse = shape.record[1]
        if len(shape.shape.parts)==1:                           # Single part polygon
            i_start = shape.shape.parts[0]
            i_end = len(shape.shape.points)
            polygon = shape.shape.points[i_start:i_end]
            path = mpltPath.Path(polygon)
            inside = path.contains_points(pts)
            SDF[:,:,c] = np.reshape(inside,lats.shape)
        else:                                                   # Deal with multiple parts to polygon one by one
            inSubTemp = np.zeros(lats.shape)
            for p in range(len(shape.shape.parts)):
                i_start = shape.shape.parts[p]
                if p<(len(shape.shape.parts)-1):
                    i_end = shape.shape.parts[p+1]
                else:
                    i_end=len(shape.shape.points)
                polygon = shape.shape.points[i_start:i_end]
                path = mpltPath.Path(polygon)
                inside = path.contains_points(pts)
                inSubTemp += np.reshape(inside,lats.shape)
            SDF[:,:,c] = inSubTemp>0
        print('Subarea '+str(c)+' '+ SDF_text[c]+' '+str(sum(sum(SDF[:,:,c])))+' points')
        if inverse != 'None': # Do inverse?
            c += 1
            SDF_text[c] = inverse
            SDF[:,:,c] = inArea * np.logical_not(SDF[:,:,c-1]) # Must be in main area, then inverse of subarea
            print('Subarea '+str(c)+' '+ SDF_text[c]+' '+str(sum(sum(SDF[:,:,c])))+' points')
    return(nSDF,SDF,SDF_text,inArea,wholeArea,sf)