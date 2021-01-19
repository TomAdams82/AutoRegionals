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
import pandas as pd
from scipy import interpolate
import matplotlib as mpl
import matplotlib.pyplot as plt
import pdb



def loadNWP(loadGrib,NWPfldr,NWPfn2,mdlrun,progDay,data_timestep,nwpfn,maplims):
    # load NWP
    runyr = int(mdlrun[:4])
    runmth = int(mdlrun[4:6])
    runday = int(mdlrun[6:8])
    runhr = int(mdlrun[8:10])
    
    runDT = datetime.datetime(runyr,runmth,runday,runhr)
    if loadGrib == False:
        prog = (progDay-1)*24 + 15 # To local time, whether 00Z or 12Z run
        progDT = runDT + datetime.timedelta(hours = prog)
        NWPfn = NWPfldr+'day'+str(progDay)+'.npy'
        TP, CC, HCC, U, V = np.load(NWPfn)
        lats, lons = np.load(NWPfn2)
        if runhr==0:
            daytxt =  (progDT+datetime.timedelta(days=1)).strftime( '%a %d %b') 
        else:        # 00Z run - will need to change if have 06Z runs too?
            daytxt =  progDT.strftime( '%a %d %b') 
    else:
        fields = ['Total precipitation','Convective precipitation','10 metre U wind component',
                  '10 metre V wind component','Low cloud cover','Medium cloud cover',
                  'High cloud cover'] 
        
        # Load a day's worth of NWP, in 3hr timesteps
        
        for i,p in enumerate(np.arange(0,27,data_timestep)): # start at 0, but record data from +3 or +6 - 15Z or 18Z
            nTsteps = (24.0/data_timestep)
            print('Loading data '+str(int(100*(i/(nTsteps+1))))+'%')
            if runhr == 12:
                prog = (progDay-1)*24 + p
            else:
                prog = (progDay-1)*24 + p + 12 # take to next complete day
            progDT = runDT + datetime.timedelta(hours = prog)
            if prog == 0: # Note there is a special analysis grib with land-sea mask etc.
                progfn = nwpfn+mdlrun+'/N1D'+runDT.strftime('%y%m%d%H0001')+progDT.strftime('%d%H011')        
            else:
                progfn = nwpfn+mdlrun+'/N1D'+runDT.strftime('%y%m%d%H0001')+progDT.strftime('%d%H001')
                       
            grbs = pygrib.open(progfn) 
            grbs.seek(0) 
            TPg = grbs.select(name=fields[0])[0]
            Ni = TPg.Ni # Make sure get same domain
            #CPg = grbs.select(name=fields[1])[0]
            Ug = grbs.select(name=fields[2],Ni=Ni)[0]
            Vg = grbs.select(name=fields[3],Ni=Ni)[0]
            LCCg = grbs.select(name=fields[4])[0]
            MCCg = grbs.select(name=fields[5])[0]
            HCCg = grbs.select(name=fields[6])[0]
                  
            TP_acc, lats, lons = TPg.data(lat1 = maplims[3], lat2 = maplims[2], lon1 = maplims[0], lon2 = maplims[1])
            #CP, lats, lons = CPg.data(lat1 = maplims[3], lat2 = maplims[2], lon1 = maplims[0], lon2 = maplims[1])
            if i>0:
                LCC, lats, lons = LCCg.data(lat1 = maplims[3], lat2 = maplims[2], lon1 = maplims[0], lon2 = maplims[1])
                MCC, lats, lons = MCCg.data(lat1 = maplims[3], lat2 = maplims[2], lon1 = maplims[0], lon2 = maplims[1]) # 0 to 1
                hCC, lats, lons = HCCg.data(lat1 = maplims[3], lat2 = maplims[2], lon1 = maplims[0], lon2 = maplims[1])
                Uraw, latsU, lonsU = Ug.data(lat1 = maplims[3], lat2 = maplims[2], lon1 = maplims[0], lon2 = maplims[1])
                Vraw, latsV, lonsV = Vg.data(lat1 = maplims[3], lat2 = maplims[2], lon1 = maplims[0], lon2 = maplims[1])
            grbs.close()
            if i == 0:
                TP = np.zeros([TP_acc.shape[0],TP_acc.shape[1],int(nTsteps)])
                CC = np.zeros([TP_acc.shape[0],TP_acc.shape[1],int(nTsteps)])
                HCC = np.zeros([TP_acc.shape[0],TP_acc.shape[1],int(nTsteps)])
                U = np.zeros([TP_acc.shape[0],TP_acc.shape[1],int(nTsteps)])
                V = np.zeros([TP_acc.shape[0],TP_acc.shape[1],int(nTsteps)])
                prevTP = TP_acc
            else:
                TP[:,:,i-1] = (TP_acc - prevTP) * (1000/data_timestep) # Convert m to mm per hour
                CC[:,:,i-1] = np.max(np.dstack((LCC,MCC)),axis=2) # take max of low of mid cloud
                HCC[:,:,i-1] = hCC 
                if Uraw.size != TP_acc.size: # need to reinterpolate winds - why does the data change resolution?
                    res1 = lons[0,1]-lons[0,0]
                    res2 = lonsU[0,1]-lonsU[0,0]
                    print('Prog '+str(prog)+' at different resolution for U - reinterpolating...')    
                    print('TP at '+str(res1)+'degree, U at '+str(res2)+'degree')
                    f = interpolate.interp2d(lonsU[0,:], latsU[:,0], Uraw, kind='linear')   
                    U[:,:,i-1] = np.flipud(f(lons[0,:],lats[:,0])) * 1.94384 # ms-1 to kt
                else:
                    U[:,:,i-1] = Uraw * 1.94384 # ms-1 to kt
                if Vraw.size != TP_acc.size: # need to reinterpolate winds - why does the data change resolution?
                    res1 = lons[0,1]-lons[0,0]
                    res2 = lonsV[0,1]-lonsV[0,0]
                    print('Prog '+str(prog)+' at different resolution for V - reinterpolating...')    
                    print('TP at '+str(res1)+'degree, V at '+str(res2)+'degree')
                    f = interpolate.interp2d(lonsV[0,:], latsV[:,0], Vraw, kind='linear')   
                    V[:,:,i-1] = np.flipud(f(lons[0,:],lats[:,0])) * 1.94384 # ms-1 to kt
                else:                
                    V[:,:,i-1] = Vraw * 1.94384 # ms-1 to kt

        daytxt =  progDT.strftime( '%a %d %b')    
    # Save each days NWP to a new file
    NWPfn = NWPfldr+'day'+str(progDay)+'.npy'
    np.save(NWPfn,[TP,CC,HCC,U,V])
    # Dont save lat longs - done at start
    # np.save(NWPfn2,[lats,lons])
    return TP,CC,HCC,U,V,daytxt
    # Define temporal density functions

def loadTDF(dataFn,data_timestep):   # 1 is X, 0 is allowing Y 
    dfm = pd.read_excel(dataFn, sheetname='Time Functions', skiprows = 1) # dataframe
    nTDF = dfm['Number of types'][0]
    TDF = np.zeros([nTDF,24/data_timestep]).astype('bool')
    TDF_text = np.zeros(nTDF).astype('U64')
    if data_timestep==3:
        dfm = pd.read_excel(dataFn, sheetname='Time Functions', skiprows = 8) # skip to data
    if data_timestep==6:        
        dfm = pd.read_excel(dataFn, sheetname='Time Functions', skiprows = (10+nTDF)) # skip to data
    for i in range(nTDF):
        d = dfm.loc[dfm['Number'] == i].values[0]
        TDF_text[i] = dfm['Text string'][i]
        if data_timestep==3:
            TDF[i,:] = d[3:11].astype('bool')    
        if data_timestep==6:
            TDF[i,:] = d[3:7].astype('bool')   
    
#    fig = plt.figure(figsize=(20, 20), facecolor='w', edgecolor='w')
#    gs = mpl.gridspec.GridSpec(3, 2)
#    x = np.arange(8)
#    utc_labels = np.concatenate([(np.arange(15,24,3)),(np.arange(0,15,3))]).astype(str)
#    for i in range(5):
#        ax=plt.subplot(gs[i])    
#        ax.plot(x,TDF[i,:],'b-')
#        ax.set_title(TDF_text[i],fontsize=14)
#        ax.set_xticklabels(utc_labels)
#        ax.set_xlabel('UTC time')
#        ax.set_yticks([0,1])
#        ax.set_ylim([-.2,1.2])
#        ax.set_yticklabels(['YYY','XXX'])
#
#    fig.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.2, hspace=0.5)
    return nTDF,TDF,TDF_text
    
    # Define phenomena density functions
def loadDFs(dataFn):    
    # Number of weather types
    dfm = pd.read_excel(dataFn, sheetname='Weather Types', skiprows = 1) # dataframe
    nDF = dfm['Number of types'][0]
    
    # Load distribution bins
    dfm = pd.read_excel(dataFn, sheetname='Weather Types', skiprows = 3) # dataframe
    d = dfm.loc[dfm['Distribution bins'] == 0].values[0]
    TP_bins = d[2:9]
    d = dfm.loc[dfm['Distribution bins'] == 1].values[0]
    CC_bins = d[2:6]/100.0    
    # Load neames, icons and distributions
    DF_unnorm_TP = np.zeros([nDF,len(TP_bins)-1])
    DF_unnorm_CC = np.zeros([nDF,len(CC_bins)-1])
    DF_unnorm_HCC = np.zeros([nDF,len(CC_bins)-1])
    DF_TP = np.zeros([nDF,len(TP_bins)-1])
    DF_CC = np.zeros([nDF,len(CC_bins)-1])
    DF_HCC = np.zeros([nDF,len(CC_bins)-1])
    DF_text = np.zeros(nDF).astype('U64')
    DF_icon = np.zeros(nDF).astype('str')
    DF_priority = np.zeros([nDF,3]) # weighting priority for differnt inputs, applied to KL-div scores
    banned_list = np.zeros([nDF,nDF]).astype('bool')
    
    dfm = pd.read_excel(dataFn, sheetname='Weather Types', skiprows = 8) # skip to data
    for i in range(nDF):
        d = dfm.loc[dfm['Number'] == i].values[0]
        DF_unnorm_TP[i,:] = d[3:9].astype('float')
        DF_text[i] = dfm['Text string'][i]
        DF_icon[i] = dfm['Icon'][i]
    # Cloud (LCC + MCC)
    dfm = pd.read_excel(dataFn, sheetname='Weather Types', skiprows = 10+nDF) # skip to data
    for i in range(nDF):
        d = dfm.loc[dfm['Number'] == i].values[0]
        DF_unnorm_CC[i,:] = d[3:6].astype('float')
    # HCC
    dfm = pd.read_excel(dataFn, sheetname='Weather Types', skiprows = 12+(2*nDF)) # skip to data
    for i in range(nDF):
        d = dfm.loc[dfm['Number'] == i].values[0]
        DF_unnorm_HCC[i,:] = d[3:6].astype('float')
    # Parameter priorities     # Set parameter priorities for TP, CC and HCC
    dfm = pd.read_excel(dataFn, sheetname='Weather Types', skiprows = 14+(3*nDF)) # skip to data
    for i in range(nDF):
        d = dfm.loc[dfm['Number'] == i].values[0]
        DF_priority[i,:] = d[3:6].astype('float')

    # Normalise
    for i in range(nDF):
        DF_TP[i,:] = DF_unnorm_TP[i,:] / np.sum(DF_unnorm_TP[i,:])
        DF_CC[i,:] = DF_unnorm_CC[i,:] / np.sum(DF_unnorm_CC[i,:])
        DF_HCC[i,:] = DF_unnorm_HCC[i,:] / np.sum(DF_unnorm_HCC[i,:])
    
    # Load banned list
    dfm = pd.read_excel(dataFn, sheetname='Banned List', skiprows = 8) # skip to data
    for i in range(nDF):
        d = dfm.loc[dfm['Number'] == i].values[0]
        banned_list[i,:] = d[3:(3+nDF)].astype('bool')

    #
#    fig = plt.figure(figsize=(20, 20), facecolor='w', edgecolor='w')
#    gs = mpl.gridspec.GridSpec(3, 5)
#    nlabels = len(CC_bins)-1  
#    precip_labels = np.zeros(nlabels).astype(str)
##    precip_labels[0]=''
##    precip_labels[-1]=''
#    for i in range(nlabels):
#        precip_labels[i] = str(TP_bins[i])+'-'+str(TP_bins[i+1])
#        
#    x = np.arange(nlabels)
#    for i in range(15):
#        ax=plt.subplot(gs[i])    
#        ax.plot(x,DF_CC[i,:],'b-')
#        ax.set_title(DF_text[i],fontsize=14)
#        ax.set_xticks(x)
#        ax.set_xticklabels(precip_labels, rotation = 'vertical')
#        ax.set_xlim([-1,nlabels])
#        ax.set_xlabel('Precip amount / mm')
#        ax.set_ylim([0,1])
#        ax.set_yticks([0,1])
#    fig.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.2, hspace=0.5)
    return nDF,DF_TP,DF_CC,DF_HCC,DF_priority,DF_text,DF_icon,banned_list,TP_bins,CC_bins


    # Define phenomena density functions
def load_windDFs(dataFn):    
    # Number of weather types
    dfm = pd.read_excel(dataFn, sheetname='Wind Types', skiprows = 1) # dataframe
    nWDF = dfm['Number of types'][0]
    
    # Load distribution bins
    dfm = pd.read_excel(dataFn, sheetname='Wind Types', skiprows = 3) # dataframe
    d = dfm.loc[dfm['Distribution bins'] == 0].values[0]
    W_bins = d[2:17]
    # Load neames, icons and distributions
    DF_unnorm_U = np.zeros([nWDF,len(W_bins)-1])
    DF_unnorm_V = np.zeros([nWDF,len(W_bins)-1])
    DF_U = np.zeros([nWDF,len(W_bins)-1])
    DF_V = np.zeros([nWDF,len(W_bins)-1])
    DF_wind_text = np.zeros(nWDF).astype('U64')
    DF_wind_icon = np.zeros(nWDF).astype('str')
    DF_wind_priority = np.zeros([nWDF,2]) # weighting priority for differnt inputs, applied to KL-div scores
    banned_list_winds = np.zeros([nWDF,nWDF]).astype('bool')
    # U
    dfm = pd.read_excel(dataFn, sheetname='Wind Types', skiprows = 8) # skip to data
    for i in range(nWDF):
        d = dfm.loc[dfm['Number'] == i].values[0]
        DF_unnorm_U[i,:] = d[3:(2+len(W_bins))].astype('float')
        DF_wind_text[i] = dfm['Text string'][i]
        DF_wind_icon[i] = dfm['Icon'][i]
    # V
    dfm = pd.read_excel(dataFn, sheetname='Wind Types', skiprows = 10+nWDF) # skip to data
    for i in range(nWDF):
        d = dfm.loc[dfm['Number'] == i].values[0]
        DF_unnorm_V[i,:] = d[3:(2+len(W_bins))].astype('float')

    # Parameter priorities     # Set parameter priorities for TP, CC and HCC
    dfm = pd.read_excel(dataFn, sheetname='Wind Types', skiprows = 12+(2*nWDF)) # skip to data
    for i in range(nWDF):
        d = dfm.loc[dfm['Number'] == i].values[0]
        DF_wind_priority[i,:] = d[3:5].astype('float')
    # Normalise
    for i in range(nWDF):
        DF_U[i,:] = DF_unnorm_U[i,:] / np.sum(DF_unnorm_U[i,:])
        DF_V[i,:] = DF_unnorm_V[i,:] / np.sum(DF_unnorm_V[i,:])
    
    # Load banned list
    dfm = pd.read_excel(dataFn, sheetname='Banned Wind List', skiprows = 8) # skip to data
    for i in range(nWDF):
        d = dfm.loc[dfm['Number'] == i].values[0]
        banned_list_winds[i,:] = d[3:(3+nWDF)].astype('bool')

#    
#    fig = plt.figure(figsize=(20, 20), facecolor='w', edgecolor='w')
#    gs = mpl.gridspec.GridSpec(7, 8)
#    nlabels = len(W_bins)-1  
#    precip_labels = np.zeros(nlabels).astype(str)
##
#    for i in range(nlabels):
#        precip_labels[i] = str(W_bins[i])+'-'+str(W_bins[i+1])
#        
#    x = np.arange(nlabels)
#    for i in range(49):
#        ax=plt.subplot(gs[i+7])    
#        ax.plot(x,DF_U[i,:],'b-')
#        ax.plot(x,DF_V[i,:],'r-')
#        ax.set_title(DF_wind_text[i],fontsize=10)
#        ax.set_xticks(x)
#        ax.set_xticklabels(precip_labels, rotation = 'vertical',fontsize=6)
#        ax.set_xlim([-1,nlabels])
#        ax.set_ylim([0,1])
#        ax.set_yticks([0,1])
#    fig.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.2, hspace=0.5)
#    fig.suptitle('U(blue), V(red) wind component kt')
##    pdb.set_trace()
#    
    return nWDF,DF_U,DF_V,DF_wind_priority,DF_wind_text,DF_wind_icon,banned_list_winds,W_bins

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

def loadLatLongs(loadGrib,NWPfn2,mdlrun,nwpfn,maplims):
    # load NWP
    if loadGrib == False:
        lats, lons = np.load(NWPfn2)
    else:
        runyr = int(mdlrun[:4])
        runmth = int(mdlrun[4:6])
        runday = int(mdlrun[6:8])
        runhr = int(mdlrun[8:10])        
        runDT = datetime.datetime(runyr,runmth,runday,runhr)
        fields = ['Total precipitation'] 
        progDay = 2# why are winds different for day 1? day 1 is  + 0 day
        # Load a single timestep of NWP for lats and longs 
        p = 0 # prog
        prog = (progDay-1)*24 + p
        progDT = runDT + datetime.timedelta(hours = prog)
        if prog == 0: # Note there is a special analysis grib with land-sea mask etc.
            progfn = nwpfn+mdlrun+'/N1D'+runDT.strftime('%y%m%d%H0001')+progDT.strftime('%d%H011')        
        else:
            progfn = nwpfn+mdlrun+'/N1D'+runDT.strftime('%y%m%d%H0001')+progDT.strftime('%d%H001')
       
        grbs = pygrib.open(progfn) 
        grbs.seek(0) 
        TPg = grbs.select(name=fields[0])[0]  
        TP_acc, lats, lons = TPg.data(lat1 = maplims[3], lat2 = maplims[2], lon1 = maplims[0], lon2 = maplims[1])        
        grbs.close()
        
        # Save lats and longs
        np.save(NWPfn2,[lats,lons])
    return lats,lons
    # Define temporal density functions
    
def loadTextChanges(dataFn):   # 1 is X, 0 is allowing Y 
    dfm = pd.read_excel(dataFn, sheetname='Text rules', skiprows = 1) # dataframe
    nRules = dfm['Number of types'][0]
    ruleCombos = np.zeros([nRules,3]).astype('int')
    ruleReplaces = np.zeros([nRules,2]).astype('U64')
    dfm = pd.read_excel(dataFn, sheetname='Text rules', skiprows = 8) # skip to data
    ruleCombos[:,0] = dfm['Type 1'].values
    ruleCombos[:,1] = dfm['Type 2'].values
    ruleCombos[:,2] = dfm['1st or 2nd occurrence'].values
    ruleReplaces[:,0] = dfm['Replace'].values
    ruleReplaces[:,1] = dfm['With'].values

    return ruleCombos, ruleReplaces 