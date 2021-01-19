#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 11:09:39 2020

@author: sparky
"""

import numpy as np
#import matplotlib as mpl
#mpl.use('AGG')
import matplotlib.pyplot as plt
import datetime
#from matplotlib import gridspec
#import matplotlib.path as mpltPath
#import matplotlib.patches as mpatches
#from mpl_toolkits.basemap import Basemap
import pdb
from plotTimeSeriesMaps_clouds import plotTimeSeriesMaps, plotTimeSeriesDFs,compareTimeSeriesDFs
import loadDF_winds as loadDF
import DFanalysis_clouds as DFanalysis
import os
import sys
import time

loadGrib = True # Only applies to first subregion anyway

rootfldr = os.path.abspath(os.path.dirname(sys.argv[0]))+'/'
shpfldr = rootfldr+'data/shapefiles/'
optfldr = rootfldr
nwpfn = '/data/models/ecmwf-prod/det/'
regionlist = ['NLD','ALD','COP','WKO','WTO','BOP','ROT','TPO','TMN','THP','GSB','HBY','WRP','TKI','WGI','MNU','KAP','WGN',
              'NSN','BLR','WLD','MRB','CH','CPL','NOTA','DUN','COTA','CLU','SLD','CNY_HI','SN_LKS']
             # Areas to do
tdy = datetime.datetime.utcnow()
if tdy.hour<9:
    tdy = tdy - datetime.timedelta(days=1)
    runhr = '12'
    maxdays = 10
else:
    if tdy.hour<21:
        runhr = '00'
        maxdays = 9 # do 9 days on 00Z forecast
    else:
        runhr = '12'
        maxdays = 10
        
mdlrun = tdy.strftime('%Y%m%d')+runhr
print('Processing Run '+mdlrun)
#mdlrun = '2021011412'
NWPfldr = rootfldr+'NWPfolder/' # presaved NWP folder
NWPfn2 = rootfldr+'NWPfolder/NWPlatlongs.npy' # presaved lats and longs for domain
optfn = rootfldr+'Forecasts/Forecast_'+mdlrun+'Z.txt'
dataFn = rootfldr+'data/DFdata.xlsx' # Where DF data saved
nParams = 3 # TP, CC, HCC
maplims = np.array([165.,180.,-33.,-48.]) # Total domain of NWP to load from
plt.close('all')
Cut_off_dont_try = 2 # If any type is above this all timesteps, don't try
startTick = time.time()
fout = open(optfn,'w') # Open file on first day 
fout.writelines('Model run '+mdlrun+'Z\n')
# Load lats and longs - SOME winds are differnet resolution
lats,lons = loadDF.loadLatLongs(loadGrib,NWPfn2,mdlrun,nwpfn,maplims) # Assume lats and longs are same through model run - res doesn't chane
# If res does change, just get lats and longs returned in loadNWP, and move spatial bit inside main time loop
# load in DFs
nDF,DF_TP,DF_CC,DF_HCC,DF_priority,DF_text,DF_icons,banned_list, TP_bins,CC_bins = loadDF.loadDFs(dataFn)
# Load winds DFs
nWDF,DF_U,DF_V,DF_winds_priority,DF_winds_text,DF_winds_icons,banned_winds_list, W_bins = loadDF.load_windDFs(dataFn)
# Load in text change rules
ruleCombos, ruleReplaces = loadDF.loadTextChanges(dataFn)
for region in regionlist:
    print('Processing '+region)
    if region != regionlist[0]: # on second region
        loadGrib = False # Don't reload NWP second time round
    fout.writelines('\n\n'+region+' forecast\n')
    shpfn = shpfldr + region+ 'subregions.shp'
    nSDF,SDF,SDF_text,inArea,wholeArea,sf = loadDF.loadSDF(lats,lons,shpfn)
    
    for progDay in range(1,maxdays+1): # 1 is day 1, 10 is day 10
        if progDay <(maxdays-3): # do to day 6 for 12Z, day 5 for 00Z in 3hr data
            data_timestep = 3 # 3 hour data
        else:
            data_timestep = 6 # 6 hour data
        TP,CC,HCC,U,V,daytxt = loadDF.loadNWP(loadGrib,NWPfldr,NWPfn2,mdlrun,progDay,data_timestep,nwpfn,maplims)
        nTDF,TDF,TDF_text = loadDF.loadTDF(dataFn,data_timestep)
        nTsteps = 24/data_timestep  
        print('Loaded data for ' + daytxt)    
        # Main loop - create XXX and YYY vs SSS vs t matrices
        XXXmatrix = np.zeros([nDF,nSDF,nTsteps])+np.inf # score of each possible XXX vs SSS at every t
        YYYmatrix = np.zeros([nDF,nSDF,nTsteps])+np.inf # same for YYY
        
        for XXXn in range(nDF):
            #print(str(XXXn)+'/'+str(nDF))
            for SSSn in range(nSDF):
                for t in range(nTsteps):
                    isWholeArea = SSSn == wholeArea
                    XXXmatrix[XXXn,SSSn,t] = DFanalysis.simpleScoreDF(False,TP[:,:,t],CC[:,:,t],HCC[:,:,t],DF_TP[XXXn,:],DF_CC[XXXn,:],DF_HCC[XXXn,:],DF_priority[XXXn,:],SDF[:,:,SSSn],isWholeArea,inArea,TP_bins,CC_bins,'X')
                    YYYmatrix[XXXn,SSSn,t] = DFanalysis.simpleScoreDF(False,TP[:,:,t],CC[:,:,t],HCC[:,:,t],DF_TP[XXXn,:],DF_CC[XXXn,:],DF_HCC[XXXn,:],DF_priority[XXXn,:],SDF[:,:,SSSn],isWholeArea,inArea,TP_bins,CC_bins,'Y')
    
        # Now loop through and add up scores according to TDF        
        weightingFnc = np.ones(nTsteps) # Set weighting for text
        if data_timestep == 3:
            weightingFnc[0]=0.5 # 15Z
            weightingFnc[-1]=0.5 # 12Z
        if data_timestep == 6:
            weightingFnc[-1]=0.5 # 12Z
            
        score_matrix = np.zeros([nDF,nDF,nTDF,nSDF])+np.inf
    
        for XXXn in range(nDF):
            #print(str(XXXn)+'/'+str(nDF))
            # Add drop out if none of NWP looks like this type, i.e. all day just this type has high min score
            min_score = np.min(XXXmatrix[XXXn,wholeArea,:])
            if min_score>Cut_off_dont_try:
                dmy =0
                #print('Omitting '+DF_text[XXXn]+' as primary due to minimum score of '+str(min_score))
            else:
                for TTTn in range(nTDF):
                    #if TTTn>0: # if TTTn == 0 there is no YYY, XXX is constant
                    for YYYn in range(nDF):
                        # Catch banned combinations e.g. showers and showers, or scattered rain and a few showers
                        if banned_list[XXXn,YYYn] == True: # 
                            for SSSn in range(nSDF):
                                score = np.zeros(nTsteps)
                                for t in range(nTsteps):
                                    if TDF[TTTn,t]==1:
                                        score[t] = weightingFnc[t] * XXXmatrix[XXXn,wholeArea,t]
                                    else:
                                        if SSSn == wholeArea:
                                            if TDF_text[TTTn] == 'XXX but YYY': # prevent Rain but showers with no spatial qualifier type scenarios
                                                score[t] = np.inf
                                            else:
                                                score[t] = weightingFnc[t] * YYYmatrix[YYYn,wholeArea,t]    
                                        else:
                                            scoreX = weightingFnc[t] * XXXmatrix[XXXn,SSSn,t]
                                            scoreY = weightingFnc[t] * YYYmatrix[YYYn,SSSn,t]
                                            score[t] = np.mean([scoreX,scoreY])
    
                                score_matrix[XXXn,YYYn,TTTn,SSSn] = np.mean(score) 
            #                    txt = getText(DF_text,TDF_text,SDF_text,XXXn,YYYn,TTTn,SSSn,ruleCombos, ruleReplaces)
            #                    print(txt+' score '+ str(np.mean(score)))                                          
                        else:
                            score_matrix[XXXn,YYYn,TTTn,:] = np.inf
    
        
        # Show best solution
        sortedScores = np.sort(score_matrix.flatten())
    #    for i in range(50):
    #        XXn,YYn,TTn,SSn = np.where(score_matrix==sortedScores[i])
    #        XXXn=XXn[0]
    #        YYYn=YYn[0]
    #        TTTn=TTn[0]
    #        SSSn=SSn[0]
    #        
    #        txt = DFanalysis.getText(DF_text,TDF_text,SDF_text,XXXn,YYYn,TTTn,SSSn,ruleCombos, ruleReplaces)
    #        print('Rank '+str(i)+': '+txt+' score '+ str(score_matrix[XXXn,YYYn,TTTn,SSSn]))
        # Take best
        XXn,YYn,TTn,SSn = np.where(score_matrix==sortedScores[0])
        XXXn=XXn[0] 
        YYYn=YYn[0] 
        TTTn=TTn[0]
        SSSn=SSn[0]
        
        # Recalculate NWP DF and guess DF as time series
        nNwp = np.sum(inArea).astype('float')
        nwpDF_XXX_TP = np.zeros([len(TP_bins)-1,nTsteps])
        nwpDF_XXX_CC = np.zeros([len(CC_bins)-1,nTsteps])
        nwpDF_XXX_HCC = np.zeros([len(CC_bins)-1,nTsteps])
        nwpDF_YYY_TP = np.zeros([len(TP_bins)-1,nTsteps])
        nwpDF_YYY_CC = np.zeros([len(CC_bins)-1,nTsteps])
        nwpDF_YYY_HCC = np.zeros([len(CC_bins)-1,nTsteps])
    
        guessDF_XXX_TP = np.zeros([len(TP_bins)-1,nTsteps])
        guessDF_XXX_CC = np.zeros([len(CC_bins)-1,nTsteps])
        guessDF_XXX_HCC = np.zeros([len(CC_bins)-1,nTsteps])
        guessDF_YYY_TP = np.zeros([len(TP_bins)-1,nTsteps])
        guessDF_YYY_CC = np.zeros([len(CC_bins)-1,nTsteps])
        guessDF_YYY_HCC = np.zeros([len(CC_bins)-1,nTsteps])
    
        scoreXXX = np.zeros([nTsteps,4])
        scoreYYY = np.zeros([nTsteps,4])
        score = np.zeros([nTsteps,4])
        isWholeArea = SSSn == wholeArea
        for t in range(nTsteps):
            if TDF[TTTn,t] == 1:
                scoreXXX[t,0],scoreXXX[t,1],scoreXXX[t,2],scoreXXX[t,3],nwpDF_XXX_TP[:,t],guessDF_XXX_TP[:,t],nwpDF_XXX_CC[:,t],guessDF_XXX_CC[:,t],nwpDF_XXX_HCC[:,t],guessDF_XXX_HCC[:,t] = DFanalysis.simpleScoreDF(True,TP[:,:,t],CC[:,:,t],HCC[:,:,t],DF_TP[XXXn,:],DF_CC[XXXn,:],DF_HCC[XXXn,:],DF_priority[XXXn,:],SDF[:,:,SSSn],isWholeArea,inArea,TP_bins,CC_bins,'X')
                score[t,:]=scoreXXX[t,:]
            else:
                scoreYYY[t,0],scoreYYY[t,1],scoreYYY[t,2],scoreYYY[t,3],nwpDF_YYY_TP[:,t],guessDF_YYY_TP[:,t],nwpDF_YYY_CC[:,t],guessDF_YYY_CC[:,t],nwpDF_YYY_HCC[:,t],guessDF_YYY_HCC[:,t] = DFanalysis.simpleScoreDF(True,TP[:,:,t],CC[:,:,t],HCC[:,:,t],DF_TP[YYYn,:],DF_CC[YYYn,:],DF_HCC[YYYn,:],DF_priority[YYYn,:],SDF[:,:,SSSn],isWholeArea,inArea,TP_bins,CC_bins,'Y')
                score[t,:]=scoreYYY[t,:]
                if SSSn != wholeArea:
                    scoreXXX[t,0],scoreXXX[t,1],scoreXXX[t,2],scoreXXX[t,3],nwpDF_XXX_TP[:,t],guessDF_XXX_TP[:,t],nwpDF_XXX_CC[:,t],guessDF_XXX_CC[:,t],nwpDF_XXX_HCC[:,t],guessDF_XXX_HCC[:,t] = DFanalysis.simpleScoreDF(True,TP[:,:,t],CC[:,:,t],HCC[:,:,t],DF_TP[XXXn,:],DF_CC[XXXn,:],DF_HCC[XXXn,:],DF_priority[XXXn,:],SDF[:,:,SSSn],isWholeArea,inArea,TP_bins,CC_bins,'X')
                    score[t]=np.mean([scoreXXX[t],scoreYYY[t]])
        
        txt = DFanalysis.getText(DF_text,TDF_text,SDF_text,XXXn,YYYn,TTTn,SSSn,ruleCombos, ruleReplaces)
        allDayIcon = DFanalysis.getIcon(XXXmatrix[:,wholeArea,:],weightingFnc,DF_icons)
        partDayIcons = np.zeros(4).astype('str')
        for i in range(4):
            if data_timestep == 3:
                tstart = i*2
                tstop = (i+1)*2
                weight = np.array([1,1])
                partDayIcons[i] = DFanalysis.getIcon(XXXmatrix[:,wholeArea,tstart:tstop],weight,DF_icons)
            if data_timestep == 6:
                tstart = i
                tstop = (i+1)
                weight = np.array([1])
                partDayIcons[i] = DFanalysis.getIcon(XXXmatrix[:,wholeArea,tstart:tstop],weight,DF_icons)                    
        
        #outpt[progDay-1] = daytxt+'; day '+str(progDay)+'; '+txt+'; Icon;'+ allDayIcon + '; PartDayIcon;' + partDayIcons[0] + ';' + partDayIcons[1] + ';' + partDayIcons[2] + ';' + partDayIcons[3] +';\n'
        #print(outpt[progDay-1])
        maplims2 = np.array([np.min(lons[inArea])-0.25,np.max(lons[inArea])+0.25,np.max(lats[inArea])+0.25,np.min(lats[inArea]-0.25)]) # Domain of region
    #    if ((progDay>1) & (progDay<9)):
    #        plotTimeSeriesDFs(0,TP_bins,nwpDF_XXX_TP,nwpDF_YYY_TP,guessDF_XXX_TP,guessDF_YYY_TP,scoreXXX[:,1],scoreYYY[:,1],data_timestep,txt)
    #        plotTimeSeriesDFs(1,CC_bins,nwpDF_XXX_CC,nwpDF_YYY_CC,guessDF_XXX_CC,guessDF_YYY_CC,scoreXXX[:,2],scoreYYY[:,2],data_timestep,txt)
    #        plotTimeSeriesDFs(2,CC_bins,nwpDF_XXX_HCC,nwpDF_YYY_HCC,guessDF_XXX_HCC,guessDF_YYY_HCC,scoreXXX[:,3],scoreYYY[:,3],data_timestep,txt)
    #        
    #        plotTimeSeriesMaps(0,data_timestep,lats,lons,TP,sf,maplims2,inArea,TDF[TTTn],SDF_text[SSSn],SDF,SSSn,txt)
    #        plotTimeSeriesMaps(1,data_timestep,lats,lons,CC,sf,maplims2,inArea,TDF[TTTn],SDF_text[SSSn],SDF,SSSn,txt)
    #        plotTimeSeriesMaps(2,data_timestep,lats,lons,HCC,sf,maplims2,inArea,TDF[TTTn],SDF_text[SSSn],SDF,SSSn,txt)
        
    
    
    
    
    
    
    
    
    
    
        # Repeat all of that for winds
        XXXmatrix = np.zeros([nWDF,nSDF,nTsteps])+np.inf # score of each possible XXX vs SSS at every t
        YYYmatrix = np.zeros([nWDF,nSDF,nTsteps])+np.inf # same for YYY
        
        for XXXn in range(nWDF):
            #print(str(XXXn)+'/'+str(nWDF))
            for SSSn in range(nSDF):
                for t in range(nTsteps):
                    isWholeArea = SSSn == wholeArea
                    
                    XXXmatrix[XXXn,SSSn,t] = DFanalysis.simpleScoreDF_winds(False,U[:,:,t],V[:,:,t],DF_U[XXXn,:],DF_V[XXXn,:],DF_winds_priority[XXXn,:],SDF[:,:,SSSn],isWholeArea,inArea,W_bins,'X')
                    YYYmatrix[XXXn,SSSn,t] = DFanalysis.simpleScoreDF_winds(False,U[:,:,t],V[:,:,t],DF_U[XXXn,:],DF_V[XXXn,:],DF_winds_priority[XXXn,:],SDF[:,:,SSSn],isWholeArea,inArea,W_bins,'Y')
    
        # Now loop through and add up scores according to TDF        
        weightingFnc = np.ones(nTsteps) # Set weighting for text
        if data_timestep == 3:
            weightingFnc[0]=0.5 # 15Z
            weightingFnc[-1]=0.5 # 12Z
        if data_timestep == 6:
            weightingFnc[-1]=0.5 # 12Z
            
        score_matrix = np.zeros([nWDF,nWDF,nTDF,nSDF])+np.inf
    
        for XXXn in range(nWDF):
            #print(str(XXXn)+'/'+str(nWDF))
            # Add drop out if none of NWP looks like this type, i.e. all day just this type has high min score
            min_score = np.min(XXXmatrix[XXXn,wholeArea,:])
            if min_score>Cut_off_dont_try:
                dmy =0
                #print('Omitting '+DF_winds_text[XXXn]+' as primary due to minimum score of '+str(min_score))
            else:
                for TTTn in range(nTDF):
                    #if TTTn>0: # if TTTn == 0 there is no YYY, XXX is constant
                    for YYYn in range(nWDF):
                        # Catch banned combinations e.g. showers and showers, or scattered rain and a few showers
                        if banned_winds_list[XXXn,YYYn] == True: # 
                            for SSSn in range(nSDF):
                                score = np.zeros(nTsteps)
                                for t in range(nTsteps):
                                    if TDF[TTTn,t]==1:
                                        score[t] = weightingFnc[t] * XXXmatrix[XXXn,wholeArea,t]
                                    else:
                                        if SSSn == wholeArea:
                                            if TDF_text[TTTn] == 'XXX but YYY': # prevent Rain but showers with no spatial qualifier type scenarios
                                                score[t] = np.inf
                                            else:
                                                score[t] = weightingFnc[t] * YYYmatrix[YYYn,wholeArea,t]    
                                        else:
                                            scoreX = weightingFnc[t] * XXXmatrix[XXXn,SSSn,t]
                                            scoreY = weightingFnc[t] * YYYmatrix[YYYn,SSSn,t]
                                            score[t] = np.mean([scoreX,scoreY])
    
                                score_matrix[XXXn,YYYn,TTTn,SSSn] = np.mean(score) 
            #                    txt = getText(DF_text,TDF_text,SDF_text,XXXn,YYYn,TTTn,SSSn,ruleCombos, ruleReplaces)
            #                    print(txt+' score '+ str(np.mean(score)))                                          
                        else:
                            score_matrix[XXXn,YYYn,TTTn,:] = np.inf
    
        
        # Show best solution
        sortedScores = np.sort(score_matrix.flatten())
    #    for i in range(50):
    #        XXn,YYn,TTn,SSn = np.where(score_matrix==sortedScores[i])
    #        XXXn=XXn[0]
    #        YYYn=YYn[0]
    #        TTTn=TTn[0]
    #        SSSn=SSn[0]
    #        
    #        txt = DFanalysis.getText(DF_winds_text,TDF_text,SDF_text,XXXn,YYYn,TTTn,SSSn,ruleCombos, ruleReplaces)
    #        print('Rank '+str(i)+': '+txt+' score '+ str(score_matrix[XXXn,YYYn,TTTn,SSSn]))
        # Take best
        XXn,YYn,TTn,SSn = np.where(score_matrix==sortedScores[0])
        XXXn=XXn[0] 
        YYYn=YYn[0] 
        TTTn=TTn[0]
        SSSn=SSn[0]
        
        # Recalculate NWP DF and guess DF as time series
        nNwp = np.sum(inArea).astype('float')
        nwpDF_XXX_U = np.zeros([len(W_bins)-1,nTsteps])
        nwpDF_XXX_V = np.zeros([len(W_bins)-1,nTsteps])
        nwpDF_YYY_U = np.zeros([len(W_bins)-1,nTsteps])
        nwpDF_YYY_V = np.zeros([len(W_bins)-1,nTsteps])
    
        guessDF_XXX_U = np.zeros([len(W_bins)-1,nTsteps])
        guessDF_XXX_V = np.zeros([len(W_bins)-1,nTsteps])
        guessDF_YYY_U = np.zeros([len(W_bins)-1,nTsteps])
        guessDF_YYY_V = np.zeros([len(W_bins)-1,nTsteps])
    
        scoreXXX = np.zeros([nTsteps,3])
        scoreYYY = np.zeros([nTsteps,3])
        score = np.zeros([nTsteps,3])
        isWholeArea = SSSn == wholeArea
        for t in range(nTsteps):
            if TDF[TTTn,t] == 1:
                scoreXXX[t,0],scoreXXX[t,1],scoreXXX[t,2],nwpDF_XXX_U[:,t],guessDF_XXX_U[:,t],nwpDF_XXX_V[:,t],guessDF_XXX_V[:,t] = DFanalysis.simpleScoreDF_winds(True,U[:,:,t],V[:,:,t],DF_U[XXXn,:],DF_V[XXXn,:],DF_winds_priority[XXXn,:],SDF[:,:,SSSn],isWholeArea,inArea,W_bins,'X')
                score[t,:]=scoreXXX[t,:]
            else:
                scoreYYY[t,0],scoreYYY[t,1],scoreYYY[t,2],nwpDF_YYY_U[:,t],guessDF_YYY_U[:,t],nwpDF_YYY_V[:,t],guessDF_YYY_V[:,t] = DFanalysis.simpleScoreDF_winds(True,U[:,:,t],V[:,:,t],DF_U[YYYn,:],DF_V[YYYn,:],DF_winds_priority[YYYn,:],SDF[:,:,SSSn],isWholeArea,inArea,W_bins,'Y')
                score[t,:]=scoreYYY[t,:]
                if SSSn != wholeArea:
                    scoreXXX[t,0],scoreXXX[t,1],scoreXXX[t,2],nwpDF_XXX_U[:,t],guessDF_XXX_U[:,t],nwpDF_XXX_V[:,t],guessDF_XXX_V[:,t] = DFanalysis.simpleScoreDF_winds(True,U[:,:,t],V[:,:,t],DF_U[XXXn,:],DF_V[XXXn,:],DF_winds_priority[XXXn,:],SDF[:,:,SSSn],isWholeArea,inArea,W_bins,'X')
                    score[t]=np.mean([scoreXXX[t],scoreYYY[t]])
        
        txt_winds = DFanalysis.getText_winds(DF_winds_text,TDF_text,SDF_text,XXXn,YYYn,TTTn,SSSn)
#        plotTimeSeriesDFs(3,W_bins,nwpDF_XXX_U,nwpDF_YYY_U,guessDF_XXX_U,guessDF_YYY_U,scoreXXX[:,1],scoreYYY[:,1],data_timestep,txt_winds)
#        plotTimeSeriesDFs(4,W_bins,nwpDF_XXX_V,nwpDF_YYY_V,guessDF_XXX_V,guessDF_YYY_V,scoreXXX[:,2],scoreYYY[:,2],data_timestep,txt_winds)     
#        chosenTimestep = 5
#        compareTimeSeriesDFs(3,W_bins,nwpDF_XXX_U[:,chosenTimestep],DF_U,XXXmatrix[:,wholeArea,chosenTimestep],DF_winds_text,XXXn,txt_winds)
#        compareTimeSeriesDFs(4,W_bins,nwpDF_XXX_V[:,chosenTimestep],DF_V,XXXmatrix[:,wholeArea,chosenTimestep],DF_winds_text,XXXn,txt_winds)
#        plotTimeSeriesMaps(3,data_timestep,lats,lons,U,sf,maplims2,inArea,TDF[TTTn],SDF_text[SSSn],SDF,SSSn,txt_winds)
#        plotTimeSeriesMaps(4,data_timestep,lats,lons,V,sf,maplims2,inArea,TDF[TTTn],SDF_text[SSSn],SDF,SSSn,txt_winds)


        allDayIcon_winds = DFanalysis.getIcon(XXXmatrix[:,wholeArea,:],weightingFnc,DF_winds_icons)
        if allDayIcon_winds != 'none':
            allDayIcon = allDayIcon_winds # Wind trups other weather
        for i in range(4):
            if data_timestep == 3:
                tstart = i*2
                tstop = (i+1)*2
                weight = np.array([1,1])
                partDayIcons_winds = DFanalysis.getIcon(XXXmatrix[:,wholeArea,tstart:tstop],weight,DF_winds_icons)
                if partDayIcons_winds != 'none':
                    partDayIcons[i] = partDayIcons_winds # Wind trumps other weather
            if data_timestep == 6:
                tstart = i
                tstop = (i+1)
                weight = np.array([1])
                partDayIcons_winds = DFanalysis.getIcon(XXXmatrix[:,wholeArea,tstart:tstop],weight,DF_winds_icons)                    
                if partDayIcons_winds != 'none':
                    partDayIcons[i] = partDayIcons_winds # Wind trumps other weather    
        outpt = daytxt+'; day '+str(progDay)+'; '+txt+'. '+txt_winds+'; Icon;'+ allDayIcon + '; PartDayIcon;' + partDayIcons[0] + ';' + partDayIcons[1] + ';' + partDayIcons[2] + ';' + partDayIcons[3] +';\n'
        print(outpt)
        
        fout.writelines(outpt)
tdy = datetime.datetime.utcnow()
totalTime = int((time.time()-startTick)/60)
outpt = 'Finished at '+tdy.strftime('%H:%M %a %d-%b-%Y')+' in '+str(totalTime)+' minutes.'
fout.writelines(outpt)
fout.close()
