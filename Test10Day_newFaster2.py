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
#from matplotlib import gridspec
#import matplotlib.path as mpltPath
#import matplotlib.patches as mpatches
#from mpl_toolkits.basemap import Basemap
import pdb
from plotTimeSeriesMaps import plotTimeSeriesMaps, plotTimeSeriesDFs
import loadDF
import DFanalysis

shpfn = '/mnt/store/tadams@met.co.nz/Autoregionals/data/TKI_subregions4.shp'
optfldr = '/mnt/store/tadams@met.co.nz/Autoregionals/'
nwpfn = '/data/models/ecmwf-prod/det/'
ROI = 'Taranaki' # region of interest - just one for now!
mdlrun = '2021011312'
NWPfn = '/mnt/store/tadams@met.co.nz/Autoregionals/data/NWP.npy' # presaved NWP
NWPfn2 = '/mnt/store/tadams@met.co.nz/Autoregionals/data/NWPlatlongs.npy' # presaved NWP

#maplims = np.array([165.,180.,-33.,-48.]) # Total domain of NWP to load from
maplims = np.array([173.5,175.,-38.5,-40.]) # Taranaki
plt.close('all')
loadGrib = True
Cut_off_dont_try = 2 # If any type is above this all timesteps, don't try
# load in DFs
lats, lons = np.load(NWPfn2)
nDF,DF,DF_text,banned_list_precip, TP_bins = loadDF.loadDFs()
nSDF,SDF,SDF_text,inArea,wholeArea,sf = loadDF.loadSDF(lats,lons,shpfn)
outpt = np.zeros(10).astype('U128')  

for progDay in range(1,11): # 1 is day 1, 10 is day 10
    if progDay <7:
        data_timestep = 3 # 3 hour data
    else:
        data_timestep = 6 # 6 hour data
    lats,lons,TP,daytxt = loadDF.loadNWP(loadGrib,NWPfn,NWPfn2,mdlrun,progDay,data_timestep,nwpfn,maplims)
    nTDF,TDF,TDF_text = loadDF.loadTDF(data_timestep)
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
                XXXmatrix[XXXn,SSSn,t] = DFanalysis.simpleScoreDF(False,TP[:,:,t],DF[XXXn,:,0],SDF[:,:,SSSn],isWholeArea,inArea,TP_bins,'X')
                YYYmatrix[XXXn,SSSn,t] = DFanalysis.simpleScoreDF(False,TP[:,:,t],DF[XXXn,:,0],SDF[:,:,SSSn],isWholeArea,inArea,TP_bins,'Y')

    # Now loop through and add up scores according to TDF        
    weightingFnc = np.ones(nTsteps) # Set weighting for text
#    if data_timestep == 3:
#        weightingFnc[0]=0.5 # 15Z
#        weightingFnc[-1]=0.5 # 12Z
#    if data_timestep == 6:
#        weightingFnc[-1]=0.5 # 12Z
        
    score_matrix = np.zeros([nDF,nDF,nTDF,nSDF])+np.inf

    for XXXn in range(nDF):
        #print(str(XXXn)+'/'+str(nDF))
        # Add drop out if none of NWP looks like this type, i.e. all day just this type has high min score
        min_score = np.min(XXXmatrix[XXXn,wholeArea,:])
        if min_score>Cut_off_dont_try:
            print('Omitting '+DF_text[XXXn]+' as primary due to minimum score of '+str(min_score))
        else:
            for TTTn in range(nTDF):
                #if TTTn>0: # if TTTn == 0 there is no YYY, XXX is constant
                for YYYn in range(nDF):
                    # Catch banned combinations e.g. showers and showers, or scattered rain and a few showers
                    if DFanalysis.DF_banned_list(0,XXXn,YYYn,banned_list_precip): # 0 as only precip so far, but this is for parameter
                        for SSSn in range(nSDF):
                            score = np.zeros(nTsteps)
                            for t in range(nTsteps):
                                if TDF[TTTn,t]==1:
                                    score[t] = weightingFnc[t] * XXXmatrix[XXXn,wholeArea,t]
                                else:
                                    if SSSn == wholeArea:
                                        score[t] = weightingFnc[t] * YYYmatrix[YYYn,wholeArea,t]    
                                    else:
                                        scoreX = weightingFnc[t] * XXXmatrix[XXXn,SSSn,t]
                                        scoreY = weightingFnc[t] * YYYmatrix[YYYn,SSSn,t]
                                        score[t] = np.mean([scoreX,scoreY])

                            score_matrix[XXXn,YYYn,TTTn,SSSn] = np.mean(score)
        #                    txt = getText(DF_text,TDF_text,SDF_text,XXXn,YYYn,TTTn,SSSn)
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
#        txt = DFanalysis.getText(DF_text,TDF_text,SDF_text,XXXn,YYYn,TTTn,SSSn)
#        print('Rank '+str(i)+': '+txt+' score '+ str(score_matrix[XXXn,YYYn,TTTn,SSSn]))
    # Take best
    XXn,YYn,TTn,SSn = np.where(score_matrix==sortedScores[0])
    XXXn=XXn[0] 
    YYYn=YYn[0] 
    TTTn=TTn[0]
    SSSn=SSn[0]
    
    # Recalculate NWP DF and guess DF as time series
    nNwp = np.sum(inArea).astype('float')
    nwpDF_XXX = np.zeros([len(TP_bins)-1,nTsteps])
    nwpDF_YYY = np.zeros([len(TP_bins)-1,nTsteps])
    guessDF_XXX = np.zeros([len(TP_bins)-1,nTsteps])
    guessDF_YYY = np.zeros([len(TP_bins)-1,nTsteps])
    scoreXXX = np.zeros([nTsteps])
    scoreYYY = np.zeros([nTsteps])
    score = np.zeros([nTsteps])
    isWholeArea = SSSn == wholeArea
    for t in range(nTsteps):
        if TDF[TTTn,t] == 1:
            scoreXXX[t],nwpDF_XXX[:,t],guessDF_XXX[:,t] = DFanalysis.simpleScoreDF(True,TP[:,:,t],DF[XXXn,:,0],SDF[:,:,SSSn],isWholeArea,inArea,TP_bins,'X')
            score[t]=scoreXXX[t]
        else:
            scoreYYY[t],nwpDF_YYY[:,t],guessDF_YYY[:,t] = DFanalysis.simpleScoreDF(True,TP[:,:,t],DF[YYYn,:,0],SDF[:,:,SSSn],isWholeArea,inArea,TP_bins,'Y')
            score[t]=scoreYYY[t]
            if SSSn != wholeArea:
                scoreXXX[t],nwpDF_XXX[:,t],guessDF_XXX[:,t] = DFanalysis.simpleScoreDF(True,TP[:,:,t],DF[XXXn,:,0],SDF[:,:,SSSn],isWholeArea,inArea,TP_bins,'X')
                score[t]=np.mean([scoreXXX[t],scoreYYY[t]])
    
    txt = DFanalysis.getText(DF_text,TDF_text,SDF_text,XXXn,YYYn,TTTn,SSSn)
    outpt[progDay-1] = daytxt+' day '+str(progDay)+' '+txt
    print(outpt[progDay-1])
    #plotTimeSeriesDFs(TP_bins,nwpDF_XXX,nwpDF_YYY,guessDF_XXX,guessDF_YYY,scoreXXX,scoreYYY,data_timestep,txt)
    
    #plotTimeSeriesMaps(data_timestep,lats,lons,TP,sf,maplims,inArea,TDF[TTTn],SDF_text[SSSn],SDF,SSSn,txt)
    #pdb.set_trace()
#    
    