# -*- coding: utf-8 -*-
"""
Created on Sun May  2 20:58:22 2021

@author: pablo

Here, we change the number of species in the simulations and see how 
it affects degradation (could have been also diversity, or growth). We plot for each propagation method, how
changing number of species (x axes) affects the response variable. We then fit a line and do a spearman correlation test.


"""

import numpy as np
import random
import matplotlib.pyplot as plt
from math import sqrt,exp,log
import copy
import pickle
from datetime import datetime
import os
import inspect
from scipy.spatial import distance
import pandas as pd
import seaborn as sns
import pingouin as pg
import scipy
import itertools
import matplotlib.ticker as mticker
f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
g = lambda x,pos : "${}$".format(f._formatSciNotation('%1.5e' % x))
fmt = mticker.FuncFormatter(g)

def abline(slope, intercept):
    """Plot a line from slope and intercept"""
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--')

# path to data
path = "/Users/pablo/Desktop/code_paper/data/210822_subsamples/"
path21 = "/Users/pablo/Desktop/code_paper/data/210822/" #data for 15 species is in a different folder (this is the standard)
# sme information about the dataset
seeds = ["22", "23", "24", "25", "26"]

conditions = [   "d3_s", "p_s", "m_s", "n_r",  "pip_s",  "mip_s"]
#renanme conditions for the plot
conditions_rename = [ "DS", "PS", "MS", "NS",  "PIS", "MIS" ] #important, same order as above

#now anyway I´ll take the last round when it went extinct (again, I might rethink about extinction condition)
rd = 49 #the round at which we want to plot (starting at 0)
samples = ["6", "9", "12", "15"] #number of species we have
num_samples = 5 #number of subsamples we do for each number of species:
                # we always keep the same 15 species but do subsamples of e.g. 9 speices, since the species
                # we chose should matter, we don´t have only 1 subsample but 5 different ones.
shapes = [".", "2", "+", "x",  "d", "p", "4", "h"] #for teh plot each subsample will have a shape
repeats = 5 #number of repeats in these simulations 

path_save = "/Users/pablo/Desktop/220628_plots_mitri/" #path to save results

                    
#%%LOAD DATA FOR EACH TUBE
#DEGRADTION, AUC AND ALPHA DIV
deg_all = {}
auc_all = {}
sp_div_all = {}
for subsample in samples:
    deg_all[subsample] = {}
    auc_all[subsample] = {}
    sp_div_all[subsample] = {}
    for i_s in range(num_samples):
        deg_all[subsample][i_s] = {}
        auc_all[subsample][i_s] = {}
        sp_div_all[subsample][i_s] = {}

        for seed in seeds:
            deg_all[subsample][i_s][seed] = {}
            auc_all[subsample][i_s][seed] = {}
            sp_div_all[subsample][i_s][seed] = {}
            
            for prop in conditions:
                deg_all[subsample][i_s][seed][prop] = []
                auc_all[subsample][i_s][seed][prop] = []
                sp_div_all[subsample][i_s][seed][prop] = []
                
                for rp in range(repeats):
                    if subsample != "15":
                        df_here = pd.read_csv(path + seed + "_" + subsample+ "/" + str(i_s) + "/" + prop + "/" + "repeat" + str(rp) + "/df_grid.csv")
                        #rd = list(df_here["round"])[-1] #omit this line if we want to exclude species extinct before
                    else: #for 15 species the sample is in a different folder
                        if i_s == 0: #here we only have one species subsets, the 15 species
                            df_here = pd.read_csv(path21 + seed + "/"  + prop + "/" + "repeat" + str(rp) + "/df_grid.csv")
                    
                    if rd in list(df_here["round"]):

                        deg_here = list(df_here.loc[df_here["round"] ==rd, "deg_score"])
                        pop_0 = list (df_here.loc[df_here["round"] ==rd, "tot_pop_0"])
                        auc_here = list(df_here.loc[df_here["round"] ==rd, "tot_auc"])
                        div_here = list(df_here.loc[df_here["round"] ==rd, "sp_div_H"])
                        
                        deg_all[subsample][i_s][seed][prop].append(deg_here)
                        auc_all[subsample][i_s][seed][prop].append([x-(y*rd) for x,y in zip(auc_here,pop_0)])
                        sp_div_all[subsample][i_s][seed][prop].append([exp(x) for x in div_here])
                        
                    else:
                        deg_all[subsample][i_s][seed][prop].append([np.nan])
                        auc_all[subsample][i_s][seed][prop].append([np.nan])
                        sp_div_all[subsample][i_s][seed][prop].append([np.nan])
                        print("problem: extinct") #if you ever get this problem take 210824_correlations_tubes as an example to how to solve it


#%%PLOT 
plt.rcParams['svg.fonttype'] = "none"#to properly import .svg plots to inkscape


#MAX
fig = plt.figure(figsize=(18,11))
i=1

for i_p,prop in enumerate(conditions):
    if prop == "pip_s":
        i+=1
    plt.subplot(3,4,i)
    #for the correlation we fit the line with all the data
    x_all = []
    y_all = []
    for s,seed in enumerate(seeds):
        color = sns.color_palette("colorblind")[s]
        
        for subsample in samples: #number of species
            for i_s in range(num_samples): #which of the samples for this number of species
                shape = shapes[i_s]
                x = []
                y = []
                for rp in range(repeats):
                    if subsample != "15": 
                        #check it is not NaN
                        if np.max(deg_all[subsample][i_s][seed][prop][rp]) == np.max(deg_all[subsample][i_s][seed][prop][rp]):
                            y.append(max(deg_all[subsample][i_s][seed][prop][rp]))
                            x.append(int(subsample))
                            y_all.append(max(deg_all[subsample][i_s][seed][prop][rp]))
                            x_all.append(int(subsample))
                    else: #for 15 species we only have one subsample
                        if i_s == 0:
                            #check it is not NaN
                            if np.max(deg_all[subsample][i_s][seed][prop][rp]) == np.max(deg_all[subsample][i_s][seed][prop][rp]):
                                y.append(max(deg_all[subsample][i_s][seed][prop][rp]))
                                x.append(int(subsample))
                                y_all.append(max(deg_all[subsample][i_s][seed][prop][rp]))
                                x_all.append(int(subsample))
                if len(x) != 0: #since will be empty for few repeats (1-4) of 15 species
                    #fig = plt.figure(figsize=(7, 5))
                    plt.scatter(x, y,
                               linewidths=1, alpha=0.7,
                               color=color, s = 60, marker = shape)
                    plt.tick_params(axis="both",labelsize=17)

    
    plt.ylim(0,1)
    #plt.xlim(-5000,225000)
    plt.xlabel("Number of species", fontsize=17)
    plt.ylabel("max(deg)", fontsize=17)
    plt.title("round: {}, prop: {}".format(int(rd),conditions_rename[i_p]),fontsize=17)
    #plt.colorbar()
    plt.tight_layout()

    m, b = np.polyfit(x_all, y_all, 1)
    plt.plot(np.array(x_all), m*np.array(x_all) + b, color = "black")
    rho, pval = scipy.stats.spearmanr(x_all,y_all)
    correct_height = 0
    if min(y_all) < 0.19:
        correct_height = 0.82
    plt.text(6, correct_height + 0.08, "y = " + str(round(m,5)) + "x + " + str(round(b,2)), style='italic', fontsize=11)#, bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 8}
    plt.text(6, correct_height + 0.02, "Spr.rho: {}, p: {}".format(round(rho,5), fmt(pval)), style='italic', fontsize=11)
   
    i+=1
    
    
plt.savefig(path_save + "maxdeg_subsamples{}_ok.pdf".format(int(rd)))
plt.savefig(path_save + "maxdeg_subsamples{}_ok.svg".format(int(rd)))
plt.show()

#%%PLOT GRID
sns.set_style("white")
data = []
headers = samples




for i_p,prop in enumerate(conditions):
    data_here = []
    for subsample in samples:
        y_all = [] #here we append a value for each N_comms value 
        for i_s in range(num_samples):
            for seed in seeds:
    
                for rp in range(repeats):
                    #check it is not Nan value (not equal to itself)
                    if subsample != "15":
                        if np.max(deg_all[subsample][i_s][seed][prop][rp]) == np.max(deg_all[subsample][i_s][seed][prop][rp]):
                            y_all.append(np.max(deg_all[subsample][i_s][seed][prop][rp]))
                    else: #for 15 species we only have one subsample
                        if i_s == 0:
                            if np.max(deg_all[subsample][i_s][seed][prop][rp]) == np.max(deg_all[subsample][i_s][seed][prop][rp]):
                                y_all.append(np.max(deg_all[subsample][i_s][seed][prop][rp]))

        data_here.append(np.median(y_all))
    data.append(data_here)


df = pd.DataFrame(data,columns=headers, index = conditions_rename)
fig, ax = plt.subplots(1,1)
img=ax.imshow(df, cmap="viridis",vmin=0.2, vmax=0.8) 
#rename the axes
ax.set_xticks(list(range(len(data[1]))))
ax.set_xticklabels(headers,fontsize=12)
ax.set_yticks(list(range(len(conditions_rename))))
ax.set_yticklabels(conditions_rename,fontsize=12)
#add the colorbar
fig.colorbar(img)
#plt.xticks(rotation=45)
plt.title("Median(max(degradation))",fontsize=14)
plt.xlabel("Number of species",fontsize=14)
plt.ylabel("propagation",fontsize=14)
plt.savefig(path_save + "grid_median_subsamples_{}.pdf".format(rd))
plt.savefig(path_save + "grid_median_subsamples_{}.svg".format(rd))