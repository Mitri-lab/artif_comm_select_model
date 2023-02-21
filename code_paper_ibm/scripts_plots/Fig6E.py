# -*- coding: utf-8 -*-
"""
Created on Sun May  2 20:58:22 2021

@author: pablo

Here, we change the dilution bottleneck between rounds in the simulations and see how 
it affects degradation (could have been also diversity, or growth). We plot for each propagation method, how
changing dilution affects the response variable. We then fit a line and do a spearman correlation test.

Since bottleneck is a percentage, except in disassembly which is a fix number of cells, we change the bottleneck
by multiplying times a factor.
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


#where the data is stored
path = "/Users/pablo/Desktop/code_paper/data/211026_bottleneck/"
path21 = "/Users/pablo/Desktop/code_paper/data/210822/" #for b_times = 1, i.e. the standars, the data is stored somewhere else
#some information about this data:
rd = 49 #the round at which we want to plot
repeats = 10 #number of repeats in these simulations 
num_tubes = 21
#seeds for the random number generator, corresponding to species sets
seeds = ["22", "23", "24", "25", "26"]
#times which factor we multiply the bottleneck
b_times = ["0.2", "0.5", "1", "2", "5"]
#propagation methods
conditions = [  "d3_s", "p_s", "m_s", "n_r",  "pip_s",  "mip_s"]
#renanme conditions for the plot
conditions_rename = [ "DS", "PS", "MS", "NS",  "PIS", "MIS" ]

path_save = "/Users/pablo/Desktop/220628_plots_mitri/" #path to save results
                    
#%%LOAD DATA FOR EACH TUBE
#DEGRADTION, AUC AND ALPHA DIV
deg_all = {}
auc_all = {}
sp_div_all = {}
for num_tube in b_times:
    deg_all[num_tube] = {}
    auc_all[num_tube] = {}
    sp_div_all[num_tube] = {}
    for seed in seeds:
        deg_all[num_tube][seed] = {}
        auc_all[num_tube][seed] = {}
        sp_div_all[num_tube][seed] = {}
        
        for prop in conditions:
            deg_all[num_tube][seed][prop] = []
            auc_all[num_tube][seed][prop] = []
            sp_div_all[num_tube][seed][prop] = []
            
            for rp in range(repeats):
                #append instead of extend, so one list for each repeat
                if num_tube != "1": 
                    df_here = pd.read_csv(path + num_tube + "/" + seed + "/" + prop + "/" + "repeat" + str(rp) + "/df_grid.csv")
                if num_tube == "1": #bottleneck times 1 is the standar used for regular simulations and data is stored in a different folder.
                    df_here = pd.read_csv(path21  + seed + "/" + prop + "/" + "repeat" + str(rp) + "/df_grid.csv")
                
                #rd = list(df_here["round"])[-1] #omit this line if we want to exclude species extinct before
                if rd in list(df_here["round"]):
                    deg_here = list(df_here.loc[df_here["round"] ==rd, "deg_score"])
                    pop_0 = list (df_here.loc[df_here["round"] ==rd, "tot_pop_0"])
                    auc_here = list(df_here.loc[df_here["round"] ==rd, "tot_auc"])
                    div_here = list(df_here.loc[df_here["round"] ==rd, "sp_div_H"])
                    if deg_here:
                        deg_all[num_tube][seed][prop].append(deg_here)
                    else:
                        deg_all[num_tube][seed][prop].append(np.nan)
                    auc_all[num_tube][seed][prop].append([x-(y*rd) for x,y in zip(auc_here,pop_0)])
                    sp_div_all[num_tube][seed][prop].append([exp(x) for x in div_here])
                else:
                    deg_all[num_tube][seed][prop].append([np.nan])
                    auc_all[num_tube][seed][prop].append([np.nan])
                    sp_div_all[num_tube][seed][prop].append([np.nan])
                    print("problem: extinct") 

#%%
"""
#PLOT MEDIAN OF EACH REPEAT, DO TESTS AND FIT LINE WITH ALL THE DATA
fig = plt.figure(figsize=(16,10))
i=1

for prop in conditions:
    plt.subplot(2,3,i)
    #for the correlation we fit the line with all the data
    x_all = []
    y_all = []
    for s,seed in enumerate(seeds):
        color = sns.color_palette("colorblind")[s]
        x = []
        y = []
        x_med = []
        y_med = []
        for num_tube in tubes:
            for rp in range(repeats):
                y.extend(deg_all[num_tube][seed][prop][rp])
                x.extend([int(num_tube) for j in deg_all[num_tube][seed][prop][rp]])
                y_med.append(np.median(deg_all[num_tube][seed][prop][rp]))
                x_med.append(int(num_tube))
        x_all.extend(x)
        y_all.extend(y)
        
        
    
        #fig = plt.figure(figsize=(7, 5))
        plt.scatter(x_med, y_med,
                   linewidths=1, alpha=0.7,
                   color=color, s = 70,
                   vmin = 0.15,
                   vmax = 0.85)
        plt.tick_params(axis="both",labelsize=17)
        
        
        plt.ylim(0,1)
        #plt.xlim(-5000,225000)
        plt.xlabel("num_tubes", fontsize=17)
        plt.ylabel("mmedian(deg)", fontsize=17)
        plt.title("Rounds: {}, prop: {}. Test w/ all tubes".format(int(rd),prop),fontsize=17, )
        #plt.colorbar()
        plt.tight_layout()
    
    m, b = np.polyfit(x_all, y_all, 1)
    plt.plot(np.array(x_all), m*np.array(x_all) + b, color = "black")
    plt.text(6, 0.9, "y = " + str(round(m,5)) + "x + " + str(round(b,2)), style='italic', bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 8}, fontsize=14)
    rho, pval = scipy.stats.spearmanr(x_all,y_all)
    plt.text(6, 0.8, "Spr. rho: {}, p-val: {}".format(round(rho,5), round(pval,9)), style='italic', fontsize=14)
    i+=1
    
    
plt.savefig(path_save + "mediandeg_tubes{}.png".format(int(rd)),dpi=500)
plt.show()
"""
#%%

#for some reason I can´t find, when importing graphs to inkscape in pdf, it does not show the text, so I will use
#istead .svg, but there I text will be an icon, unless I set the folowing:
plt.rcParams['svg.fonttype'] = "none"
#MAX DEGRADATION. TEST AND FIT INCLUDING ONLY THE BEST TUBE OF EACH REPEAT
fig = plt.figure(figsize=(16,10))
i=1

#importnat here I have x axis in log2 scale, since for example I double or reduce to half teh bottelneck
for i_p,prop in enumerate(conditions):
    if prop == "pip_s":
        i+=1
    plt.subplot(3,4,i)
    #for the correlation we fit the line with all the data
    x_all = []
    y_all = []
    for s,seed in enumerate(seeds):
        color = sns.color_palette("colorblind")[s]
        x_max = []
        y_max = []
        for num_tube in b_times:
            for rp in range(repeats):
                #check that is not np.nan
                if np.max(deg_all[num_tube][seed][prop][rp]) == np.max(deg_all[num_tube][seed][prop][rp]):
                    y_max.append(np.max(deg_all[num_tube][seed][prop][rp]))
                    x_max.append(np.log2(float(num_tube)))
                    y_all.append(np.max(deg_all[num_tube][seed][prop][rp]))
                    x_all.append(np.log2(float(num_tube)))
            
    
        #fig = plt.figure(figsize=(7, 5))
        plt.scatter(x_max, y_max,
                   linewidths=1, alpha=0.7,
                   color=color, s = 60)
        plt.tick_params(axis="both",labelsize=17)
        
        
        plt.ylim(0,1)
        #plt.xlim(-5000,225000)
        plt.xlabel("2-fold change of bottleneck", fontsize=17)
        plt.ylabel("max(deg)", fontsize=17)
        plt.title("round: {}, prop: {}".format(int(rd),conditions_rename[i_p]),fontsize=17 )
        #plt.colorbar()
        plt.tight_layout()
    
    m, b = np.polyfit(x_all, y_all, 1)
    plt.plot(np.array(x_all), m*np.array(x_all) + b, color = "black")
    rho, pval = scipy.stats.spearmanr(x_all,y_all)
    correct_height = 0
    if min(y_all) < 0.19:
        correct_height = 0.82
    plt.text(-2.4, correct_height+ 0.08, "y = " + str(round(m,5)) + "x + " + str(round(b,2)),  fontsize=10) #bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 8},
    plt.text(-2.4, correct_height+0.02, "Spr rho: {}, p-val: {}".format(round(rho,5), fmt(pval)),  fontsize=10)
    plt.tight_layout()
    i+=1
    
    
plt.savefig(path_save + "maxdeg_bottleneck{}.pdf".format(int(rd)))
plt.savefig(path_save + "maxdeg_bottleneck{}.svg".format(int(rd)))
plt.show()

#%%PLOT GRID
sns.set_style("white")
data = []
headers = b_times




for i_p,prop in enumerate(conditions):
    data_here = []
    for num_tube in b_times:
        y_all = [] #here we append a value for each N_comms value 
        for s,seed in enumerate(seeds):
            for rp in range(repeats):
                #check it is not Nan value (not equal to itself)
                if np.max(deg_all[num_tube][seed][prop][rp]) == np.max(deg_all[num_tube][seed][prop][rp]):
                    y_all.append(np.max(deg_all[num_tube][seed][prop][rp]))
                    
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
plt.xlabel("Factor multiplied by bottleneck",fontsize=14)
plt.ylabel("propagation",fontsize=14)
plt.savefig(path_save + "grid_median_bottleneck_{}.pdf".format(rd))
plt.savefig(path_save + "grid_median_bottleneck_{}.svg".format(rd))

    
    



