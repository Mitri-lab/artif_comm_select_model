# -*- coding: utf-8 -*-
"""
Created on Sun May  2 20:58:22 2021

@author: pablo

Here, we change the number of tubes or communities where invassion occurs (f_inv) in the simulations and see how 
it affects degradation (could have been also diversity, or growth). We plot for each propagation method, how
changing f_invasion (x axes) affects the response variable. We then fit a line and do a spearman correlation test.

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
#used to round pval
f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
g = lambda x,pos : "${}$".format(f._formatSciNotation('%1.5e' % x))
fmt = mticker.FuncFormatter(g)

import matplotlib.ticker as mticker
def abline(slope, intercept):
    """Plot a line from slope and intercept"""
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--')

#the two funcitons bellow are taken from https://stackoverflow.com/questions/8671808/matplotlib-avoiding-overlapping-datapoints-in-a-scatter-dot-beeswarm-plot
"""
def rand_jitter(arr):
    stdev = .01 * (max(arr) - min(arr))
    return arr + np.random.randn(len(arr)) * stdev

def jitter(x, y, s=20, c='b', marker='o', cmap=None, norm=None, vmin=None, vmax=None, alpha=None, linewidths=None, verts=None, hold=None, **kwargs):
    return plt.scatter(rand_jitter(x), rand_jitter(y), s=s, c=c, marker=marker, cmap=cmap, norm=norm, vmin=vmin, vmax=vmax, alpha=alpha, linewidths=linewidths, **kwargs)
"""


#path to data
path = "/Users/pablo/Desktop/code_paper/data/210823_invasion/"
#for the condition of f_inv = 0.25, data is stored somewere else
path21 = "/Users/pablo/Desktop/code_paper/data/210822/"
#info about the dataset
rd = 49 #the round at which we want to plot (starting to count from 0)
repeats = 10 #number of repeats in these simulations 
num_tubes = 21
num_toxs = 10
#seeds we use for the random number generator 
seeds = [ "22","23", "24", "25", "26"]
#differentr values we use for f_inv, well the names of the folder, a bit more exact values are in integer tube
tubes = ["0_05", "0_15", "0_25", "0_35", "0_45", "0_55", "0_62", "0_715", "0_81", "0_905", "1"]

#for 0_25 species the sample is in a different folder
#I correct a bit the frequencies to which they correspond, since either way I just take an exact number of tubes
integer_tube = [0.05, 0.145, 0.24, 0.335, 0.43, 0.525, 0.62, 0.715, 0.81, 0.905, 1]
conditions = [   "d3_s", "pip_s", "mip_s"]
conditions_rename = [   "DS", "PIS", "MIS"]



path_save = "/Users/pablo/Desktop/220628_plots_mitri/" #path to save results
                    
#%%LOAD DATA 
#DEGRADTION, AUC AND ALPHA DIV
deg_all = {}
auc_all = {}
sp_div_all = {}
for num_tube in tubes:
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
                if num_tube != "0_25": 
                    df_here = pd.read_csv(path + num_tube + "/" + seed + "/" + prop + "/" + "repeat" + str(rp) + "/df_grid.csv")
                if num_tube == "0_25": #f_inv=0.25 is the standard used for regular simulations and data is stored in a different folder.
                    df_here = pd.read_csv(path21  + seed + "/" + prop + "/" + "repeat" + str(rp) + "/df_grid.csv")

                #rd = list(df_here["round"])[-1] #omit this line if we want to exclude species extinct before
                if rd in list(df_here["round"]):
                    deg_here = list(df_here.loc[df_here["round"] ==rd, "deg_score"])
                    pop_0 = list (df_here.loc[df_here["round"] ==rd, "tot_pop_0"])
                    auc_here = list(df_here.loc[df_here["round"] ==rd, "tot_auc"])
                    div_here = list(df_here.loc[df_here["round"] ==rd, "sp_div_H"])
                    
                    deg_all[num_tube][seed][prop].append(deg_here)
                    auc_all[num_tube][seed][prop].append([x-(y*rd) for x,y in zip(auc_here,pop_0)])
                    sp_div_all[num_tube][seed][prop].append([exp(x) for x in div_here])
                else:
                    deg_all[num_tube][seed][prop].append([np.nan])
                    auc_all[num_tube][seed][prop].append([np.nan])
                    sp_div_all[num_tube][seed][prop].append([np.nan])
                    print("problem: extinct") #if you ever get this problem take 210824_correlations_tubes as an example to how to solve it




#%%MAX DEGRADATION
plt.rcParams['svg.fonttype'] = "none"#to properly import .svg plots to inkscape


fig = plt.figure(figsize=(21,9))
i=1

for i_p,prop in enumerate(conditions):
    plt.subplot(2,4,i)
    #for the correlation we fit the line with all the data
    x_all = []
    y_all = []
    for s,seed in enumerate(seeds):
        color = sns.color_palette("colorblind")[s]
        x_max = []
        y_max = []
        for i_t, num_tube in enumerate(tubes):
            for rp in range(repeats):
                #check that is not np.nan
                if np.max(deg_all[num_tube][seed][prop][rp]) == np.max(deg_all[num_tube][seed][prop][rp]):
                    y_max.append(max(deg_all[num_tube][seed][prop][rp]))
                    x_max.append(integer_tube[i_t])
                    y_all.append(max(deg_all[num_tube][seed][prop][rp]))
                    x_all.append(integer_tube[i_t])
        
        
    
        #fig = plt.figure(figsize=(7, 5))
        plt.scatter(x_max, y_max,
                   linewidths=1, alpha=0.7,
                   color=color, s = 60)
        plt.tick_params(axis="both",labelsize=17)
        
        
        plt.ylim(0,1)
        #plt.xlim(-5000,225000)
        plt.xlabel("fraction of tubes with invasion", fontsize=17)
        plt.ylabel("max(deg)", fontsize=17)
        plt.title("round: {}, prop: {}".format(int(rd),conditions_rename[i_p]),fontsize=17, )
        #plt.colorbar()
        plt.tight_layout()
    
    m, b = np.polyfit(x_all, y_all, 1)
    plt.plot(np.array(x_all), m*np.array(x_all) + b, color = "black")
    rho, pval = scipy.stats.spearmanr(x_all,y_all)
    correct_height = 0
    if min(y_all) < 0.19:
        correct_height = 0.82
    plt.text(0.05, correct_height + 0.08, "y = " + str(round(m,5)) + "x + " + str(round(b,2)), style='italic', fontsize=12)#bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 8}
    plt.text(0.05, correct_height + 0.02, "Spr rho: {}, p-val: {}".format(round(rho,5), fmt(pval)), style='italic', fontsize=12)
    i+=1
    
    
plt.savefig(path_save + "maxdeg_invasion{}.pdf".format(int(rd)))
plt.savefig(path_save + "maxdeg_invasion{}.svg".format(int(rd)))

plt.show()

#%%PLOT GRID
sns.set_style("white")
data = []
headers = [str(x) for x in integer_tube]

for i_p,prop in enumerate(conditions):
    data_here = []
    for num_tube in tubes:
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
ax.set_xticklabels(["1","3","5", "7","9", "11", "13", "15", "17", "19", "21"],fontsize=12)
ax.set_yticks(list(range(len(conditions_rename))))
ax.set_yticklabels(conditions_rename,fontsize=12)
#add the colorbar
fig.colorbar(img)
plt.xticks(rotation=0)
plt.title("Median(max(degradation))",fontsize=14)
plt.xlabel("Number of tubes with invasion",fontsize=14)
plt.ylabel("propagation",fontsize=14)
plt.savefig(path_save + "grid_median_invasion_{}.pdf".format(rd))
plt.savefig(path_save + "grid_median_invasion_{}.svg".format(rd))






