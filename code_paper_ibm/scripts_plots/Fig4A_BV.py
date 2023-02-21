# -*- coding: utf-8 -*-
"""
Created on Sat Jun 26 19:35:24 2021

@author: pablo

Script to check and plot how many combinations each method has explored (cumulative over rounds).
Here we show the average over repeats and seeds to condense all the information in one plot.
To do so, we consider each repeat of each seed independent, so we have 50 data points (5seeds*10repeats),
and we do the average of those 50. 

added np.nanmean and np.nanstd, to avoid those "extinct" repeats when doing the mean
"""

#let´s try to plot the number of combinations explored

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
#import pingouin as pg

import itertools

#path to the data
#path = "/Users/pablo/Desktop/data_mitri/210822/"
path = "/Users/pablo/Desktop/code_paper/data/210822/" 
#some info about these simulations
seeds = ["22", "23", "24", "25", "26"] #random seeds will be renamed as 1-5 just for the plot
#propagation methods we are using.
conditions = ["d3_s", "p_s", "m_s", "pip_s","mip_s", "n_r", "d3_r", "p_r", "m_r", "pip_r" , "mip_r"] 
#IMPORTANT: in the same order as above
conditions_rename = ["DS", "PS", "MS", "PIS", "MIS", "NS", "DR", "PR", "MR", "PIR", "MIR" ]

repeats = 10
rounds = 50
num_tubes = 21
#were to save the plots that we generate
#path_save = "/Users/pablo/Desktop/211111_results/"
path_save = "/Users/pablo/Desktop/220628_plots_mitri_b/"

#in this version of the plot we consider all the 50 repeats (10 of each of the 5 random seeds),
#independent, we could have done the mean for each seed and then the mean over seeds, just smthing to consider.

#%% load data
number_combs = {}
for seed in seeds:
    number_combs[seed] = {}
    for prop in conditions:
        number_combs[seed][prop] = []
        for repeat in range(repeats):
            num_combs = []
            combs = []
            df = pd.read_csv(path + seed + "/" + prop + "/" + "repeat" + str(repeat) + "/df_sp.csv")
            for rd in range(rounds):
                df_here = df.loc[df["round"] == rd]
                for tube in range(num_tubes):
                    comb = list(df_here.loc[df_here["tube"]==tube, "sp"])
                    combs.append(tuple(comb))
                num_combs.append(len(set(combs)))
                combs = list(set(combs)) #I guess this will save a bit of ram
            number_combs[seed][prop].append(num_combs)
            
            
#reorganise the information
number_combs_all50 = {} #add each repeat for each seed independently
for prop in conditions:
    number_combs_all50[prop] = []
    for seed in seeds:
        number_combs_all50[prop].extend(number_combs[seed][prop])


#%%PLOT
plt.rcParams['svg.fonttype'] = "none" #to later properly import the text in .svg into inkscape

colours = sns.color_palette("Set2")

t=range(rounds)
sns.set_style("darkgrid")

dict_of_2 = {}#to assing a colour I´ll take advantage of the 2 first caracters that are unique by method
linestyles = {"s":"-","r":"--"}
i_color = 0
for i_p,prop in enumerate(conditions):
    
    if prop[:2] in dict_of_2: #if we have already appended the one under selection we assing same color to random
        color_here = dict_of_2[prop[:2]]
    else: 
        color_here = colours[i_color]
        dict_of_2[prop[:2]] = color_here
        i_color +=1
    line_here = linestyles[prop[-1]]
    
    mean = [np.nanmean(x) for x in zip(*number_combs_all50[prop])]
    std = [np.nanstd(x) for x in zip(*number_combs_all50[prop])]

    mean = np.array(mean)
    std = np.array(std)
    plt.plot(t,mean,label=conditions_rename[i_p], color=color_here, linestyle= line_here, linewidth=3,alpha=0.8)
    #plt.fill_between(t, mean - std, mean + std, color=color_here, alpha=0.1, edgecolor='none')
    plt.legend(ncol=4, fontsize=12.5,framealpha=0.5, title = "propagation")
    plt.ylim(0, 500)
    plt.title("10 repeats, 5 species sets", fontsize=16)
    plt.xlabel("Rounds",fontsize=16)
    plt.ylabel("Explored combinations",fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
plt.savefig(path_save + "explored_combs_one_plot.pdf", bbox_inches='tight')
plt.savefig(path_save + "explored_combs_one_plot.svg", bbox_inches='tight')

plt.show()


