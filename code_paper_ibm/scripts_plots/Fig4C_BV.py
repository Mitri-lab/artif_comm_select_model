# -*- coding: utf-8 -*-
"""
Created on Sat Jun 26 19:35:24 2021

@author: pablo

Script to calculate the beta diversity, here calculated as the mean Bray Curtis distance
in the 210 possible pair of tubes. Tubes are describe as a vector of 15 numbers, each of them
indicating the final popualtion for each of the 15 species. If a species is not present, just 0.

Here we combine everything into 1 plot, considering all the 10 repeats for the 5 sp sets as
independent. Then I can do the mean and just show each propagation with a different color.

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
# import pingouin as pg
import itertools



# Let´s load the data
#path = "/Users/pablo/Desktop/data_mitri/210822/"
path = "/Users/pablo/Desktop/code_paper/data/210822/"  

#some info about this dataset
num_tubes = 21
rounds = 50
num_tubes = 21
repeats = 10
#seeds used for the random number generator, renamed as 1-5 for the plot
seeds = ["22", "23", "24", "25", "26"]
#propagation methods
conditions = ["d3_s", "p_s", "m_s", "pip_s","mip_s", "n_r", "d3_r", "p_r", "m_r", "pip_r", "mip_r"] #atention: in the draft paper the last item in this list was mistakenly mip_s (instead of the right one mip_r, that matches the rename MIR) this should be corrected for the last version
#IMPORTANT: in the same order as above
conditions_rename = ["DS", "PS", "MS", "PIS", "MIS","NS", "DR", "PR", "MR", "PIR", "MIR"]

#where to save the plots
#path_save = "/Users/pablo/Desktop/211111_results/"
path_save = "/Users/pablo/Desktop/220628_plots_mitri_b/"


#%% load data

pairs = list(itertools.combinations(list(range(21)),2))

# I calculate the bray curtis distance between each pair of tubes, and then I do the mean of the 210 combinations.
beta_div = {} #diversity in the metacommunity
for seed in seeds:
    beta_div [seed] = {}
    for prop in conditions:
        beta_div[seed][prop] = []
        for repeat in range(repeats):
            div_o_r = []
            df = pd.read_csv(path + seed + "/" + prop + "/" + "repeat" + str(repeat) + "/df_sp.csv")
            for rd in range(rounds):
                df_here = df.loc[df["round"] == rd]
                tubes = []
                for tube in range(num_tubes):
                    this_tube = [0 for x in range(15)]
                    sps = list(df_here.loc[df_here["tube"]==tube, "sp"])
                    pops = list(df_here.loc[df_here["tube"]==tube, "final_pop"])
                    for i in range(len(sps)):
                        this_tube[int(sps[i])] = pops[i]
                    tubes.append(this_tube)
                #calculate diversity, i add little constants to avoid math problems
                betas = []
                for pair in pairs:
                    beta = distance.braycurtis(tubes[pair[0]], tubes[pair[1]])
                    betas.append(beta)
                div_o_r.append(np.mean(betas))
            beta_div[seed][prop].append(div_o_r)
            

#reorganise the information
beta_div_all50 = {} #add each repeat for each seed independentluy

for prop in conditions:

    beta_div_all50[prop] = []
    for seed in seeds:

        beta_div_all50[prop].extend(beta_div[seed][prop])

#%%now the plot
plt.rcParams['svg.fonttype'] = "none"

#colors for the plots, I was looking for something in colorbrewere, but there are no categorical palettes with more than 8 colors and for colorblind peple,
#I liked set2,which is colorblind friendly, but only has 8 colors, in case we want to use it:
    #sns.color_palette("Set2")
#Then I found a blog where they recommended the following list.
#https://gist.github.com/thriveth/8560036
CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']
#this only has 9 colors, if we want more change it
#UPDATE, I could instead use the Set 2 and distinguish between selection (solid) and random (dashed)
colours = sns.color_palette("Set2")

t=range(rounds)
#plt.grid()
sns.set_style("darkgrid")
dict_of_2 = {} #to assing a colour I´ll take advantage of the 2 first caracters that are unique by method
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
    mean = [np.nanmean(x) for x in zip(*beta_div_all50[prop])]
    std = [np.nanstd(x) for x in zip(*beta_div_all50[prop])]

    mean = np.array(mean)
    std = np.array(std)
    plt.plot(t,mean,label=conditions_rename[i_p], color=color_here, linestyle= line_here, linewidth=3,alpha=0.8)
    #plt.fill_between(t, mean - std, mean + std, color=color_here, alpha=0.1, edgecolor='none')
    #i += 1
    plt.legend(ncol=4, fontsize=12.5,framealpha=0.5, title = "propagation")
    plt.ylim(0, 1.35)
    plt.title("10 repeats, 5 species sets", fontsize=16)
    plt.xlabel("Rounds",fontsize=16)
    plt.ylabel("Mean (final beta div B.C.)",fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
plt.savefig(path_save + "beta_div_bray_curtis_one_plot.pdf", bbox_inches='tight')
plt.savefig(path_save + "beta_div_bray_curtis_one_plot.svg", bbox_inches='tight')

plt.show()
