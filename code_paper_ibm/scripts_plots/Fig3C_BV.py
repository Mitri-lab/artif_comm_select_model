# -*- coding: utf-8 -*-
"""
Created on Sat Jun 26 19:35:24 2021

@author: pablo


Script to plot diversity, caculated as exp(shanon index). We plot the three diversities.
    - alpha = mean diversity inside each tube
    - gamma = diversity in the metacommunity (if we don´t consider tubes but mix all the cells together)
    - beta = variance between tubes (gamma/alpha)
In this version I combine in one plot all the information, I do the mean including all repeats and all seeds.
(for now just including all 50 datasets, like, considering each repeat of each seed an independent datapoint,
 like this I can show the sd if I want, although I guess it will be pretty messy and at the end I will not display it)

Notes: 
    added np.nanmean and np.nanstd when plotting, this is, if some values are missing for some 
    repeat still do the mean and std with the remaining values

211201 for alpha diversity, to take the mean of teh tubes, instead of exp(mean(H)) that I was doing before,
    now I changed to mean(exp(H))
    Explanation by Björn: Let’s say that the species i=1,2,34 each have a relative abundance pi. Then the Shannon diversity is H=-sum_i(pi log(pi)) as always, and the true diversity in this community is D = exp(H).
                        The alpha diversity is the average diversity between communities, which could be either mean(H) or mean(D). But, in order to calculate beta = gamma/alpha we want to use the true diversities, which should be mean( exp(H) )
                        Averaging the Shannon first, exp( mean(H) ) forms a geometric (instead of arithmetic) average between communities for the alpha diversity, which could make a difference in rare cases but I don’t think it makes a lot of difference in practice.
    (I´ll add _ok at the end of the name of plots with this correction, gamma shouldn´t change but I also rename it just to check everything is ok)

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
# import pingouin as pg

import itertools


#path to data
#path = "/Users/pablo/Desktop/data_mitri/210822/"
path = "/Users/pablo/Desktop/code_paper/data/210822/"  
#some info about this dataset
rounds = 50
num_tubes = 21
repeats = 10
#seeds I use for the random number generator, for the plot they are renamed as 1-5
seeds = ["22", "23", "24", "25", "26"]
#propagation methods we use
conditions = ["d3_s", "p_s", "m_s", "pip_s","mip_s", "n_r", "d3_r", "p_r", "m_r", "pip_r" , "mip_r"] 
#IMPORTANT: in the same order as above
conditions_rename = ["DS", "PS", "MS", "PIS", "MIS", "NS", "DR", "PR", "MR", "PIR", "MIR" ]


#where to save the plots
#path_save = "/Users/pablo/Desktop/data_mitri/211111_results/"
path_save = "/Users/pablo/Desktop/220628_plots_mitri_b/"

#%% load data
gamma_div = {} #diversity in the metacommunity
alpha_div = {}

for seed in seeds:
    gamma_div[seed] = {}
    alpha_div[seed] = {}
    
    for i_p,prop in enumerate(conditions):
        gamma_div[seed][prop] = []
        alpha_div[seed][prop] = []
        for repeat in range(repeats):
            gamma_div_o_r = []
            alpha_div_o_r = []
            df_sp = pd.read_csv(path + seed + "/" + prop + "/" + "repeat" + str(repeat) + "/df_sp.csv")
            df_grid = pd.read_csv(path + seed + "/" + prop + "/" + "repeat" + str(repeat) + "/df_grid.csv")
            for rd in range(rounds):
                div = list(df_grid.loc[df_grid["round"] ==rd, "sp_div_H"])
                alpha_div_o_r.append(np.mean([np.exp(x) for x in div]))
                df_here = df_sp.loc[df_sp["round"] == rd]
                all_sp = {}
                
                for tube in range(num_tubes):
                    sps = list(df_here.loc[df_here["tube"]==tube, "sp"])
                    pops = list(df_here.loc[df_here["tube"]==tube, "final_pop"])
                    for i in range(len(sps)):
                        try:
                            all_sp[str(sps[i])] += pops[i]
                        except:
                            all_sp[str(sps[i])] = pops[i]
                #calculate diversity, i add little constants to avoid math problems
                H_sp = -1*sum([v/(sum(all_sp.values())+0.001)*log((v+0.001)/(sum(all_sp.values()) + 0.001)) for v in all_sp.values()])
                gamma_div_o_r.append(np.exp(H_sp))
            gamma_div[seed][prop].append(gamma_div_o_r)
            alpha_div[seed][prop].append(alpha_div_o_r)

#calculate beta diversity (I have another script where I do it by bray curtis, just in case)
# beta = gamma/alpha (expressing each of them already as effective number of species)
beta_div = {}
for seed in seeds:
    beta_div[seed] = {}
    for prop in conditions:
        beta_div[seed][prop] = []
        for rp in range(repeats):
            beta_div[seed][prop].append(np.divide(gamma_div[seed][prop][rp], alpha_div[seed][prop][rp]))
            
#%%Reorganise the information for this plot, remove the distinction between seeds
gamma_div_all50 = {} #add each repeat for each seed independentluy
alpha_div_all50 = {}
beta_div_all50 = {}


for prop in conditions:
    gamma_div_all50[prop] = []
    alpha_div_all50[prop] = []
    beta_div_all50[prop] = []
    for seed in seeds:
        gamma_div_all50[prop].extend(gamma_div[seed][prop])
        alpha_div_all50[prop].extend(alpha_div[seed][prop])
        beta_div_all50[prop].extend(beta_div[seed][prop])


#%%PLOT
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
#GAMMA DIVERSITY

t=range(rounds)
#plt.grid()
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
    
    mean = [np.nanmean(x)for x in zip(*gamma_div_all50[prop])]
    std = [np.nanstd(x)for x in zip(*gamma_div_all50[prop])]

    mean = np.array(mean)
    std = np.array(std)
    plt.plot(t,mean,label=conditions_rename[i_p], color=color_here, linestyle=line_here, linewidth=3,alpha=0.8)
    #plt.fill_between(t, mean - std, mean + std, color=color_here, alpha=0.1, edgecolor='none')
    i += 1
    plt.legend(ncol=4, fontsize=12.5,framealpha=0.5, title = "Propagation")
    plt.ylim(1, 11)
    plt.title("10 repeats, 5 species sets", fontsize=16)
    plt.xlabel("Rounds",fontsize=16)
    plt.ylabel("Mean (final gamma div)",fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    
plt.savefig(path_save + "gamma_div_one_plot_ok.pdf",  bbox_inches='tight')
plt.savefig(path_save + "gamma_div_one_plot_ok.svg",  bbox_inches='tight')
plt.show()

#ALPHA DIVERSITY

t=range(rounds)
#plt.grid()
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
    mean = [np.nanmean(x)for x in zip(*alpha_div_all50[prop])]
    std = [np.nanstd(x)for x in zip(*alpha_div_all50[prop])]

    mean = np.array(mean)
    std = np.array(std)
    plt.plot(t,mean,label=conditions_rename[i_p], color=color_here, linestyle=line_here, linewidth=3,alpha=0.8)
    #plt.fill_between(t, mean - std, mean + std, color=color_here, alpha=0.1, edgecolor='none')
    i += 1
    plt.legend(ncol=4, fontsize=12.5,framealpha=0.5, title = "Propagation")
    plt.ylim(1, 5.2)
    plt.title("10 repeats, 5 species sets", fontsize=16)
    plt.xlabel("Rounds",fontsize=16)
    plt.ylabel("Mean (final alpha div)",fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
plt.savefig(path_save +"alpha_div_one_plot_ok.pdf", bbox_inches='tight')
plt.savefig(path_save +"alpha_div_one_plot_ok.svg", bbox_inches='tight')
plt.show()

#BETA DIVERSITY

#This one is already the exponent
t=range(rounds)
#plt.grid()
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
    mean = [np.nanmean(x) for x in zip(*beta_div_all50[prop])]
    std = [np.nanstd(x) for x in zip(*beta_div_all50[prop])]

    mean = np.array(mean)
    std = np.array(std)
    plt.plot(t,mean,label=conditions_rename[i_p], color=color_here, linestyle = line_here, linewidth=3,alpha=0.8)
    #plt.fill_between(t, mean - std, mean + std, color=color_here, alpha=0.1, edgecolor='none')
    i += 1
    plt.legend(ncol=4, fontsize=12.5,framealpha=0.5, title = "Propagation")
    plt.ylim(0, 7)
    plt.title("10 repeats, 5 species sets", fontsize=16)
    plt.xlabel("Rounds",fontsize=16)
    plt.ylabel("Mean (final beta div)",fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
plt.savefig(path_save + "beta_div_one_plot_ok.pdf", bbox_inches='tight')
plt.savefig(path_save + "beta_div_one_plot_ok.svg", bbox_inches='tight')
plt.show()
            
