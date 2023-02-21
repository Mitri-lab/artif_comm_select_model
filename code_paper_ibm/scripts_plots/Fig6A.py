# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 11:15:01 2021

@author: pablo

We took the best community after 50 rounds with each propagation method and then
transfer them for 25 rounds without selection (NS propagation), here we plot the 
degradation, diversity or AUC over those 25 rounds. Each plot corresponds to a species set.

#update 20221011, corrected mistakes in one_plot version where we average all the species sets
"""

import pickle
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import os 
from math import sqrt,exp,log
from itertools import zip_longest

#where to find the data
path = "/Users/pablo/Desktop/code_paper/data/210822_2_ok/"

#some info about these results
prior_rounds = 50 #how many rounds have been done in the previous experiment
rounds = 25 #how many furter rounds we run without selection
repeats = 10
#All the propagations we are using disassembly (d), propagule(p), migrant_pool(m), no_select(n),
#either under selection (s) or random (r) treatment
#conditions = ["d2_s", "d2_s_anc", "p_s", "p_s_anc", "m_s", "m_s_anc", "pip_s", "pip_s_anc", "mip_s", "mip_s_anc", "n_r", "n_r_anc", "d_s", "d_s_anc" ]
conditions = ["d3_s", "d3_s_anc", "p_s", "p_s_anc", "m_s", "m_s_anc", "pip_s", "pip_s_anc", "mip_s", "mip_s_anc", "n_r", "n_r_anc"]
#to rename for the plot, important in the same order as above
#conditions_rename = ["DS", "DS_anc", "PS", "PS_anc", "MS", "MS_anc", "PIS", "PIS_anc", "MIS", "MIS_anc", "NS", "NS_anc", "D4S", "D4S_anc"]
conditions_rename = [ "DS", "DS_anc", "PS", "PS_anc", "MS", "MS_anc", "PIS", "PIS_anc", "MIS", "MIS_anc", "NS", "NS_anc"]
seeds = ["22", "23", "24", "25", "26"] #for the plot they will be renamed as 1-5

#where to stoer the results
path_save = "/Users/pablo/Desktop/220628_plots_mitri/"
#%% GENERATE DATA STRUCTURES
median_deg = {} 
max_deg = {}
mean_auc = {}
mean_div = {}

for seed in seeds:
    median_deg[seed] = {} #the keys correspond to the propagation method, inside they have a list, and
                      #this list has more lists inside, each one for one repeat. Each position of those lists
                      #shows the median degradation (over tubes) at the end of a round for each of the rounds.
    max_deg[seed] = {}
    mean_auc[seed] = {}
    mean_div[seed] = {}
    
    for prop in conditions:
        median_deg[seed][prop] = []
        max_deg[seed][prop] = []
        mean_auc[seed][prop] = []
        mean_div[seed][prop] = []
    

        for rp in range(repeats):
            
            df_here = pd.read_csv(path + seed+"/{}/repeat{}/df_grid.csv".format(prop,rp))
            deg_here = []
            max_here = []
            auc_here = []
            sp_div_here = []
            
            for rd in range (rounds):
                if rd in list(df_here["round"]):
                    deg = list(df_here.loc[df_here["round"] == rd, "deg_score"])
                    pop_0 = list(df_here.loc[df_here["round"] == rd, "tot_pop_0"])
                    auc = list(df_here.loc[df_here["round"] == rd, "tot_auc"])
                    div = list(df_here.loc[df_here["round"] == rd, "sp_div_H"])
                    deg_here.append(np.median(deg))
                    max_here.append(max(deg))
                    auc_here.append(np.mean(auc))#.append(np.mean([x/(y*101) for x,y in zip(auc,pop_0)]))
                    sp_div_here.append(np.mean([exp(x) for x in div]))
                #else:
                   # deg_here.append(np.nan)
            median_deg[seed][prop].append(deg_here)
            max_deg[seed][prop].append(max_here)
            mean_auc[seed][prop].append(auc_here)
            mean_div[seed][prop].append(sp_div_here)
            


#%%
#And now we also plot the summary plot (for med degradation)
plt.rcParams['svg.fonttype'] = "none" #needed to properly recognise letters of .svg plots in inkscape
colours = sns.color_palette("Set2")
sns.set_style("darkgrid")
def plot(dictionary, what, ylim0, ylim1):
    """

    Parameters
    ----------
    dictionary : the dictionary of what we want to plot
    what : a string of what we are pliting, just for labelling
    ylim0: lower limit of y axis
    ylim1: upper limit of y axis
    
    Returns
    -------
    None.

    """
    fig = plt.figure(figsize=(20,10))
    
    for i_s,seed in enumerate(seeds):
        plt.subplot(2, 3, i_s+1)
        #plt.grid()
        t=range(prior_rounds, prior_rounds + rounds)
        i=0
        #colors = {"d2":"#e09112", "p_":"#12b3b3", "m_":"#963b79", "n_":"gray", "d_": "black",  "mi": "#630f39", "pi":"#09613f" }
        dict_of_2={}
        i_color = 0
        linestyles = {"s":"-","r":"--"}
        for k,v in dictionary[seed].items():
            #if k == "d_s":
                #i_color +=1 #Otherwise D4S would have the color before corresponding to NS
            if k[:2] in dict_of_2: #if we have already appended the one under selection we assing same color to random
                color_here = dict_of_2[k[:2]]
            else: 
                color_here = colours[i_color]
                dict_of_2[k[:2]] = color_here
                i_color +=1
            #color_here = colors[k]
            style = "-"
            if "anc" in k:
                style = "dotted"

            #If some repeat is missing data here, go to previous versions of this code I solve that there. 
            med = [np.mean(x)for x in zip(*v)]
            std = [np.std(x)for x in zip(*v)]
            while len(med)<rounds:
                med.append(np.nan)
                std.append(np.nan)
            med = np.array(med)
            std = np.array(std)
            plt.plot(t,med,label=conditions_rename[i], linewidth=3, color = color_here, ls = style, alpha = 0.7)
            #plt.fill_between(t, med - std, med + std, alpha=0.1, edgecolor='none', color = color_here)
            i += 1
        if i_s == 4:
            plt.legend(ncol = 1, fontsize=14,framealpha=0.8, title = "Propagation", loc='center left', bbox_to_anchor=(1, 0.5))
        plt.ylim(ylim0, ylim1) 
        plt.xlabel("Further rounds",fontsize=18)
        plt.ylabel("Mean ({})".format(what),fontsize=18)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.title("Sp. set: " + str(i_s), fontsize =18)
        plt.tight_layout()

    plt.savefig( path_save + "transfers_of_best_com_{}_all.pdf".format(what), bbox_inches='tight')
    plt.savefig( path_save + "transfers_of_best_com_{}_all.svg".format(what), bbox_inches='tight')
    fig.show()
        
plot(max_deg, "max_deg", 0.2, 0.9)
plot(mean_auc, "mean_auc",20000, 180000)
plot(mean_div, "mean_div",1,5)
#%%PLOT AVERAGING ALL TEH 5 SPECIES SETS

#reorganise the information
median_deg_all = {}
max_deg_all = {}
mean_auc_all = {}
mean_div_all = {}
for prop in conditions:
    median_deg_all[prop] = []
    max_deg_all[prop] = []
    mean_auc_all[prop] = []
    mean_div_all [prop] = []
    for seed in seeds:
        median_deg_all[prop].extend(median_deg[seed][prop])
        max_deg_all[prop].extend(max_deg[seed][prop])
        mean_auc_all[prop].extend(mean_auc[seed][prop])
        mean_div_all[prop].extend(mean_div[seed][prop])


colours = sns.color_palette("Set2")
sns.set_style("darkgrid")
def plot(dictionary, what, ylim0, ylim1):
    """

    Parameters
    ----------
    dictionary : the dictionary of what we want to plot
    what : a string of what we are pliting, just for labelling
    ylim0: lower limit of y axis
    ylim1: upper limit of y axis
    
    Returns
    -------
    None.

    """
    fig = plt.figure(figsize=(8,5))
    all_seeds = {}
    t=range(prior_rounds, prior_rounds + rounds)
    dict_of_2={}
    i_color = 0
    linestyles = {"s":"-","r":"--"}
    i = 0
        #colors = {"d2":"#e09112", "p_":"#12b3b3", "m_":"#963b79", "n_":"gray", "d_": "black",  "mi": "#630f39", "pi":"#09613f" }

    for k,v in dictionary.items():
        
        
        #if k == "d_s":
            #i_color +=1 #Otherwise D4S would have the color before corresponding to NS
        if k[:2] in dict_of_2: #if we have already appended the one under selection we assing same color to random
            color_here = dict_of_2[k[:2]]
        else: 
            color_here = colours[i_color]
            dict_of_2[k[:2]] = color_here
            i_color +=1
        #color_here = colors[k]
        style = "-"
        if "anc" in k:
            style = "dotted"

        #If some repeat is missing data here, go to previous versions of this code I solve that there. 
        med = [np.mean(x)for x in zip(*v)]
        std = [np.std(x)for x in zip(*v)]
        while len(med)<rounds:
            med.append(np.nan)
            std.append(np.nan)
        med = np.array(med)
        std = np.array(std)
        plt.plot(t,med,label=conditions_rename[i], linewidth=3, color = color_here, ls = style, alpha = 0.7)
        #plt.fill_between(t, med - std, med + std, alpha=0.1, edgecolor='none', color = color_here)
        i += 1

    plt.legend(ncol = 1, fontsize=14,framealpha=0.8, title = "Propagation", loc='center left', bbox_to_anchor=(1, 0.5))
    plt.ylim(ylim0, ylim1) 
    plt.xlabel("Further rounds",fontsize=18)
    plt.ylabel("Mean ({})".format(what),fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.title("All Sp. sets, 10 repeats each" , fontsize =18)
    plt.tight_layout()

    plt.savefig( path_save + "transfers_of_best_com_{}_all_oneplot.pdf".format(what), bbox_inches='tight')
    plt.savefig( path_save + "transfers_of_best_com_{}_all_oneplot.svg".format(what), bbox_inches='tight')

        
plot(max_deg_all, "max_deg", 0.2, 0.9)
plot(mean_auc_all, "mean_auc",20000, 180000)
plot(mean_div_all, "mean_div",1,5)