# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 11:15:01 2021

@author: pablo

Script to plot results of the simulations as boxplots for the 10 points of each repeat, separating by species sets.
In particular we show the difference in performance between a given round (usually last) and the first one.
Usually we plot degradation, but here we also retrieve information to plot growth (AUC), or alpha diversity(shanon)

"""

import pickle
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import os 
from math import sqrt,exp,log
import itertools
from scipy.stats import wilcoxon
from scipy.stats import shapiro



#path where the data we are using is stored

path = "/Users/pablo/Desktop/data_mitri/210822/"   

#information about these data
num_reps = 10 #number of repeats
rounds = 50 #manually specify all the rounds we have to store data
rd_plot = 49 #At which round to plot; -1 for the last round regardless of which it is (including extinct repeats)
#seeds we used (names of the folders)
seeds = ["22", "23", "24", "25", "26"]
#in the plots seeds will be renamed from 1-len(seeds), change seed_name below to avoid this

conditions = [ "d3_s", "p_s", "m_s",  "d3_r", "p_r", "m_r", "n_r", "pip_s", "pip_r", "mip_s", "mip_r"]
#for the final plots we will give different names than those from the model(here IN THE SAME ORDER!)
conditions_rename = ["DS", "PS", "MS",   "DR", "PR", "MR",  "NS","PIS", "PIR", "MIS", "MIR"]
#next we specify which conditions to include in the plot (already renamed)
conditions_to_plot =  ["DS", "PS", "MS", "DR", "PR", "MR" , "NS"]
#conditions_to_plot =  [ "PIS", "MIS"]

#path where we want to save the plots and df with test
# path_save = "/Users/pablo/Desktop/211111_results/"
path_save = "/Users/pablo/Desktop/220628_plots_mitri_b/"




#%% LOAD THE DATA
#create a data frame containing the following information
columns_df= ["seed","prop","repeat","max_deg_dif","med_deg_dif", "auc_dif",  "div_dif"]
df = pd.DataFrame( columns = columns_df)
#rename the seeds 1-5
seed_name = 1
for seed in seeds:
    #All the propagations we are using disassembly (d), propagule(p), migrant_pool(m), no_select(n),
    #either under selection (s) or random (r) treatment
    

    median_deg = {} #the keys correspond to the propagation method, inside they have a list, and
                          #this list has more least inside, each one for one repeat. Each of those least
                          #show the median degradation (over tubes) at the end of a round for each of
                          #the rounds.
    max_deg = {}
    mean_auc = {}
    mean_div = {}
    for prop in conditions:
        median_deg[prop] = []
        max_deg[prop] = []
        mean_auc[prop] = []
        mean_div[prop] = []
    #Now we loop trough the files
    for prop in conditions:
        for rp in range(num_reps):
            df_here = pd.read_csv(path + seed + "/" + prop +  "/repeat" + str(rp) + "/df_grid.csv")
            #to store info over rounds
            deg_here = []
            max_here = []
            auc_here = []
            sp_div_here = []
            for rd in range (rounds): #could just check the rounds in df_here["round"], but anyway
                if rd in list(df_here["round"]):
                    #create some list with info for the 21 tubes in this repeat and round
                    deg = list(df_here.loc[df_here["round"] == rd, "deg_score"])
                    pop_0 = list(df_here.loc[df_here["round"] == rd, "tot_pop_0"])
                    auc = list(df_here.loc[df_here["round"] == rd, "tot_auc"])
                    div = list(df_here.loc[df_here["round"] == rd, "sp_div_H"])
                    #add info to over rounds
                    deg_here.append(np.median(deg))
                    max_here.append(max(deg))
                    auc_here.append(np.mean([x-(y*rd_plot) for x,y in zip(auc,pop_0)])) #in auc we discount area initial pop
                    #auc_here.append(np.mean([x/(y*101 +0.000001) for x,y in zip(auc,pop_0)]))
                    sp_div_here.append(np.mean([exp(x) for x in div]))
    
            median_deg[prop].append(deg_here)
            max_deg[prop].append(max_here)
            mean_auc[prop].append(auc_here)
            mean_div[prop].append(sp_div_here)

    #create the dataframe with the difference between round_plot and initial round
    for i_p, prop in enumerate(conditions):
        for repeat in range(num_reps):
            df_now = pd.DataFrame([[str(seed_name),
                                     conditions_rename[i_p],
                                     repeat,
                                     max_deg[prop][repeat][rd_plot]-max_deg[prop][repeat][0],
                                     median_deg[prop][repeat][rd_plot]-median_deg[prop][repeat][0],
                                     mean_auc[prop][repeat][rd_plot]-mean_auc[prop][repeat][0],
                                     mean_div[prop][repeat][rd_plot]-mean_div[prop][repeat][0]]], 
                                    columns = columns_df)
                                     
    
            df = pd.concat([df,df_now])
    seed_name += 1

#%%PLOT

#for some reason I can´t find, when importing graphs to inkscape in pdf, it does not show the text, so I will use
#istead .svg, but there I text will be an icon, unless I set the folowing:
plt.rcParams['svg.fonttype'] = "none"
#Make a data frame including just the conditions to plot
df_plot = df[df['prop'].isin(conditions_to_plot)]

def my_boxplot(x, y):
    #make boxplot
    #sns.set(style="ticks")
    sns.set_style("darkgrid")
    sns.set_context( rc={"grid.linewidth": 2})
    fig=plt.figure(figsize=(len(conditions_to_plot)*1.2, 5))
    #if only 2 conditions *1.8, otherwise *1.2
    plt.axhline(y=0, color='black', linestyle='dotted', alpha = 0.5, lw=2, zorder=1) 

    bp= sns.boxplot(x, y=y, 
                     data=df_plot, 
                     palette="colorblind", 
                      hue='seed', 
                      boxprops={"zorder":1},
                      whiskerprops={"zorder":10},
                      linewidth=2,zorder=2) #zorder does not work here though
    #showfliers= False #in case we add the stripplot bellow
    #bp= sns.stripplot(y=y, x=x, data=df_plot, jitter=True,dodge=True, marker='o', alpha=0.6,hue='seed',color='black')
    
    # get legend information from the plot object
    handles, labels = bp.get_legend_handles_labels()
    plt.ylim(-0.9, 0.7)
    # specify just one legend
    plt.legend(handles[0:5], labels[0:5], title="sp_set", fontsize=16, ncol=5, loc="lower left", framealpha=0.5)
    #plt.legend("")
    plt.xticks(fontsize= 28)
    plt.yticks(fontsize=28)
    plt.title("Diff. in 50 rounds", fontsize = 28)
    plt.ylabel("Max. degradation", fontsize = 28)
    plt.xlabel("Propagation", fontsize = 28)
    plt.tight_layout()
    plt.savefig(path_save + "boxplot{}_.svg".format(y))
    plt.savefig(path_save + "boxplot{}_.pdf".format(y))
    
#remember to modify the y axes label according to what you want to plot.
def my_swarmplot(x, y):
    # Make swarm plot
    sns.set_style("darkgrid")
    sns.set_context( rc={"grid.linewidth": 2})
    fig=plt.figure(figsize=(len(conditions_to_plot)*1.2, 5))
    #if only 2 conditions *1.8, otherwise *1.2
    plt.axhline(y=0, color='black', linestyle='dotted', alpha = 0.5, lw=4, zorder=1) 

    bp= sns.swarmplot(x=x, y=y, hue='seed', 
                      data=df_plot, 
                      palette="colorblind", 
                      dodge=False,size=6) 
    
    # get legend information from the plot object
    handles, labels = bp.get_legend_handles_labels()
    plt.ylim(-0.9, 0.7)
    # specify just one legend
    plt.legend(handles[0:5], labels[0:5], title="sp_set", fontsize=16, ncol=5, loc="lower left", framealpha=0.5)
    #plt.legend("")
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    # plt.title("Diff. in 50 rounds", fontsize = 28)
    plt.ylabel("Difference in max. degradation", fontsize = 18)
    plt.xlabel("Method", fontsize = 18)
    plt.tight_layout()
    plt.savefig(path_save + "swarmplot{}_.svg".format(y))
    plt.savefig(path_save + "swarmplot{}_.pdf".format(y))

what = "max_deg_dif"
my_swarmplot("prop", what)
#my_swarmplot("prop", "max_deg_dif")
#my_swarmplot("prop", "med_deg_dif")
#my_boxplot("prop", "max_deg_dif")
#my_boxplot("prop", "med_deg_dif")
#my_boxplot("prop", "auc_dif")
#my_boxplot("prop", "div_dif")


#%%STATISTICAL TESTING
for prop in conditions_rename:
    x = df.loc[df["prop"] == prop,what]
    print(shapiro(x))
#most of them con be normal, but just to be safe I´ll do wilcoxon


#EVERYONE VS NS
df_sw = pd.DataFrame( columns = ["prop1", "prop2","sw", "p_val"] )
#chose what to plot from the df
for prop in conditions_rename:
    if prop == "NS":
        continue
    one = df.loc[df["prop"] == prop, what]
    two = df.loc[df["prop"] == "NS",what]
    d = np.subtract(one,two)
    sw, p = wilcoxon(d)
    df_now = pd.DataFrame([[prop, "n_r", sw, p]], columns = ["prop1", "prop2","sw", "p_val"] )
    df_sw = pd.concat([df_sw,df_now])
df_sw.to_csv(path_save + "dfw_vs_NS_"+what+".csv")

#EVERYONE VS ROUND 0
df_sw = pd.DataFrame( columns = ["prop", "sw", "p_val"] )
#chose what to plot from the df
for prop in conditions_rename:
    one = df.loc[df["prop"] == prop, what]
    sw, p = wilcoxon(one,alternative='greater')
    df_now = pd.DataFrame([[prop, sw, p]], columns = ["prop", "sw", "p_val"] )
    df_sw = pd.concat([df_sw,df_now])
df_sw.to_csv(path_save + "dfw_vs_round0_"+what+".csv")

#COMPARE EACH PROP SELECTION AGAINST CONTROL
props = ["D", "P", "M", "MI", "PI"]
df_sw = pd.DataFrame( columns = ["prop1", "prop2","sw", "p_val"] )
for prop in props:
    one = df.loc[df["prop"] == prop +"S", what]
    two = df.loc[df["prop"] == prop + "R",what]
    d = np.subtract(one,two)
    sw, p = wilcoxon(d)
    df_now = pd.DataFrame([[prop +"S", prop + "R", sw, p]], columns = ["prop1", "prop2","sw", "p_val"] )
    df_sw = pd.concat([df_sw,df_now])
df_sw.to_csv(path_save + "dfw_S_vs_R_"+what+".csv")

#COMPARE PROPAGATIONS UNDER SELECTION
props = ["DS", "PS", "MS", "PIS", "MIS"]
combs = list(itertools.combinations(props, 2))
df_sw = pd.DataFrame( columns = ["prop1", "prop2","sw", "p_val"] )
for comb in combs:
    one = df.loc[df["prop"] == comb[0],what]
    two = df.loc[df["prop"] == comb[1],what]
    d = np.subtract(one,two)
    sw, p = wilcoxon(d)
    df_now = pd.DataFrame([[comb[0], comb[1], sw, p]], columns = ["prop1", "prop2","sw", "p_val"] )
    df_sw = pd.concat([df_sw,df_now])
df_sw.to_csv(path_save + "dfw_vs_S_"+what+".csv")

#%%TESTING MED
"""

for prop in conditions_rename:
    x = df.loc[df["prop"] == prop,"med_deg_dif"]
    print(shapiro(x))
#most of them con be normal, but just to be safe I´ll do wilcoxon

#EVERYONE VS NS
df_sw = pd.DataFrame( columns = ["prop1", "prop2","sw", "p_val"] )
#chose from the dataframe which column to plot
what = "max_deg_dif"
for prop in conditions_rename:
    if prop == "NS":
        continue
    one = df.loc[df["prop"] == prop, what]
    two = df.loc[df["prop"] == "NS",what]
    d = np.subtract(one,two)
    sw, p = wilcoxon(d)
    df_now = pd.DataFrame([[prop, "NS", sw, p]], columns = ["prop1", "prop2","sw", "p_val"] )
    df_sw = df_sw.append(df_now)
df_sw.to_csv(path_save + "df_w_vs_NS_{}.csv".format(what))

#COMPARE EACH PROP SELECTION AGAINST CONTROL
props = ["D", "P", "M", "PI", "MI"]
df_sw = pd.DataFrame( columns = ["prop1", "prop2","sw", "p_val"] )
what = "max_deg_dif"
for prop in props:
    one = df.loc[df["prop"] == prop +"S", what]
    two = df.loc[df["prop"] == prop + "R",what]
    d = np.subtract(one,two)
    sw, p = wilcoxon(d)
    df_now = pd.DataFrame([[prop +"S", prop + "R", sw, p]], columns = ["prop1", "prop2","sw", "p_val"] )
    df_sw = df_sw.append(df_now)
df_sw.to_csv(path_save + "df_w_S_vs_R_{}.csv".format(what))

#COMPARE PROPAGATIONS UNDER SELECTION
props = ["DS", "PS", "MS", "PIS", "MIS"]
what = "max_deg_dif"
combs = list(itertools.combinations(props, 2))
df_sw = pd.DataFrame( columns = ["prop1", "prop2","sw", "p_val"] )
for comb in combs:
    one = df.loc[df["prop"] == comb[0], what]
    two = df.loc[df["prop"] == comb[1],what]
    d = np.subtract(one,two)
    sw, p = wilcoxon(d)
    df_now = pd.DataFrame([[comb[0], comb[1], sw, p]], columns = ["prop1", "prop2","sw", "p_val"] )
    df_sw = df_sw.append(df_now)
df_sw.to_csv(path_save + "df_w_vs_S_{}.csv".format(what))

"""
    
