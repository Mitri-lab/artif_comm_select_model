# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 08:20:07 2021

@author: pablo

Script to take the general results of the simulations, see which is the predomiant community at the end 
and check how this community will perform, regarding degradation, in the ranking of all the possible ancestral communities.
"""

import pickle
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import os 
from collections import Counter
from ast import literal_eval
import itertools
from scipy.stats import wilcoxon
from scipy.stats import shapiro

#path with the results of the evolution experiment
# path = "/Users/pablo/Desktop/data_mitri/210822/"
path = "/Users/pablo/Desktop/code_paper/data/210822/" 
#path for the ranking of ancestral communities (simulations of all the possible 32767 combinations of 15 species)
# path_rank = "/Users/pablo/Desktop/data_mitri/211005_all_communities/"
path_rank = "/Users/pablo/Desktop/code_paper/data/211005_all_communities/"

#some information about these results (210822)
num_species=15
num_repeats=10
num_tubes = 21
rd = 49
#The seeeds we have used for the rng generator and that name the directory to store resutls
seeds = ["22", "23", "24", "25", "26"] #renamed as 1-5 just for the plot
#Propagations to load
conditions = [ "d3_s", "p_s", "m_s",  "d3_r", "p_r", "m_r", "n_r", "pip_s", "pip_r", "mip_s", "mip_r"]
#for the final plots we will give different names than those from the model(here IN THE SAME ORDER!)
conditions_rename = [ "DS", "PS", "MS",   "DR", "PR", "MR",  "NS","PIS", "PIR", "MIS", "MIR"]
#next we specify which conditions to include in the plot (already renamed)
#conditions_to_plot =  [ "DS", "PS", "MS", "DR", "PR", "MR" , "NS"]
conditions_to_plot =  [ "PIS", "MIS"]



#path where to store the plots
# path_save = "/Users/pablo/Desktop/211111_results/"
path_save = "/Users/pablo/Desktop/220628_plots_mitri_b/"


maxrank = 2**num_species-1
#%%UPDATE ON 210712, NOW THE PREDOMINANT COMMUNITY IS JUST THE MOST ABUNDANT COMMUNITY
#I found the predominant combinations in the last round of each 
def most_frequent(List):
    occurence_count = Counter(List)
    return occurence_count.most_common(1)[0][0]


#LOAD PREDOMINANT COMMUNITIES
dict_combs = {}
for seed in seeds:
    dict_combs[seed] = {}
    for prop in conditions:
        list_combs = []
        for rp in range(num_repeats):
            comb_tubes=[]
            df = pd.read_csv(path + seed + "/" + prop + "/repeat{}/df_sp.csv".format(rp))
            df = df.loc[df["round"] ==rd] #like this we wont include those repeats that went extinct 
            #df = df.loc[df["round"] == list(df["round"])[-1]] #like that we always take the communities in the last round before extinction
            for tube in range(num_tubes):
                tube_here = df.loc[df["tube"] ==tube, "sp"]
                comb_tubes.append(tuple([int(x) for x in tube_here]))

            list_combs.append(most_frequent(comb_tubes))
        dict_combs[seed][prop] = list_combs
            #plt.hist(df_here.loc[df_here["round"] ==49, "sp"])
            #plt.show()


    #and let´s save dict_combs
    open_file = open(path + "dict_combs", "wb")
    pickle.dump(dict_combs, open_file)
    open_file.close()
    


#%%CHECK THE POSITION OF THE SELECTED COMBINATIONS IN THE ANCESTRAL RANKING

#load all possible combinations
positions = {}
for seed in seeds:
    positions[seed] = {}
    df_all = pd.read_csv(path_rank + seed + "/df.csv")
    df_all['com'] = df_all.com.apply(lambda x: literal_eval(str(x)))
    #only quartets 
    #df_all = df_all[df_all['com'].map(len)<5]
    #order communities by descending degradation score
    df_all = df_all.sort_values("deg_score", ascending = False) 
    df_all = df_all.reset_index() 
    
    for prop in conditions:
        positions[seed][prop] = []
        for comb in dict_combs[seed][prop]:
            try:
                row = df_all.loc[df_all["com"] == tuple([int(x) for x in comb])]
                positions[seed][prop].append(row.index[0])
            except:
                print(seed,prop,comb)
                positions[seed][prop].append(maxrank)
                # positions[seed][prop].append(np.nan)
                



#Accomodate the information in a dataframe as it will be easier to plot with sns
#when taking the log for the plot, I avoid problems by adding +1
columns_here= ["seed","prop","repeat","position", "log2_position"]
df = pd.DataFrame( columns = columns_here)
seed_name = 1 #we rename seeds so they appear in the plot starting at 1.
for seed in seeds:
    for i_p, prop in enumerate(conditions):
        for repeat in range(num_repeats):
            df_now = pd.DataFrame([[str(seed_name),
                                     conditions_rename[i_p],
                                     repeat,
                                     positions[seed][prop][repeat],
                                     np.log2(positions[seed][prop][repeat]+1)]], 
                                    columns = columns_here)#.fillna(0,inplace=True)
            df = pd.concat([df,df_now],ignore_index=True)
    seed_name += 1


#%% BOXPLOT

#for some reason I can´t find, when importing graphs to inkscape in pdf, it does not show the text, so I will use
#istead .svg, but there I text will be an icon, unless I set the folowing:
plt.rcParams['svg.fonttype'] = "none"
#Make a data frame including just the conditions to plot
df_plot = df[df['prop'].isin(conditions_to_plot)]
df_plot.reset_index()

def my_boxplot(x_, y_):
    #make boxplot
    fig=plt.figure(figsize=(len(conditions_to_plot)*1.2, 5))
    #for two conditions *1.8 (and no legend), otherwise *1.2
    
    
    sns.set(style="ticks")
    sns.set_style("darkgrid")
    sns.set_context( rc={"grid.linewidth": 2})
    bp= sns.boxplot(x=x_, y=y_, 
                     data=df_plot,palette="colorblind", 
                      hue='seed', 
                      boxprops={"zorder":1},
                      whiskerprops={"zorder":10},
                      linewidth=2,zorder=2)
    # make grouped stripplot
    #bp= sns.stripplot(y=y_, x=x_, data=df_plot,jitter=True,dodge=True, marker='o', alpha=0.5,hue='seed', color='black')
    
    # get legend information from the plot object (to not repeat if we include stripplot)
    handles, labels = bp.get_legend_handles_labels()
    #plt.ylim(-1, 21)
    # specify just one legend
    #plt.legend(handles[0:5], labels[0:5], title="sp_seed", fontsize=16, ncol=5, loc="lower right", framealpha=0.5)
    #plt.legend("")
    plt.xticks(fontsize= 28)
    plt.yticks(fontsize=28)
    plt.title("P.C. at round 50", fontsize = 28)
    plt.ylabel("$log_{2}$(rank)", fontsize = 28)
    plt.xlabel("Propagation", fontsize = 28)
    plt.tight_layout()
    plt.savefig(path_save + "boxplot{}.svg".format(y_))
    plt.savefig(path_save + "boxplot{}.pdf".format(y_))
    
def my_swarmplot(x_, y_):
    #make swarmplot
    fig=plt.figure(figsize=(len(conditions_to_plot)*1.8, 5))
    #for two conditions *1.8 (and no legend), otherwise *1.2
    
    
    sns.set(style="ticks")
    sns.set_style("darkgrid")
    sns.set_context( rc={"grid.linewidth": 2})
    bp= sns.swarmplot(x=x_, y=y_, 
                      data=df_plot,palette="colorblind", 
                      hue='seed', 
                      dodge=False,size=6)
    # make grouped stripplot
    #bp= sns.stripplot(y=y_, x=x_, data=df_plot,jitter=True,dodge=True, marker='o', alpha=0.5,hue='seed', color='black')
    
    # get legend information from the plot object (to not repeat if we include stripplot)
    handles, labels = bp.get_legend_handles_labels()
    plt.ylim(-1, 16)
    # specify just one legend
    #plt.legend(handles[0:5], labels[0:5], title="sp_seed", fontsize=16, ncol=5, loc="lower right", framealpha=0.5)
    #plt.legend("")
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.title("P.C. at round 50", fontsize = 18)
    plt.ylabel("$log_{2}$(rank)", fontsize = 18)
    plt.xlabel("Propagation", fontsize = 18)
    plt.tight_layout()
    plt.savefig(path_save + "boxplot{}_PISMIS.svg".format(y_))
    plt.savefig(path_save + "boxplot{}_PISMIS.pdf".format(y_))

# my_boxplot("prop", "log2_position")
my_swarmplot("prop", "log2_position")

### Print min rank for each condition
dfmin = df.groupby(["seed","prop","repeat"]).min().reset_index()
dfmin.loc[:,("seed","prop","repeat","log2_position")].to_csv(path_save + "dfmin_ranks.csv")

#%%STATISTICAL TESTING
#check for normality
for prop in conditions_rename:
    x = df.loc[df["prop"] == prop,"position"]
    print(shapiro(x))
    
#since most of them are not normally distributed I´ll use wilkinson test to compare

#COMPARE EACH PROP SELECTION AGAINST CONTROL
props = conditions_rename
what = "position"
combs = list(itertools.combinations(props, 2))
df_sw = pd.DataFrame( columns = ["prop1", "prop2","sw", "p_val"] )
for comb in combs:
    one = df.loc[df["prop"] == comb[0], what]
    two = df.loc[df["prop"] == comb[1], what]
    # d = np.subtract(one,two)
    sw, p = wilcoxon(x=one,y=two)
    df_now = pd.DataFrame([[comb[0], comb[1], sw, p]], columns = ["prop1", "prop2","sw", "p_val"] )
    df_sw = pd.concat([df_sw,df_now],ignore_index=True)
df_sw.to_csv(path_save + "dfw_rankinf_vs_all2.csv")
