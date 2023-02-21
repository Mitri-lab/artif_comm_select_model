# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 08:20:07 2021

@author: pablo



Boxplots of mean sum(f) of the species in the community, distienguishing by sp.set
when calculating the sumf I do not weight by the species population. I do the mean sumf per species
(here yes weighting by population), and then I do the mean of all the species as long as they are present in the tube at the last round.

Here we add as a red dotted line the mean sum(f) per community obtained if we do random combinations of communities.
(I´t is really close to 0.5, which actually makes sense)
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
#path = "/Users/pablo/Desktop/data_mitri/210822/"
path = "/Users/pablo/Desktop/code_paper/data/210822/" 

#information about these simulations
num_repeats=10
num_tubes = 21
last_round = 49
#The seeeds we have
seeds = ["22", "23", "24", "25", "26"] #when ploting they are renamed starting from 1
#Propagations to load
#Propagations to load
conditions = [ "d3_s", "p_s", "m_s", "d3_r",  "p_r", "m_r", "n_r", "pip_s", "pip_r", "mip_s", "mip_r"]
#for the final plots we will give different names than those from the model(here IN THE SAME ORDER!)
conditions_rename = [ "DS", "PS", "MS",  "DR", "PR", "MR",  "NS","PIS", "PIR", "MIS", "MIR"]
#next we specify which conditions to include in the plot (already renamed)
conditions_to_plot =  [ "DS", "PS", "MS", "DR", "PR", "MR" , "NS"]
#conditions_to_plot =  [ "PIS", "MIS"]



#path where to store the plots
#path_save = "/Users/pablo/Desktop/211111_results/"
path_save = "/Users/pablo/Desktop/220628_plots_mitri_b/"

#%%Calculate the random value of sumf
#we take the combination of cells at round 0(could also just be a random combination of species, for now I leave this)
#and time 0 so that we don´t have to worry about weighting by strain population or not.
random_sumf=[]
prop = "d3_s" #the propagation does not really matter, all the propagations will have the same

for seed in seeds:
    #load the species for this seed
    open_file = open(path + "{}/{}/Species".format(seed,prop), "rb")
    sp_here = pickle.load(open_file)
    open_file.close()
    
    #load the species combinatios for this seed
    open_file = open(path + "{}/{}/index_repeats".format(seed,prop), "rb")
    i_rp_here = pickle.load(open_file)
    open_file.close()
    
    #calculate the mean
    
    for rp in i_rp_here:
        for tube in rp:
            for sp in tube:
                random_sumf.append(sum(sp_here[sp]["fs"]))
                
        


#%%LOAD SUM IN THE LAST ROUND
dict_sumf = {} #here it will have, for each seed and propagation, it will have a value for each tube, but this value
#is the result of doing a simple mean over the species in the tube, just by presence without considering their population.
#we agreed however on taking a mean value for each species so they don´t have more weight if they have more strians,
#and there, I do take into account the population of each strain (for the mean per species)
for seed in seeds:
    dict_sumf[seed] = {}
    for prop in conditions:
        dict_sumf[seed][prop] = [] #inside there will be one list for each repeat
        for rp in range(num_repeats):
            sumf_rep = [] #it will have one value for each tube in a repeat(later on, we take the mean of all the tubes in a repeat)
            df=pd.read_csv(path + seed + "/" + prop + "/repeat{}/df_st.csv".format(rp))
            rd =  list(set(df["round"]))[-1]
            if rd != last_round:
                print(seed, prop, rd) #we need to know if the simulations ended earlier at a given propagation
            df = df.loc[df["round"]== rd]
            #turn the column of fs values into tupple if I haven´t done that in the main scriot already
            df['st'] = df.st.apply(lambda x: literal_eval(str(x)))
            #obtain sumf and weigting (population)
            for tube in range(num_tubes):
                df_tube = df.loc[df["tube"] == tube]
                if len(df_tube)>0:
                    tube_mean_sumf = []
                    for sp in set(df_tube["sp"]): #inside the species we do weight by the population of each strain
                        df_tube_here = df_tube.loc[df_tube["sp"]==sp]
                        sumf = [np.sum(fs) for fs in df_tube_here["st"]]
                        populations = list(df_tube_here["final_pop"])
                        tube_mean_sumf.append(np.average(sumf, weights= populations))
                    sumf_rep.append(np.mean(tube_mean_sumf))
                else:
                    sumf_rep.append(np.nan)
            dict_sumf[seed][prop].append(sumf_rep) 
            


#%%CREATE DF

#Accomodate the information in a dataframe as it will be easier to plot with sns
#also for sumf we take the mean of the 21 tubes in each repeat (ignoring nan values: if a tube went extinct)
columns_here= ["seed","prop","repeat","mean_sumf"]
df = pd.DataFrame( columns = columns_here)
seed_name = 1 #we rename seeds so they appear in the plot starting at 1.
for seed in seeds:
    for i_p, prop in enumerate(conditions):
        for repeat in range(num_repeats):
            df_now = pd.DataFrame([[str(seed_name),
                                     conditions_rename[i_p],
                                     repeat,
                                     np.nanmean(dict_sumf[seed][prop][repeat])]], 
                                    columns = columns_here)
            df = pd.concat([df,df_now])
    seed_name += 1


#%% BOXPLOT

#for some reason I can´t find, when importing graphs to inkscape in pdf, it does not show the text, so I will use
#istead .svg, but there I text will be an icon, unless I set the folowing:
plt.rcParams['svg.fonttype'] = "none"
#Make a data frame including just the conditions to plot
df_plot = df[df['prop'].isin(conditions_to_plot)]

def my_boxplot(x_, y_):
    #make boxplot
    fig=plt.figure(figsize=(len(conditions_to_plot)*1.2, 5))
    #normally *1.2
    #to plot just two conditions *1.8, and then readapt the figure in inkscape..., no legend in this case
    
    sns.set(style="ticks")
    sns.set_style("darkgrid")
    sns.set_context( rc={"grid.linewidth": 2})
    bp= sns.boxplot(x=x_, y=y_, 
                     data=df_plot,palette="colorblind", 
                      hue='seed', 
                      boxprops={"zorder":1},
                      whiskerprops={"zorder":10},
                      linewidth=2,zorder=2)
    plt.axhline(y=np.mean(random_sumf), color='red', linestyle='dotted', alpha = 0.8, lw=2, zorder=1) 
    # make grouped stripplot
    #bp= sns.stripplot(y=y_, x=x_, data=df_plot,jitter=True,dodge=True, marker='o', alpha=0.5,hue='seed', color='black')
    
    # get legend information from the plot object (to not repeat if we include stripplot)
    handles, labels = bp.get_legend_handles_labels()
    plt.ylim(0, 1)
    # specify just one legend
    plt.legend(handles[0:5], labels[0:5], title="Sp. set:", fontsize=13, ncol=1, loc="upper right", framealpha=0.5)
    #plt.legend("")
    plt.xticks(fontsize= 28)
    plt.yticks(fontsize=28)
    #plt.title("round 50", fontsize = 28)
    plt.ylabel("Average total investment", fontsize = 28)
    plt.xlabel("Method", fontsize = 28)
    plt.tight_layout()
    plt.savefig(path_save + "boxplot{}_sumf_not_weighted.svg".format(y_))
    plt.savefig(path_save + "boxplot{}_sumf_not_weighted.pdf".format(y_))

def my_swarmplot(x_, y_):
    # Make swarm plot
    sns.set_style("darkgrid")
    sns.set_context( rc={"grid.linewidth": 2})
    fig=plt.figure(figsize=(len(conditions_to_plot)*1.2, 5))
    #if only 2 conditions *1.8, otherwise *1.2
    plt.axhline(y=np.mean(random_sumf), color='red', linestyle='dotted', alpha = 0.5, lw=2, zorder=1) 

    bp= sns.swarmplot(x=x_, y=y_, hue='seed', 
                      data=df_plot, 
                      palette="colorblind", 
                      dodge=False,size=6) 
    
    # get legend information from the plot object
    handles, labels = bp.get_legend_handles_labels()
    plt.ylim(0, 1)
    # specify just one legend
    plt.legend(handles[0:5], labels[0:5], title="Sp. set", fontsize=16, ncol=5, loc="lower left", framealpha=0.5)
    #plt.legend("")
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    # plt.title("Diff. in 50 rounds", fontsize = 28)
    plt.ylabel("Average total investment", fontsize = 18)
    plt.xlabel("Method", fontsize = 18)
    plt.tight_layout()
    plt.savefig(path_save + "swarmplot{}_sumf_not_weighted.svg".format(y_))
    plt.savefig(path_save + "swarmplot{}_sumf_not_weighted.pdf".format(y_))

my_swarmplot("prop", "mean_sumf")

#%%STATISTICAL TESTING
#check for normality
for prop in conditions_rename:
    x = df.loc[df["prop"] == prop,"mean_sumf"]
    print(shapiro(x))
    
#since most of them are not normally distributed I´ll use wilkinson test to compare

#COMPARE AGAINST NO-SELECTION CONTROL
what = "mean_sumf"
df_sw = pd.DataFrame( columns = ["prop1", "prop2","sw", "p_val"] )
for prop in conditions_rename:
    if prop == "NS":
        continue
    one = df.loc[df["prop"] == prop, what]
    two = df.loc[df["prop"] == "NS",what]
    d = np.subtract(one,two)
    sw, p = wilcoxon(d)
    df_now = pd.DataFrame([[prop, "n_r", sw, p]], columns = ["prop1", "prop2","sw", "p_val"] )
    df_sw = pd.concat([df_sw,df_now])
df_sw.to_csv(path_save + "dfw_vs_NS.csv")

#EVERYONE VS ROUND 0
df_sw = pd.DataFrame( columns = ["prop", "sw", "p_val"] )
#chose what to plot from the df
for prop in conditions_rename:
    one = df.loc[df["prop"] == prop, what]
    d = np.subtract(one,0.5)
    sw, p = wilcoxon(d,alternative='greater')
    df_now = pd.DataFrame([[prop, sw, p]], columns = ["prop", "sw", "p_val"] )
    df_sw = pd.concat([df_sw,df_now])
df_sw.to_csv(path_save + "dfw_vs_round0.csv")

#COMPARE EACH PROP SELECTION AGAINST CONTROL
props = ["D", "P", "M", "MI", "PI"]
df_sw = pd.DataFrame( columns = ["prop1", "prop2","sw", "p_val"] )
for prop in props:
    one = df.loc[df["prop"] == prop +"S", what]
    two = df.loc[df["prop"] == prop +"R",what]
    d = np.subtract(one,two)
    sw, p = wilcoxon(d)
    df_now = pd.DataFrame([[prop +"S", prop + "R", sw, p]], columns = ["prop1", "prop2","sw", "p_val"] )
    df_sw = pd.concat([df_sw,df_now])
df_sw.to_csv(path_save + "dfw_S_vs_R.csv")

#COMPARE PROPAGATIONS UNDER SELECTION
props = ["DS", "PS", "MS", "PIS", "MIS"]
combs = list(itertools.combinations(props, 2))
df_sw = pd.DataFrame( columns = ["prop1", "prop2","sw", "p_val"] )
for comb in combs:
    one = df.loc[df["prop"] == comb[0], what]
    two = df.loc[df["prop"] == comb[1],what]
    d = np.subtract(one,two)
    sw, p = wilcoxon(d)
    df_now = pd.DataFrame([[comb[0], comb[1], sw, p]], columns = ["prop1", "prop2","sw", "p_val"] )
    df_sw = pd.concat([df_sw,df_now])
df_sw.to_csv(path_save + "dfw_vs_S.csv")


#COMPARE ALL (JUST IN CASE)
props = conditions_rename
what = "mean_sumf"
combs = list(itertools.combinations(props, 2))
df_sw = pd.DataFrame( columns = ["prop1", "prop2","sw", "p_val"] )
for comb in combs:
    one = df.loc[df["prop"] == comb[0], what]
    two = df.loc[df["prop"] == comb[1], what]
    d = np.subtract(one,two)
    sw, p = wilcoxon(d)
    df_now = pd.DataFrame([[comb[0], comb[1], sw, p]], columns = ["prop1", "prop2","sw", "p_val"] )
    df_sw = df_sw.append(df_now)
df_sw.to_csv(path_save + "dfw_sumf_vs_all.csv")
