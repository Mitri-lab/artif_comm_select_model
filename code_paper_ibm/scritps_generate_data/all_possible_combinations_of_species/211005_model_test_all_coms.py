# -*- coding: utf-8 -*-
"""
Spyder Editor

Script to simulate all the possible combinations of species. We do it in triplicates, like, for each possible combination we seed 3 identical tubes 
and let them grow. This script is made in python, to be run in a python interpreter (we should run it 5 times, one for each of the 5 species sets that
correspond to seeds 22-26).

"""
#%% Section 1
#----------------------#
#    IMPORT PACKAGES   #
#----------------------#
import numpy as np
import pandas as pd
#import random
from collections import Counter
from collections.abc import Iterable
from math import sqrt,exp,log
import matplotlib.pyplot as plt
from datetime import datetime
import os
import pickle
import copy
import itertools
#import time as ttime



#%% Section 2
#----------------------#
#   PARAMETER CHOICE   #
#----------------------#

path="/home/pguridi/211005/" #where we will store things

#SET RANDOM NUMBER GENERATOR WITH RANDOM SEED
seed = 26
rng = np.random.default_rng(seed)


#GENERAL PARAMETERS:
total_sp = 15
N_comms = 3 #number of different communities in the simulation (would be like replicates we do in each case, like having 3 tubes each with the same starting conditions)
# In this version of the model each tube contains different species
round_lim = 1 #number of rounds
time_lim = 80 #time steps within each round

num_repeats = 1 #each repeat is a simualtion with the number of tubes (here we don´t need)

#PROPAGATION #we won´t use any of this here (so I delete the parameters related to that).

#INSIDE EACH TUBE:

in_cells0 = 10 #initial number of cells of each of the species                    
N_nutr = 4 #number of nutrients (edit in_nuts below as well if you modify this)
in_nuts = [2000 for i in range(N_nutr)] #amount of each of the nutrients
N_tox = 10 #number of toxins (edit in_toxs and below as well if you modify this)
in_toxs = [700 for i in range( N_tox)] #amounts of each of the toxins

#mutation #we won´t mutate anything, but just in case I keep it
par_mut = "fs" #Parameters to mutate. Mutation function designed to keep range [0,1], modify it 
              #to mutate parameters with a different range. In this model we can only mutate 1 
              #parameter at a time, modify mutation() and count() funtcions to include more.
mu_mut = 0 #Mutation probability
sd_mut = 0.4 #Standard deviation of the distribution to mutate. We multiply the value times a 
             #sample of a normal distribution centred around 1, with this sd. 
#death function
K_d = (sum(in_toxs))/10 #this is the K in the Hill function
n_d = 2 #Hill coeficient for the toxin
#division

scn = [] #I think it si more convinient as a list rather than as a dict as above:
for j in range (N_nutr):
    scn.append(1/in_nuts[j])

    


#And we save all these choices to a .txt file:
date = datetime.now()#YY_mm_dd_H_M_S
now = date.strftime("%Y_%m_%d_%H%M%S") + "_" + str(seed) #we add the date and time to each folder so they are unique
os.mkdir(path + now)
os.mkdir(path + now + "/plots")
f=open(path + now + '/INFO.txt','w')
f.write("Random seed: "+ str(seed) +'\n')
f.write("Number of species: "+ str(total_sp) +'\n')
f.write("Initial omout of each species per tube: "+ str(in_cells0) +'\n')
f.write("Initial nutrients: "+ str(in_nuts) +'\n')
f.write("Initial toxins: "+ str(in_toxs) +'\n')
f.write("Total tubes: "+ str(N_comms) +'\n')


f.write("Rounds: "+ str(round_lim) +'\n')
f.write("Time steps: "+ str(time_lim) +'\n')
f.write("Mutate: "+ str(mu_mut) + str(par_mut) + str(sd_mut) +'\n')
f.write("Death: "+ str(K_d) + "^" + str(n_d) + '\n')
f.write("Scaling for nutrients in division: " + str(scn) + '\n')
#for item in species: 
    #f.write("%s\n" % item)
f.close()




#%% Section 3
#----------------------#
#       FUNCTIONS      #
#----------------------#
def generate_species(n):
    """
    Function to randomly generate species (i.e. combination of parameter values).
    
    Parameters of the function
    --------------------------
    n : number of species to generate

    Returns
    -------
    A list with the randomly generated species, each of them as dictionaries, where keys are 
    parameter names and values the parameter values.

    """
    species_list = []
    for i in range(n):
        c={}
        c["t"] = i #species identity, it is called "s" in the methods of the paper
        c["a"] = rng.beta(2,2)    #activation probability [0,1]
        c["r"] = rng.beta(2,2)    #replication probability [0,1]
        c["v"] = 1 #for now I consider we don´t have it
        #now we create the n values for the corresponding number of nutrients we normalise so 
        #"n" values are stored in a list called "ns"
        ns_0 = rng.uniform (0,1,N_nutr)
        sparsity_n = rng.binomial(1,p=0.5,size=N_nutr)
        if sum(sparsity_n) == 0: #they should consume at least one of the nutrients
            correct = rng.choice(range(N_nutr))
            sparsity_n[correct] = 1
        ns_1 = np.multiply(ns_0,sparsity_n)
        c["ns"] = list(ns_1/sum(ns_1))   #so now they do not sum 1, but any number in [0,1]
        #investment for the degradation of each toxin: "f" values stored in "fs"
        fs_0 = rng.uniform (0,1, N_tox)
        sparsity_f = rng.binomial(1,p=0.5,size= N_tox) #they won´t degrade some toxins
        if sum(sparsity_f) == 0: #they should degrade at least one of the toxins
            correct = rng.choice(range( N_tox))
            sparsity_f[correct] = 1
        fs_1 = np.multiply(fs_0,sparsity_f)
        sum_f = rng.uniform(0,1) # the f values will sum to a number up to a number between 0 and 1.
        c["fs"] = [sum_f*x/sum(fs_1) for x in fs_1] #here the sum all all fs will be a number between
        # 0 and 1, after considering the sparsity. If we want we can always mutate this sum of fs. If we use u!=1, then
        #it should be considered in this scaling
        #Toxic effect: "m" values stored in "ms"
        ms_0 = rng.uniform (0.001,0.02, N_tox) 
        sparsity_m = rng.binomial(1,p=0.5,size= N_tox)
        if sum(sparsity_m) == 0: #at least one of the toxins should be toxic for the species
            correct = rng.choice(range( N_tox))
            sparsity_m[correct] = 1
        c["ms"] = list(np.multiply(ms_0,sparsity_m))
        c["p0"] = 0 #population of inactivated cells
        c["p1"] = 0 #population of activated cells.
        #For the methods of the per we consider S_i = p0 + p1
    
        species_list.append(c)
    return(species_list)



def species_combinations(num_species, sp_per_tube, num_tubes, num_repeats):
    """
    Function to generate species combinations for all the tubes and for all the repeats. 

    Parameters
    ----------
    num_species : number
        number of total species
    sp_per_tube : number
        number of species per community initially seeded
    num_tubes : number
        number of tubes or communities
    num_repeats : number
        number of repeats

    Returns
    -------
    List of lists: each repeat has a list and whithin this repeat lies there is a list corresponding to the  
    species combinatio at round 0.

    """
    list_sp = list(range(num_species))
    combinations = []
    for rp in range(num_repeats): 
        check_all_sp = [] 
        control = 0 #to avoid an endless loop if it is not possible or little likely to have all the species
        while len(set(check_all_sp))< num_species: #Check all the species are in the metacommunity at least once
            check_all_sp = [] 
            index_sp = []
            for t in range(num_tubes):
                tube = rng.choice(list_sp,sp_per_tube,replace = False) #all the species in a tube are different
                index_sp.append(tube.copy())
                check_all_sp.extend(tube.copy())
            control += 1
            if control>1000:
                break
        combinations.append(index_sp)
    return(combinations)


def mutation(item, sd):
    """
    Function to mutate (change the value) of one of the parameters in a cell, 
    This function is used inside the replication function (assumptions that mutations occur at the division step)


    Parameters
    ----------
    item : value of the particular parameter we want to mutate or list if there are more than 1
    sd : in mutation we multiply the old value times a sample from normal distribution centred 
    around 1, this is the standard deviation of that distribution. (number)

    Returns
    -------
    The mutated value or when the parameter is a list, a list including the mutated value(s) 

    """
    
    if isinstance(item, (int, float)): #if it is a number
        new_item = item
        new_item *= rng.lognormal(mean=0, sigma=sd, size=None) 
        #we keep all the parameters in the range [0,1]
        if new_item > 1: 
            new_item = 1
        elif new_item < 0:
            new_item = 0          
    else: #if it is a list, any element can be mutated
        new_item = item.copy()
        if sum(item) > 0: #check it is not empty
            #Mutate any of the items in the list with the same prob
            indexes = [i for i, x in enumerate(new_item) if x] #those different from 0
            mutate_index = rng.binomial(1,p=1/len(indexes),size= len(indexes)) #elements in the list that will mutate
            if sum(mutate_index) == 0: #at least one element in the list should mutate
                correct = rng.choice(range(len(mutate_index)))
                mutate_index[correct] = 1
            for index in [i for i,x in enumerate(mutate_index) if x]:
                new_item[indexes[index]] *=  rng.lognormal(mean=0, sigma=sd, size=None) 
            #the sum of all the elements what can´t be higher than 1
            if sum(new_item)>1:
                new_item = [x/sum(new_item) for x in new_item]
            #due to numerical problems with floats in python, to make sure it is not higher than 1, repeat the loop:
            if sum(new_item)>1:
                new_item = [x/(sum(new_item)+1e-10) for x in new_item]



    return (new_item)
            


def count(tube):
    """
    Function to count the population of species and strains. 
    
    Parameters
    ----------
    tube : any list of cells to count (usually will be grid[t], where t is the number of a particular tube)
 

    Returns
    -------
    It will return a tupple in which:
        - First element: dictionary (output of Counter function) where keys are species and values their population. 
        - Second element: dictionary where a key indicates the species, and the value is a dictionary inside
            which the keys are the mutable parmater (here fs) of a strain and the value the strain population.
        - Third element: the diversity at the species level (sp_div), calculated as the shanon index of the species H. 
        - Forth element: the diversity at the strain level (st_div) for all the species originally present in the 
            tube, so it is a list. #If we want just one value for the strain diversity we could do the average.

    """
    sp_t = {} #species in this tube
    st_t = {} #strains in this tube 
    
    for cell in tube:
        #add population to the counter of each species (create the counter if it does not exist)
        #population = activated (p1) + inactivated cells (p0); in the methods of the paper is called S_i
        sp_t[cell["t"]] = sp_t.get(cell["t"], 0) + cell["p0"] + cell["p1"]
        #st_t for each species, a counter of the strains
        if cell["t"] in st_t:
            st_t[cell["t"]][tuple(cell[par_mut])] = st_t[cell["t"]].get(tuple(cell[par_mut]), 0) + cell["p0"] + cell["p1"]
        else:
            st_t[cell["t"]] = {tuple(cell[par_mut]):  cell["p0"] + cell["p1"]}

    sp_count = dict(sorted(sp_t.items())) #having species sorted helps a bit later on
    st_count = dict(sorted(st_t.items()))
    
    #species diversity
    H_sp = -1*sum([v/sum(sp_count.values())*log(v/sum(sp_count.values())) for v in sp_count.values()])
    #strain diversity
    H_st = []
    for sp_counter in st_count.values():
        H_st.append(-1*sum([(v/sum(sp_counter.values()))*log(v/sum(sp_counter.values())) for v in sp_counter.values()]))
    return (sp_count, st_count, H_sp, H_st) #np.mean([exp(x) for x in H_st])
    

    
def simpler_count(tube):
    """
    Function to count the population of each species. We use this simpler count function to count
    the number of bacteria of each species at each time step. 

    Parameters
    ----------
    tube : any list of cells to count (usually will be grid[t], where t is the number of a particular tube)
 

    Returns
    -------
    It will return a dictionary (output of Counter function, but sorted by keyvalue, i.e. species)
    where keys are species and values their population. 

    """
    sp_t = {} #species in this tube
    
    #First we get the species and strains we have in this tube
    for cell in tube:
        #add to the counter of each species (create it if not present)
        sp_t[cell["t"]] = sp_t.get(cell["t"], 0) + cell["p0"] + cell["p1"]
        

    return (dict(sorted(sp_t.items()))) 


def count_population(tube):
    """
    Function to obtain the total population of a tube, adding all the species and strains.

    Parameters
    ----------
    tube : list of cells
    Returns
    -------
    It returns the number of total cells in the tube, without caring about species,activated state or anything

    """
    list_count=[]
    for cell in tube:
        list_count.append(cell["p0"])
        list_count.append(cell["p1"])
    return(sum(list_count))



def consumption(final, initial):
    """
    Function to assess degradation in a particular tube at a particular moment.
    It can also work for the nutrient consumption 

    Parameters
    ----------
    final : list with the final amount of nutrients or toxins
    initial: list with the initial amount of nutrients or toxins

    Returns
    -------
    The degradation score or nutrient uptake score in that tube, where 0 means 
    no consumption or nutrients or toxins, and 1 means consumption of all the nutreints
    or toxins
    """
    
    #just the definition in the paper methods of the above function, the result should be the same
    score = 1 - sqrt(1/len(final) * sum((a/b)**2 for a,b in zip(final,initial)))
    
    return(score)


def consumption_old(final, initial):
    """
    Originally using this funciton, it is just the same thing as the above (well for numerical issues with floats in 
                                                                            python maybe like the 15th decimal number 
                                                                            changes, but it is the same)
    """
    
    score = 1 - sqrt(sum([a**2 for a in final]))/sqrt(sum([b**2 for b in initial]))

    return(score)



def dilution(community, ratio):#implement a function to count the cells and select only 1%
    """
    Function to make a random dilution from a tube in the propagation step (used in propagule, migrant pool and no_sel)
    CAREFUL-> It is a "real dilution": it also deletes de selected cells from the original tube. This could lead to a 
    cell with p0=0, p1=0 in the original tube. Those would be eliminated later on in the death step, although anyway 
    the original tube is not used in the way teh model is defined.
    It also deactivates activated cells in the dilution process.
    
    Parameters
    ----------
    community: list of cells that we want to dilute
    ratio:  the fraction (0,1] of cells from the community that we want to take at random

    Returns
    -------
    a list of selected cells. Deletes from community the selected cells
    """
    
    selected_cells = []
    if count_population(community) > 0: #if not empty 
        control = 0
        while len(selected_cells) == 0: #we will return at least one cell if the community is not empty
            for i,cell in enumerate(community):
                pop = cell["p0"] + cell["p1"]
                select = rng.poisson(pop*ratio)
                if select:
                    #we cannot select more cells than the ones we have
                    if select > pop:
                        select = pop
                    cell["p1"] = 0 #deactivate also the cells in the community
                    cell["p0"] = pop - select #we remove the selected cells so we cannot select later
                    cell_here = copy.deepcopy(cell)
                    cell_here["p0"]= select
                    selected_cells.append(cell_here)
            control += 1
            if control > 1000:
                break

    return (selected_cells)



def sp_finder(tube, sp):
    """
    Funciotn to return all the cells of a particular species in a tube

    Parameters
    ----------
    tube : list of cells to find a species

    sp : which species (number)

    Returns
    -------
    A list with the cells of that species in that tube


    """

    cells = copy.deepcopy([cell for cell in tube if cell["t"] == sp])
    
    return (cells)
    
    

def sample_or_amplify(community, size):
    """
    Function to be used in disassembly propagation.   
    We assume we can amplify cells for this reason, they are not deleted from the original tube when chosen.
    
    Parameters
    ----------
    community: list of cells that we want to dilute
    size: number of cells to get approximatelly in the smapling/amplification

    Returns
    -------
    From the list you provide(community) it will return another list of size approximatelly "size". Since the process 
    is based in poisson sampling, it can actually change a bit. If community is bigger than size it will do it by 
    random drawing. It community is smaller than size, it will randomly amplify the items in community.

    """
    new_list = []
    tot_pop = count_population(community)
    if tot_pop>0: 
        sample_here = copy.deepcopy(community) #copy the list of cells to not modify the original
        control = 0
        while len(new_list) == 0: #It sheldomly happens, but if the list is empty the process is repeated
            for cell in sample_here:
                pop = cell["p0"] + cell["p1"]
                num_cells = rng.poisson(size*pop/tot_pop)
                if num_cells:
                    cell["p0"] = num_cells
                    cell["p1"] = 0
                    new_list.append(cell) #we have copied above so it should work
            control += 1
            if control > 1000:
                break
    
    return(new_list)



def get_cmap(n, name='hsv'):
    """Returns a function that maps each index in 0, 1, ..., n-1 to a distinct RGB color; the keyword argument name 
    must be a standard mpl colormap name.
    This generates enough colors depending on the total number of species that will be used to distinguish species when 
    ploting within a round.
    #https://stackoverflow.com/questions/14720331/how-to-generate-random-colors-in-matplotlib
    """
    return plt.cm.get_cmap(name, n)

#I already create the colormap I need
multiple_colors = get_cmap(total_sp + int(0.2*total_sp)) #Generate few more colors bcause the extremes are quite similar


#%% Section 4
#----------------------#
#  WORKING STRUCTURES  #
#----------------------#

#LOAD OR GENERATE SPECIES

#generate

species = generate_species(total_sp)
#Save the species as a .txt
f=open(path + now + '/species.txt','w')
f.write("The species in this simulation are:")
for item in species: 
    f.write("%s\n" % item)
f.close()
#Save the species in pickle format (to reload as dict and use if needed)
open_file = open(path + now + "/Species" , "wb")
pickle.dump(species, open_file)
open_file.close()


#Whenever we need to load them:
"""
open_file = open("/Users/pablo/Desktop/data_mitri/210530/species", "rb")
species = pickle.load(open_file)
open_file.close()
"""



#%% SECTION 5:
#-------------------------------#
#  LET´S CREATE THE STRUCTURE   #
#-------------------------------#


columns_here= ["com",  "final_pop", "max_pop", "fold_change", "auc", "deg_score", "nut_score", "sp_div_H"]

 
df = pd.DataFrame( columns = columns_here)


#Generate all the posible combinations with 15 species in communities of 4



comms = list(itertools.combinations(list(range(15)),15)) + \
    list(itertools.combinations(list(range(15)),14)) + \
        list(itertools.combinations(list(range(15)),13)) +\
            list(itertools.combinations(list(range(15)),12)) + \
                list(itertools.combinations(list(range(15)),11)) +\
                    list(itertools.combinations(list(range(15)),10)) +\
                        list(itertools.combinations(list(range(15)),9)) +\
                            list(itertools.combinations(list(range(15)),8)) +\
                                list(itertools.combinations(list(range(15)),7)) +\
                                    list(itertools.combinations(list(range(15)),6)) +\
                                        list(itertools.combinations(list(range(15)),5)) +\
                                            list(itertools.combinations(list(range(15)),4)) +\
                                                list(itertools.combinations(list(range(15)),3)) +\
                                                    list(itertools.combinations(list(range(15)),2)) +\
                                                        list(itertools.combinations(list(range(15)),1))

#quartets = pd.read_csv("/Users/pablo/Desktop/210526/Species_quartets_20210526-120138.csv", header=0)

counteeer = 0 #to save the progress not all the time but after a number of simulations
for tup in comms:
    
    #CREATE DATA STRUCTURE (GRID OF TUBES)
    grid = [[] for i in range(N_comms)] #Grid to store cells 
    nuts = [] #Grid to store nutrients
    toxs = [] #Grid to store toxins
    #In grid we will have the tubes, in each of them I store all the cells
    #Then we have a similar grid for nutrients(nuts) and another one for toxins (toxs)
    #Each tube in nuts or toxs corresponds to the tube in the same position in the grid
    
    for t in range(N_comms):
        nuts.append([*in_nuts]) 
        toxs.append([*in_toxs])

        for i in range(len(tup)):
            cell_here = copy.deepcopy(species[tup[i]])
            cell_here["p0"]=in_cells0
            grid[t].append(cell_here)
            
        rng.shuffle(grid[t]) 
        # This would go here, not anymore: grid[t].insert(0,[t, *in_nuts, *in_toxs])
    
    
    #%% Section 6
    #----------------------#
    #  MODEL WITHIN ROUNDS #
    #----------------------#
    
    #LET´S GO
    
    #here I am not using rounds, but just for the plot to indicate it
    rd = 0
    
    #STRUCTURES FOR DATA STORAGE FOR EACH TIME STEP WITHIN A ROUND
    #Each of them is a dictionary, each tube as a key, and inside a list of lists, one list for
    #each bacterial species, each nutrient or each toxins respectively. In each of the list we 
    #append one value at each time step. The dicts will be filled later on:
    dict_grid = {} 
    dict_nuts = {} 
    dict_toxs = {} 
    #Next we store the population at the begining and the species we have in this particular tube:
    pop_0 = []
    sp_tube = [] #I need to save the species in each tube for the count later on
    #and finally some results at the end of the round
    deg_list = [] #list to store what each of the tubes or communities has degraded
    nut_list = [] #list to strore the nutrient uptake (same formula as degradation)
    div_list = [] #list to store the species diversity of each tube at the end of a round
    pop_list = [] #list to store the final population 
    auc_list = []
    max_pop_list = []
    fc_list = []

   
    for t in range(N_comms):
        pop_0.append (simpler_count(grid[t]))
        sp_tube.append(list(pop_0[t].keys()))
        dict_grid["t{0}".format(t)] = [[y] for y in pop_0[t].values()] 
        #we replenish the initial nutrients and toxins at every round.
        dict_nuts["t{0}".format(t)] = [[y] for y in in_nuts] 
        dict_toxs["t{0}".format(t)] = [[y] for y in in_toxs]
    
    
    tm = 0
    while tm < time_lim: # For each time step
                        
        #DEGRADATION-ACTIVATION
        #---------------------#
        for t,tube in enumerate(grid):
            for cell in tube:
                
                
                #Rename some variables for simplicity, not really needed but makes it look a bit nicer
                ns = cell["ns"] #all the nutrient consumption values
                fs = cell["fs"] #Same as above for toxind deg. 
                v = cell["v"]
                p0 = cell["p0"] #population of non-activated cells
                p1 = cell["p1"] #population of activated cells
                pt = p0+p1 #total population (named as S_i in the code of the paper)
                
                #Check if all required nutrients are present
                enough = [x>=n for x,n in zip(nuts[t],ns)] # array of booleans, equal to 0,1 when multiplying
                if sum(enough) >= 1: #equivalent to " if any(enough):" but this seems to be slightly more efficient
                    #if they have nutrients they can degrade and/or activate then scale when they don´t have all 
                    #the nutrients, so they degrade less, and activate with less probability, although they will 
                    #still take 1 unit * investment in total of the remaining nutrients to activate and degrade.
                    scaling = sum([a*b for a,b in zip(ns,enough)])  #If all the required nutrients are present, 
                                                                    #scaling = 1, so it has no effect.
                    #new_ns scales the ns to still consume 1 unit in total if we are missing some of the nutrients 
                    #new_ns will be equal to ns if we have all the nutrients. Note that when there is few nutirents,
                    #this could make that those nutrients that were enough above are not anymore, but this is not
                    #as the population that can afford to consume nutrients is calculated below.
                    #I add the +1e-10 to not divide by 0 now that we have a sparse model
                    new_ns = [a*b for a,b in zip(ns,enough)]; new_ns = [x /(sum(new_ns)+1e-10) for x in new_ns]
                
                    #check the populaiton size we can afford to activate and degrade based on the nutrients
                    afford = [(nuts[t][j])/(new_ns[j]) for j,x in enumerate(new_ns) if x]
                    if afford:
                        
                        max_afford = int(min(afford)) #maximal number (integer) of cells that we can afford, based on the 
                                                      #limiting nutrient (then they can change metabolism and 
                                                      #actually continue dividing and degrading)
                    else:
                        max_afford == 0

                    #DEGRADATION
                    if max_afford > 0:
                        if max_afford < pt:
                            pt = max_afford
                        for k in range(N_tox):
                            if fs[k] != 0: #not really needed, but reduces time cost with sparse degradation
                                if toxs[t][k] >= (fs[k]**v)*scaling*pt: 
                                    #We have more toxin than what the strain degrades per time step (so it can degrade 
                                    #and consume the nutrients it requires)
                                    toxs[t][k] -= (fs[k]**v)*scaling*pt #Degrades de corresponding toxin
                                    for j in range(N_nutr):
                                        nuts[t][j] -= pt*new_ns[j]*fs[k]**v #Consumes nutrients proportionaly to what 
                                                                            #is statedin the genome, the nutrients we 
                                                                            #have, and the toxin we have degraded.                                   
                                
                                elif toxs[t][k] > 0: #if we dont meet the above but still have some toxins
                                    #if desired, matematically the investment (fs[k]**u)**v) can be removed in this chunk
                                    pt_here = toxs[t][k]/((fs[k]**v)*scaling)
                                    toxs[t][k] = 0
                                    for j in range(N_nutr):
                                        nuts[t][j] -= pt_here*new_ns[j]*fs[k]**v
                       
                    #ACTIVATION   
                    if max_afford > 0:
                        
                        #costly to remain activated
                        i_div = (1-sum(fs))**v #investment in division; rename it for simplicity
                        if max_afford >= p1: 
                            max_afford_here = max_afford - p1
                            for j in range(len(ns)):
                                nuts[t][j] -= p1*new_ns[j]*i_div
                        else:
                            cell["p0"] += (p1 - max_afford) #cells that can´t afford to stay activated deactivate
                            cell["p1"] = max_afford
                            max_afford_here = 0
                            for j in range(len(ns)):
                                nuts[t][j] -= max_afford*new_ns[j]*i_div
                        
                        #now let´s check to activate
                        requirement = sum([a*b*c for a,b,c in zip(new_ns,scn,nuts[t])]) #they activate proportionally 
                                                                                        #to the nutrients they have
                        #number of cells that will activate:
                        cells_activate = rng.poisson(cell["a"]*scaling*requirement*p0*i_div)
                        if cells_activate > p0:
                                cells_activate = p0
                        if cells_activate > max_afford_here:
                            cells_activate = max_afford_here  
                        cell["p0"] -= cells_activate
                        cell["p1"] += cells_activate
                        for j in range(len(ns)):
                            nuts[t][j] -= cells_activate*new_ns[j]*i_div
                            #Consume of each nutrient what is stated times investment in div division
                            
                    elif max_afford == 0: #since it is costly to stay activated, deactivate if there are no nuts
                        cell["p0"] += cell["p1"]
                        cell["p1"] = 0
                
                elif enough == 0: 
                    cell["p0"] += cell["p1"]
                    cell["p1"] = 0
                    
                    
        #REPLICATION
        #----------#
        for t,tube in enumerate(grid): #repeating the for loop seems to help a bit in efficiency
            i = 0
            original_len = len(tube)
            while i < original_len: #while loop to not loop through newly generated cells
                cell = tube[i] #local variables are more efficient
                i_div = (1-sum(cell["fs"]))**cell["v"] #investemnt, rename for simplicity
                num_new_cells = rng.poisson(cell["r"]*i_div*cell["p1"])
                if num_new_cells > cell["p1"]:
                    num_new_cells = cell["p1"]
                num_mutants = rng.poisson(num_new_cells*mu_mut)
                if num_mutants > num_new_cells:
                    num_mutants = num_new_cells
                cell["p1"] -= num_new_cells
                cell["p0"] += (num_new_cells*2 - num_mutants)
                for mutant in range(num_mutants):
                    new_cell = copy.deepcopy(cell)
                    new_cell["p0"] = 1
                    new_cell["p1"] = 0
                    new_cell[par_mut] =  mutation(new_cell[par_mut], sd_mut).copy()
                    tube.append(new_cell) 
            
                i += 1
        
        #DEATH
        for t in range(N_comms): 
            for c in reversed(range(len(grid[t]))): #reversed to not interfere the looping when deleting a cell
                ms =  grid[t][c]["ms"] #values for the toxicity of each nutrient
                p0 = grid[t][c]["p0"]
                p1 = grid[t][c]["p1"]
                #Now we multiply the hill function of each toxin times the toxicity of that toxin (m0,m1...)
                #both activated or unactivated cells can die
                num_death_p0 = rng.poisson (p0*sum([m * (tox**n_d)/(tox**n_d + K_d**n_d) for m,tox in zip(ms, toxs[t])]))
                if num_death_p0 >p0:
                   num_death_p0 = p0
                num_death_p1 = rng.poisson (p1*sum([m * (tox**n_d)/(tox**n_d + K_d**n_d) for m,tox in zip(ms, toxs[t])]))
                if num_death_p1 >p1:
                   num_death_p1 = p1
                grid[t][c]["p0"] -= num_death_p0
                grid[t][c]["p1"] -= num_death_p1
                
                if grid[t][c]["p0"] + grid[t][c]["p1"] == 0: #when there are no cells of this strain
                    grid[t].pop(c) 
                    #grid[t].remove(grid[t][c]) #; this was before, .pop() is slightly better
            
            #After division and death we randomize the order of cells to not give preferences in the next round
            rng.shuffle(grid[t])
    
        #STORAGE AT EACH TM WITHIN THE ROUND
        #We store information at every tm in the round
        for t in range( N_comms):
            pop_tm = simpler_count(grid[t])
            for i in range(len(sp_tube[t])):
                #we have to check that the species did not go extint
                try:
                    dict_grid["t{0}".format(t)][i].append(pop_tm[sp_tube[t][i]])
                except:
                    dict_grid["t{0}".format(t)][i].append(0)
            for i in range(len(in_nuts)):
                dict_nuts["t{0}".format(t)][i].append(nuts[t][i])
            for i in range(len(in_toxs)):
                dict_toxs["t{0}".format(t)][i].append(toxs[t][i]) 
        tm += 1
        
    #STORE SOME INFO AT THE END OF THE ROUND
    for t in range(N_comms):
        #Degradation per tube
        deg_list.append(consumption(toxs[t], in_toxs))
        #Nutrient uptake per tube:
        nut_list.append(consumption(nuts[t], in_nuts))
        #get the population count at the last time step of the round
        pop = count(grid[t]) #See the count function to understand the structure
        #Calculate population fold change for each bacterium:
        fc = []
        for k,v in pop_0[t].items():
            try:
                final = pop[0][k]
            except:
                final = 0
            fc.append(final/(v+0.00001)) #log with base 10

        fc_list.append(np.mean(fc))


        #And as we have stored info per time step we can calulate the AUC
        auc = []
        maxim = []
        for x in dict_grid["t{0}".format(t)]:
            auc.append(sum(x))
            maxim.append(max(x))
        #store general info
        pop_list.append(sum(pop[0].values()))
        div_list.append(pop[2])
        auc_list.append(sum(auc))
        max_pop_list.append(np.mean(maxim))
                
                
    
    #Now we plot what we have stored
    
    fig=plt.figure(figsize=(15, 4))
    plt.suptitle("Round = " + str(rd), fontsize=16)
    columns = 3
    rows = 1
    for t in range(N_comms):
        plt.subplot(rows, columns, t+1)
        top = 0 #just to get the highest population at some point in this tube (aesthetics)
        for i,j in enumerate(dict_grid["t{0}".format(t)]):
            top_here = max(j)
            if top_here > top:
                top = top_here
            plt.plot(range(time_lim+1), j, linewidth=2)
        plt.ylim (0, top* 1.3)
        #Contiue here
        plt.text(1,top*1.24, "Sp: " + str(sp_tube[t]))
        plt.text(1,top*1.08, "D:" + str(round(deg_list[t],3)))
        plt.text(1,top*1.16, "N:" + str([round(x) for x in nuts[t]]))
        
        plt.title("tube"+str(t))
        plt.xlabel("time steps", fontsize=12)
        plt.ylabel("population", fontsize=12)
        #for ax in fig.get_axes():
            #ax.label_outer()
    #Save the plot
    plt.savefig( path + now + "/plots"+ "/com "+ str(tup) + ".pdf")
    plt.close(fig)
    #Save a list with dict_grid, dict_nuts and dict_toxs in this round (pickle format)
    #open_file = open(path + now + "/param_" + str(parameter)  + "/" + str(number), "wb")
    #pickle.dump([dict_grid, dict_nuts, dict_toxs], open_file)
    #open_file.close()
    #Save the species_tube in pickle format (to reload as dict and use if needed)
    #open_file = open(path + now + "/param_"+ str(parameter)  + "/" + str(number), "wb")
    #pickle.dump(sp_tube, open_file)
    #open_file.close()
    
        

    
    df_now = pd.DataFrame([[ tup,
                             np.mean(pop_list),
                             np.mean(max_pop_list),
                             np.mean(fc_list),
                             np.mean(auc_list),
                             np.mean(deg_list),
                             np.mean(nut_list),
                             np.mean(div_list)]], 
                            columns = columns_here)
                             
    
    df = df.append(df_now)
    
    
    if counteeer%500 == 0:
        df.to_csv(path + now + "/df.csv" ,index=False)
    
    counteeer += 1
          
df.to_csv(path + now + "/df.csv" ,index=False)



        
