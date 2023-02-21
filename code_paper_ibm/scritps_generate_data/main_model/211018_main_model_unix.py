# -*- coding: utf-8 -*-
"""
Spyder Editor

On top of the changes in 210817_model_strains_poisson_ok, here:
    Remove "u" parameter which was used as an exponent for each of the f values.
    
At each repeat I reinitialise the repeat, but I rather do it with the seed we feed at the beginning,
shouldn´t matter, but just in case somebody could argue a weird linking between repeats of different seeds.

Save plot as pdf instead of png

210829
solve an issue in the count, just in case we divide by 0, which should not happen, but added the correction just in case

211018
Correct a bug that Björn found, in the propagation when adding invasion in "pip" and "mip" I was wrongly updating 
sp_to_include, in the second sp_to_include, I had before sp_in_play. The right one is:
    sp_to_include = list(set(sp_to_include)-{invader}) #The added species will not be added again

I also removed all the propagations that we are not using

Here I run all the propagation methods I want in the same script in order to need less cores and be able to use Marc´s 
computer in the lab.

Here I just want to try the propagation d3_s with the actual community.

220601
#generate species combinations only once, like the species, so they are the same for all propagations.
"""


#%% Section 1
#----------------------#
#    IMPORT PACKAGES   #
#----------------------#
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import pickle
import copy
import time as timer
from datetime import datetime
from math import sqrt,exp,log
import sys
#from collections import Counter
#import random



#%% Section 2
#----------------------#
#   PARAMETER CHOICE   #
#----------------------#

#SET RANDOM NUMBER GENERATOR WITH ANY SEED FOR REPRODUCIBLE RESULTS
#seed = int(sys.argv[1]) # when running the script in linux and feeding the seed as external variable
seed = int(sys.argv[1]) # when running the script in linux and feeding the seed as external variable
rng = np.random.default_rng(seed)
#it should be ok to generate the seed here, I generate species only for the first propagation methods,
#and then the seed is restarted for each repeat, so it should be reproducible.

#Another rng will be initialised within each repeat. It allows to have identical results when runing e.g. repeats 0:4 
#and then 4:8, as when running direclty 0:8.

#CREATE THE DIRECTORY TO STORE THE RESULTS
enter_path = "/home/pguridi/220601/" #specify path to store results of the simulation
path = enter_path + str(seed) + "/" #create a folder named as the seed (to seed we could append date+time for instance)
if os.path.exists(path):
    pass
else:
    try: #still had to add this because of some trouble when looping on bash to run this script several times
        os.mkdir(path)
    except:
        pass
#GENERAL PARAMETERS:
total_sp = 15 #total number of species in the simulation (i.e. per species set)
N_comms = 21 #number of different communities or tubes in the simulation 
beta = 1/3 #fraction of communities to be transferred
N_beta = int(N_comms*beta) #number of tubes to be transfered or propagated at the end of a round
round_lim = 50 #number of rounds
time_lim = 80 #time steps within each round (between t_0 and t_end)
list_repeats = range(0,10) #list with the repeats to run; 
                           #a repeat is a replicate of the process with different starting communities
N_repeats = 10 #regardless of which repeats we run, we generate the indexes of all of them for reproducibility

#INSIDE EACH TUBE:
N_spc0 = 4 #number of different species initially alocated per communitie (i.e. per tube)
in_cells = [10 for x in range(N_spc0)] #initial number of cells (population) of each of the species
in_cells0 = in_cells[0] #generic number of initial cells for a strain, used later for the st_repo
N_nutr = 4 #number of nutrients 
in_nuts = [2000 for i in range(N_nutr)] #amount of each of the nutrients supplied at the beginning of a round
N_tox = 10 #number of toxins 
in_toxs = [700 for i in range( N_tox)] #amounts of each of the toxins supplied at the beginning of a round

#PROPAGATION
#choose the propagations to run (one after another)
list_prop = ["d3", "d3", "mip", "mip", "pip"]
#(disassembly = "d3", propagule = "p", migrant pool = "m", propagule with invasion = "pip", migrant pool with invasion = "mip", no selection = "n")
#choose the treatments (selection ("s") or random("r"), for each of the propagations above, for "no selection" it actually does not matter what we chose here but something must be chosen)
list_treat = ["s", "r", "s", "r", "r"]
#we loop over the list of propagations selected above
for index_prop in range(len(list_prop)):
    #The current propagation method and treatment
    treat = list_treat[index_prop]
    prop = list_prop[index_prop]
    
    #dilution ratios, they must be higher than 0
    d_p = 0.05 #dilution ratio in propagule
    d_m = 0.05 #dilution ratio in migrant pool
    d_n = 0.05 #dilution ratio in no selection
    f_migr = 0.25 #fraction of tubes where migration occurs in disasembly (later implemented in other propagations)
    epsilon = N_comms/N_beta #threshold to consider a community extinct: when the population is lower than this threshold.
    s_d = 1 #sensitivity disassembly, minimun number of cells of one species detected in the disassembly process
    b_d_1 = 10 #bottleneck disassembly 1
    b_d_2 = b_d_1 #bottleneck disassembly 2 (not used, it used to simulate a second amplification due to colony growht)
    
    #DIVISION
    scn = [] #scaling for nutrient amounts playign a role in activation (first step of division)
    for j in range (N_nutr):
        scn.append(1/in_nuts[j])
    #mutation, which occurs in replication (second step of division)
    par_mut = "fs" #Parameter(s) to mutate (if more than one parameter is going to mutate use a list, and check mutation function)
    mu_mut = 0.01 #Mutation probability
    sd_mut = 0.4 #Standard deviation of the distribution (lognormal) to mutate the original value by multiplying times
                 #a random sample from this distribution. 
    
    #DEATH FUNCTION
    K_d = (sum(in_toxs))/10 #this is the K in the Hill function
    n_d = 2 #Hill coeficient for the toxin
    
    #We save plots within a round, the grid at time 0, and some lists with information with some periodicity. It can be 
    #done at every round, but to save a bit of time and storage, the periodicity to store can be chosen. By default all the
    #rounds 0-5 and the last round will be saved (to change that go further down where "weplothere" is defined). 
    #In between, chose the periodicity in number of rounds:
    save_periodicity = 10
    
    
    #SAVE ALL THESE CHOICES TO A .txt FILE:
    date = datetime.now()#YY_mm_dd_H_M_S
    #date.strftime("%Y_%m_%d_%H%M%S") 
    now = str(prop) +"_" + str(treat) # now is indeed the path to the directory, originally called now because it included
                                      # date and time to always be unique
    
    
    

    
    #Create a folder for the propagation method, inside the seed folder
    if os.path.exists(path + now):
        pass
    else:
        try: #still had to add this because of some trouble when looping on bash to run this script several times
            os.mkdir(path + now)
        except:
            pass
    f=open(path + now + '/INFO.txt','w')
    f.write("Date: " + date.strftime("%Y_%m_%d_%H%M%S") + '\n')
    f.write("Random seed: "+ str(seed) +'\n')
    f.write("Number of species: "+ str(total_sp) +'\n')
    f.write("Initial omout of each species per tube: "+ str(in_cells) +'\n')
    f.write("Initial nutrients: "+ str(in_nuts) +'\n')
    f.write("Initial toxins: "+ str(in_toxs) +'\n')
    f.write("Total tubes: "+ str(N_comms) +'\n')
    f.write("Tubes to transfer: "+ str(N_beta) +'\n')
    f.write("Propagation: "+ str(prop) + "_" + str(treat) + '\n')
    if prop == "d":
        f.write("bd1: "+ str(b_d_1) +'\n')
        f.write("bd2: "+ str(b_d_2) +'\n')
        f.write("p_inv_d: "+ str(f_migr) +'\n')
    if prop == "m":
        f.write("d_m: "+ str(d_m) +'\n')
    if prop == "p":
        f.write("d_p: "+ str(d_p) +'\n')
    f.write("Rounds: "+ str(round_lim) +'\n')
    f.write("Time steps: "+ str(time_lim) +'\n')
    f.write("Mutate: "+ str(mu_mut) + str(par_mut) + str(sd_mut) +'\n')
    f.write("Death: "+ str(K_d) + "^" + str(n_d) + '\n')
    f.write("Scaling for nutrients in division: " + str(scn) + '\n')
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
        try:
            H_sp = -1*sum([(v/(sum(sp_count.values())))*log(v/sum(sp_count.values())) for v in sp_count.values()])
        except:
            H_sp=0
        #strain diversity
        H_st = []
        for sp_counter in st_count.values():
            try:
                H_st.append(-1*sum([(v/sum(sp_counter.values()))*log(v/sum(sp_counter.values())) for v in sp_counter.values()]))
            except:
                H_st.append(0) 
        
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
    #------------------------------#
    #  SPCIES AND ITS COMBINATION  #
    #------------------------------#
    
    #LOAD OR GENERATE SPECIES
    
    #generate species
    if treat == list_treat[0] and prop == list_prop[0]: #generate species only once so they are shared for all the propagations
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
    #load species stored as pickle format:
    """
    open_file = open(path + "210421_species", "rb")
    species = pickle.load(open_file)
    open_file.close()
    
    #and we also save them
    open_file = open(path + now + "/species" , "wb")
    pickle.dump(species, open_file)
    open_file.close()
    """
    
    #SPECIES COMBINATION
    if treat == list_treat[0] and prop == list_prop[0]: #generate species combinations only once, like the species, so they are the same for all propagations.
        INDEX_REP = species_combinations(total_sp, N_spc0, N_comms, N_repeats)
    #and save them
    open_file = open(path + now +"/index_repeats", "wb")
    pickle.dump(INDEX_REP, open_file)
    open_file.close()
    
    #or load the order in which to combine the species in each repeat:
    """
    open_file = open( path + "210421_index_repeats", "rb")
    INDEX_REP = pickle.load(open_file)
    open_file.close()
    """
    
    
    #%% SECTION 5:
    #-------------------------#
    #  LOOPING OVER REPEATS   #
    #-------------------------#
    
    t1 = timer.time() #Needed to print the time it took at the end (cumulative for each repeat)
    #There we go
    for rp in list_repeats: 
        #rp = int(sys.argv[3]) #Not used anymore, used once from linux (if so, list_repeats would be just 1 item).
        rng = np.random.default_rng(seed) #restart the random number generator in each repeat 
                                        #(anyway with different species sets will lead to different outcomes)
        os.mkdir(path + now + "/repeat" + str(rp)) #make a folder inside the seed folder with the name of the repeat
        os.mkdir(path + now + "/repeat" + str(rp) + "/per_round") #folder inside repeat to store results at some round
        #Document to keep track of the extinctions within a repeat: when a selected tube does not have enough population
        f=open(path + now + "/repeat" + str(rp)+ '/extinctions.txt','w')
        f.close()
        
        #CREATE DATA STRUCTURE (GRID OF TUBES)
        nuts = [] #Grid to store nutrients
        toxs = [] #Grid to store toxins
        grid = [[] for i in range(N_comms)] #list to store the tubes or communities. Whithin a tube, each item(dictinary) 
                                            #corresponds to a population or strain characterised by all the parameters 
                                            #(theta) assigned to a species (see generate_species fuction)
        
        #Each tube in nuts or toxs corresponds to the tube in the same position in the grid
        
        #POPULATE THE DATA STRUCTURES
        sp_in_play = [] #stores how many times each species appears in the metacommunity, used in disassembly propagation 
        for t in range(N_comms):
            nuts.append([*in_nuts]) 
            toxs.append([*in_toxs])
            for i, sp in enumerate(INDEX_REP[rp][t]): #Here we add the desired species
                cell = copy.deepcopy(species[sp]) 
                sp_in_play.append(cell["t"])
                cell["p0"] += in_cells[i]
                grid[t].append(cell)
            rng.shuffle(grid[t]) # We do a random shuffling of each tube with cells
    
        
        st_repo = {}
        for sp in set(sp_in_play):
            cell_here = copy.deepcopy(species[sp])
            cell_here["p0"] = in_cells0
            st_repo[sp] = [cell_here]
        
    
        #---------------------------------#
        #     DATA STORAGE STRUCTURES     #
        #---------------------------------#
        #TO STORE AT THE END OF EACH ROUND (to see data storage within a round see model section)
        #Create the following dataframes
        df_nuts = pd.DataFrame( columns = ["round", "tube", "nut", "amount"])
        df_toxs = pd.DataFrame( columns = ["round", "tube", "tox", "amount"])
        df_grid = pd.DataFrame( columns = ["round", "tube", "sp_div_H","tot_pop_0","tot_final_pop", "tot_auc", "mean_fc", \
                                           "deg_score", "nut_score"])
        df_sp = pd.DataFrame( columns = ["round", "tube", "sp", "pop_0", "final_pop", "st_div_H", "auc", "fc"])
        df_st = pd.DataFrame( columns = ["round", "tube", "sp", "st", "final_pop"])
        #to track when mutants appear from which ancestor they are comming
        df_mutants = pd.DataFrame( columns = ["round", "tube", "sp", "ancestor", "mutant"]) 
        #Store the selected tubes in each round
        col_df_tubes = ["round"]
        for t in range(N_beta):
            col_df_tubes.append("tube" + str(t))
        df_tubes = pd.DataFrame (columns = col_df_tubes)
        #Store the ancestry: at each round, which tube is comming from which other in a preivous round, at round 0 all of them are original
        col_df_ancestry = ["round"]
        for t in range(N_comms):
            col_df_ancestry.append("tube" + str(t))
        df_ancestry = pd.DataFrame (columns = col_df_ancestry)
        df_ancestry.loc[0] = [0]+[np.nan for x in range(N_comms)]
        #Later on, at each round the content to fill in all these data frame will be added.
    
        
        #----------------------#
        #  MODEL WITHIN ROUNDS #
        #----------------------#
    
        rd = 0  #round
        while rd < round_lim: #For each round
            
            #STRUCTURES FOR DATA STORAGE FOR EACH TIME STEP WITHIN A ROUND
            #Dictionaries, each tube as a key, and inside a list of lists, with one list for each bacterial species (same 
            #order as in "sp_tube"), each nutrient or each toxins respectively. In each of the list we append one value at 
            #each time step. 
            dict_nuts = {} 
            dict_toxs = {} 
            dict_grid = {} 
            #Few lists, to store one value per tube:
            pop_0 = [] #list of dicts. They store population of each tube at the begining of the round
            sp_tube = [] #list of lists. They store species in each tube at the begining of the round
            deg_list = [] #store degradation score of each community at the of of the round
            nut_list = [] #store nutrient consumption score of each community at the end of the round
            #List to store population over time steps:
            auc_tot = [] #population of each tube at each time step (one list for each tube, each entry is one time step)
            
            #Let´s fill in data at round 0
            for t in range(N_comms): 
                pop_0.append (simpler_count(grid[t]))
                sp_tube.append(list(pop_0[t].keys()))
                auc_tot.append([sum(pop_0[t].values())]) #auc tot includes the initial population without growing
                dict_grid["t{0}".format(t)] = [[x] for x in pop_0[t].values()] 
                #replenish the initial nutrients and toxins at every round.
                dict_nuts["t{0}".format(t)] = [[x] for x in in_nuts] 
                dict_toxs["t{0}".format(t)] = [[x] for x in in_toxs]
            
            #START THE ROUND
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
                            df_mutants.loc[len(df_mutants)] = [rd, t, cell["t"], tuple(cell[par_mut]), tuple(new_cell[par_mut])]
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
                    
                    #add the total population per tube to later be able to calculate auc
                    auc_tot[t].append(count_population(tube))
      
                
                #STORAGE WITHIN A ROUND (AT EACH TIME STEP)
                #Store information at every tm at some rd, we can decide the periodicity:
                weplothere = 0 #just to save a loop with two conditions later on
                if rd < 5 or (rd+1)%save_periodicity == 0 or rd == (round_lim - 1):
                    weplothere = 1
                    for t in range(N_comms):
                        pop_tm = simpler_count(grid[t])
                        for i in range(len(sp_tube[t])):
                            #we have to check that the species did not go extint
                            try:
                                dict_grid["t{0}".format(t)][i].append(pop_tm[sp_tube[t][i]])
                            except:
                                dict_grid["t{0}".format(t)][i].append(0)
                        for j in range(len(in_nuts)):
                            dict_nuts["t{0}".format(t)][j].append(nuts[t][j])
                        for k in range(len(in_toxs)):
                            dict_toxs["t{0}".format(t)][k].append(toxs[t][k]) 
    
                #Increase the time step
                tm += 1
            
            #STORAGE AT THE END OF A ROUND (LAST TIME STEP)
            for t in range(N_comms): 
                #Degradation per tube:
                deg_list.append(consumption(toxs[t], in_toxs))
                #Nutrient uptake per tube:
                nut_list.append(consumption(nuts[t], in_nuts))
                #Store nutrients at the end of the round
                for j in range(len(nuts[t])):
                    df_nuts.loc[len(df_nuts)] = [rd, t, j, nuts[t][j]]
                #Store toxins at the end of the round
                for k in range(len(toxs[t])):
                    df_toxs.loc[len(df_toxs)] = [rd, t, k, toxs[t][k]]
                
                #get the population count at the last time step of the round:
                pop = count(grid[t]) #See the count function to understand the structure
                #Calculate population fold change for each species
                fc = []
                for sp,initial in pop_0[t].items():
                    try:
                        final = pop[0][sp] #pop is the final population
                    except:
                        final = 0 #if that species went extinct
                    fc.append(final/(initial + 0.00001)) #should not happen, but the 0.00001 is added to avoid division by 0.
                auc = [sum(auc_tot[t])] 
                in_pop = [auc_tot[t][0]]
                if weplothere:
                #If we have stored info per time step we can calulate the AUC per species
                    auc = []
                    in_pop = []
                    for x in dict_grid["t{0}".format(t)]:
                        auc.append(sum(x))
                        in_pop.append(x[0])
                
                #store general info
                df_grid.loc[len(df_grid)] = [rd, t, pop[2], sum(in_pop), sum(pop[0].values()), sum(auc), np.mean(fc), deg_list[t], nut_list[t]]
        
                #store species
                for index, k_sp in enumerate(sp_tube[t]):
                    auc_here = np.nan
                    in_pop_here = np.nan
                    if weplothere:
                        #they are calculated from dict_grid, so their index should match the one in sp_tube
                        auc_here = auc[index] 
                        in_pop_here = in_pop[index]
                    try:
                        final_pop = pop[0][k_sp]
                        #st_div = pop[3][index] #if some other species is missing index could not match, corrected below
                        st_div = pop[3][list(pop[0]).index(k_sp)]
                        
                    except:
                        final_pop = 0
                        st_div = np.nan
                    
                    df_sp.loc[len(df_sp)] = [rd, t, k_sp, in_pop_here, final_pop, st_div , auc_here, fc[index]]
                    
                #Store strains for each species
                for k_sp in pop[0].keys():
                    for k_st,v_st in pop[1][k_sp].items():
                        df_st.loc[len(df_st)] = [rd, t, k_sp, k_st, v_st]
    
    
            #---------------------#
            # PLOT WITHIN A ROUND #
            #---------------------#
            
            #PLOT ALL THE TUBES AT SOME ROUNDS
            columns = 7
            if weplothere: #the round periodicity we decided above
                fig=plt.figure(figsize=(38, int(5* int((N_comms+columns-1)/columns))))
                plt.suptitle("Round = " + str(rd), fontsize=16)
                rows = int((N_comms+columns-1)/columns)
                for t in range(N_comms):
                    plt.subplot(rows, columns, t+1)
                    top = 0 #just to get the highest population at some point in this tube (aesthetics)
                    for i,j in enumerate(dict_grid["t{0}".format(t)]):
                        top_here = max(j)
                        if top_here > top:
                            top = top_here
                        plt.plot(range(time_lim+1), j, linewidth=2,color=multiple_colors(sp_tube[t][i]))
                    plt.ylim (0, top* 1.3 + 1)
                    #Add in text some information of each tube
                    plt.text(0.3,top*1.24 + 0.9, "Sp: " + str(sp_tube[t])) #species present at time step = 0
                    plt.text(0.3,top*1.16 + 0.8, "N:" + str([round(x) for x in nuts[t]])) #nutrients remaining
                    plt.text(0.3,top*1.08 + 0.7, "D:" + str(round(deg_list[t],3))) #degradation score
                    
                    plt.title("tube"+str(t))
                    plt.xlabel("time steps", fontsize=12)
                    plt.ylabel("population", fontsize=12)
                    #for ax in fig.get_axes():
                        #ax.label_outer()
                #Save the plot
                plt.savefig( path + now + "/repeat" + str(rp)  + "/per_round" + "/round" + str(rd)+".pdf", dpi=80)
                plt.close(fig)
                #Save a list with dict_grid, dict_nuts and dict_toxs in this round (pickle format)
                open_file = open(path + now + "/repeat" + str(rp) + "/per_round" + "/round" + str(rd), "wb")
                pickle.dump([dict_grid, dict_nuts, dict_toxs], open_file)
                open_file.close()
                #Save the species_tube in pickle format (to reload as dict and use if needed)
                open_file = open(path + now + "/repeat" + str(rp) + "/sp_tube" + str(rd), "wb")
                pickle.dump(sp_tube, open_file)
                open_file.close()
    
    
            #---------------#
            #  PROPAGATION  #
            #---------------#
            if rd < (round_lim - 1):
                #CHOOSE TUBES
                #create a list with tuples including (deg_value , index) of the selected tubes to transfer
                if treat == "s": #selection treatment
                    selected_tubes = sorted(((value, index) for index, value in enumerate(deg_list)), reverse=True)[:N_beta]
                elif treat == "r": #random treatment
                    #position 0 of the tupple is 1 since degradation here does not matter-> equal weights in disassembly
                    selected_tubes = [(1,tube) for tube in rng.choice(N_comms, N_beta, replace = False)]
                #We save the selected tubes
                df_tubes.loc[len(df_tubes)] = [rd]+[x[1] for x in selected_tubes]
                
                #Check if all the selected tubes have enough cells (epsilon amount of cells), otherwise chose another tube 
                #and report it
                if prop == "n":
                    pass #in previous versions, the run was stopping if some tube went extinct, but that is not a problem
                    #for t in range(N_comms):
                        #if len(grid[t]) == 0:
                            #rd = round_lim
                            #break
                else: #propagule, migrant pool, or disassembly
                    extinction_checking = 0
                    for i_t, tup in enumerate(selected_tubes):
                        if count_population(grid[tup[1]]) < epsilon: 
                            extinction_checking +=1 
                            #when a tube does not have enough cells, select other tube, but keep track that it is something
                            #anecdotic, since it could introduce some bias in the random treatment.
                            #Under selection, take the following tube in the ranking
                            if treat == "s":
                                while (N_comms-N_beta) > extinction_checking:
                                    new_selected_tube = sorted(((value, index) for index, value in enumerate(deg_list)), \
                                                               reverse=True)[N_beta - 1 + extinction_checking : N_beta  + extinction_checking]
                                    if count_population(grid[new_selected_tube[0][1]]) < epsilon: 
                                        extinction_checking +=1
                                    else:
                                        break
        
                            #Under random treatment, the new tube is randomly chosen
                            if treat == "r":
                                while (N_comms-N_beta) > extinction_checking:
                                    available_tubes = set(range(N_comms)) - set([x[1] for x in selected_tubes])
                                    new_selected_tube =  [(1,tube) for tube in rng.choice(list(available_tubes), 1, replace = False)]
                                    available_tubes -= {new_selected_tube[0][1]}
                                    if count_population(grid[new_selected_tube[0][1]]) < epsilon: 
                                        extinction_checking +=1
                                    else:
                                        break
                            selected_tubes[i_t] = new_selected_tube[0]
                    #report if there was a extinction
                    if extinction_checking:
                        f=open(path + now + "/repeat" + str(rp)+ '/extinctions.txt','a')
                        f.write("Repeat: {}; Round: {}; Extinctions: {} \n".format(rp,rd,extinction_checking))
                        f.close()
                    if extinction_checking == (N_comms-N_beta):
                        #print("Ups no more cells")
                        rd = round_lim #so I break the loop and the repeat gets extinct 
                #*********#
                # RESTART #
                #*********#
                #For the next round, refresh nutrients and toxins, then the cell grid will be restarted in the propagation 
                #method:
                nuts = [] #Grid to store nutrients
                toxs = [] #Grid to store toxins
                for t in range(N_comms):
                    nuts.append([*in_nuts]) 
                    toxs.append([*in_toxs])
                
                #-----------------#
                #   NO SELECTION  #
                #-----------------#
                #this is just a control where we just dilute all the tubes, thus w/o any selection
                if prop == "n":
                    #we reseed the cells
                    grid2 = [] #new grid to store cells gor the next round
                    ancestors = [] #track which tube is comming from which tube
                    for t in range(N_comms):
                        grid2.append(dilution(grid[t],d_n)) #dilution deactivates those activated cells
                        ancestors.append(t)
                    grid = copy.deepcopy(grid2)
                    df_ancestry.loc[len(df_ancestry)] = [rd+1]+[t for t in ancestors]
                
                #-------------#
                #  PROPAGULE  #
                #-------------#
                elif prop == "p":
                    #we reseed the cells
                    grid2 = [] #Grid to store cells 
                    group_size = int(N_comms/N_beta) 
                    reminder = N_comms - group_size*N_beta #this is just to adjust in case number of tubes/ tubes to next 
                                                            #is not exact
                    ancestors = []                                                    
                    for prior in range(N_beta):
                        if prior == 0:
                            fix = 0
                        else:
                            fix = 1
                        for t in range(prior*group_size + reminder*fix, (prior+1)*group_size + reminder): 
                            grid2.append(dilution(grid[selected_tubes[prior][1]],d_p)) #we deactivate those activated cells
                            ancestors.append(selected_tubes[prior][1])
                    grid = copy.deepcopy(grid2)
                    df_ancestry.loc[len(df_ancestry)] = [rd+1]+[t for t in ancestors]
                    
    
                
                #---------------------------------#
                #  PROPAGULE  +  INVASION POISSON #
                #---------------------------------#
                elif prop == "pip":
                    #we reseed the cells
                    grid2 = [] #Grid to store cells 
                    group_size = int(N_comms/N_beta) 
                    reminder = N_comms- group_size*N_beta #this is just to adjust in case number
                    ancestors = []                        #of tubes/ tubes to next is not exact
                    for prior in range(N_beta):
                        if prior == 0:
                            fix = 0
                        else:
                            fix = 1
                        for t in range(prior*group_size + reminder*fix, (prior+1)*group_size + reminder): 
                            grid2.append(dilution(grid[selected_tubes[prior][1]],d_p)) #we deactivate those activated cells
                            ancestors.append(selected_tubes[prior][1])
                            
                    grid = copy.deepcopy(grid2)
                    df_ancestry.loc[len(df_ancestry)] = [rd+1]+[t for t in ancestors]
                    
                    
                    
                    #INVASION
                    tubes_to_migrate = rng.choice (N_comms,int(N_comms*f_migr), replace = False)
                    #for now I use the same parameter (prob) as in disassembly, but considering that we do
                    #not know about the species that are there and any species can arrive
                    for t in tubes_to_migrate:
                        num_invaders = 1 + rng.poisson(0.5)
                        sp_to_include = list(set(sp_in_play))
                        for i_i in range(num_invaders):
                            invader = rng.choice(sp_to_include) #invader is the species number
                            grid[t].extend(sample_or_amplify(st_repo[invader], b_d_1))
                            sp_to_include = list(set(sp_to_include)-{invader}) #The added species will not be added again
                        
                        rng.shuffle(grid[t]) 
    
    
                
                #----------------#
                #  MIGRANT POOL  #
                #----------------#
                #ancestry in migrant pool would be just mixing all the tubes in df tubes
                elif prop == "m":
                    #we store the cells to transfer
                    pool = []
                    for tup in selected_tubes:
                        pool.extend(grid[tup[1]])
                    grid = [] #I am not using the grid anymore, so I can erase it
                    for t in range(N_comms):
                        grid.append(dilution(pool,d_m/N_beta)) #here percentages is scaled by the number of tubes we have mixed.
                    
    
                
                #-----------------------------------#
                #  MIGRANT POOL + INVASSION POISSON #
                #-----------------------------------#
                elif prop == "mip":
                    #we store the cells to transfer
                    pool = []
                    for tup in selected_tubes:
                        pool.extend(grid[tup[1]])
                    grid = [] #I am not using the grid anymore, so I can erase it
                    for t in range(N_comms):
                        grid.append(dilution(pool,d_m/N_beta)) #here percentages is scaled by the 
                                                                   #number of tubes we have mixed.
            
                    #INVASION
                    tubes_to_migrate = rng.choice (N_comms,int(N_comms*f_migr), replace = False)
                    #for now I use the same parameter (prob) as in disassembly, but considering that we do
                    #not know about the species that are there and any species can arrive
                    for t in tubes_to_migrate:
                        num_invaders = 1 + rng.poisson(0.5)
                        sp_to_include = list(set(sp_in_play))
                        for i_i in range(num_invaders):
                            invader = rng.choice(sp_to_include) #invader is the species number
                            grid[t].extend(sample_or_amplify(st_repo[invader], b_d_1))
                            sp_to_include = list(set(sp_to_include)-{invader}) #added species are not reintroduced
                        rng.shuffle(grid[t]) 
    
                        
                

                
                #-----------------#
                #  DISASSEMBLY_3  #
                #-----------------#
                #Current version of disassembly we are using, main changes vs previous version:
                #here we have all the same species per tube, but we would like to seed tubes so we maintain the strains in each selected tube tube.
                #otherwise in the implemented d_2, we would just take the strain of each species coming from the best performing tube.
                #we still keep the st_repo for migration.
                #minor change: I also updated so there is no invasion if the 15 species are present.
                
                elif prop=="d3": 
                    N_spc = [len(sp_tube[x[1]]) for x in selected_tubes] #number of species ger tube
                    E = [0 for i in N_spc] #counter of extinction for each of the selected tubes
    
                    selected_sp = [[] for i in range(N_beta)] #species present in the selected tubes at the beginning of the round
                    
                    #UPDATE THE FOSSIL RECORD OF STRAINS (although we later take all the strains, here we just save one representative of each species)
                    control_update = [0 for sp in range(total_sp)]
                    for i,tup in enumerate(selected_tubes): 
                        for sp in sp_tube[tup[1]]: #species in the tube at the beginning of the round
                            selected_sp[i].append(sp)
                            list_now = sp_finder(grid[tup[1]],sp)
                            if count_population(list_now) >= s_d:
                                if control_update[sp] == 0:
                                    st_repo[sp] = copy.deepcopy(list_now)
                                    control_update[sp] = 1 #like that we only update with the strains from the best possible tube
                            else:
                                E[i] += 1
    
                    #RESCALE DEGRADATION SCORE (D_hat)
                    D_hat = [] #we score times penalty
                    for i,tup in enumerate(selected_tubes):
                        #D_hat.append(tup[0] * (N_spc[i] - E[i])) 
                        D_hat.append(tup[0] * ((N_spc[i] - E[i])/N_spc[i]) )
                    #normalise so the scores sum up to 1
                    norm_D_hat = [x/sum(D_hat) for x in D_hat] 
                    
                    #CREATE THE PUTATIVE TUBES (weighting by D_hat how many communities each selected tube will generate)
                    i_pt = rng.choice(N_beta, N_comms, p=norm_D_hat) #index of the putative tubes
                    putative_tubes = []
                    sp_we_have = [] #species in the metacommunity, needed for the migration step, to prioritise those not present species
                    ancestors = [selected_tubes[x][1] for x in i_pt]
                    df_ancestry.loc[len(df_ancestry)] = [rd+1]+[t for t in ancestors]
                    for t in i_pt:
                        sp_we_have.extend(selected_sp[t])
                        pt_here = []
                        for sp in selected_sp[t]:
                            #THE FOLLOWING LINE IS WHAT I CHANGE IN d_3
                            #even though we have the strain re
                            if sp_finder(grid[selected_tubes[t][1]],sp):
                                pt_here.append(sample_or_amplify(sp_finder(grid[selected_tubes[t][1]],sp),b_d_1))
                            else: #could happen that the sp went stinct so we take it from the repo
                                pt_here.append(sample_or_amplify(st_repo[sp],b_d_1))
                        putative_tubes.append(pt_here)
                        
                    #IMMIGRATION-EMIGRATION STEP
                    #choose the tubes from where species will be removed (EMIGRANT)
                    tubes_to_emigrate = rng.choice (N_comms,int(N_comms*f_migr), replace = False)
                    for t in tubes_to_emigrate:
                        num_emigrants = 1 + rng.poisson(0.5)
                        for i_i in range(num_emigrants):
                            if len(putative_tubes[t])>1: #we don´t want to end up w/o species
                                # which of the species already present in the tube will be replaced (we try to not remove 
                                #those unique species in the metacommunity)
                                emigrant = rng.choice(len(putative_tubes[t])-1) #emigrant is a position (1st, 2nd...)
                                #but let´s check that it is not unique in the metacomunity
                                options = list(range(len(putative_tubes[t]))); rng.shuffle(options)
                                for x in options:
                                    #we check to not remove a unique species in the metacomunity
                                    if sp_we_have.count(putative_tubes[t][x][0]["t"]) > 1:
                                        emigrant = x
                                        break
                                    #if all are unique we will still make the replacement with a random emigrant species selected above
            
                                #we update the list and populate the putative tube
                                #emigrant is a position in the tube
                                #invader is the number of th species
                                sp_we_have.remove(putative_tubes[t][emigrant][0]["t"])
                                del putative_tubes[t][emigrant]
                                
                    #choose the fraction of the tubes to get a IMMIGRANT (or invader), in this case we still get it from the st_repo
                    tubes_to_immigrate = rng.choice (N_comms,int(N_comms*f_migr) , replace=False)
                    #select the species that will arrive to the tube, prioritisig those not present in the metacommunity
                    for t in tubes_to_immigrate:
                        num_invaders = 1 + rng.poisson(0.5)
                        for i_i in range(num_invaders):
                            sp_to_include = list(set(sp_in_play) - set(sp_we_have))
                            if len(sp_to_include) > 0:
                                invader = rng.choice(sp_to_include) #invader is the species number
                            else:
                                #if all the species are in the metacommunity we randomly select (except for those species
                                #already present in this particular tube)
                                sp_in_tube = [x[0]["t"] for x in putative_tubes[t]]
                                if list(set(sp_in_play) - set(sp_in_tube)):
                                    invader = rng.choice(list(set(sp_in_play) - set(sp_in_tube)))
                                #UPDATE: I set as comments the folowing two lines since we don´t want to invade if the species is already there.
                                #else: #if all the species are in the tube, this would be like not invading
                                    #invader = rng.choice(list(set(sp_in_play)))
    
                            putative_tubes[t].append(sample_or_amplify(st_repo[invader], b_d_1))
                            #and update the list traking the species that we have in the metacomunity
                            sp_we_have.append(invader)
                    
                    #RESTART THE MATRIX (GRID)
                    grid = [[] for x in range(N_comms)]
                    for t in range(N_comms):
                        for sp_list in putative_tubes[t]:
                            #we amplify to b_d_2 (now I removed it os b_d_2 = b_d_1, but I leave this 
                            #flexibility for the future).
                            if b_d_1 == b_d_2:
                                grid[t].extend(sp_list)
                            else:
                                grid[t].extend(sample_or_amplify(sp_list, b_d_2)) 
                        rng.shuffle(grid[t]) 
                        
        
            #INCREASE THE ROUND COUNTER
            rd += 1
            
            #-------------------------------------------------
            #save initial grid
            if rd == (round_lim - 1) or (rd+1)%save_periodicity == 0: #note that rd has been increased just above
                #store the intial grid (round periodicity can be modified)
                open_file = open(path + now + "/repeat" + str(rp) + "/grid0_round_" + str(rd), "wb")
                pickle.dump(grid, open_file)
                open_file.close()
        
        #STORE IN THE LAST ROUND
        df_grid.to_csv(path + now + "/repeat" + str(rp) + "/df_grid.csv")
        df_sp.to_csv(path + now + "/repeat" + str(rp) + "/df_sp.csv")
        df_st.to_csv(path + now + "/repeat" + str(rp) +"/df_st.csv")
        df_mutants.to_csv(path + now + "/repeat" + str(rp) +"/df_mutants.csv")
        df_nuts.to_csv(path + now + "/repeat" + str(rp) +"/df_nuts.csv")
        df_toxs.to_csv(path + now + "/repeat" + str(rp) +"/df_toxs.csv")
        df_tubes.to_csv(path + now + "/repeat" + str(rp) +"/df_tubes.csv")
        df_ancestry.to_csv(path + now + "/repeat" + str(rp) +"/df_ancestry.csv")
        #store the grid just in case
        open_file = open(path + now + "/repeat" + str(rp) + "/grid100_round_" + str(rd-1), "wb")
        pickle.dump(grid, open_file)
        open_file.close()
        
        t2 = timer.time() 
    
    
    
        f=open(path + now + '/timer' + str(rp) + '.txt','w')
        f.write("Time it took: "+ str(t2-t1) +'\n')
        f.close()
    
    
    
    
    #%% Final section
