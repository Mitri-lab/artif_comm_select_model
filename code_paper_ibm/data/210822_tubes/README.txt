Simulations changing the number of communities (tubes) in the metacommunity for the experimental set-up. 
Those with 21 tubes are comming from 210822 simulations
Each folder is named as the number of communities the simulations inside have.

The script to run this simulations is in this folder "211018_model_strains_poisson_tubes_labcomp"
Within the script we just need to chose "N_comms" (here as a list and the simulation for each metacommunity size will be run one after another)
We also should choos which propagations to run with "list_prop" and "list_treat"
The model is designed to be run from unix, and the seed for the random number generator (and that generates the species set) is feed externally.
You can check how to run in python in the file "runloop_tubes_labcomp", or adapt and put a seed as a number to directly run with a python IDE. 

Inside each folder corresponds to a species set (here named 22-26 as the seeds we used for the random number generator in
the script). Inside, we have a folder corresponding to each propagation method. The fist part of the name indicates the method:
- d3 : disassembly 
- p : propagule
- m : migrant pool
- pip : propagule with invasion
- mip : migrant pool with invasion
- n_r : no selection
The last part of the name indicates whether the treatment was selection (_s) or random control(_r)

Inside the folder for each propagation method, we have one folder for each repeat (a repeat is a different arrangement 
of the species in community at the beginning of the experiment, but always keeping the same 15 species for the 10 repeats).

You can see inside the .csv files with all the information (you can see a description in the word document "ibm_guide") and
a folder per round where we include plots at some particular rounds. Well, although here I to save space I have removed basically
every file and only keep "df_grid" which is the one we are using for the plot.

