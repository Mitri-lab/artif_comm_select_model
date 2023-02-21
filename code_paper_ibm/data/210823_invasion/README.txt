
Results when changing the invasion fraction (i.e. number of tubes where the invasion occurs at each propagation step),
to run for all the methods having invasion (i.e. diassembly "d3", mifrant pool with invasion "mip" and propagule with invasion "pip"). 
The script to generate the data is present in this folder (211018_model_strains_poisson_invasion_marc_computer), it is written in python,
but was runned in unix to paralelise one run per seed (=species set), thus we use 5 cores to run it.In "runloop" txt you can see how we run it.
Within the script we should choose the fraction of tubes with migration ("f_migr"), right now we input it as a list and each condition is run one
after another. We should also choose which propagations to run in list_prop and list_treat, for each condition, they will be run one after another. 
We do 10 repeats for each seed. Same conditions should be the standard ones we used for the simulations. If you prefer to simple run for instance
with a python IDE you can input the seed you want (we used 22-26) and run one at a time.

Each folder is named as teh corresponding f_migr value (fraction of tubes with migration = int(f_mirg*21)), changing the "." by "_". 
	0_15
	0_25
	0_35
	0_45
	0_55
	0_62
	0_715
	0_81
	0_905
	1

Anyway, the right thing to do might be just to indicate the number of tubes that will get invasion, respectively: 1, 3, 5, 7, ..., 21.

The results of the folder 0_25 are stored in 210822 simulation.

Inside the folder for each propagation method, we have one folder for each repeat (a repeat is a different arrangement 
of the species in community at the beginning of the experiment, but always keeping the same 15 species for the 10 repeats).

You can see inside the .csv files with all the information (you can see a description in the word document "ibm_guide") and
a folder per round where we include plots at some particular rounds. Well, although here I to save space I have removed basically
every file and only keep "df_grid" which is the one we are using for the plot.




