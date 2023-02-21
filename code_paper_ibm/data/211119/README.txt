Here we run the evolved communities in the last round with ancestral species. This is, we take the communties
in the last round, and put the same species in the same initial populations but replacing evolved strains by
ancestral species. Then we let this community grow for 1 round of 80 time steps.

In order to run it you should choose within the code the seed (species set to run, here we used 22-26; the seed 
will work as starting poing for the random number generators to first generate the species and then run the model).
You should also choose and propagations you want to run ("list_prop" and "list_treat")

Each folder corresponds to a species set (here named 22-26 as the seeds we used for the random number generator in
the script). Inside, we have a folder corresponding to each propagation method (from which we take the last community
in this case). The fist part of the name indicates the mthod:
- d3 : disassembly (the number of species can change due to immigration and emigration)
- p : propagule
- m : migrant pool
- pip : propagule with invasion
- mip : migrant pool with invasion
- n_r : no selection
The last part of the name indicates whether the treatment was selection (_s) or random control(_r)

Inside the folder for each propagation method, we have one folder for each repeat (a repeat is a different arrangement 
of the species in community at the beginning of the experiment, but always keeping the same 15 species for the 10 repeats).

You can see inside the .csv files with all the information (you can see a description in the word document "ibm_guide")
From the current directory I have removed to save space (within each repeat of each propagation method):
	- "per_round" folder. This folder included a growth plot and some structure storing information.
	-  "grid0" and grid100 files, which are structures storing the whole system at the beginning and the end of the simulation, could be useful	+
		to check something but not needed for the plots.

