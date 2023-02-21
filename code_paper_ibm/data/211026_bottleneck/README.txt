In this folder I run the part of the "sensitivity analysis" where I cahnge the bottleneck to 0.2, 0.5, 2, or 5 (or 20) times the original. 
To run this script I´ll use 5 cores, one per seed. The script (211025_model_strains_poisson_bottleneck_labcomp) is written in python 
but designed to be run from unix, like that we can externally call each seed (see "runloop" file), and paralelise the run of the 5 species sets.
If you rader you can just input the seed in python and run it from a python IDE.
Within the python code you should choose "d_factor" times which we will multiply the tilution factor.
Also you should choose the propagation methods in "list_prop" and "list_treat".
We do 10 repeats for each seed. Same conditions, should be the standard ones we used for the simulations.


The bottleneck is 10 cells (approx) in disassembly, and 5% of the tube in the other methods, in order to change the
bottleneck we multiply it times a factor of 0.2, 0.5, 2, or 5. For each condition the whole analysis is done.

We have folders named as teh factor times which we multiply the bottleneck. (1 would be the standard which is 210822 folder)
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
a folder per round where we include plots at some particular rounds. Well, here I have deleted almost every file to save space
and just keep df_grid which is what we use to generate the plots.

