Here we use the results 210822, and take the community that renders the best results at the end of round 50 (round 49 when starting at 0),
in particular at time step 0 (initial bacteria present there). Then we run 10 exact replicates of this initial grid during 25 rounds under no selection treatment. 
We do this with two versions:
- In the regular one it´s the exact community with the evolved species that is grown for 25 rounds without selection.
- In the directories named at the end "_anc" (generated with the script named at end "_ancestral" ), we take this initial community round 50, 
	and within the community we replace each instance of a species (that could have mutations caused by evolution) by ancestral species. 

The cody is writen in python and can be run for instance with a python IDE (I was using spyder). Before runing you should choose the seed (species set, here 22-26)
to run and also the propagation methods (as a list in "props_grid") to run (will be run one after another). Note that seed and propagation method should match the results from where we want
to take the best community and run without slection (NS propagation)

Results will be stored in folders named as the propagation method. Inside each folder corresponds to a species set (here named 22-26 as the seeds we used for the random number generator in
the script). Inside, we have a folder corresponding to each propagation method. The fist part of the name indicates the method
(from which we take the final grid in this case):
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
a folder per round where we include plots at some particular rounds.



