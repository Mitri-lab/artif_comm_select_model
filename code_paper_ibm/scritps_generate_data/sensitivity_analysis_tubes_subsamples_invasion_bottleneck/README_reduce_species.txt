Here I include results when reducing the number od species. 
Still keep the same sets of 15 species, but I subsample a number of species of the size indicated in the folder as seed_sizesubsample
Since the species that happen to be in the subsample matter, I take 6 subsamples in each case (and for each of them I run 5 repeats).
The script is prepared to run from unix (so we can paralelise, for instance by species sets), see runloop_reduce_species file to see how we run.
The seed (that generates the species sets) should be chosen externally in unix; then within the python script we should choos which propagation
methods to run "list_prop" and "list_treat", we can also choose how many species has a subsamle "subsample_size", or how many subsamples of each 
size we take "subsample_index", then we can also choose how many repeated runs we do for each subsample "list_repeats".

In the results we´ll have folders named as the species set followed by the number of species in teh subsample. Inside we´ll have a folder corresponding
to the particular subsample with tha number of species. Inside, we´ll have a folder corresponding to each propagation method. The fist part of the name indicates the method:
- d3 : disassembly 
- p : propagule
- m : migrant pool
- pip : propagule with invasion
- mip : migrant pool with invasion
- n_r : no selection
The last part of the name indicates whether the treatment was selection (_s) or random control(_r)

Inside the folder for each propagation method, we´ll have one folder for each repeat (a repeat is a different arrangement 
of the species in community at the beginning of the experiment, but always keeping the same 15 species for the 10 repeats).
You´ll see inside the .csv files with all the information (you can see a description in the word document "ibm_guide") and
a folder per_round where we include plots at some particular rounds.

