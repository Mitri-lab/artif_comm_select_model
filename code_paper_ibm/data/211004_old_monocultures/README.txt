Here we have ancestral monocultures corresponding to the evolved communities. For each of the communities in the 
last round, we take separatelly each of the species, and run separatelly each of the species present (the evolved strains)
and run a simulation with them for a run. Like, we take the initial grid at round number 50 (49 if we start to count from 0), and 
run monocultures of the species that populate each community.
Each monoculture is run in triplicates and then we take the average.

The model is written in python, so can be run for instance in a python IDE. 

In the script we should feed the last grids from the main data (210822) and also chose which propagtions we want to run (they 
will be run one after another, they should be indicated in a list named "props_grid". If we choose a lot of propagations here, 
at least my computer runs out of memory and stops at some point, so at least in a regular computer better to do like 3-4 per run.
