Here I save the code for the IBM of the actual version of the paper, note that we discovered some errors there, we simulated corrected them and the trends are maintained,
but we want to correct them when the reviewers ask for some extra analysis (in the actual paper the mistakes remain). Mistakes/bugs are the following ones:

1- Dilution without replacement for NS, propagule and migrant pool. since it is without replacement I was removing the cells from the original tube, 
then the population is shrinking and so thus the probability to put cells in the further new tubes comming from that selected tube. thus in propagule 
every 3 tubes, tube 1 will have the highest amounts of cells and tube 3 the lowest. Then in migrant pool every new tube will have less cells that the 
previous one. This can be seen in the presentation 20221024_meeting

2- Error that could lead to stop growth and degradation of a strain in some cases when a nutrient is almost over, and if no other
strain consume that little bit of that nutrient, it will remain like that untill the end of the run. Probably does not have a big impact,
but still. More detailed explanation:  When a strain consumes nutrient j and this nutrietn is in an amount n_ij> N_j > new_n_ij 
(where new_n_ij is the rescaled version when some other nutrients are over), this strain will stop growing and dividing in this round unless another strain
consumes a bit of this nutrient for it to be consider depleted by this strain, to solve this when a strain checks if a nutrient is enough, 
instead of chekcing if N_j > n_ij it´s better to just do N_j>1, like that we won´t consume nutrients if they are bellow 1 unit, but we solve the problem.

3- Minor mistake I found in disassembly emigration: when we randomly choose the species to be removed from a tube (well is not totally random as we check 
it is not unique in the metacommunity), I was not including one of the species in the tube (the last one due to range that does not inlcude the last number).


