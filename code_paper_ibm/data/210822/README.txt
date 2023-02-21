THIS FOLDER CONTAINS THE MAIN RESULTS WITH MY FINAL MODEL , AND THUS THE DATA TO INCLUDE IN THE PAPER(AT LEAST UP TO DATE)
Data generated with the very last version of the data (here the new thing is that the seed is restarted at each repeat,
so I can paralelise at the repeat level if I want; due to stochasticity this might change a bit the results). 

Anyway, here I have the conditions I have been using:
5 species sets of 15 species (seeds 22-26, not really a reason for that)
	10 repeats
		50 rounds
21 tubes, and I repeat all that for all the propagation methods under selection(_s) or random(_r) tratment:
	disassembly old(d_), 
	disassembly decoupling immigration-emigration(d2_),
	disassembly current version (d3_)
	propagule(p_), 
	migrant_pool(m_), 
	propagule with invasion (pip_), 
	migrant pool with invasion(pip_)
and no selection (n_r)

NOTE IN THIS VERSION TO UPLOAD ONLINE I REMOVE SOME FILES TO SAVE SPACE:
Within each propagation: 
	- per_round directory (which has per species growth plots and data structures with information for each round)
	- grid0 round 10-40 (the initial grid, I only save the one at round 49)
	- df_mutants (a data frame tracking all teh mutations that appear and their ancestry, but that we are not using in the end)