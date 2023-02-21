In this folder we store the data from the IBM required to generate the plots in the actual version of the paper
Directories are named as referred in the plotscrips to avoid confusion, what they include:

- 210822: main results of the simulations (needed for all figures from the IBM)

- 210822_2_ok: transfer the best communties that evolved in each propagation method,
	just diluting them without selection. Needed for Fig6A

- 210822_subsamples: part of the sensitivity analysis, run the experiment with a lower number of 
	species in the metacommunity. Needed for Fig6C.

- 210822_tubes: part of the sensitivity analysis, run the experiment with a different number of 
	communities (tubes) in the metacommunity. Needed for Fig6D.

- 210823_invasion: part of the sensitivity analysis, change the number of tubes with invasion 
	(for the propagation methods with invasion) at each transfer event. Needed for Fig6B.

- 211004_old_monocultures: we run monocultures of those species that evolved in the final communities.
	Needed for Fig3G.

- 211005_all_communities: Run all the possible combinations of 15 species (the 32762 possible ones) 
	Needed for Fig2E, Fig2F, and Fig3G.

- 211026_bottleneck: part of the sensitivity analysis, change the dilution bottleneck at each transfer event.
	Needed for Fig6E.

- 211119: take the communities in the last round of the experiment, replace the evolved species by
	the corresponding ancestral ones and run one round of growth. Needed for Fig5B and Fig5C.



In each of the cases we include the python script to generate these data, as well as a README file
with some extra explanations. To understand the structure of the model you can check the script
that generates our main results ("210822" directory), and the explanations we give in the ibm_guide.
(all the other versions are based on this main one). In ibm_guide you can also find an explanation
of what contains all the data that the model generates, although most of it has been removed from
these directories to save space.

Directories can be named as different species sets we use (22-26, named as the random seed that generated them),
the different propagation methods we use that can be under selection (_s) or random (_r) treatment. Sometimes
they also refer to the repeat (use the same conditions and species but with a different random initial combination 
of species). The rest might be specific to each particular simulation, check the README in teh directory or the ibm_guide.
 