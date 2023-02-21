Here we have the code to generate data for the "sensitivity analysis" we show in Fig6B-E.

Note that these scripts generate a lot of data, so can become quite heavy to run. Basically
we change the conditions to:

- invasion: we change the number of tubes with invasion (among those propagation methods having
	invasion) at the propagation step.

- tubes: different numeber of tubes (communities) in the metacommunity for the evolutionary rounds.

- reduce_species: the metacommunity will contain not 15 species but a lower number.

- bottleneck: change the dilution bottleneck at the propagation step.

All these scripts where adapted to run from Unix with the corresponding "runloop", to see
the specificities of each script how to run and how data is stored, check the corresponding
"README" file in each case.