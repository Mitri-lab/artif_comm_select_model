##### PARAMETER SETTINGS FOR:
###

#### EXPERIMENTAL and COMMUNITY PARAMS
Ntubes = 21                   # Number of test tubes <==> number of communities
Nspec = 4                     # Number of species, originally and per community, excluding mutants
Nnutr = 4                     # Number of nutrient substrates
Ntox  = 10                    # Number of toxic substrates
Ntransfers = 50               # Number of transfers to simulate
Nreps = 10
Nspecpool = 15

Ntsf_stability = 25

y0_spec = 1.0e2                # Initial abundance of cells per species
y0_nutr = 1.0e2                # Initial (and inflow!) concentration of nutrients
y0_tox  = 1.0e2                # Initial (and inflow!) concentration of toxins

dilution   = 1.0e3             # Dilution rate at transfer
frac_bneck = 0.33              # Selection bottleneck - which fraction of communities are transferred
mixrate = 0.22
rate_migr_poisson = 0.5

#### COMPUTATIONAL PARAMS FOR BATCH CULTURE
Ntstep = 1000
t0, tstop = (0.0, 100)

#### MUTATION RATE and EFFECT SIZE
mutdev   = 0.4     # Mutation drawn from logNormal(1.0, mutdev)
mutfac   = 2.0     # When mutations drawn from Beta(mutfac,mutfac)
mutrate  = 5.0e-2  # Mutation rate per transfer event
