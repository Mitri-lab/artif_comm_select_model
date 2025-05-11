
##### SCRIPT TO INITIALISE SPECIES POOL FOR MULTIPLE REPEATED RUNS OF A SELECTION EXPERIMENT
## Using
##   selsims_coexmut_hilldegr.py
##   selection_ODEmodels_candidates.py
##   set_cases.py
##   growthpars.py
## Written by bvessman@gmail.com
## 2021.04.15
## Cleaned on 2025.05.11
###

### Import standard modules
import sys
import copy
import scipy
import argparse
import itertools
import numpy as np
from pathlib import Path, PurePath
from datetime import date, datetime

## Collect time stamp
datestr = datetime.now().strftime("%Y%m%d")
timestr = datetime.now().strftime("%Y%m%d-%H%M%S")

## Set random number generator
parser = argparse.ArgumentParser(description='Set a random seed.')
parser.add_argument('seed', metavar='N', type=int, nargs='?',
                    help='an integer to set a random seed')

args = parser.parse_args()
if args.seed:
    rg_seed = args.seed
else:
    rg_seed = 0000
rng = np.random.default_rng(seed=rg_seed)

### Import own libraries
import casepars
import growthpars
import AST_printfuns # as printfs
import AST_SelRepeats as repfuns
import seltheory_functions as selfuns
import seltheory_ODEmodels as smodels
from selsimsettings import *

###################################################################
###################################################################
######## Initialisation
###### Step 0: Define simul params

#### EXPERIMENTAL and COMMUNITY PARAMS
Nreg = Nnutr +Ntox
Nplaces = 2*Nreg

path_basedir = PurePath(casepars.basedir)
sizedir = str(Nspec)+'x'+str(Nnutr)+'x'+str(Ntox)

#### OUTPUT PARAMS
pyversion = sys.version_info[1] # To decide how to format print

#### Step 1a: BASE DIRECTORY for OUTPUT
##   Check if directory exists, otherwise create it
datedir = datestr+'_'+str(rg_seed)
path_outdir = path_basedir / sizedir / datedir
Path(path_outdir).mkdir(parents=True, exist_ok=True)

#### Step 1b: HEADERS for OUTPUT FILES
headtup_main  = ('Tsf', 'Tube')
## States
headtup_pops  = headtup_main+('Species', 'Abundance')
headtup_AUCs  = headtup_main+('Species', 'AUC')
headtup_nutrs = headtup_main+('Nutrient', 'Concentration')
headtup_toxs  = headtup_main+('Toxin', 'Concentration')

## Community and species parameters for all communities from species set
headtup_comms   = ('Tube','Tdegr','Ndegr','Div','AUC','Popsize')
headtup_species = ('Tube','Species','Popsize','AUC','sum_f')

## Growth parameters, fixed and mutable
headtup_mutpars = headtup_main+('Spec','Substr','Value')
#headtup_fixpars = ('Tube','Spec','Substr','Value') # Growth, death, investment
#                                                   # r, m, f are matrices
#                                                   # of size [Nspec x Nsubstr]

## Species ID and ancestry - tube m originated in previous tube n
headtup_specID  = headtup_main+('Tube_prev', 'Species_index', 'Species_ID')

## Fossil record - for each species, save the last occurence
#headtup_fossil  = ('Tsf','Species','Tube_origin')

## Collect all headers
headertups = (headtup_pops, headtup_AUCs,
              headtup_nutrs, headtup_toxs,
              headtup_mutpars,
              #headtup_fixpars,
              headtup_specID)

headtups_allcomms = (headtup_comms, headtup_species)

#### Step 2: PREPARE PARAMS AND DRAW SPECIES POOL
exppars_specpool = {'size':(Nspecpool, Nnutr, Ntox),       ### Nx, Ny, Nz defined in selsimsettingsXYZ.py
                    'sparse':(growthpars.sparseflag,
                              growthpars.fracsparse)} ### sparseflag and fracsparse defined in growthparsXYZ.py

grpars  = {'r0':growthpars.r0,
           'm0':growthpars.m0,
           'f0':growthpars.f0(Ntox),
           'Y0':growthpars.Y0,
           'd0':growthpars.d0,
           'KN0':growthpars.KN0,
           'KT0':growthpars.KT0,
           'nHill':growthpars.cHill,
           'fpow':growthpars.fpow}

pars_specpool = selfuns.draw_CRparams_longform(rng,exppars_specpool,grpars)
pars_specpool['specid'] = list(range(Nspecpool))
pars_specpool['nHill']  = growthpars.cHill
pars_specpool['fpow']   = growthpars.fpow
#pars_specpool = {
#    'specid': list(range(Nspecpool)),
#    'finv': fmat, # [Nspec, Ntox]
#    'm': mmat,    # [Nspec, Ntox]
#    'r': rmat,    # [Nspec, Nnutr]
#    'KN': KNvec,  # [Nspec, ] OR SCALAR
#    'KT': KTvec,  # [Nspec, ] OR SCALAR
#    'Y': Yvec,    # [Nspec, ]
#    'd': dvec}    # [Nspec, ]

#### Step 2: PRINT SPECIES POOL
#print('Printing common parameters to ... ')
#print(path_outdir)
#print('Timestamp '+timestr+' printed to file Timestamp')
AST_printfuns.print_specpool(path_outdir, timestr, pars_specpool)

#### Step 3: ASSEMBLE PARAM SETS
exppars =  {'y0substrate':(y0_nutr, y0_tox),
            'y0spec':y0_spec,
            'Nsize': (Nspec, Nnutr, Ntox, Ntubes),
            'size': (Ntubes, Nspec, Nnutr, Ntox),
            'Ntransf':Ntransfers,
            'Ntsf_stability':Ntsf_stability,
            'bneckdil':(dilution, frac_bneck),
            'mixrate':(mixrate, rate_migr_poisson),
            'RNG':rng
            }

## NOTE: 'ftup' might need to be changed
computepars_ODE = {'ftup':casepars.fgrowth,
                   'timeframe': (Ntstep, t0, tstop),
                   'mutpars': (mutrate, mutdev, mutfac)}

## In case, run all communities
runstability = casepars.runstability
runall = casepars.runall
if runall:
    #### Step 3: DRAW COMMUNITIES
    ##   3a: Define how many co-culture communities there are in total
    Ncombs = 0
    for n in range(1,Nspecpool+1):
        Ncombs += int(scipy.special.comb(N=Nspecpool, k=n, exact=False))

    ## Pre-allocate the number of elements in the community list
    commlist_all = [tuple(),]*Ncombs

    #### Draw and put the species combinations in a list
    i0 = 0
    for n in range(1,Nspecpool+1):
        listtmp   = list(itertools.combinations(list(range(Nspecpool)),n))
        Ncombstmp = len(listtmp)
        commlist_all[i0:(i0+Ncombstmp)] = copy.deepcopy(listtmp)
        i0 = i0+Ncombstmp

    #### Step 4: FOR EACH REPEAT, SHUFFLE
    str_repprint = 'Starting \'all methods\' ... '
    with open(path_outdir / 'Outfile.out', 'a') as f:
        f.write(str_repprint+'\n')

    ## 4.1 Shuffle species pool into Ncomms x Nspecpercomm
    commpars_all = selfuns.assemble_allcomms_fromlong(commlist_all, pars_specpool, (Nnutr, Ntox))

    ## 4.2 Prepare for and run each selection method
    #### Step 1a: PREPARE DIRECTORIES for OUTPUT
    path_allout  = path_outdir / 'Allcomms'
    Path(path_allout).mkdir(parents=True, exist_ok=True)

    #### Step 5a: PREPARE OUTPUT FILES
    fname_comms = 'Community_quantities_allcomms_seed_'+str(rg_seed)+'_'+timestr
    fname_specs = 'Species_quantities_allcomms_seed_'+str(rg_seed)+'_'+timestr

    ## And the corresponding paths
    path_comms = path_allout/fname_comms
    path_specs = path_allout/fname_specs

    ### Make files, print headers
    fnametup = (path_comms, path_specs)

    #### Step 5b: PRINT HEADERS
    for fname, headers in zip(fnametup,headtups_allcomms):
        with open(fname,'w') as f:
            f.write(''.join(f'{el:15s}' for el in headers)+'\n')

    ##### Step 5c: Define filenames to print output
    printpars = {'fnames':{'fname_comms': path_comms,
                           'fname_specs': path_specs}}

    ###### Step 6: RUN, FORREST! RUN!
    y0tup = (y0_spec, y0_nutr, y0_tox)
    repfuns.compute_oneround_allcoms(commpars_all,y0tup,printpars,computepars_ODE)


#### PREPARE TO RUN ALL METHODS
printdict = {'paramheaders':headertups,
             'basedir':path_outdir,
             'datetup':(datestr, timestr, 0)}

#### Step 3: FOR EACH REPEAT, SHUFFLE
list_fmethods = [selfuns.compute_dilutions_lineages,
                 selfuns.compute_dilutions_topN_nosel,
                 selfuns.compute_dilutions_propagule,
                 selfuns.compute_dilutions_Dis_nosel,
		         selfuns.compute_dilutions_disassembly,
                 selfuns.compute_propagule_with_immigration,
                 selfuns.compute_propagule_withim_nosel
                 ]
list_methodnames = ['Lineages',
                    'PropNosel',
                    'Propagule',
                    'DisNosel',
                    'Disassembly',
                    'PropWithIm',
                    'PwINosel'
                    ]
for rep in range(Nreps):
    #print('Starting selection methods, repeat '+str(rep)+' ... ')
    str_repprint = 'Starting selection methods, repeat '+str(rep)
    with open(path_outdir / 'Outfile.out', 'a') as f:
        f.write(str_repprint+'\n')
    printdict['datetup'] = (datestr, timestr, rep)

    # 3.0 Initialise the fossil record
    record = selfuns.init_fossilrecord( pars_specpool, y0_spec )

    ## 3.1 Shuffle species pool into Ncomms x Nspecpercomm
    commlist, pars_comm = selfuns.assemble_CRcomms_fromlong(rng,pars_specpool, exppars['size'])
    # commlist, pars_comm = selfuns.assemble_init_comms_fromlong(rng,pars_specpool, exppars['size'])
    pars_comm['nHill']  = growthpars.cHill
    pars_comm['fpow']   = growthpars.fpow

    ## 3.2 Prepare for and run each selection method
    # HOWTO RUN SELECTION METHODS IN PARALLEL
    for f_method,str_method in zip(list_fmethods,list_methodnames):
        ## Copy common parameters to avoid modifying in-place
        record_method = copy.deepcopy(record)
        pars_method  = copy.deepcopy(pars_comm) # Mutations occur in place
        ## 3.21 Run each selection method
        casedict = {'methods': (f_method, str_method),
                    'cases': (casepars.fgrowth, 'Full')}
        pars_out, y_outs = repfuns.AST_Run_Seltransfer(pars_method,record_method,printdict, casedict, exppars, computepars_ODE)

        if runstability:
            #### Rerun no-selection case on selected communities
            ## Step 1: Update print function with new headers
            headertups_stability = (headtup_pops, headtup_AUCs,
                                    headtup_nutrs, headtup_toxs)
            print_stability = {'paramheaders':headertups_stability,
                               'basedir':path_outdir,
                               'datetup':(datestr, timestr, rep)}
            ## Step 2: RUN
            repfuns.AST_Run_Stability(pars_out, y_outs, print_stability, casedict, exppars, computepars_ODE)
