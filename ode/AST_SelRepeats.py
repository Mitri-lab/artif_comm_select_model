
##### SCRIPT TO SIMULATE MULTIPLE SELECTION RUNS WITH TOP-N (aka Propagule) SELECTION METHOD
## Using
##   selsims_coexmut_hilldegr.py
## Written by Bjorn.Vessman@unil.ch
## 2020.04.16
###

import sys
import copy
import numpy as np
import seltheory_functions as selfuns
from pathlib import Path, PurePath

###################################################################
###################################################################

def compute_oneround_allcoms(list_commpars,y0tup,printpars,computepars):
    #### One round of growth of all n-species communities

    #### Read experimental parameters
    y0, y0_nutr, y0_tox = y0tup
    Nspec1, Nnutr = list_commpars[0]['r'].shape
    Nspec1, Ntox  = list_commpars[0]['m'].shape

    #### Output files organised in dictionary
    fnames_dict = printpars['fnames']
    fn_comms = fnames_dict['fname_comms']
    fn_specs = fnames_dict['fname_specs']

    ### Arguments for tube community growth in step 1
    Ntstep, t0, T = computepars['timeframe']
    ftup = computepars['ftup']

    #### Compute each communities
    for ic, commpars in enumerate(list_commpars):
        IDs_comm = commpars['specid']
        Nspec = len(IDs_comm)

        #### Step 0: Prepare parameters
        y0tmp = np.hstack((np.repeat(y0, Nspec),
                           np.repeat(y0_nutr, Nnutr),
                           np.repeat(y0_tox, Ntox)))
        comppars_co = (Ntstep, t0, T, y0tmp)

        #### Grow populations and compute degradation
        tco, youts = selfuns.computeODE(ftup, comppars_co, commpars)
        yout_spec  = youts[-1,:Nspec]
        yout_nutr  = youts[-1,Nspec:(Nspec+Nnutr)]
        yout_tox   = youts[-1,(Nspec+Nnutr):]

        #### Compute quantities to print
        I_AUC = np.trapz(np.subtract(youts[:,:Nspec],y0), dx=(T-t0)/(Ntstep+1), axis=0)
        Poptot = np.sum(yout_spec)
        AUCtot = np.sum(I_AUC)
        relpop = np.divide(np.add(yout_spec,1.0),np.sum(yout_spec)+Nspec)
        H = -1.0*np.sum(np.multiply(relpop,np.log(relpop)))
        divShannon = np.exp(H)

        #### Calculate norm of toxin concentrations
        Trel  = np.power(np.divide(yout_tox, y0_tox),2.0)
        Tmean = np.divide(np.sum(Trel), Ntox)
        Tdegr = np.subtract(1.0,np.sqrt(Tmean))

        #### Calculate norm of nutrient concentrations
        Nrel  = np.power(np.divide(yout_nutr, y0_nutr),2.0)
        Nmean = np.divide(np.sum(Nrel), Nnutr)
        Ndegr = np.subtract(1.0,np.sqrt(Nmean))

        #### Print 1: community pars
        # Output file per comm: Comm, Tdegr, Ndegr, Diversity, AUC_tot, Pop_tot
        with open(fn_comms,'a') as f:
            f.write(''.join(['{0:12d}'.format(i) for i in [ic,]])+'\t')
            arr_outputs = [Tdegr, Ndegr, divShannon, AUCtot, Poptot]
            # print(arr_outputs)
            f.write(''.join(['{0:12.2e}'.format(a) for a in arr_outputs])+'\n')

        #### Print 2: species pars
        # Output file per species: Comm, Spec, Popsize, AUC, sum_f
        with open(fn_specs,'a') as f:
            if Nspec<2:
                sumf = np.sum(commpars['finv'])
                f.write(''.join(['{0:12d}'.format(i) for i in [ic,IDs_comm[0]]])+'\t')
                arr_outputs = [yout_spec[0],I_AUC[0], sumf]
                f.write(''.join(['{0:12.2e}'.format(a) for a in arr_outputs])+'\n')
            else:
                for idxs, (ID,y_end,AUC_ps,f_perspec) in enumerate(zip(IDs_comm,yout_spec,I_AUC,commpars['finv'])):
                    sumf = np.sum(f_perspec)
                    f.write(''.join(['{0:12d}'.format(i) for i in [ic,ID]])+'\t')
                    arr_outputs = [y_end, AUC_ps, sumf]
                    # print(arr_outputs)
                    f.write(''.join(['{0:12.2e}'.format(a) for a in arr_outputs])+'\n')

    return()

#### Growth function to use for the ecological interactions
def AST_Run_Seltransfer(commpars,record, printdict, casedict, exppars, computepars_ODE):
    #### Function version of Simsel_Dis_202103.py
    # INPUT:
    # commpars
    # record
    # printdict
    # casedict
    # exppars
    # computepars_ODE
    #
    # OUTPUT: None

    #### STEP 0: Read input
    headertups   = printdict['paramheaders']
    path_basedir = printdict['basedir']
    datestr, timestr, rep = printdict['datetup']
    (Ntubes, Nspec, Nnutr, Ntox) = exppars['size']
    y0_spec = exppars['y0spec']

    f_selmethod, str_selmethod = casedict['methods']
    f_growth, str_cases = casedict['cases']

    ## Derived from input
    Nreg = Nnutr +Ntox
    Nplaces = 2*Nreg


    #### Step 1a: PREPARE DIRECTORIES for OUTPUT
    # basedir contains a base directory, date and size of community
    # Add directories for selection method (Disassembly, Propagule ... )
    # path_basedir = PurePath(basedir)
    path_outdir  = path_basedir / str_selmethod / str_cases
    Path(path_outdir).mkdir(parents=True, exist_ok=True)

    #### Step 1b: PREPARE OUTPUT FILES
    fname_states_pop  = 'States_popul_'+timestr+'_rep_'+str(rep)
    fname_states_AUC  = 'States_AUC_'+timestr+'_rep_'+str(rep)
    fname_states_nutr = 'States_nutrient_'+timestr+'_rep_'+str(rep)
    fname_states_tox  = 'States_toxin_'+timestr+'_rep_'+str(rep)
    fname_invest = 'Invest_'+timestr+'_rep_'+str(rep)
    fname_specID = 'SpeciesID_'+timestr+'_rep_'+str(rep)
    ## In case we mutate growth, mortality, yield and/or degradation
    # fname_growth = 'Growth_'+timestr+'_rep_'+str(rep)

    ## And the corresponding paths
    path_states_pop  = path_outdir/fname_states_pop
    path_states_AUC  = path_outdir/fname_states_AUC
    path_states_nutr = path_outdir/fname_states_nutr
    path_states_tox  = path_outdir/fname_states_tox
    path_invest  = path_outdir/fname_invest
    path_specID  = path_outdir/fname_specID
    # path_growth = path_outdir/fname_growth

    ### Make files, print headers
    fnametup = (path_states_pop, path_states_AUC,
                path_states_nutr, path_states_tox,
                path_invest,
                # path_growth,
                path_specID)

    #### Step 1c: PRINT HEADERS
    for fname, headers in zip(fnametup,headertups):
        with open(fname,'w') as f:
            f.write(''.join(f'{el:15s}' for el in headers)+'\n')


    ##### Step 2a: Define initial state
    y0_pertub = np.hstack(([y0_spec]*Nspec, np.zeros(Nplaces-Nspec)))
    y0tubes   = np.reshape(np.tile(y0_pertub,Ntubes),(Ntubes,Nplaces))

    ##### Step 2b: Define filenames to print output
    printpars = {'fnames':{'fname_states': (path_states_pop,path_states_nutr,path_states_tox),
                           'fname_AUC': path_states_AUC,
                           'fname_mutpars': path_invest,
                           #'fname_fixpars': (path_growth,),
                           'fname_specID': path_specID}}

    ##### Step 2c: Define parameters for ODE
    # Eventually, we might want to re-run for several growth functions
    # computepars_ODE['ftup'] = f_growth

    ###### Step 3: RUN, FORREST! RUN!
    retpars, yreturn = f_selmethod(commpars,record,y0tubes,exppars,printpars,computepars_ODE)
    ## f_selmethod is alias for e.g.
    #  selfuns.compute_dilutions_Dis_nosel(commpars,y0tubes,exppars,printpars,computepars_ODE)
    return( retpars, yreturn )


#### Repeat the no-selection case to estimate stability of selected communities
def AST_Run_Stability(commpars, y0tubes, printdict, casedict, exppars, computepars):
    ####
    # INPUT:
    # commpars     dict of community growth parameters
    #                passed from last round of selection
    # y0tubes      array [Ntubes x Nplaces] of population sizes
    #                passed from last round of selection
    # printdict    dict: headers for output files, path to base directory, date
    # casedict     dict: growth function
    # exppars
    # computepars
    #
    # OUTPUT: None

    #### STEP 0: Read input
    headertups   = printdict['paramheaders']
    path_basedir = printdict['basedir']
    datestr, timestr, rep = printdict['datetup']
    (Ntubes, Nspec, Nnutr, Ntox) = exppars['size']
    # y0_spec = exppars['y0spec']

    f_selmethod, str_selmethod = casedict['methods']
    f_growth, str_cases = casedict['cases']

    ## Derived from input
    Nreg = Nnutr +Ntox
    Nplaces = 2*Nreg

    #### Step 1a: PREPARE DIRECTORIES for OUTPUT
    # basedir contains a base directory, date and size of community
    # Add directories for selection method (Disassembly, Propagule ... )
    # path_basedir = PurePath(basedir)
    path_outdir  = path_basedir / str_selmethod / str_cases / 'Stability'
    Path(path_outdir).mkdir(parents=True, exist_ok=True)

    #### Step 1b: PREPARE OUTPUT FILES
    fname_states_pop  = 'Stability_popul_'+timestr+'_rep_'+str(rep)
    fname_states_AUC  = 'Stability_AUC_'+timestr+'_rep_'+str(rep)
    fname_states_nutr = 'Stability_nutrient_'+timestr+'_rep_'+str(rep)
    fname_states_tox  = 'Stability_toxin_'+timestr+'_rep_'+str(rep)
    ## ID, finv less relevant
    # fname_specID = 'SpeciesID_'+timestr+'_rep_'+str(rep)
    # fname_invest      = 'Stability_'+timestr+'_rep_'+str(rep)

    ## And the corresponding paths
    path_states_pop  = path_outdir/fname_states_pop
    path_states_AUC  = path_outdir/fname_states_AUC
    path_states_nutr = path_outdir/fname_states_nutr
    path_states_tox  = path_outdir/fname_states_tox
    # path_invest   = path_outdir/fname_invest
    # path_specID = path_outdir/fname_specID

    ### Make files, print headers
    fnametup = (path_states_pop, path_states_AUC,
                path_states_nutr, path_states_tox,
                # path_invest, path_specID
                )

    #### Step 1c: PRINT HEADERS
    for fname, headers in zip(fnametup,headertups):
        with open(fname,'w') as f:
            f.write(''.join(f'{el:15s}' for el in headers)+'\n')

    ##### Step 2a: Initial state
    # y0tubes is now of size [Ntubes x Nplaces], saved from rounds

    ##### Step 2b: Define filenames to print output
    printpars = {'fnames':{'fname_states': (path_states_pop,path_states_nutr,path_states_tox),
                           'fname_AUC': path_states_AUC,
                           #'fname_mutpars': path_invest,
                           #'fname_fixpars': (path_growth,),
                           #'fname_specID': path_specID
                           }}

    ###### Step 3: RUN, FORREST! RUN!
    # selfuns.compute_dilutions_lineages(commpars,record,y0tubes,exppars,printpars,computepars)
    # Make sure that in the calling script:
    #   exppars is updated so that Ntransfers is correct
    #   y0tubes, commpars is passed from the last round of selection
    #   headertups is only N, T, S and AUC
    selfuns.compute_dilutions_stability(commpars,y0tubes,exppars,printpars,computepars)
    return( )
