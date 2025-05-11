## FUNCTIONS TO PRINT STUFF FOR THE ARTIFICIAL SELECTION THEORY (AST) PROJECT
##
## Written by Bjorn Vessman (UNIL DMF, bjorn.vessman@unil.ch)
## From the 15th of April 2021 and onwards
## Log
##  2021-04-15   First edition. Basic print_specpool()


def print_specpool(outpath, timestr, pars):
    # Having drawn a pool of species (a dict of parameter arrays for ODE model)
    # we print the parameters to files for future purposes.
    # INPUT
    # outdir      A pathlib.PurePath to the directory where

    ###### HEADERS for PARAMETER FILES
    headtup_main = ('Species', )
    headtup_pars_mats = headtup_main+('Substr','Value') # Growth, death, investment
                                                        # r, m, f are matrices
                                                        # of size [Nspec x Nsubstr]
    headtup_pars_vec  = headtup_main+('Value',) # Yield and degradation
                                                # Y, d are vectors
                                                # of length [Nspec]
    ## Collect all headers
    headertups = (headtup_pars_mats, headtup_pars_mats, headtup_pars_mats,
                  headtup_pars_vec, headtup_pars_vec)


    ###### PREPARE PARAMETER FILES and PRINT HEADERS
    fnout_growth = 'Growth_'+timestr
    fnout_mort   = 'Mortality_'+timestr
    fnout_invest = 'Invest_'+timestr
    fnout_yield  = 'Yield_'+timestr
    fnout_degr   = 'Degr_'+timestr

    path_growth = outpath / fnout_growth
    path_mort   = outpath / fnout_mort
    path_invest = outpath / fnout_invest
    path_yield  = outpath / fnout_yield
    path_degr   = outpath / fnout_degr

    ### Make files, print headers
    fnametup = (path_invest, path_growth, path_mort,
                path_yield, path_degr)
    for fname, headers in zip(fnametup,headertups):
        with open(fname,'w') as f:
            f.write(''.join(f'{el:12s}' for el in headers)+'\n')

    fntime = outpath / 'Timestamp'
    with open(fntime,'a') as f:
        f.write(''.join(f'{timestr:16s}')+'\n')
        # f.write(''.join('{0:16s}'.format(timestr))+'\n')

    #### Step 3c: PRINT PARAMETERS
    Nspecpool = len(pars['Y'])
    ### Print matrix parameters for investment, growth, mortality
    fntup = (path_invest, path_growth, path_mort)
    for trait, fntrt in zip(('finv', 'r', 'm'),fntup):
        with open(fntrt,'a') as f:
            for spec in range(Nspecpool):
                traitarr = pars[trait][spec,:]
                for tr,val in enumerate(traitarr):
                    f.write(''.join(f'{el:12d}' for el in [spec,tr])+'\t')
                    f.write(''.join(f'{val:12.2e}')+'\n')
                    #f.write(''.join('{0:12d}'.format(el) for el in [spec,tr])+'\t')
                    #f.write(''.join('{0:12.2e}'.format(val))+'\n')

    ### Print fixed parameters for yield and degradation
    for trait, fntrt in zip(('Y','d'),(path_yield, path_degr)):
        with open(fntrt,'a') as f:
            for spec in range(Nspecpool):
                trval = pars[trait][spec]
                f.write(''.join(f'{spec:12d}')+'\t')
                f.write(''.join(f'{trval:12.2e}')+'\n')
                #f.write(''.join('{0:12d}'.format(spec))+'\t')
                #f.write(''.join('{0:12.2e}'.format(trval))+'\n')

    return( )
