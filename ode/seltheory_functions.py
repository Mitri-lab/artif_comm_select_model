import sys
import copy
import scipy
import random
# import pandas as pd
import numpy as np
import scipy.integrate as integr
# from numpy import linalg as LA
from datetime import date, datetime

"""
Created on Mon 27.12.2021 by bjorn (bjorn.vessman@unil.ch, bvessman@gmail.com)

Functions to simulate the selection methods (propagule, disassembly and no-selection) on consumer-resource models.
The functions have been ported from previous iterations and reflect the current (Dec 2021) implementation and understanding of the methods.

"""

##########################################################################
##########################################################################

def draw_CRparams_longform(rng,exppars,grpars):
    #### Function to initialise a set of tubes with communities of certain size
    ##   Long-form version that outputs a species pool as a list. The species can
    ##    then be combined in different communities
    ### INPUT ###
    ## exppars    Experimental parameters regarding
    ##             size: number of tubes, species,  nutrients and toxins
    ##             initial conditions: inoculum size, substrate concentrations and flow rate
    ## grpars     A set of pre-defined parameter values for growth, death etc.

    ### OUTPUT ###
    ## pars       Growth parameters of each species

    #### Collect experimental parameters
    Nspecpool, Nnutr, Ntox = exppars['size']
    flag_sparse, frac_sparse = exppars['sparse']
    spflag_r, spflag_m, spflag_f = flag_sparse
    spfrac_r, spfrac_m, spfrac_f = frac_sparse

    #### Basic parameters
    r0low, r0high = grpars['r0'] # Basic growth rate
    m0low, m0high = grpars['m0'] # Basic death rate

    Ymean, Ydev = grpars['Y0'] # Drawn log-normally
    dmean, ddev = grpars['d0'] # Drawn log-normally

    #### Sparse case, where species only degrade some nutrients
    if spflag_r:
        spmat_r = rng.binomial(n=1,p=1.0-spfrac_r,size=(Nspecpool,Nnutr))

    if spflag_m:
        spmat_m = rng.binomial(n=1,p=1.0-spfrac_m,size=(Nspecpool,Ntox))

    if spflag_f:
        spmat_f = rng.binomial(n=1,p=1.0-spfrac_f,size=(Nspecpool,Ntox))

    for idxc,(spr, spm, spf) in enumerate(zip(spmat_r,spmat_m,spmat_f)):
        if (spr==0).all():
            spmat_r[idxc,rng.integers(low=0,high=Nnutr)] = 1
        if (spm==0).all():
            spmat_m[idxc,rng.integers(low=0,high=Ntox)] = 1
        if (spf==0).all():
            spmat_f[idxc,rng.integers(low=0,high=Ntox)] = 1

    ####### Draw community parameters at random, by the above basic parameters
    Yvec = rng.lognormal(mean=Ymean, sigma=Ydev, size=Nspecpool)
    dvec = rng.lognormal(mean=dmean, sigma=ddev, size=Nspecpool)
    #Ylogs = rng.uniform(low=Ymean, high=Ydev, size=Nspecpool)
    #dlogs = rng.uniform(low=dmean, high=ddev, size=Nspecpool)
    #Yvec = np.power(10.0,Ylogs)
    #dvec = np.power(10.0,dlogs)

    #### Draw and sparsify
    rmat = rng.uniform(low=r0low, high=r0high, size=(Nspecpool,Nnutr))
    mmat = rng.uniform(low=m0low, high=m0high, size=(Nspecpool,Ntox))
    fmat = rng.uniform(low=0.0, high=1.0, size=(Nspecpool,Ntox))
    rmat *= spmat_r
    mmat *= spmat_m
    fmat *= spmat_f

    #### Re-scale the investment 'f' to sum to d~Uni(0,1)
    for idxs,ftmp in enumerate(fmat):
        fmat[idxs,:] /= np.sum( ftmp )
        fmat[idxs,:] *= rng.uniform(low=0,high=1,size=1)

    ### Draw substrate parameters
    if type(grpars['KN0']) is tuple:
        KN0low, KN0hi = grpars['KN0']
        KT0low, KT0hi = grpars['KT0']
        KNvec = rng.uniform(low=KN0low, high=KN0hi, size=Nspecpool)
        KTvec = rng.uniform(low=KT0low, high=KT0hi, size=Nspecpool)
    else:
        KTvec = grpars['KT0']
        KNvec = grpars['KN0']

    pars = {'finv': fmat,
            'm': mmat,
            'r': rmat,
            'KN': KNvec,
            'KT': KTvec,
            'Y': Yvec,
            'd': dvec}

    return( pars )


def assemble_allcomms_fromlong(commlist, grpars_long, sizetup):
    ### Functions to assemble communities for CR model
    ##  INPUT
    #     commlist       list() list of communities to assemble. In case of 'all' communities, these have
    #     grpars_long    dict() of growth parameters, structure as output of draw_CRparams_longform()
    #     sizetup        tuple of community sizes
    ##  OUTPUT
    #     pars           dict() of parameters, each of shape [Ncomms, Nspec [, Ntrait]]

    #### Collect experimental parameters
    Nnutr, Ntox = sizetup
    Ncommstot = len(commlist)

    ### Pre-allocate list of growth parameters
    dict_allcomms = [dict(),]*Ncommstot
    #mlist = [np.array(),]*Ncommstot
    #flist = [np.array(),]*Ncommstot
    #Yvec  = [np.array(),]*Ncommstot
    #dvec  = [np.array(),]*Ncommstot

    ### Given species indices, collect parameters
    for idxcomm, ids_specs in enumerate(commlist):
        Nspec = len(ids_specs)
        rmattmp = np.zeros((Nspec,Nnutr))
        mmattmp = np.zeros((Nspec,Ntox))
        fmattmp = np.zeros((Nspec,Ntox))
        Yvectmp = np.zeros(Nspec)
        dvectmp = np.zeros(Nspec)
        for idx_in_comm, specID in enumerate(ids_specs):
            ## idxs_comm is the position in the new community:  idxs_comm = 0, 1, 2, ...
            ## specID is the species identity: specID = 0,1,2,3,4, ...
            rmattmp[idx_in_comm,:] = grpars_long['r'][specID,:]
            mmattmp[idx_in_comm,:] = grpars_long['m'][specID,:]
            fmattmp[idx_in_comm,:] = grpars_long['finv'][specID,:]
            Yvectmp[idx_in_comm]   = grpars_long['Y'][specID]
            dvectmp[idx_in_comm]   = grpars_long['d'][specID]

        pars_thiscomm = {'size':(Nspec,Nnutr,Ntox),
                         'specid':list(ids_specs),
                         'finv':fmattmp,
                         'm':mmattmp,
                         'r':rmattmp,
                         'KN': grpars_long['KN'],
                         'KT': grpars_long['KT'],
                         'Y':Yvectmp,
                         'd':dvectmp,
                         'fpow':grpars_long['fpow'],
                         'nHill':grpars_long['nHill'],
                         }
        dict_allcomms[idxcomm] = copy.deepcopy(pars_thiscomm)

    return( dict_allcomms )

def assemble_CRcomms_fromlong(rng,grpars_long, sizetup):
    ### Functions to assemble communities for CR model
    ##  INPUT
    #     grpars_long    dict() of growth parameters, structure as output of draw_CRparams_longform()
    #     sizetup        tuple of community sizes
    ##  OUTPUT
    #     pars           dict() of parameters, each of shape [Ncomms, Nspec [, Ntrait]]

    #### Collect experimental parameters
    Ncomms, Nspec, Nnutr, Ntox = sizetup

    Nreg    = Nnutr+Ntox
    Nplaces = 2*Nreg
    Nmut    = Nplaces-Nspec

    ### Pre-allocate matrices for each parameters
    rmat = np.zeros((Ncomms,Nplaces,Nnutr))
    mmat = np.zeros((Ncomms,Nplaces,Ntox))
    fmat = np.zeros((Ncomms,Nplaces,Ntox))

    IDmat = np.multiply(-1.0,np.ones((Ncomms,Nplaces)))
    Yvec = np.zeros((Ncomms,Nplaces))
    dvec = np.zeros((Ncomms,Nplaces))

    ### Shuffle species indices to randomize community assembly
    if 'specid' in grpars_long.keys():
        speclist = copy.deepcopy(grpars_long['specid'])
        Nspecpool = len(speclist)
    else:
        Nspecpool = len(grpars_long['r'])
        speclist = list(range(Nspecpool))
    #rng.shuffle(speclist)

    # Assign species to communities, at random, without repeating species in comms
    idx_tip = Nspecpool%Ncomms # Tipping point comm index, below which we put Nspecmax_grtd+1 guaranteed
    Nspecmax_grtd = Nspecpool//Ncomms+1
    commlist = [[-2]*Nspec]*Ncomms
    for ic in range(Ncomms):
        specs_tmp = [-1]*Nspec
        if Nspecpool<Ncomms:
            ### In most realistic cases, we have less species than the number of comms
            if ic<Nspecpool:
                spec_tointroduce = speclist[ic]
                specs_tmp[0] = spec_tointroduce

                specs_notpres = [s for s in speclist if s!=spec_tointroduce]
                specs_tmp[1:] = rng.choice(specs_notpres,size=Nspec-1,replace=False)
            else:
                specs_tmp = list(rng.choice(speclist,size=Nspec,replace=False))
        else:
            ### In case we have _more_ species than the number of comms
            if ic<idx_tip:
                ## Up to this index, we put 'Nspecmax_grtd' species per community
                idxs = ic*Nspecmax_grtd
                idxe = (ic+1)*Nspecmax_grtd
            else:
                ## From there, we put Nspecmax_grtd-1
                #idxs = idx_tip*Nspecmax_grtd+(ic-idx_tip)*(Nspecmax_grtd-1)
                idxs = ic*(Nspecmax_grtd-1)+idx_tip
                idxe = idxs+Nspecmax_grtd-1

            # Calculate number of non-guaranteed spots in this case?
            Nguaranteed = idxe-idxs
            Nnonguarant = Nspec-Nguaranteed

            # Introduce guaranteed and nonguaranteed
            spec_tointroduce = speclist[idxs:idxe]
            specs_tmp[:Nguaranteed] = spec_tointroduce

            specs_notpres = list(set(speclist)-set(spec_tointroduce))
            specs_tmp[Nguaranteed:] = list(rng.choice(specs_notpres,size=Nnonguarant,replace=False))

        commlist[ic] = specs_tmp

    ### Given species indices, collect parameters
    for idxcomm, ids_specs in enumerate(commlist):
        for idx_in_comm, specID in enumerate(ids_specs):
            ## idxs_comm is the position in the new community:  idxs_comm = 0, 1, 2, ...
            ## specID is the species identity: specID = 0,1,2,3,4, ...

            rmattmp = grpars_long['r'][specID,:]
            mmattmp = grpars_long['m'][specID,:]
            fmattmp = grpars_long['finv'][specID,:]

            rmat[idxcomm,idx_in_comm,:] = rmattmp
            mmat[idxcomm,idx_in_comm,:] = mmattmp
            fmat[idxcomm,idx_in_comm,:] = fmattmp

            IDmat[idxcomm,idx_in_comm] = speclist[specID]
            Yvec[idxcomm,idx_in_comm]  = grpars_long['Y'][specID]
            dvec[idxcomm,idx_in_comm]  = grpars_long['d'][specID]

    if type(grpars_long['KN']) is np.ndarray:
        KN = np.zeros((Ncomms,Nplaces))
        KT = np.zeros((Ncomms,Nplaces))
        for idxcomm, ids_specs in enumerate(commlist):
            for idx_in_comm, specID in enumerate(ids_specs):
                KN[idxcomm,idx_in_comm] = grpars_long['KN'][specID]
                KT[idxcomm,idx_in_comm] = grpars_long['KT'][specID]
    else:
        KN = grpars_long['KN']
        KT = grpars_long['KT']

    pars = {'size':(Nplaces,Nnutr,Ntox),
            'specid':IDmat,
            'finv':fmat,
            'm':mmat,
            'r':rmat,
            'KN': KN,
            'KT': KT,
            #'nHill': grpars_long['nHill'],
            'Y':Yvec,
            'd':dvec,
            #'fpow':grpars_long['fpow']
           }

    return( commlist, pars )


def init_fossilrecord( pars_longform, y0 ):
    #### INPUT
    # pars_longform   dict with keys 'd','Y','r','m','finv'
    #     the corresponding values are the parameter values, of sizes
    #     [d] = Nspec_pool
    #     [Y] = Nspec_pool
    #     [r] = Nspec_pool x Nnutr
    #     [m] = Nspec_pool x Ntox
    #     [finv] = Nspec_pool x Ntox
    # y0              float or array of initial population size
    #### OUTPUT
    # dict_record    a dict() of length Nspecpool (the number of species in the initial pool)
    #                The dictionary is keyed by species id. Each entry is in itself
    #                a dictionary of growth parameters r, m, f, Y, d in addition to
    #                metadata: ID of species, isolated from which tube in which round
    #                and the number of strains at present (useful within simulations)

    ### Collect problem size
    Nspec_pool, Nnutr = pars_longform['r'].shape
    Nspec_pool, Ntox  = pars_longform['m'].shape

    ### Define the per-species record
    # Record what species, where and when. How many strains and the (initial) pop size
    # 'specid' is by default the index 0,1,2..., BUT
    #  still useful to have a 'specid' parameter if species are denoted 1, 114, 1w ...
    record_perspec = {'specid':-1,
                      'last_tube':0,
                      'round':0,
                      'Nstrains':1,
                      'y0':y0}

    for key in pars_longform.keys():
        values = pars_longform[key]
        if type(values) is np.ndarray: # d, Y, r, m, f
            if values.ndim>1:
                # r, m, f are matrices of size Nspec_pool x Ntrait
                record_perspec[key] = np.zeros(len(values[0]))
            else:
                # Y, d are arrays of length Nspec_pool
                record_perspec[key] = 0.0
        #else:
        # KN, KT are scalars. Not needed in the record at present
        # record_perspec[key] = 0.0
        # continue

    ### Repeat the per-species record and initialise with values
    dict_record = dict()
    for idxs in range(Nspec_pool):
        dict_tmp = copy.deepcopy(record_perspec)
        dict_tmp['specid'] = idxs
        # dict_tmp['Nstrains'] = 1
        # dict_tmp['last_tube'] = 0
        # dict_tmp['transfer'] = 0

        for key in pars_longform.keys():
            values = pars_longform[key]
            if type(values) is np.ndarray: # d, Y, r, m, f
                dict_tmp[key] = values[idxs]
        dict_record[idxs] = dict_tmp

    return(dict_record)


###########################################################################
############### Functions for growth and mutations
###########################################################################


def computeODE(fun, comppars, funpars, jac=None):
    # f, jac = fjac
    x      = integr.ode(fun,jac).set_integrator('dopri5')

    N, t0, T, y0 = comppars
    Nstates = len(y0)
    xout  = np.zeros([N+1,Nstates])
    xout[0,:] = y0

    x.set_initial_value(y0,t0).set_f_params(funpars).set_jac_params(funpars)
    dt = (T-t0)/N
    t  = np.linspace(t0,T,N+1)
    idx_temp = 0

    while x.successful() and idx_temp<N:
        idx_temp += 1
        x_temp = x.integrate(x.t+dt)
        xout[idx_temp] = x_temp

    return( t, xout )


def mutate_fixf_random(rng,specpars,idxtup,mutargs):
    #### For a given tube of index 'itb', mutate a trait 'itr'
    ### INPUT
    ## specpars = {'ftup','mtup','rtup' ... },
    ##             each trait is a list of length Ntubes,
    ##             of n-dimensional arrays of size (Nspec+1)x(trait_len)
    ## idxtup   = (i_tube, i_spec, i_tr)
    ##             Indices marking tube, species and trait to mutate
    ## mutargs  = (muttr, mutfac, idxmut)
    ##             Trait, which deviation for drawIndices marking which tube, species and trait to mutate

    itb, isp, itrs = idxtup
    mutdev, mutfac, idx_mutant = mutargs

    ### Change trait value --- Note that 'finv' is only non-cheater trait to change
    fmat = specpars['finv']
    Ntubes, Nspec, Ntox = fmat.shape

    ## Mutations modify 'f' by a value from a lognormal(0, mutdev) distribution
    trcp = copy.deepcopy( fmat[itb,isp,:] )
    # mutmult = 0.5+np.random.beta(mutfac,mutfac) # +0.5 makes 'mutmult' centered on 1.0

    mutmult = rng.lognormal(mean=0.0, sigma=mutdev, size=len(itrs))
    for itr, mm in zip(itrs,mutmult):
        newval = np.multiply(trcp[itr], mm)
        if newval<0:  # This should not be possible, but just in case
            newval = 0.0
        
        ### Introduce the modified trait
        trcp[itr] = newval
    
    #### Old version of the code 
    ##   a species can gain functions with a small probability
    ##   i.e. if trcp[itr]==0, then sample a new non-zero value here
    #if trcp[itr]>0.0: # In case the strain is already able to degrade
    #    newval = np.multiply(trcp[itr], mutmult)
    #    if newval<0:  # This should not be possible, but just in case
    #        newval = 0.0
    #else: # Then trcp[itr]=0.0 and the strain will acquire a new ability
    #    # mutadd = np.abs(np.subtract(mutmult, 1.0)) # Something small between 0 and 1
    #    mutadd = rng.uniform(low=0,high=1,size=1)    # Something small between 0 and 1
    #    newval = np.add(trcp[itr], mutadd)

    ### Check that sum(f^u) <= 1
    fpowu, fpowv = specpars['fpow']
    fpowsum = np.sum(np.power(trcp,fpowu))
    if fpowsum>1.0:
        trcp /= fpowsum

    fmat[itb,idx_mutant,:] = copy.deepcopy( trcp )

    ## Collect all (other) growth params
    for tr in ('r','m'):
        mattmp = copy.deepcopy( specpars[tr] )
        mattmp[itb,idx_mutant,:] = mattmp[itb,isp,:]
        specpars[tr] = copy.deepcopy( mattmp )

    ## Degradation capability and biomass yield
    dtmp = copy.deepcopy( specpars['d'][itb] )
    Ytmp = copy.deepcopy( specpars['Y'][itb] )

    if dtmp.ndim>1:
        dtmp[idx_mutant,:] = dtmp[isp,:]
        specpars['d'][itb,:,:] = dtmp
    else:
        dtmp[idx_mutant] = dtmp[isp]
        specpars['d'][itb,:] = dtmp

    Ytmp[idx_mutant]   = Ytmp[isp]
    specpars['Y'][itb,:] = Ytmp

    ## Finally, assign resident ID to mutant
    matID = specpars['specid']
    matID[itb,idx_mutant] = matID[itb,isp]
    return()


def mutate_tubes_sequential(rng,y0tubes,Nsize,mutpars,pars):
    ##### For each tube, each species ...
    ##    if mutrate > r~Uni(0,1) then mutate that species
    ##
    ####### Input
    ## y0tubes   A list or array, length Ntubes, of initial abundances
    ##             [y0_spec1, ..., y0specN, y0mut1, ... y0mutN]
    ## Nsize     tuple (N_species, N_nutrients, N_toxic, N_tubes)
    ## mutpars   tuple (mutprob, flag_mutrm, mutdev)
    ## pars_r    trait parameters of resident

    Nspec, Nnutr, Ntox, Ntubes = Nsize
    mutprob, mutdev, mutfac = mutpars

    ## Draw mutations at random for each tube ...
    for itb, y_tube in enumerate(y0tubes):
        idxfree = [i for i,y in enumerate(y_tube) if y<1.0] # Find all free spots
        idxsurv = [i for i,y in enumerate(y_tube) if y>0.0] # And the indices of surviving species

        Nfree = len(idxfree)
        Nsurv = len(idxsurv)

        ##### Check that there exists ...
        #  a) a free spot for the mutants and
        #  b) a species to mutate
        if Nfree>0 and Nsurv>0:
            ### There is a fixed number of slots for populations at present
            if Nfree>Nsurv:
                ## We can draw as many mutants as surviving strains
                Nmuts = rng.binomial(n=Nsurv,p=mutprob)
            else:
                ## The number of establishing mutants constrained by the number of free spots
                Nmuts = rng.binomial(n=Nfree,p=mutprob)

            ### If we have drawn any mutations
            if Nmuts>0:
                mutating_strains = rng.choice(idxsurv,size=Nmuts,replace=False)
                for rank_mut, mutidx in enumerate(mutating_strains):
                    idxfr = idxfree[rank_mut]
                    ##  Find one or more non-zero traits to mutate
                    itr_nonzero = [i for i,val in enumerate(pars['finv'][itb,mutidx]) if val>0]
                    itr_chosen  = rng.binomial(n=1,p=1.0/len(itr_nonzero),size=len(itr_nonzero))
                    if sum(itr_chosen)<1:
                        itrs = rng.choice(itr_nonzero,size=1)
                    else:
                        itrs = [itr_nonzero[i] for i,v in enumerate(itr_chosen) if v>0]
                    mutargs = (mutdev, mutfac, idxfr)
                    mutate_fixf_random(rng,pars,(itb,mutidx,itrs),mutargs)
                    y0tubes[itb,idxfr] = 1.0e2

    return()


def mutate_and_compete_fixf(pars,y0,exppars,comppars):
    ##### Given a set of tubes in a particular transfer round
    ###    perform one growth period, competing between variants, and dilute

    y0_nutr, y0_tox = exppars['y0substrate']
    dilution        = exppars['dilution']
    Nsize           = exppars['Nsize']
    N, t0, T = comppars['timeframe']
    ftup     = comppars['ftup']
    mutpars  = comppars['mutpars']
    rng      = comppars['RNG']

    Nspec, Nnutr, Ntox, Ntubes = Nsize
    Nreg  = Nnutr+Ntox
    Nplaces = 2*Nreg
    Nmut  = Nplaces-Nspec

    ### For each tube
    argsco = copy.deepcopy(pars)
    xouts = np.zeros((Ntubes,Nplaces+Nreg))
    I_AUC = np.zeros((Ntubes,Nplaces))
    for tube in range(Ntubes):
        #### Select corresponding species parameters
        argsco['d'] = pars['d'][tube].copy()
        argsco['Y'] = pars['Y'][tube].copy()
        #argsco['KN'] = pars['KN'][tube].copy()
        #argsco['KT'] = pars['KT'][tube].copy()

        #for trait in ('ftup','mtup','rtup'):
        #    fmat, Kf, hf = pars[trait]
        #    argsco[trait] = (fmat[tube].copy(), Kf[tube].copy(), hf[tube].copy())
        for trait in ('finv','m','r'):
            fmat = pars[trait]
            argsco[trait] = fmat[tube].copy()

        #### Find inoculation size in y0, define substrate concentrations
        y0tmp = np.hstack((y0[tube], np.repeat(y0_nutr, Nnutr), np.repeat(y0_tox, Ntox)))

        comppars_co = (N, t0, T, y0tmp)
        # tco, xco = selsim.computeODE(ftup, comppars_co, argsco)
        tco, xco = computeODE(ftup, comppars_co, argsco)
        xouts[tube,:] = xco[-1,:]
        I_AUC[tube,:] = np.trapz(np.subtract(xco[:,:Nplaces],xco[0,:Nplaces]),
                                 dx=(T-t0)/(N+1), axis=0)

        #### Re-dilute from output states
        y0spec = np.divide(xco[-1,:Nplaces],dilution)
        # y0[tube] = [y for y in y0spec if y>=1.0 else 0.0]
        for idxy, y in enumerate(y0spec):
            if y<1.0: # Then remove identity marker
                y0[tube,idxy] = 0.0

                #### Remove species ID and growth parameters
                pars['specid'][tube,idxy] = -1
                pars['Y'][tube,idxy] = 0.0
                for tr,trlen in zip(('r','m','finv'),(Nnutr,Ntox,Ntox)):
                    # rmat, Kr, rn = pars[tr]
                    rmat = pars[tr]
                    rmat[tube,idxy,:] = np.zeros(trlen)

                if pars['d'].ndim>2:
                    pars['d'][tube,idxy,:] = np.zeros(Ntox)
                else:
                    pars['d'][tube,idxy] = 0.0
                #fmat = pars['finv']
                #fmat[tube,idxy,:] = np.zeros(Ntox)
            else:
                y0[tube,idxy] = y

    ### Mutate if there are any free spots
    spaceformutant = y0<1.0
    if spaceformutant.any(): # Shorthand - don't bother mutating if no space for mutants
        # mutate_tubes_onlyf(y0,Nsize,mutpars,pars)
        mutate_tubes_sequential(rng,y0,Nsize,mutpars,pars)

    return((xouts,I_AUC))


def compete_nomutation(pars,y0,exppars,comppars):
    ##### Given a set of tubes in a particular transfer round
    ###    perform one growth period, competing between variants, and dilute

    y0_nutr, y0_tox = exppars['y0substrate']
    dilution        = exppars['dilution']
    Nsize           = exppars['Nsize']
    N, t0, T = comppars['timeframe']
    ftup     = comppars['ftup']

    Nspec, Nnutr, Ntox, Ntubes = Nsize
    Nreg  = Nnutr+Ntox
    Nplaces = 2*Nreg
    Nmut  = Nplaces-Nspec

    ### For each tube
    argsco = copy.deepcopy(pars)
    xouts = np.zeros((Ntubes,Nplaces+Nreg))
    I_AUC = np.zeros((Ntubes,Nplaces))
    for tube in range(Ntubes):
        #### Select corresponding species parameters
        argsco['d'] = pars['d'][tube].copy()
        argsco['Y'] = pars['Y'][tube].copy()
        #argsco['KN'] = pars['KN'][tube].copy()
        #argsco['KT'] = pars['KT'][tube].copy()

        for trait in ('finv','m','r'):
            fmat = pars[trait]
            argsco[trait] = fmat[tube].copy()

        #### Find inoculation size in y0, define substrate concentrations
        y0tmp = np.hstack((y0[tube], np.repeat(y0_nutr, Nnutr), np.repeat(y0_tox, Ntox)))

        comppars_co = (N, t0, T, y0tmp)
        # tco, xco = selsim.computeODE(ftup, comppars_co, argsco)
        tco, xco = computeODE(ftup, comppars_co, argsco)
        xouts[tube,:] = xco[-1,:]
        I_AUC[tube,:] = np.trapz(np.subtract(xco[:,:Nplaces],xco[0,:Nplaces]),
                                 dx=(T-t0)/(N+1), axis=0)

        #### Re-dilute from output states
        y0spec = np.divide(xco[-1,:Nplaces],dilution)
        for idxy, y in enumerate(y0spec):
            if y<1.0: # Then remove identity marker
                y0[tube,idxy] = 0.0

                #### Remove species ID and growth parameters
                pars['specid'][tube,idxy] = -1
                pars['Y'][tube,idxy] = 0.0
                for tr,trlen in zip(('r','m','finv'),(Nnutr,Ntox,Ntox)):
                    rmat = pars[tr]
                    rmat[tube,idxy,:] = np.zeros(trlen)

                if pars['d'].ndim>2:
                    pars['d'][tube,idxy,:] = np.zeros(Ntox)
                else:
                    pars['d'][tube,idxy] = 0.0
            else:
                y0[tube,idxy] = y

    return((xouts,I_AUC))


def mutate_tubes_sequential_disassembly(rng,y0tubes,Nsize,mutpars,pars):
    ##### For each tube, each species ...
    ##    if mutrate > r~Uni(0,1) then mutate that species
    ##
    ####### Input
    ## y0tubes   A list or array, length Ntubes, of initial abundances
    ##             [y0_spec1, ..., y0specN, y0mut1, ... y0mutN]
    ## Nsize     tuple (N_species, N_nutrients, N_toxic, N_tubes)
    ## mutpars   tuple (mutprob, flag_mutrm, mutdev)
    ## pars_r    trait parameters of resident

    Nspec, Nnutr, Ntox, Ntubes = Nsize
    mutprob, mutdev, mutfac = mutpars

    ## Draw mutations at random for each tube ...
    for itb, y_tube in enumerate(y0tubes):
        idxfree = [i for i,y in enumerate(y_tube) if y<1.0] # Find all free spots
        idxsurv = [i for i,y in enumerate(y_tube) if y>0.0] # And the indices of surviving species

        Nfree = len(idxfree)
        Nsurv = len(idxsurv)

        ##### Check that there exists ...
        #  a) a free spot for the mutants and
        #  b) a species to mutate
        if Nfree>0 and Nsurv>0:
            ### There is a fixed number of slots for populations at present
            if Nfree>Nsurv:
                ## We can draw as many mutants as surviving strains
                Nmuts = rng.binomial(n=Nsurv,p=mutprob)
            else:
                ## The number of establishing mutants constrained by the number of free spots
                Nmuts = rng.binomial(n=Nfree,p=mutprob)

            ### If we have drawn any mutations
            if Nmuts>0:
                mutating_strains = rng.choice(idxsurv,size=Nmuts,replace=False)
                for rank_mut, mutidx in enumerate(mutating_strains):
                    idxfr = idxfree[rank_mut]
                    ##  Find one or more non-zero traits to mutate
                    itr_nonzero = [i for i,val in enumerate(pars['finv'][itb,mutidx]) if val>0]
                    itr_chosen = rng.binomial(n=1,p=1.0/len(itr_nonzero),size=len(itr_nonzero))
                    if sum(itr_chosen)<1:
                        itrs = rng.choice(itr_nonzero,size=1)
                    else:
                        itrs = [itr_nonzero[i] for i,v in enumerate(itr_chosen) if v>0]
                    mutargs = (mutdev, mutfac, idxfr)
                    mutate_fixf_random(rng,pars,(itb,mutidx,itrs),mutargs)

                    ### Introduce mutant at 10% of the ancestral pop size
                    #y0tubes[itb,idxfr] = np.divide(ytot_prev,10.0)
                    y0tubes[itb,idxfr] = np.divide(y0tubes[itb,mutidx],10.0)

                    ### In disassembly we re-scale the population size
                    #IDs_tube = pars['specid'][itb]
                    #sID_ancestral = IDs_tube[mutidx]
                    ## All strains of this species
                    #idx_thisspec = [i for i,s in enumerate(IDs_tube) if s==sID_ancestral]

                    ## Re-scale the total pop size to sum to 1
                    #ytot_post = np.sum([y0tubes[itb,i] for i in idx_thisspec])
                    #for i in idx_thisspec:
                    #    y0tubes[itb,i] = y0tubes[itb,i]/ytot_post

    return()


def mutate_and_compete_dis(pars,y0,exppars,comppars):
    ##### Given a set of tubes in a particular transfer round
    ###    perform one growth period, competing between variants, and dilute

    y0_nutr, y0_tox = exppars['y0substrate']
    y0_spec         = exppars['y0spec']
    Nsize           = exppars['Nsize']
    N, t0, T = comppars['timeframe']
    ftup     = comppars['ftup']
    mutpars  = comppars['mutpars']
    rng      = comppars['RNG']

    Nspec, Nnutr, Ntox, Ntubes = Nsize
    Nreg  = Nnutr+Ntox
    Nplaces = 2*Nreg
    Nmut  = Nplaces-Nspec

    ### For each tube
    xouts = np.zeros((Ntubes,Nplaces+Nreg))
    I_AUC = np.zeros((Ntubes,Nplaces))
    xextinct = [tuple()]*Ntubes

    ### Parameters per tube
    argsco = copy.deepcopy(pars)
    parsizetup = (np.zeros((Nplaces,Ntox)),np.zeros((Nplaces,Ntox)),
                  np.zeros((Nplaces,Nnutr)),np.zeros(Nplaces),
                  np.zeros(Nplaces),np.zeros(Nplaces))
    for k,sz in zip(('finv', 'm', 'r', 'Y', 'd','specid'),parsizetup):
        #### Initialize each parameter with zero-matrices of the correct shape
        argsco[k] = sz

    specIDs_pertube = []
    for tube in range(Ntubes):
        #### Select corresponding species parameters
        np.copyto(argsco['specid'][:], pars['specid'][tube,:])
        np.copyto(argsco['d'][:], pars['d'][tube,:])
        np.copyto(argsco['Y'][:], pars['Y'][tube,:])
        # argsco['KN'] = pars['KN'][tube].copy()
        # argsco['KT'] = pars['KT'][tube].copy()

        for trait in ('finv','m','r'):
            np.copyto(argsco[trait][:], pars[trait][tube,:])

        #### Re-normalise inoculation size in y0, define substrate concentrations
        y0_thistube = y0[tube]
        IDs_thistube = pars['specid'][tube]
        unique_specids = [sid for sid in np.unique(IDs_thistube) if sid>=0]
        specIDs_pertube.append(unique_specids)

        # Step 1: calculate strain frequencies
        for sid in unique_specids:
            ids_thisspec = [i for i,s in enumerate(IDs_thistube) if s==sid]
            y0_thistube[ids_thisspec]/=np.sum(y0_thistube[ids_thisspec])
        # Step 2: Re-normalise to have y0_spec abundance per species
        y0_input = np.multiply(y0_thistube, y0_spec)
        y0tmp = np.hstack((y0_input, np.repeat(y0_nutr, Nnutr), np.repeat(y0_tox, Ntox)))

        comppars_co = (N, t0, T, y0tmp)
        # tco, xco = selsim.computeODE(ftup, comppars_co, argsco)
        tco, xco = computeODE(ftup, comppars_co, argsco)
        np.copyto(xouts[tube,:], xco[-1,:])
        np.copyto(I_AUC[tube,:],
                  np.trapz(np.subtract(xco[:,:Nplaces],xco[0,:Nplaces]),
                            dx=(T-t0)/(N+1), axis=0))

        ### For species in this tube, rescale, remove or mark extinction
        y0out = xco[-1,:Nplaces]
        # specids = pars['specid'][tube]
        # specs_unique = [s for s in np.unique(specids) if s>=0]

        ### For each species in tube ...
        for s in unique_specids:
            idces_strains = [idxs for idxs,spec in enumerate(IDs_thistube) if spec==s]

            ## Mark which strains survived and which vanished
            idx_surv = [i for i in idces_strains if y0out[i]>=1.0]
            idx_vanish = [i for i in idces_strains if y0out[i]<1.0]

            ## If any survivors, store their pop sizes
            if len(idx_surv)>0:
                for isv in idx_surv:
                    y0[tube,isv] = y0out[isv]
                ## Old version: Normalise population size to defined limits
                ## y0sum = np.sum([y0out[i] for i in idx_surv])
                #for i in idx_surv:
                #    y0[tube,i] = y0out[i]*y0_spec/y0sum
            else:
                ## If the list is empty, mark extinction.
                ## We re-introduce the species in the calling function.
                xextinct[tube] += (s,)

            ## Remove vanished strains
            for i in idx_vanish:
                ## Set pop size to zero
                y0[tube,i] = 0.0
                ## Remove species ID and growth parameters
                pars['specid'][tube,i] = -1
                pars['Y'][tube,i] = 0.0
                for tr,trlen in zip(('r','m','finv'),(Nnutr,Ntox,Ntox)):
                    # rmat, Kr, rn = pars[tr]
                    rmat = pars[tr]
                    rmat[tube,i,:] = np.zeros(trlen)

                if pars['d'].ndim>2:
                    pars['d'][tube,i,:] = np.zeros(Ntox)
                else:
                    pars['d'][tube,i] = 0.0

    ### Mutate if there are any free spots
    spaceformutant = y0<1.0
    if spaceformutant.any(): # Shorthand - don't bother mutating if no space for mutants
        # mutate_tubes_disass(y0,Nsize,mutpars,pars)
        mutate_tubes_sequential_disassembly(rng,y0,Nsize,mutpars,pars)

    return((xouts,I_AUC,xextinct,specIDs_pertube))


###########################################################################
############### Functions for dilute and transfer
###########################################################################

def compute_dilutions_stability(pars,y0,exppars,printpars,computepars):
    #### Transfer selected communities by dilution to find stability properties
    ###  INPUT
    #  'pars'     dict of community growth parameters,
    #             passed from last round of selection
    #  'y0'       [Ntubes x Nplaces] array of population sizes,
    #             passed from last round of selection
    #  'exppars'  dict of experimental parameters. Ntransf is updated to reflect
    #             the number of transfers here

    #### Read experimental parameters
    y0_nutr, y0_tox = exppars['y0substrate']
    Nspec, Nnutr, Ntox, Ntubes  = exppars['Nsize']
    Ntransfers                  = exppars['Ntsf_stability']
    dilution, frac_bneck        = exppars['bneckdil']

    #### Output files organised in dictionary
    fnames_dict  = printpars['fnames']
    fntup_state  = fnames_dict['fname_states']   # Steady-state values of state params S1,S2,...,N1,...,T1
    fn_AUCs      = fnames_dict['fname_AUC']

    ## State and mutant trait filenames
    fn_popstates,fn_nutrstates,fn_toxstates = fntup_state

    ### Arguments for tube community growth in step 1
    Nreg  = Nnutr+Ntox
    Nplaces = 2*Nreg

    tsfpars  = {'y0substrate':(y0_nutr, y0_tox),
                'dilution':dilution,
                'Nsize':exppars['Nsize']}
    comppars = {'ftup':computepars['ftup'],
                'timeframe':computepars['timeframe'], # Tuple: timeframe = (Nstep, t0, tstop)
                'RNG':exppars['RNG']}
    traittup = ('finv', )
    traitlenarr = (Ntox, )

    #### Transfers start here!
    for transfer in range(Ntransfers):
        #### Step 1: Grow, mutate and compute degradation
        youts, I_AUC = compete_nomutation(pars,y0,tsfpars,comppars)

        #### Intermission 1.1: Print end state concentrations
        with open(fn_popstates,'a') as f:
            for tube, yvectube in enumerate(youts[:,:Nplaces]):
                for spec, abund in enumerate(yvectube):
                    idxprint = [transfer,tube,spec]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(abund))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{abund:12.2e}')+' \n')

        with open(fn_nutrstates,'a') as f:
            for tube, yvectube in enumerate(youts[:,Nplaces:(Nplaces+Nnutr)]):
                for nutr, conc in enumerate(yvectube):
                    idxprint = [transfer,tube,nutr]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(conc))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{conc:12.2e}')+' \n')

        with open(fn_toxstates,'a') as f:
            for tube, yvectube in enumerate(youts[:,(Nplaces+Nnutr):]):
                for tox, conc in enumerate(yvectube):
                    idxprint = [transfer,tube,tox]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(conc))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{conc:12.2e}')+' \n')

        #### Intermission 1.4: Print AUC
        with open(fn_AUCs,'a') as f:
            for tube, AUC_pertube in enumerate(I_AUC):
                for spec, auc in enumerate(AUC_pertube):
                    idxprint = [transfer,tube,spec]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(auc))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{auc:12.2e}')+' \n')

    return()

def compute_dilutions_lineages(pars,record,y0,exppars,printpars,computepars):
    #### Transfer a given set of communities by the top-N (propagule) method
    ###  INPUT
    #  'pars' is a dict of all needed parameters
    #   Nspec, Nnutr, Ntox = commpars['size']
    #     rmat = commpars['r']     # array, shape Nspec*Nnutr
    #     mmat = commpars['m']     # array, shape Nspec*Ntox
    #     fmat = commpars['finv']  # array, shape Nspec*Ntox
    #     Yvec = commpars['Y']     # array, len Nspec
    #     dvec = commpars['d']     # array, len Nspec

    #### Read experimental parameters
    y0_nutr, y0_tox    = exppars['y0substrate']
    Nspec, Nnutr, Ntox, Ntubes  = exppars['Nsize']
    Ntransfers                  = exppars['Ntransf']
    dilution, frac_bneck        = exppars['bneckdil']
    rng                         = exppars['RNG']

    #### Output files organised in dictionary
    fnames_dict  = printpars['fnames']
    fntup_state  = fnames_dict['fname_states']   # Steady-state values of state params S1,S2,...,N1,...,T1
    fn_mutpars   = fnames_dict['fname_mutpars']  # Parameters to mutate: only (f,) at present
    # fntup_fixpars = fnames_dict['fname_fixpars']  # Fixed parameters, printed in beginning
    fn_specid    = fnames_dict['fname_specID']
    fn_AUCs      = fnames_dict['fname_AUC']

    ## State and mutant trait filenames
    fn_popstates,fn_nutrstates,fn_toxstates = fntup_state
    # fnout_growth,fnout_mort,fnout_yield,fnout_degr = fntup_fixpars

    ### Arguments for tube community growth in step 1
    Nreg  = Nnutr+Ntox
    Nplaces = 2*Nreg

    tsfpars  = {'y0substrate':(y0_nutr, y0_tox),
                'dilution':dilution,
                'Nsize':exppars['Nsize']}
    comppars = {'ftup':computepars['ftup'],
                'timeframe':computepars['timeframe'], # Tuple: timeframe = (Nstep, t0, tstop)
                'mutpars':computepars['mutpars'],
                'RNG':exppars['RNG']}
    traittup = ('finv', )
    traitlenarr = (Ntox, )

    #### Transfers start here!
    for transfer in range(Ntransfers):
        #### Step 0: print present community composition and species ID
        with open(fn_specid,'a') as f:
            for tube,matID in enumerate(pars['specid']):
                for sp,id in enumerate(matID):
                    IDprint = [transfer,tube,tube,sp,round(id)]
                    f.write(''.join('{0:12d}'.format(el) for el in IDprint)+'\n')
                    # f.write(''.join('{0:12.2e}'.format(id))+'\n')
                    # f.write(''.join(f'{el:12.2d}' for el in IDprint)+' \n')

        #### Step 0.1: Print parameter values
        #for trait, fntrt in zip(traittup,fntup_mutpars):
        #    with open(fntrt,'a') as f:
        #        for tube in range(Ntubes):
        #            parprint = np.append([transfer,tube],pars[trait][tube].flatten())
        #            f.write(''.join('{0:12.2e}'.format(el) for el in parprint)+'\n')
        #            # f.write(''.join(f'{el:12.2e}' for el in parprint)+' \n')
        #for trait, trlen, fntrt in zip(traittup,traitlenarr,fn_mutpars):
        with open(fn_mutpars,'a') as f:
            for tube,fmat in enumerate(pars['finv']):
                for sp,fvec in enumerate(fmat):
                    for tr,val in enumerate(fvec):
                        idxprint = [transfer,tube,sp,tr]
                        f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                        f.write(''.join('{0:12.2e}'.format(val))+'\n')
                        # f.write(''.join(f'{el:12d}' for el in printarr)+' \t')
                        # f.write(''.join(f'{val:12.2e}')+' \n')

        #### Step 1: Grow, mutate and compute degradation
        youts, I_AUC = mutate_and_compete_fixf(pars,y0,tsfpars,comppars)

        #### Intermission 1.1: Print end state concentrations
        with open(fn_popstates,'a') as f:
            for tube, yvectube in enumerate(youts[:,:Nplaces]):
                for spec, abund in enumerate(yvectube):
                    idxprint = [transfer,tube,spec]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(abund))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{abund:12.2e}')+' \n')

        with open(fn_nutrstates,'a') as f:
            for tube, yvectube in enumerate(youts[:,Nplaces:(Nplaces+Nnutr)]):
                for nutr, conc in enumerate(yvectube):
                    idxprint = [transfer,tube,nutr]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(conc))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{conc:12.2e}')+' \n')

        with open(fn_toxstates,'a') as f:
            for tube, yvectube in enumerate(youts[:,(Nplaces+Nnutr):]):
                for tox, conc in enumerate(yvectube):
                    idxprint = [transfer,tube,tox]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(conc))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{conc:12.2e}')+' \n')

        #with open(fn_states,'a') as f:
        #    for idxt, y_out in enumerate(youts):
        #        stateprint = np.append([transfer,idxt],y_out)
        #        f.write(''.join('{0:12.2e}'.format(el) for el in stateprint)+'\n')
        #        # f.write(''.join(f'{el:12.2e}' for el in stateprint)+' \n')

        #### Intermission 1.4: Print AUC
        with open(fn_AUCs,'a') as f:
            for tube, AUC_pertube in enumerate(I_AUC):
                for spec, auc in enumerate(AUC_pertube):
                    idxprint = [transfer,tube,spec]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(auc))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{auc:12.2e}')+' \n')

        #### Step 2: Rank tubes wrt their degradation in increasing order
        Tdegr = [np.linalg.norm(y[-Ntox:]) for y in youts] # 2-norm emphasises abundant substrate
        # Tdegr = [np.max(y[-Ntox:]) for y in youts] # Max of toxins could be numerically unstable
        # Tdegr = [np.sum(y[-Ntox:]) for y in youts] # Sum does not emphasise as much
    return(pars,y0)


def compute_dilutions_topN_nosel(pars,record,y0,exppars,printpars,computepars):
    #### Transfer a given set of communities by a NO-SELECTION VERSION of the top-N (propagule) method
    ###  INPUT
    #  'pars' is a dict of all needed parameters
    #   Nspec, Nnutr, Ntox = commpars['size']
    #     rmat = commpars['r']     # array, shape Nspec*Nnutr
    #     mmat = commpars['m']     # array, shape Nspec*Ntox
    #     fmat = commpars['finv']  # array, shape Nspec*Ntox
    #     Yvec = commpars['Y']     # array, len Nspec
    #     dvec = commpars['d']     # array, len Nspec

    #### Read experimental parameters
    y0_nutr, y0_tox    = exppars['y0substrate']
    Nspec, Nnutr, Ntox, Ntubes  = exppars['Nsize']
    Ntransfers                  = exppars['Ntransf']
    dilution, frac_bneck        = exppars['bneckdil']
    rng                         = exppars['RNG']

    #### Output files organised in dictionary
    fnames_dict   = printpars['fnames']
    fntup_state   = fnames_dict['fname_states']   # Steady-state values of state params S1,S2,...,N1,...,T1
    fn_mutpars    = fnames_dict['fname_mutpars']  # Parameters to mutate: only (f,) at present
    # fntup_fixpars = fnames_dict['fname_fixpars']  # Fixed parameters, printed in beginning
    fn_specid     = fnames_dict['fname_specID']
    fn_AUCs      = fnames_dict['fname_AUC']

    ## State and mutant trait filenames
    fn_popstates,fn_nutrstates,fn_toxstates = fntup_state
    # fnout_growth,fnout_mort,fnout_yield,fnout_degr = fntup_fixpars

    ### Arguments for tube community growth in step 1
    Nreg  = Nnutr+Ntox
    Nplaces = 2*Nreg

    tsfpars  = {'y0substrate':(y0_nutr, y0_tox),
                'dilution':dilution,
                'Nsize':exppars['Nsize']}
    comppars = {'ftup':computepars['ftup'],
                'timeframe':computepars['timeframe'], # Tuple: timeframe = (Nstep, t0, tstop)
                'mutpars':computepars['mutpars'],
                'RNG':exppars['RNG']}
    traittup = ('finv', )
    traitlenarr = (Ntox, )

    ### Arguments for selection in step 2
    Nbottleneck = round(Ntubes*frac_bneck)
    #Ncopiesofoldtubes = int(1.0/frac_bneck)
    Ncopiesofoldtubes = Ntubes//Nbottleneck

    #### Transfers start here!
    for transfer in range(Ntransfers):
        #### Print community composition and species ancestral ID
        with open(fn_specid,'a') as f:
            for tube,matID in enumerate(pars['specid']):
                for sp,id in enumerate(matID):
                    IDprint = [transfer,tube,tube,sp,round(id)]
                    f.write(''.join('{0:12d}'.format(el) for el in IDprint)+'\n')
                    # f.write(''.join(f'{el:12.2e}' for el in IDprint)+' \n')

        #### Print parameter values
        #for trait, fntrt in zip(traittup,fntup_mutpars):
        #    with open(fntrt,'a') as f:
        #        for tube in range(Ntubes):
        #            parprint = np.append([transfer,tube],pars[trait][tube].flatten())
        #            f.write(''.join('{0:12.2e}'.format(el) for el in parprint)+'\n')
        #            # f.write(''.join(f'{el:12.2e}' for el in parprint)+' \n')
        #for trait, trlen, fntrt in zip(traittup,traitlenarr,fn_mutpars):
        with open(fn_mutpars,'a') as f:
            for tube,fmat in enumerate(pars['finv']):
                for sp,fvec in enumerate(fmat):
                    for tr,val in enumerate(fvec):
                        idxprint = [transfer,tube,sp,tr]
                        f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                        f.write(''.join('{0:12.2e}'.format(val))+'\n')
                        # f.write(''.join(f'{el:12d}' for el in printarr)+' \t')
                        # f.write(''.join(f'{val:12.2e}')+' \n')

        #### Step 1: Grow, mutate and compute degradation
        youts, I_AUC = mutate_and_compete_fixf(pars,y0,tsfpars,comppars)

        #### Intermission 1.1: Print end state concentrations
        with open(fn_popstates,'a') as f:
            for tube, yvectube in enumerate(youts[:,:Nplaces]):
                for spec, abund in enumerate(yvectube):
                    idxprint = [transfer,tube,spec]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(abund))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{abund:12.2e}')+' \n')

        with open(fn_nutrstates,'a') as f:
            for tube, yvectube in enumerate(youts[:,Nplaces:(Nplaces+Nnutr)]):
                for nutr, conc in enumerate(yvectube):
                    idxprint = [transfer,tube,nutr]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(conc))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{conc:12.2e}')+' \n')

        with open(fn_toxstates,'a') as f:
            for tube, yvectube in enumerate(youts[:,(Nplaces+Nnutr):]):
                for tox, conc in enumerate(yvectube):
                    idxprint = [transfer,tube,tox]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(conc))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{conc:12.2e}')+' \n')

        #with open(fn_states,'a') as f:
        #    for idxt, y_out in enumerate(youts):
        #        stateprint = np.append([transfer,idxt],y_out)
        #        f.write(''.join('{0:12.2e}'.format(el) for el in stateprint)+'\n')
        #        # f.write(''.join(f'{el:12.2e}' for el in stateprint)+' \n')

        #### Intermission 1.3: Print AUC
        with open(fn_AUCs,'a') as f:
            for tube, AUC_pertube in enumerate(I_AUC):
                for spec, auc in enumerate(AUC_pertube):
                    idxprint = [transfer,tube,spec]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(auc))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{auc:12.2e}')+' \n')

        #### Step 2, old version: Rank tubes wrt their degradation in increasing order
        ### Previously, we would perform same step as in propagule selection,
        ##  but on a shuffled array
        # Tdegr = [np.linalg.norm(y[-Ntox:]) for y in youts] # 2-norm emphasises abundant substrate
        # np.random.shuffle(Tdegr)
        # ranks = np.argsort(Tdegr)
        # toptubes_ranked = ranks[:Nbottleneck] # Choose from shuffled top ranks

        ### Step 2: New tubes by uniformly random draws
        toptubes_ranked = rng.choice(Ntubes,size=Nbottleneck,replace=False)

        #### Step 3: Re-populate parameters
        # Old version where we drew new tubes randomly (uniform)
        # Removed since I don't think any study implemented stochasticity here
        # idx_besttubes = np.random.choice(toptubes_ranked,size=Ntubes)
        idx_besttubes = np.repeat(toptubes_ranked, Ncopiesofoldtubes)
        for tube_to, tube_from in enumerate(idx_besttubes):
            y0[tube_to] = np.array(y0[tube_from])
            pars['specid'][tube_to] = np.array(pars['specid'][tube_from])
            for trait in ('finv','r','m','Y','d'):
                pars[trait][tube_to] = copy.deepcopy(pars[trait][tube_from])

    return(pars,y0)


def compute_dilutions_propagule(pars,record,y0,exppars,printpars,computepars):
    #### Transfer a given set of communities by the top-N (propagule) method
    ###  INPUT
    #  'pars' is a dict of all needed parameters
    #   Nspec, Nnutr, Ntox = commpars['size']
    #     rmat = commpars['r']     # array, shape Nspec*Nnutr
    #     mmat = commpars['m']     # array, shape Nspec*Ntox
    #     fmat = commpars['finv']  # array, shape Nspec*Ntox
    #     Yvec = commpars['Y']     # array, len Nspec
    #     dvec = commpars['d']     # array, len Nspec
    #  'y0' is a list (length Ntubes) of inoculum sizes for each population
    #  'exppars' is a dictionary of experimental parameters
    #  'printpars' is a dictionary of paths and filenames for the output
    #  'computepars' is a dictionary for ODE parameters


    #### Read experimental parameters
    y0_nutr, y0_tox    = exppars['y0substrate']
    Nspec, Nnutr, Ntox, Ntubes  = exppars['Nsize']
    Ntransfers                  = exppars['Ntransf']
    dilution, frac_bneck        = exppars['bneckdil']
    # rng                         = exppars['RNG']

    #### Output files organised in dictionary
    fnames_dict   = printpars['fnames']
    fntup_state   = fnames_dict['fname_states']   # Steady-state values of state params S1,S2,...,N1,...,T1
    fn_mutpars    = fnames_dict['fname_mutpars']  # Parameters to mutate: only (f,) at present
    # fntup_fixpars = fnames_dict['fname_fixpars']  # Fixed parameters, printed in beginning
    fn_specid     = fnames_dict['fname_specID']
    fn_AUCs      = fnames_dict['fname_AUC']

    ## State and mutant trait filenames
    fn_popstates,fn_nutrstates,fn_toxstates = fntup_state
    # fnout_growth,fnout_mort,fnout_yield,fnout_degr = fntup_fixpars

    ### Arguments for tube community growth in step 1
    Nreg  = Nnutr+Ntox
    Nplaces = 2*Nreg

    tsfpars  = {'y0substrate':(y0_nutr, y0_tox),
                'dilution':dilution,
                'Nsize':exppars['Nsize']}
    comppars = {'ftup':computepars['ftup'],
                'timeframe':computepars['timeframe'], # Tuple: timeframe = (Nstep, t0, tstop)
                'mutpars':computepars['mutpars'],
                'RNG':exppars['RNG']}

    ### Arguments for selection in step 2
    Nbottleneck = round(Ntubes*frac_bneck)
    #Ncopiesofoldtubes = int(1.0/frac_bneck)
    Ncopiesofoldtubes = Ntubes//Nbottleneck

    #### Transfers start here!
    for transfer in range(Ntransfers):
        #### Print community composition and species ancestral ID
        with open(fn_specid,'a') as f:
            for tube,matID in enumerate(pars['specid']):
                for sp,id in enumerate(matID):
                    IDprint = [transfer,tube,tube,sp,round(id)]
                    f.write(''.join('{0:12d}'.format(el) for el in IDprint)+'\n')
                    # f.write(''.join(f'{el:12.2e}' for el in IDprint)+' \n')

        #### Intermission 1.2: Print parameter values
        #for trait, fntrt in zip(traittup,fntup_mutpars):
        #    with open(fntrt,'a') as f:
        #        for tube in range(Ntubes):
        #            parprint = np.append([transfer,tube],pars[trait][tube].flatten())
        #            f.write(''.join('{0:12.2e}'.format(el) for el in parprint)+'\n')
        #            # f.write(''.join(f'{el:12.2e}' for el in parprint)+' \n')
        #for trait, trlen, fntrt in zip(traittup,traitlenarr,fn_mutpars):
        with open(fn_mutpars,'a') as f:
            for tube,fmat in enumerate(pars['finv']):
                for sp,fvec in enumerate(fmat):
                    for tr,val in enumerate(fvec):
                        idxprint = [transfer,tube,sp,tr]
                        f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                        f.write(''.join('{0:12.2e}'.format(val))+'\n')
                        # f.write(''.join(f'{el:12d}' for el in printarr)+' \t')
                        # f.write(''.join(f'{val:12.2e}')+' \n')

        #### Step 1: Grow, mutate and compute degradation
        youts, I_AUC = mutate_and_compete_fixf(pars,y0, tsfpars,comppars)

        #### Intermission 1.1: Print end state concentrations
        with open(fn_popstates,'a') as f:
            for tube, yvectube in enumerate(youts[:,:Nplaces]):
                for spec, abund in enumerate(yvectube):
                    idxprint = [transfer,tube,spec]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(abund))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{abund:12.2e}')+' \n')

        with open(fn_nutrstates,'a') as f:
            for tube, yvectube in enumerate(youts[:,Nplaces:(Nplaces+Nnutr)]):
                for nutr, conc in enumerate(yvectube):
                    idxprint = [transfer,tube,nutr]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(conc))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{conc:12.2e}')+' \n')

        with open(fn_toxstates,'a') as f:
            for tube, yvectube in enumerate(youts[:,(Nplaces+Nnutr):]):
                for tox, conc in enumerate(yvectube):
                    idxprint = [transfer,tube,tox]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(conc))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{conc:12.2e}')+' \n')

        #with open(fn_states,'a') as f:
        #    for idxt, y_out in enumerate(youts):
        #        stateprint = np.append([transfer,idxt],y_out)
        #        f.write(''.join('{0:12.2e}'.format(el) for el in stateprint)+'\n')
        #        # f.write(''.join(f'{el:12.2e}' for el in stateprint)+' \n')

        #### Intermission 1.3: Print AUC
        with open(fn_AUCs,'a') as f:
            for tube, AUC_pertube in enumerate(I_AUC):
                for spec, auc in enumerate(AUC_pertube):
                    idxprint = [transfer,tube,spec]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(auc))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{auc:12.2e}')+' \n')

        #### Step 2: Rank tubes wrt their degradation in increasing order
        Tdegr = [np.linalg.norm(y[-Ntox:]) for y in youts] # 2-norm emphasises abundant substrate
        # Tdegr = [np.max(y[-Ntox:]) for y in youts] # Max of toxins could be numerically unstable
        # Tdegr = [np.sum(y[-Ntox:]) for y in youts] # Sum does not emphasise as much

        ranks = np.argsort(Tdegr)
        toptubes_ranked = ranks[:Nbottleneck]

        #### Step 3: Re-populate parameters
        # Old version where we draw new tubes randomly (uniform)
        # Don't think any study implemented stochasticity here
        # idx_besttubes = rng.choice(toptubes_ranked,size=Ntubes)
        idx_besttubes = np.repeat(toptubes_ranked,Ncopiesofoldtubes)
        for tube_to, tube_from in enumerate(idx_besttubes):
            y0[tube_to] = np.array(y0[tube_from])
            pars['specid'][tube_to] = np.array(pars['specid'][tube_from])
            for trait in ('finv','r','m','Y','d'):
                pars[trait][tube_to] = copy.deepcopy(pars[trait][tube_from])

    # return((copy.deepcopy(pars),youts))
    return(pars,y0)


def update_fossilrecord( pars, record, specdict_tubes, transfer, y0,
                         traittup=('d','Y','r','m','finv') ):
    #### INPUT
    # pars     dict with ODE parameters 'd','Y','r','m','finv' as output by selsim.assemble_CRcomms_fromlong()
    # record   dict (keyed by species ID) of dicts with ODE parameters and
    #          metadata as output by selsim.init_fossilrecord()
    # specdict_tubes    dict of (some subset of) species indices - which tube do we draw the new record from?
    # transfer          the present transfer
    # y0                population size of communities, size [Ntubes x Nspec]
    # traittup          tuple of trait parameters to record
    #### OUTPUT
    # record         list of length Nspecpool (the number of species in the initial pool)
    #                where each entry is a dictionary of growth parameters r, m, f, Y, d
    #                in addition to metadata: ID of species, isolated from which tube in which round
    #                and the number of strains at present (useful within simulations)

    ### For each species, collect data from 'pars' and put in 'record'
    for specid, tube_from in zip(specdict_tubes.keys(),specdict_tubes.values()):
        # Overall idea: copy the traits in 'pars' of each strains from the correct tube
        specid_int = round(specid)
        idx_strains_intube = [idx for idx, spec in enumerate(pars['specid'][tube_from]) if round(spec)==specid_int]
        Nstrains = len(idx_strains_intube)
        # idx_record = [i for i,d in enumerate(record) if d['specid']==specid_int][0]

        ### Collect metadata
        Nstrains = len(idx_strains_intube)
        dict_tmp = copy.deepcopy(record[specid_int])
        dict_tmp['specid'] = specid_int
        dict_tmp['last_tube'] = tube_from
        dict_tmp['round'] = transfer
        dict_tmp['Nstrains'] = Nstrains
        if Nstrains<2:
            idx_strains_intube=idx_strains_intube[0]

        dict_tmp['y0'] = copy.deepcopy(y0[tube_from,idx_strains_intube])
        ### Collect parameter arrays
        for trait in ('d','Y'): # d, Y
            pars_fromtube = pars[trait][tube_from]
            dict_tmp[trait] = copy.deepcopy(pars_fromtube[idx_strains_intube])

        for trait in ('r', 'm', 'finv'): # r, m, f
            pars_fromtube = pars[trait][tube_from]
            dict_tmp[trait] = copy.deepcopy(pars_fromtube[idx_strains_intube,:])

        record[specid_int] = copy.deepcopy(dict_tmp)

    return(record)


def get_species_fromrecord( pars, record, dict_extincts, y0,
                            traittup=('d','Y','r','m','finv') ):
    #### INPUT
    # pars     dict with ODE parameters 'd','Y','r','m','finv' as output by selsim.assemble_CRcomms_fromlong()
    # record   dict (keyed by species ID) of dicts with ODE parameters and
    #          metadata as output by selsim.init_fossilrecord()
    # dict_extincts  dict {species_ID: list_of_tubes} - this species is extinct in these tubes
    # y0             population size of communities, size [Ntubes x Nspec]
    # traittup       tuple of trait parameters to record
    #### OUTPUT


    ### For each species, take data from 'record' and put (back) in 'pars' and 'y0'
    for specid, tubelist in zip(dict_extincts.keys(),dict_extincts.values()):
        specid_int = round(specid)

        specrecord = record[specid_int]
        Nstrains   = specrecord['Nstrains']

        for tubeidx in tubelist:
            y0_replace = specrecord['y0']
            y0tube_tmp = y0[tubeidx,:]
            IDstube_tmp = pars['specid'][tubeidx,:]

            ### Find which and how many free indices we have in the tube
            # Future alternative: each community has a custom size, modified by resizing pars[] and y0[]
            idxfree = [i for i,spid in enumerate(IDstube_tmp) if spid<0]
            Nfree   = len(idxfree)
            if Nfree<1:
                continue
                #print('Species '+str(specid)+', Tube '+str(tubeidx))
                #print('Nstrains '+str(Nstrains)+', Nfree '+str(Nfree))
                #print('IDs in community: ', IDstube_tmp)

            if Nstrains==1: # Need a special case here
                if Nfree>0:
                    idxf = idxfree[0]
                    pars['specid'][tubeidx,idxf] = specid_int
                    y0[tubeidx,idxf] = y0_replace

                    for tr in traittup: # ('d','Y','r','m','finv')
                        pars[tr][tubeidx][idxf] = copy.deepcopy(specrecord[tr])
            else:
                if Nstrains >= Nfree:
                    # Temporary solution: take only most abundant strains from record
                    IDs_record = np.argsort(y0_replace)[-Nfree:]
                    # y0_replace = [y0_replace[sid] for sid in IDs_record]
                    ## Re-normalise since we only took a portion
                    #y0rec_tot = np.sum([y0_replace[round(sid)] for sid in IDs_record])
                    #for sid in IDs_record:
                    #    y0_replace[round(sid)]/=y0rec_tot
                else:
                    # No problem, go ahead and reintroduce the species
                    IDs_record = list(range(Nstrains))

                ### Introduce the species
                for i,i_rec in enumerate(IDs_record):
                    idxf = idxfree[i]
                    y0tube_tmp[idxf]  = y0_replace[i_rec]
                    IDstube_tmp[idxf] = specid_int
                pars['specid'][tubeidx,:] = copy.deepcopy(IDstube_tmp)
                y0[tubeidx,:] = copy.deepcopy(y0tube_tmp)

                for tr in traittup: # ('d','Y','r','m','finv')
                    partmp = copy.deepcopy(pars[tr][tubeidx])
                    for i,i_rec in enumerate(IDs_record):
                        idxf = idxfree[i]
                        partmp[idxf] = specrecord[tr][round(i_rec)]
                    pars[tr][tubeidx] = copy.deepcopy(partmp)

    return()


def emigrate(rng,pars,rate_emigr,y0,sizetup,tubes_from):
    #### Function to emigrate species from tubes, used by
    #### used by compute_dilutions_disassembly() defined below

    Nspec, Nnutr, Ntox, Ntubes = sizetup

    #### For each tube, the communities to remove one or more species
    for idxt in tubes_from:
        #### Step 1: Find how many species there are in the community
        local_IDs_set = set([round(s) for s in np.unique(pars['specid'][idxt,:]) if s>=0])
        Nlocals = len(local_IDs_set)
        if Nlocals<2:
            continue
            ### Skip cases where there is only a single species in the tube
            ##    since migration would mean complete exchange of species

        #### Step 2: Do not emigrate species that are unique with respect to metacommunity
        idx_other_tbs = [t for t in range(Ntubes) if t!=idxt]
        local_IDs_metacomm = [round(s) for s in np.unique(pars['specid'][idx_other_tbs]) if s>=0]

        # 'Local' metacomm of all tubes sans 'idxt', i.e. this one
        unique_IDs_tb_set = local_IDs_set-set(local_IDs_metacomm)
        if len(unique_IDs_tb_set)>0:
            # Do not emigrate the species that do not feature in other communities
            emigrant_pool = list(local_IDs_set-unique_IDs_tb_set)
            if len(emigrant_pool)<1:
                continue
        else:
            # All species in this tube feature also in other communities
            emigrant_pool = list(local_IDs_set)

        #### Step 3: Emigrate at least one, at most Nlocals-1
        N_emigrable = np.min([len(emigrant_pool), Nlocals-1]) # In case emigrant_pool==local_IDs_set
        N_emigr = np.min([1+rng.poisson(rate_emigr), N_emigrable])
        # Debug: fix N_emigrants = 1
        # N_emigr = 1

        #### Step 4: Find the emigrants, the species to be removed
        sp_emigr = rng.choice(emigrant_pool,replace=False,size=N_emigr)

        #### Step 5: Remove species parameters
        for spid in sp_emigr:
            idces_em = [i for i,s in enumerate(pars['specid'][idxt,:]) if s==spid]
            for ie in idces_em:
                y0[idxt][ie] = 0.0
                pars['specid'][idxt,ie] = -1 # Modify-in-place?
                for tr,trlen in zip(('r','m','finv'),[Nnutr,Ntox,Ntox]):
                    pars[tr][idxt,ie,:] = np.zeros(trlen)
                pars['d'][idxt,ie] = 0.0
                pars['Y'][idxt,ie] = 0.0
                if type(pars['KN']) is np.ndarray:
                    pars['KN'][idxt,ie] = 0.0
                    pars['KT'][idxt,ie] = 0.0

    return()


def draw_tubes_for_immigration(rng,pars,Ntubes,Nremix):
    #### Draw tubes for immigration, with the guarantee that
    #    there is space for one or more immigrants
    #    This would otherwise need to be checked if we used rng.choice()

    tubes_im = [-1]*Nremix
    tubelist_shuffled = rng.choice(Ntubes,replace=False,size=Ntubes)

    for i,t in enumerate(tubelist_shuffled):
        if i==Nremix:
            break

        specids_thistube = pars['specid'][t]
        ### Check if there is a free spot
        if (-1) in specids_thistube:
            tubes_im[i] = t

    # In case we have less than Nremix
    tubes_im = [s for s in tubes_im if s>=0]
    return(tubes_im)


def immigrate_fossilrecord(rng,pars,rate_immigr,y0,sizetup,tubes_to,record):
    #### Immigrate speces from the fossil record to the tubes 'tubes_to'
    ####   using the already defined get_species_fromrecord() that is
    ####   used in the disassembly method to re-introduce extinct species
    ####
    ##   INPUT
    ####   rng    the random number generator
    ####   pars   dict of growth parameters.
    ####          Keys ['size', 'specid', 'finv', 'm', 'r',
    ####                'KN', 'KT', 'Y', 'd', 'nHill', 'fpow']
    ####   rate_immigr    float, parameter for Poisson(rate) to draw
    ####                  the number of introduced species per tube
    ####   y0     Array of population sizes
    ####   sizetup        tuple (Nspec, Nnutr, Ntox, Ntubes)
    ####   tubes_to       Array of tubes to receive a migrant species
    ####   record         The fossil record, defined by init_fossilrecord()

    Nspec, Nnutr, Ntox, Ntubes = sizetup
    dict_tubes_immigr = {}
    dict_species_by_tubes = {}

    #### Step 0a: Define meta-community - the collection of species present in at least one tube
    ## Note that the metacomm can expand through the loop as species are introduced
    specIDs_metacomm = set([s for s in np.unique(pars['specid']) if s>=0])
    # At present, species IDs are 0, 1, 2, 3 ... while -1 is reserved for empty space
    spec_IDs_record_set = set(record.keys()) # All species in the simulation
    #### Step 0b: Immigrate all species that are not present in the metacommunity
    specIDs_notpres_set = spec_IDs_record_set-specIDs_metacomm
    for s in list(specIDs_notpres_set):
        tubechoice = rng.choice(tubes_to)
        dict_tubes_immigr[s] = [tubechoice,]
        if tubechoice in dict_species_by_tubes.keys():
            dict_species_by_tubes[tubechoice] = np.append(dict_species_by_tubes[tubechoice],[s,])
        else:
            dict_species_by_tubes[tubechoice] = [s,]

    #### Step 1: Find a (maximum) number of species to immigrate and free spots for the strains
    N_immigr_vec = rng.poisson(rate_immigr,size=len(tubes_to))
    for idxt,t in enumerate(tubes_to):
        if t not in dict_species_by_tubes.keys():
            # I.e. if there is no not-present-in-metacomm species in this tube
            N_immigr_vec[idxt] += 1 # Immigrate at least one species per tube

    #### For each 'tube-to', the communities to receive one or more migrating species
    for t,N_immigr in zip(tubes_to,N_immigr_vec):
        #### Step 2a: Find the 'local' species presence in this tube
        local_IDs_set = set([s for s in np.unique(pars['specid'][t,:]) if s>=0])
        # local_IDs_set = set(pars['specid'][idxt,:])-{-1,} # Equivalent
        pool_immigrants = spec_IDs_record_set-local_IDs_set

        #### Step 2b: Also do not repeat species that we have already marked for immigration
        if idxt in dict_species_by_tubes.keys():
            pool_immigrants -= set(dict_species_by_tubes[idxt])

        #### Step 2c: Choose immigrants
        N_immigr  = min(N_immigr, len(pool_immigrants)) # In case we have few potential immigrant species left
        sp_immigr = rng.choice(list(pool_immigrants), replace=False, size=N_immigr)
        ######### Debug: Draw only one immigrant
        #sp_immigr = rng.choice(list(pool_immigrants), replace=False, size=1)

        #### Step 3: Note all immigrants
        # Q: Also update dict_species_by_tubes?
        for spec in sp_immigr:
            if spec in dict_tubes_immigr.keys():
                dict_tubes_immigr[spec] = np.append(dict_tubes_immigr[spec],[t,])
            else:
                dict_tubes_immigr[spec] = [t,]

    #### Step 4: Introduce parameters, pop size etc for immigrants
    get_species_fromrecord( pars, record, dict_tubes_immigr, y0 )

    return()

def compute_dilutions_disassembly(pars,record,y0,exppars,printpars,computepars):
    #### Transfer a given set of communities by the disassembly method
    ###  INPUT
    #  'pars' is a dict of all needed parameters
    #   Nspec, Nnutr, Ntox = commpars['size']
    #     rmat = commpars['r']     # array, shape Nspec*Nnutr
    #     mmat = commpars['m']     # array, shape Nspec*Ntox
    #     fmat = commpars['finv']  # array, shape Nspec*Ntox
    #     Yvec = commpars['Y']     # array, len Nspec
    #     dvec = commpars['d']     # array, len Nspec
    #  'y0' is a list (length Ntubes) of inoculum sizes for each population
    #  'exppars' is a dictionary of experimental parameters
    #  'printpars' is a dictionary of paths and filenames for the output
    #  'computepars' is a dictionary for ODE parameters

    #### Read experimental parameters
    y0_spec            = exppars['y0spec']
    y0_nutr, y0_tox    = exppars['y0substrate']
    Nspec, Nnutr, Ntox, Ntubes  = exppars['Nsize']
    Ntransfers                  = exppars['Ntransf']
    dilution, frac_bneck        = exppars['bneckdil']
    mixrate, rate_migr_poisson  = exppars['mixrate']
    rng                         = exppars['RNG']

    #### Output files organised in dictionary
    fnames_dict   = printpars['fnames']
    fntup_state   = fnames_dict['fname_states']   # Steady-state values of state params S1,S2,...,N1,...,T1
    fn_mutpars    = fnames_dict['fname_mutpars']  # Parameters to mutate: only (f,) at present
    # fntup_fixpars = fnames_dict['fname_fixpars']  # Fixed parameters, printed in beginning
    fn_specid     = fnames_dict['fname_specID']
    fn_AUCs      = fnames_dict['fname_AUC']

    ## State and mutant trait filenames
    fn_popstates,fn_nutrstates,fn_toxstates = fntup_state
    # fnout_growth,fnout_mort,fnout_yield,fnout_degr = fntup_fixpars

    ### Arguments for tube community growth in step 1
    Nreg  = Nnutr+Ntox
    Nplaces = 2*Nreg

    Nsize = exppars['Nsize']
    tsfpars  = {'y0substrate':(y0_nutr, y0_tox),
                'y0spec':y0_spec,
                'dilution':dilution,
                'Nsize':Nsize}
    comppars = {'ftup':computepars['ftup'],
                'timeframe':computepars['timeframe'], # Tuple: timeframe = (Nstep, t0, tstop)
                'mutpars':computepars['mutpars'],
                'RNG':exppars['RNG']}

    ### Arguments for selection in step 2
    ytoxmax = np.sqrt(Ntox)*y0_tox
    Nremix  = round(mixrate*Ntubes)
    Nbottleneck = round(Ntubes*frac_bneck)

    #### Transfers start here!
    for transfer in range(Ntransfers):
        #### Print community composition and species ancestral ID
        with open(fn_specid,'a') as f:
            for tube,matID in enumerate(pars['specid']):
                for sp,id in enumerate(matID):
                    IDprint = [transfer,tube,tube,sp,round(id)]
                    f.write(''.join('{0:12d}'.format(el) for el in IDprint)+'\n')
                    # f.write(''.join(f'{el:12.2e}' for el in IDprint)+' \n')

        #### Intermission 1.2: Print parameter values
        #for trait, fntrt in zip(traittup,fntup_mutpars):
        #    with open(fntrt,'a') as f:
        #        for tube in range(Ntubes):
        #            parprint = np.append([transfer,tube],pars[trait][tube].flatten())
        #            f.write(''.join('{0:12.2e}'.format(el) for el in parprint)+'\n')
        #            # f.write(''.join(f'{el:12.2e}' for el in parprint)+' \n')
        #for trait, trlen, fntrt in zip(traittup,traitlenarr,fn_mutpars):
        with open(fn_mutpars,'a') as f:
            for tube,fmat in enumerate(pars['finv']):
                for sp,fvec in enumerate(fmat):
                    for tr,val in enumerate(fvec):
                        idxprint = [transfer,tube,sp,tr]
                        f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                        f.write(''.join('{0:12.2e}'.format(val))+'\n')
                        # f.write(''.join(f'{el:12d}' for el in printarr)+' \t')
                        # f.write(''.join(f'{val:12.2e}')+' \n')

        #### Step 1: Grow, mutate and compute degradation
        youts, I_AUC, y_exts, IDpertube = mutate_and_compete_dis(pars,y0,tsfpars,comppars)
        ### Make sure that y0 is normalised before competition

        #### Intermission 1.1: Print end state concentrations
        with open(fn_popstates,'a') as f:
            for tube, yvectube in enumerate(youts[:,:Nplaces]):
                for spec, abund in enumerate(yvectube):
                    idxprint = [transfer,tube,spec]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(abund))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{abund:12.2e}')+' \n')

        with open(fn_nutrstates,'a') as f:
            for tube, yvectube in enumerate(youts[:,Nplaces:(Nplaces+Nnutr)]):
                for nutr, conc in enumerate(yvectube):
                    idxprint = [transfer,tube,nutr]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(conc))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{conc:12.2e}')+' \n')

        with open(fn_toxstates,'a') as f:
            for tube, yvectube in enumerate(youts[:,(Nplaces+Nnutr):]):
                for tox, conc in enumerate(yvectube):
                    idxprint = [transfer,tube,tox]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(conc))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{conc:12.2e}')+' \n')

        #with open(fn_states,'a') as f:
        #    for idxt, y_out in enumerate(youts):
        #        stateprint = np.append([transfer,idxt],y_out)
        #        f.write(''.join('{0:12.2e}'.format(el) for el in stateprint)+'\n')
        #        # f.write(''.join(f'{el:12.2e}' for el in stateprint)+' \n')

        #### Intermission 1.4: Print AUC
        with open(fn_AUCs,'a') as f:
            for tube, AUC_pertube in enumerate(I_AUC):
                for spec, auc in enumerate(AUC_pertube):
                    idxprint = [transfer,tube,spec]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(auc))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{auc:12.2e}')+' \n')

        #### Step 2a: Rank tubes wrt their degradation in increasing order
        Tdegr = [np.linalg.norm(y[-Ntox:]) for y in youts] # 2-norm emphasises abundant substrate
        # Tdegr = [np.max(y[-Ntox:]) for y in youts] # Max of toxins could be numerically unstable
        # Tdegr = [np.sum(y[-Ntox:]) for y in youts] # Sum does not emphasise as much
        ranks = np.argsort(Tdegr)
        ## Note that sorting T concentrations in ascending order
        ## is equivalent to sorting degradation in decreasing
        toptubes_ranked = ranks[:Nbottleneck]
        topdegr_prel  = [np.divide(ytoxmax-Tdegr[r],ytoxmax) for r in toptubes_ranked]
        # topdegr_prel = [np.divide(1.0,Tdegr[r]) for r in toptubes_ranked] # To further emphasise differences
        survival = np.array([np.divide(len(IDpertube[r])-len(y_exts[r]),len(IDpertube[r])) for r in toptubes_ranked])

        #### 2b: Update the species record
        toptubes_spec = pars['specid'][toptubes_ranked]
        spec_fromtube = {}
        for tr,tb,ye in zip(toptubes_ranked, toptubes_spec, y_exts):
            specs_tube = tb[tb>=0] # Using -1 as a marker for empty slots
            specs_new  = [s for s in specs_tube if s not in spec_fromtube.keys()]
            # Double-check that the species is not extinct
            specs_new  = [s for s in specs_new if s not in ye]
            if len(specs_new)>0:
                for s in specs_new:
                    spec_fromtube[s]=tr

        record = update_fossilrecord(pars, record, spec_fromtube, transfer, y0)

        #### Step 3: Re-introduce extinct species from the (updated) fossil record
        ## First: find all species in need of replacement
        dict_extincts = {}
        for idxt, exttup in enumerate(y_exts):
            if len(exttup)>0:
                for spec in exttup:
                    if spec in dict_extincts.keys():
                        dict_extincts[spec] = np.append(dict_extincts[spec],[idxt,])
                    else:
                        dict_extincts[spec] = [idxt,]

        ## Then, if there are any such tubes, replace
        if bool(dict_extincts): # Boolean value of empty dict is False
            get_species_fromrecord(pars,record,dict_extincts,y0)

        #### Step 4: Populate new tubes by the selected communities
        survivalsum = np.sum(survival)
        if survivalsum>0:
            # I.e. if at least one community has survived
            topdegr = np.multiply(topdegr_prel,survival)
            probprop = np.divide(topdegr,np.sum(topdegr))
        else:
            # In the unlikely case where all selected tubes are extinct, choose new tubes uniformly
            probprop = np.divide(np.ones(Nbottleneck), Nbottleneck)

        idx_besttubes = rng.choice(toptubes_ranked,size=Ntubes,p=probprop)
        for tube_to, tube_from in enumerate(idx_besttubes):
            np.copyto(pars['specid'][tube_to],pars['specid'][tube_from]) # np.array()
            np.copyto(y0[tube_to,:],y0[tube_from,:])
            for trait in ('finv','r','m','Y','d'):
                # tube_from = idx_besttubes[tube_to]
                np.copyto(pars[trait][tube_to],pars[trait][tube_from])

        #### Step 5a: Emigrate populations
        tubes_em = rng.choice(Ntubes, replace=False, size=Nremix)
        emigrate(rng,pars,rate_migr_poisson,y0,Nsize,tubes_em)

        #### Step 5b: Immigrate from fossil record
        #tubes_withfreespots = [t for t,vals in enumerate(pars['specid']) if (-1) in vals]
        #Nwithfree = len(tubes_withfreespots)
        #if Nwithfree>=Nremix:
        #    tubes_im = rng.choice(tubes_withfreespots, replace=False, size=Nremix)
        #else:
        #    tubes_im = tubes_withfreespots

        tubes_im = draw_tubes_for_immigration(rng,pars,Ntubes,Nremix)
        if len(tubes_im)>0:
            immigrate_fossilrecord(rng,pars,rate_migr_poisson,y0,Nsize,tubes_im,record)

    return(pars,y0)

def compute_dilutions_Dis_nosel(pars,record,y0,exppars,printpars,computepars):
    #### Transfer a given set of communities by the disassembly method
    ###  INPUT
    #  'pars' is a dict of all needed parameters
    #   Nspec, Nnutr, Ntox = commpars['size']
    #     rmat = commpars['r']     # array, shape Nspec*Nnutr
    #     mmat = commpars['m']     # array, shape Nspec*Ntox
    #     fmat = commpars['finv']  # array, shape Nspec*Ntox
    #     Yvec = commpars['Y']     # array, len Nspec
    #     dvec = commpars['d']     # array, len Nspec
    #  'y0' is a list (length Ntubes) of inoculum sizes for each population
    #  'exppars' is a dictionary of experimental parameters
    #  'printpars' is a dictionary of paths and filenames for the output
    #  'computepars' is a dictionary for ODE parameters

    #### Read experimental parameters
    y0_spec            = exppars['y0spec']
    y0_nutr, y0_tox    = exppars['y0substrate']
    Nspec, Nnutr, Ntox, Ntubes  = exppars['Nsize']
    Ntransfers            = exppars['Ntransf']
    dilution, frac_bneck  = exppars['bneckdil']
    mixrate, rate_migr_poisson = exppars['mixrate']
    rng                   = exppars['RNG']

    #### Output files organised in dictionary
    fnames_dict   = printpars['fnames']
    fntup_state   = fnames_dict['fname_states']   # Steady-state values of state params S1,S2,...,N1,...,T1
    fn_mutpars    = fnames_dict['fname_mutpars']  # Parameters to mutate: only (f,) at present
    # fntup_fixpars = fnames_dict['fname_fixpars']  # Fixed parameters, printed in beginning
    fn_specid     = fnames_dict['fname_specID']
    fn_AUCs      = fnames_dict['fname_AUC']

    ## State and mutant trait filenames
    fn_popstates,fn_nutrstates,fn_toxstates = fntup_state
    # fnout_growth,fnout_mort,fnout_yield,fnout_degr = fntup_fixpars

    ### Arguments for tube community growth in step 1
    Nreg  = Nnutr+Ntox
    Nplaces = 2*Nreg

    Nsize = exppars['Nsize']
    tsfpars  = {'y0substrate':(y0_nutr, y0_tox),
                'y0spec':y0_spec,
                'dilution':dilution,
                'Nsize':Nsize}
    comppars = {'ftup':computepars['ftup'],
                'timeframe':computepars['timeframe'], # Tuple: timeframe = (Nstep, t0, tstop)
                'mutpars':computepars['mutpars'],
                'RNG':exppars['RNG']}

    ### Arguments for selection in step 2
    ytoxmax = np.sqrt(Ntox)*y0_tox
    Nbottleneck = round(Ntubes*frac_bneck)
    Nremix   = round(mixrate*Ntubes)

    #### Transfers start here!
    for transfer in range(Ntransfers):
        #### 0.0 Print community composition and species ancestral ID
        with open(fn_specid,'a') as f:
            for tube,matID in enumerate(pars['specid']):
                for sp,id in enumerate(matID):
                    IDprint = [transfer,tube,tube,sp,round(id)]
                    f.write(''.join('{0:12d}'.format(el) for el in IDprint)+'\n')
                    # f.write(''.join('{0:12.2e}'.format(id))+'\n')
                    # f.write(''.join(f'{el:12.2e}' for el in IDprint)+' \n')

        #### 0.1: Print parameter values
        #for trait, fntrt in zip(traittup,fntup_mutpars):
        #    with open(fntrt,'a') as f:
        #        for tube in range(Ntubes):
        #            parprint = np.append([transfer,tube],pars[trait][tube].flatten())
        #            f.write(''.join('{0:12.2e}'.format(el) for el in parprint)+'\n')
        #            # f.write(''.join(f'{el:12.2e}' for el in parprint)+' \n')
        #for trait, trlen, fntrt in zip(traittup,traitlenarr,fn_mutpars):
        with open(fn_mutpars,'a') as f:
            for tube,fmat in enumerate(pars['finv']):
                for sp,fvec in enumerate(fmat):
                    for tr,val in enumerate(fvec):
                        idxprint = [transfer,tube,sp,tr]
                        f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                        f.write(''.join('{0:12.2e}'.format(val))+'\n')
                        # f.write(''.join(f'{el:12d}' for el in printarr)+' \t')
                        # f.write(''.join(f'{val:12.2e}')+' \n')

        #### Step 1: Grow, mutate and compute degradation
        youts, I_AUC, y_exts, IDpertube = mutate_and_compete_dis(pars,y0,tsfpars,comppars)

        #### Intermission 1.1: Print end state concentrations
        with open(fn_popstates,'a') as f:
            for tube, yvectube in enumerate(youts[:,:Nplaces]):
                for spec, abund in enumerate(yvectube):
                    idxprint = [transfer,tube,spec]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(abund))+'\n')
                    # f.write(''.join(f'{el:12d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{abund:12.2e}')+' \n')

        with open(fn_nutrstates,'a') as f:
            for tube, yvectube in enumerate(youts[:,Nplaces:(Nplaces+Nnutr)]):
                for nutr, conc in enumerate(yvectube):
                    idxprint = [transfer,tube,nutr]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(conc))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{conc:12.2e}')+' \n')

        with open(fn_toxstates,'a') as f:
            for tube, yvectube in enumerate(youts[:,(Nplaces+Nnutr):]):
                for tox, conc in enumerate(yvectube):
                    idxprint = [transfer,tube,tox]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(conc))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{conc:12.2e}')+' \n')

        #with open(fn_states,'a') as f:
        #    for idxt, y_out in enumerate(youts):
        #        stateprint = np.append([transfer,idxt],y_out)
        #        f.write(''.join('{0:12.2e}'.format(el) for el in stateprint)+'\n')
        #        # f.write(''.join(f'{el:12.2e}' for el in stateprint)+' \n')

        #### Intermission 1.4: Print AUC
        with open(fn_AUCs,'a') as f:
            for tube, AUC_pertube in enumerate(I_AUC):
                for spec, auc in enumerate(AUC_pertube):
                    idxprint = [transfer,tube,spec]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(auc))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{auc:12.2e}')+' \n')

        #### Step 2, Old version where we follow the same steps as in selection
        ###          but on a shuffled array
        #Tdegr = [np.linalg.norm(y[-Ntox:]) for y in youts] # 2-norm emphasises abundant substrate
        #np.random.shuffle(Tdegr)
        #ranks = np.argsort(Tdegr)
        #toptubes_ranked = ranks[:Nbottleneck]

        #### Step 2: take random subset of the tubes
        toptubes_ranked = rng.choice(Ntubes,size=Nbottleneck,replace=False)
        survival = np.array([np.divide(len(IDpertube[r])-len(y_exts[r]),len(IDpertube[r])) for r in toptubes_ranked])

        #### 2a: Update the species record
        toptubes_spec = pars['specid'][toptubes_ranked]
        spec_fromtube = {}
        for tr,tb in zip(toptubes_ranked, toptubes_spec):
            specs_tube = tb[tb>=0] # Using -1 as a marker for empty slots
            specs_new  = [s for s in specs_tube if s not in spec_fromtube.keys()]
            if len(specs_new)>0:
                for s in specs_new:
                    spec_fromtube[s]=tr

        record = update_fossilrecord(pars, record, spec_fromtube, transfer, y0)

        #### Step 3: Re-introduce extinct species from the (updated) fossil record
        ## First: arrange the list of tubes with extinct species into a dictionary
        dict_extincts = {}
        for idxt, exttup in enumerate(y_exts):
            if len(exttup)>0:
                for spec in exttup:
                    if spec in dict_extincts.keys():
                        dict_extincts[spec] = np.append(dict_extincts[spec],[idxt,])
                    else:
                        dict_extincts[spec] = [idxt,]

        ## Then, if there are any such tubes, replace
        if bool(dict_extincts): # Boolean value of empty dict is False
            get_species_fromrecord(pars,record,dict_extincts,y0)

        #### Step 4: Populate new tubes by the selected communities
        survivalsum = np.sum(survival)
        if survivalsum>0:
            probprop = np.divide(survival,survivalsum)
        else:
            probprop = np.divide(np.ones(Nbottleneck), Nbottleneck)

        idx_besttubes = rng.choice(toptubes_ranked,size=Ntubes,p=probprop)
        for tube_to, tube_from in enumerate(idx_besttubes):
            #pars['specid'][tube_to,:] = pars['specid'][tube_from,:]
            #y0[tube_to,:] = y0[tube_from,:]
            np.copyto(pars['specid'][tube_to,:], pars['specid'][tube_from,:])
            np.copyto(y0[tube_to,:], y0[tube_from,:])
            for trait in ('finv','r','m','Y','d'):
                # tube_from = idx_besttubes[tube_to]
                # pars[trait][tube_to,:] = pars[trait][tube_from,:]
                np.copyto(pars[trait][tube_to,:], pars[trait][tube_from,:])

        #### Step 5: Migrate populations
        #### Step 5a: Emigrate populations
        tubes_em = rng.choice(Ntubes, replace=False, size=Nremix)
        emigrate(rng,pars,rate_migr_poisson,y0,Nsize,tubes_em)

        #### Step 5b: Immigrate from fossil record
        #tubes_withfreespots = [t for t,vals in enumerate(pars['specid']) if (-1) in vals]
        #Nwithfree = len(tubes_withfreespots)
        #if Nwithfree>=Nremix:
        #    tubes_im = rng.choice(tubes_withfreespots, replace=False, size=Nremix)
        #else:
        #    tubes_im = tubes_withfreespots

        tubes_im = draw_tubes_for_immigration(rng,pars,Ntubes,Nremix)
        if len(tubes_im)>0:
            immigrate_fossilrecord(rng,pars,rate_migr_poisson,y0,Nsize,tubes_im,record)

    return(pars,y0)


def immigrate_fossilrecord_blind(rng,pars,rate_immigr,y0,sizetup,tubes_to,record):
    #### Immigrate speces from the fossil record to the tubes 'tubes_to'
    ####   using the already defined get_species_fromrecord() that is
    ####   used in the disassembly method to re-introduce extinct species
    ####   This version is 'blind' in the sense that it introduces any species, 
    ####   without taking into account if it is already in the community or not
    ####
    ##   INPUT
    ####   rng    the random number generator
    ####   pars   dict of growth parameters.
    ####          Keys ['size', 'specid', 'finv', 'm', 'r',
    ####                'KN', 'KT', 'Y', 'd', 'nHill', 'fpow']
    ####   rate_immigr    float, parameter for Poisson(rate) to draw
    ####                  the number of introduced species per tube
    ####   y0     Array of population sizes
    ####   sizetup        tuple (Nspec, Nnutr, Ntox, Ntubes)
    ####   tubes_to       Array of tubes to receive a migrant species
    ####   record         The fossil record, defined by init_fossilrecord()

    Nspec, Nnutr, Ntox, Ntubes = sizetup
    dict_tubes_immigr = {}

    #### Step 0: Define meta-community - the collection of species present in at least one tube
    # At present, species IDs are 0, 1, 2, 3 ... while -1 is reserved for empty space
    pool_immigrants = set(record.keys()) # All species in the simulation

    #### Step 1: Add at least one species per tube
    N_immigr_vec = np.add(1,rng.poisson(rate_immigr,size=len(tubes_to)))

    #### For each 'tube-to', the communities to receive one or more migrating species
    for t,N_immigr in zip(tubes_to,N_immigr_vec):
        ## In the unlikely event that we attempt to immigrate more species than available
        N_immigr = min(N_immigr, len(pool_immigrants)) 
        
        #### Step 2: Choose immigrants
        sp_immigr = rng.choice(list(pool_immigrants), replace=False, size=N_immigr)
        ######### Debug: Draw only one immigrant
        #sp_immigr = rng.choice(list(pool_immigrants), replace=False, size=1)

        #### Step 3: Note which tubes to put immigrants in
        # Q: Also update dict_species_by_tubes?
        for spec in sp_immigr:
            if spec in dict_tubes_immigr.keys():
                dict_tubes_immigr[spec] = np.append(dict_tubes_immigr[spec],[t,])
            else:
                dict_tubes_immigr[spec] = [t,]

    #### Step 4: Introduce parameters, pop size etc for immigrants
    get_species_fromrecord( pars, record, dict_tubes_immigr, y0 )

    return()


def compute_propagule_with_immigration(pars,record,y0,exppars,printpars,computepars):
    #### Modification of the propagule method where we
    ###  immigrate a few species to some selected communities
    ###  INPUT
    #  'pars' is a dict of all needed parameters
    #   Nspec, Nnutr, Ntox = commpars['size']
    #     rmat = commpars['r']     # array, shape Nspec*Nnutr
    #     mmat = commpars['m']     # array, shape Nspec*Ntox
    #     fmat = commpars['finv']  # array, shape Nspec*Ntox
    #     Yvec = commpars['Y']     # array, len Nspec
    #     dvec = commpars['d']     # array, len Nspec
    #  'y0' is a list (length Ntubes) of inoculum sizes for each population
    #  'exppars' is a dictionary of experimental parameters
    #  'printpars' is a dictionary of paths and filenames for the output
    #  'computepars' is a dictionary for ODE parameters


    #### Read experimental parameters
    y0_nutr, y0_tox    = exppars['y0substrate']
    Nspec, Nnutr, Ntox, Ntubes  = exppars['Nsize']
    Ntransfers                  = exppars['Ntransf']
    dilution, frac_bneck        = exppars['bneckdil']
    mixrate, rate_migr_poisson  = exppars['mixrate']
    rng                         = exppars['RNG']

    #### Output files organised in dictionary
    fnames_dict   = printpars['fnames']
    fntup_state   = fnames_dict['fname_states']   # Steady-state values of state params S1,S2,...,N1,...,T1
    fn_mutpars    = fnames_dict['fname_mutpars']  # Parameters to mutate: only (f,) at present
    # fntup_fixpars = fnames_dict['fname_fixpars']  # Fixed parameters, printed in beginning
    fn_specid     = fnames_dict['fname_specID']
    fn_AUCs      = fnames_dict['fname_AUC']

    ## State and mutant trait filenames
    fn_popstates,fn_nutrstates,fn_toxstates = fntup_state
    # fnout_growth,fnout_mort,fnout_yield,fnout_degr = fntup_fixpars

    ### Arguments for tube community growth in step 1
    Nreg  = Nnutr+Ntox
    Nplaces = 2*Nreg

    tsfpars  = {'y0substrate':(y0_nutr, y0_tox),
                'dilution':dilution,
                'Nsize':exppars['Nsize']}
    comppars = {'ftup':computepars['ftup'],
                'timeframe':computepars['timeframe'], # Tuple: timeframe = (Nstep, t0, tstop)
                'mutpars':computepars['mutpars'],
                'RNG':exppars['RNG']}

    ### Arguments for selection in step 2
    Nsize = exppars['Nsize']
    Nremix   = round(mixrate*Ntubes)
    Nbottleneck = round(Ntubes*frac_bneck)
    #Ncopiesofoldtubes = int(1.0/frac_bneck)
    Ncopiesofoldtubes = Ntubes//Nbottleneck


    #### Transfers start here!
    for transfer in range(Ntransfers):
        #### Print community composition and species ancestral ID
        with open(fn_specid,'a') as f:
            for tube,matID in enumerate(pars['specid']):
                for sp,id in enumerate(matID):
                    IDprint = [transfer,tube,tube,sp,round(id)]
                    f.write(''.join('{0:12d}'.format(el) for el in IDprint)+'\n')
                    # f.write(''.join(f'{el:12.2e}' for el in IDprint)+' \n')

        #### Intermission 1.2: Print parameter values
        #for trait, fntrt in zip(traittup,fntup_mutpars):
        #    with open(fntrt,'a') as f:
        #        for tube in range(Ntubes):
        #            parprint = np.append([transfer,tube],pars[trait][tube].flatten())
        #            f.write(''.join('{0:12.2e}'.format(el) for el in parprint)+'\n')
        #            # f.write(''.join(f'{el:12.2e}' for el in parprint)+' \n')
        #for trait, trlen, fntrt in zip(traittup,traitlenarr,fn_mutpars):
        with open(fn_mutpars,'a') as f:
            for tube,fmat in enumerate(pars['finv']):
                for sp,fvec in enumerate(fmat):
                    for tr,val in enumerate(fvec):
                        idxprint = [transfer,tube,sp,tr]
                        f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                        f.write(''.join('{0:12.2e}'.format(val))+'\n')
                        # f.write(''.join(f'{el:12d}' for el in printarr)+' \t')
                        # f.write(''.join(f'{val:12.2e}')+' \n')

        #### Step 1: Grow, mutate and compute degradation
        youts, I_AUC = mutate_and_compete_fixf(pars,y0, tsfpars,comppars)

        #### Intermission 1.1: Print end state concentrations
        with open(fn_popstates,'a') as f:
            for tube, yvectube in enumerate(youts[:,:Nplaces]):
                for spec, abund in enumerate(yvectube):
                    idxprint = [transfer,tube,spec]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(abund))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{abund:12.2e}')+' \n')

        with open(fn_nutrstates,'a') as f:
            for tube, yvectube in enumerate(youts[:,Nplaces:(Nplaces+Nnutr)]):
                for nutr, conc in enumerate(yvectube):
                    idxprint = [transfer,tube,nutr]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(conc))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{conc:12.2e}')+' \n')

        with open(fn_toxstates,'a') as f:
            for tube, yvectube in enumerate(youts[:,(Nplaces+Nnutr):]):
                for tox, conc in enumerate(yvectube):
                    idxprint = [transfer,tube,tox]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(conc))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{conc:12.2e}')+' \n')

        #with open(fn_states,'a') as f:
        #    for idxt, y_out in enumerate(youts):
        #        stateprint = np.append([transfer,idxt],y_out)
        #        f.write(''.join('{0:12.2e}'.format(el) for el in stateprint)+'\n')
        #        # f.write(''.join(f'{el:12.2e}' for el in stateprint)+' \n')

        #### Intermission 1.3: Print AUC
        with open(fn_AUCs,'a') as f:
            for tube, AUC_pertube in enumerate(I_AUC):
                for spec, auc in enumerate(AUC_pertube):
                    idxprint = [transfer,tube,spec]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(auc))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{auc:12.2e}')+' \n')

        #### Step 2: Rank tubes wrt their degradation in increasing order
        Tdegr = [np.linalg.norm(y[-Ntox:]) for y in youts] # 2-norm emphasises abundant substrate
        # Tdegr = [np.max(y[-Ntox:]) for y in youts] # Max of toxins could be numerically unstable
        # Tdegr = [np.sum(y[-Ntox:]) for y in youts] # Sum does not emphasise as much

        ranks = np.argsort(Tdegr)
        toptubes_ranked = ranks[:Nbottleneck]

        #### Step 3: Re-populate parameters
        # Old version where we draw new tubes randomly (uniform)
        # Don't think any study implemented stochasticity here
        # idx_besttubes = rng.choice(toptubes_ranked,size=Ntubes)
        idx_besttubes = np.repeat(toptubes_ranked,Ncopiesofoldtubes)
        for tube_to, tube_from in enumerate(idx_besttubes):
            y0[tube_to] = np.array(y0[tube_from])
            pars['specid'][tube_to] = np.array(pars['specid'][tube_from])
            for trait in ('finv','r','m','Y','d'):
                pars[trait][tube_to] = copy.deepcopy(pars[trait][tube_from])

        #### Step 4: Immigrate some species to some (randomly chosen) communities
        tubes_im = draw_tubes_for_immigration(rng,pars,Ntubes,Nremix)
        if len(tubes_im)>0:
            immigrate_fossilrecord_blind(rng,pars,rate_migr_poisson,y0,Nsize,tubes_im,record)

    return(pars,y0)


def compute_propagule_withim_nosel(pars,record,y0,exppars,printpars,computepars):
    #### Transfer a given set of communities by a NO-SELECTION VERSION of the top-N (propagule) method
    ###  INPUT
    #  'pars' is a dict of all needed parameters
    #   Nspec, Nnutr, Ntox = commpars['size']
    #     rmat = commpars['r']     # array, shape Nspec*Nnutr
    #     mmat = commpars['m']     # array, shape Nspec*Ntox
    #     fmat = commpars['finv']  # array, shape Nspec*Ntox
    #     Yvec = commpars['Y']     # array, len Nspec
    #     dvec = commpars['d']     # array, len Nspec

    #### Read experimental parameters
    y0_nutr, y0_tox    = exppars['y0substrate']
    Nspec, Nnutr, Ntox, Ntubes  = exppars['Nsize']
    Ntransfers                  = exppars['Ntransf']
    dilution, frac_bneck        = exppars['bneckdil']
    mixrate, rate_migr_poisson  = exppars['mixrate']
    rng                         = exppars['RNG']

    #### Output files organised in dictionary
    fnames_dict   = printpars['fnames']
    fntup_state   = fnames_dict['fname_states']   # Steady-state values of state params S1,S2,...,N1,...,T1
    fn_mutpars    = fnames_dict['fname_mutpars']  # Parameters to mutate: only (f,) at present
    # fntup_fixpars = fnames_dict['fname_fixpars']  # Fixed parameters, printed in beginning
    fn_specid     = fnames_dict['fname_specID']
    fn_AUCs      = fnames_dict['fname_AUC']

    ## State and mutant trait filenames
    fn_popstates,fn_nutrstates,fn_toxstates = fntup_state
    # fnout_growth,fnout_mort,fnout_yield,fnout_degr = fntup_fixpars

    ### Arguments for tube community growth in step 1
    Nreg  = Nnutr+Ntox
    Nplaces = 2*Nreg

    tsfpars  = {'y0substrate':(y0_nutr, y0_tox),
                'dilution':dilution,
                'Nsize':exppars['Nsize']}
    comppars = {'ftup':computepars['ftup'],
                'timeframe':computepars['timeframe'], # Tuple: timeframe = (Nstep, t0, tstop)
                'mutpars':computepars['mutpars'],
                'RNG':exppars['RNG']}
    traittup = ('finv', )
    traitlenarr = (Ntox, )

    ### Arguments for selection in step 2
    Nsize = exppars['Nsize']
    Nremix   = round(mixrate*Ntubes)
    Nbottleneck = round(Ntubes*frac_bneck)
    #Ncopiesofoldtubes = int(1.0/frac_bneck)
    Ncopiesofoldtubes = Ntubes//Nbottleneck

    #### Transfers start here!
    for transfer in range(Ntransfers):
        #### Print community composition and species ancestral ID
        with open(fn_specid,'a') as f:
            for tube,matID in enumerate(pars['specid']):
                for sp,id in enumerate(matID):
                    IDprint = [transfer,tube,tube,sp,round(id)]
                    f.write(''.join('{0:12d}'.format(el) for el in IDprint)+'\n')
                    # f.write(''.join(f'{el:12.2e}' for el in IDprint)+' \n')

        #### Print parameter values
        #for trait, fntrt in zip(traittup,fntup_mutpars):
        #    with open(fntrt,'a') as f:
        #        for tube in range(Ntubes):
        #            parprint = np.append([transfer,tube],pars[trait][tube].flatten())
        #            f.write(''.join('{0:12.2e}'.format(el) for el in parprint)+'\n')
        #            # f.write(''.join(f'{el:12.2e}' for el in parprint)+' \n')
        #for trait, trlen, fntrt in zip(traittup,traitlenarr,fn_mutpars):
        with open(fn_mutpars,'a') as f:
            for tube,fmat in enumerate(pars['finv']):
                for sp,fvec in enumerate(fmat):
                    for tr,val in enumerate(fvec):
                        idxprint = [transfer,tube,sp,tr]
                        f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                        f.write(''.join('{0:12.2e}'.format(val))+'\n')
                        # f.write(''.join(f'{el:12d}' for el in printarr)+' \t')
                        # f.write(''.join(f'{val:12.2e}')+' \n')

        #### Step 1: Grow, mutate and compute degradation
        youts, I_AUC = mutate_and_compete_fixf(pars,y0,tsfpars,comppars)

        #### Intermission 1.1: Print end state concentrations
        with open(fn_popstates,'a') as f:
            for tube, yvectube in enumerate(youts[:,:Nplaces]):
                for spec, abund in enumerate(yvectube):
                    idxprint = [transfer,tube,spec]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(abund))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{abund:12.2e}')+' \n')

        with open(fn_nutrstates,'a') as f:
            for tube, yvectube in enumerate(youts[:,Nplaces:(Nplaces+Nnutr)]):
                for nutr, conc in enumerate(yvectube):
                    idxprint = [transfer,tube,nutr]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(conc))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{conc:12.2e}')+' \n')

        with open(fn_toxstates,'a') as f:
            for tube, yvectube in enumerate(youts[:,(Nplaces+Nnutr):]):
                for tox, conc in enumerate(yvectube):
                    idxprint = [transfer,tube,tox]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(conc))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{conc:12.2e}')+' \n')

        #with open(fn_states,'a') as f:
        #    for idxt, y_out in enumerate(youts):
        #        stateprint = np.append([transfer,idxt],y_out)
        #        f.write(''.join('{0:12.2e}'.format(el) for el in stateprint)+'\n')
        #        # f.write(''.join(f'{el:12.2e}' for el in stateprint)+' \n')

        #### Intermission 1.3: Print AUC
        with open(fn_AUCs,'a') as f:
            for tube, AUC_pertube in enumerate(I_AUC):
                for spec, auc in enumerate(AUC_pertube):
                    idxprint = [transfer,tube,spec]
                    f.write(''.join('{0:12d}'.format(el) for el in idxprint)+'\t')
                    f.write(''.join('{0:12.2e}'.format(auc))+'\n')
                    # f.write(''.join(f'{el:d}' for el in idxprint)+' \t')
                    # f.write(''.join(f'{auc:12.2e}')+' \n')

        #### Step 2, old version: Rank tubes wrt their degradation in increasing order
        ### Previously, we would perform same step as in propagule selection,
        ##  but on a shuffled array
        # Tdegr = [np.linalg.norm(y[-Ntox:]) for y in youts] # 2-norm emphasises abundant substrate
        # np.random.shuffle(Tdegr)
        # ranks = np.argsort(Tdegr)
        # toptubes_ranked = ranks[:Nbottleneck] # Choose from shuffled top ranks

        ### Step 2: New tubes by uniformly random draws
        toptubes_ranked = rng.choice(Ntubes,size=Nbottleneck,replace=False)

        #### Step 3: Re-populate parameters
        # Old version where we drew new tubes randomly (uniform)
        # Removed since I don't think any study implemented stochasticity here
        # idx_besttubes = np.random.choice(toptubes_ranked,size=Ntubes)
        idx_besttubes = np.repeat(toptubes_ranked, Ncopiesofoldtubes)
        for tube_to, tube_from in enumerate(idx_besttubes):
            y0[tube_to] = np.array(y0[tube_from])
            pars['specid'][tube_to] = np.array(pars['specid'][tube_from])
            for trait in ('finv','r','m','Y','d'):
                pars[trait][tube_to] = copy.deepcopy(pars[trait][tube_from])

        #### Step 4: Immigrate some species to some (randomly chosen) communities
        tubes_im = draw_tubes_for_immigration(rng,pars,Ntubes,Nremix)
        if len(tubes_im)>0:
            immigrate_fossilrecord_blind(rng,pars,rate_migr_poisson,y0,Nsize,tubes_im,record)

    return(pars,y0)
