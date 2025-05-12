import sys
import copy
import scipy
# import pandas as pd
import numpy as np

##########################################################################
#### Hill functions for batch and chemostat
##########################################################################

def f_Hillf(x,args):
    #### Hill function sum( f*x^n/(x^n +K^n) )
    fmat  = args["fmax"]   # Max substrate effect (max growth rate or max death rate)
    Khalf = args["Khalf"]  # Half-max coefficient for substrate
    cHill = args["n"]      # Hill coeff for toxic compound(s)

    f_num   = np.power(x,cHill)
    f_denom = np.add(f_num,np.power(Khalf,cHill))
    f_div   = np.divide(f_num, f_denom)

    fout    = np.dot(fmat,f_div)
    return( fout )

def f_Hillf_mono(x,args):
    #### Alias not used anymore
    return( f_Hillf(x,args) )

def f_Hillf_1Nutr(x,args):
    #### Alias not used anymore
    return( f_Hillf(x,args) )

def f_HillfNotDot(x,args):
    #### The per-substrate version of the Hill function, with no dot product
    ##    used to calculate per-compound effect on  with degradation
    fmat  = args["fmax"]   # Max substrate effect (max growth rate or max death rate)
    Khalf = args["Khalf"]  # Half-max coefficient for substrate
    cHill = args["n"]      # Hill coeff for toxic compound(s)

    f_num   = np.power(x,cHill)
    f_denom = np.add(f_num,np.power(Khalf,cHill))
    f_div   = np.divide(f_num, f_denom)
    return( np.multiply(fmat,f_div) )

def f_HillfNotDot_mono(x,args):
    #### Old alias, not used at present
    return( f_HillfNotDot(x,args) )

def f_Monodf(x,args):
    #### Hill function sum( f*x^n/(x^n +K^n) )
    fmat  = args["fmax"]   # Max substrate effect (max growth rate or max death rate)
    Khalf = args["Khalf"]  # Half-max coefficient for substrate

    f_denom = np.add(x, Khalf)
    f_div   = np.divide(x, f_denom)
    fout    = np.dot(fmat,f_div)
    return( fout )

def dfHilldx(x,args):
    #### Differentiation with respect to substrate x, of Hill function sum( f*x^n/(x^n +K^n) )
    fmat  = args["fmax"]   # Max substrate effect (max growth rate or max death rate)
    Khalf = args["Khalf"]  # Half-max coefficient for substrate
    cHill = args["n"]      # Hill coeff for toxic compound(s)

    Nspec, Nchem = fmat.shape # Output will be of same shape

    KpowHill = np.power(Khalf, cHill)
    XpowHill = np.power(x, cHill)
    dHill_denom = np.add(KpowHill, XpowHill)
    dHill_ratio = np.divide(KpowHill, np.power(dHill_denom,2.0))          # Length Nchem
    nXHill      = np.multiply(cHill, np.power(x, np.subtract(cHill,1.0))) # Length Nchem
    dHill_tot   = np.multiply(dHill_ratio,nXHill)                         # Length Nchem

    dHill_repS = np.reshape(np.tile(dHill_tot, Nspec),(Nspec, Nchem))
    return( np.multiply(fmat,dHill_repS) )

##########################################################################
#### Per-capita growth functions
##########################################################################

#### Lotka--Volterra growth rate
def f_grpc_gLV(t, y, args):
    #### Input
    # y       nx1 array
    # args    dictionary with keywords
    #          'rvec',  nx1 array of mono-culture growth
    #          'alpha', nxn interaction matrix
    #### Output
    # Lotka-Volterra
    #   dy/dt = (r+Ay)*y
    A_mat = args["Amat"]
    rvec  = args["rvec"]

    gr_percap = np.add(rvec, np.dot(A_mat, y))
    return( gr_percap )

#### Consumer--resource growth rate
def f_grpc_CR(t, x, args):

    ### Max growth rate 'r' and half-max constant 'K_N'
    rmat = args["r"]    # Growth effect of nutrient, matrix of size [Npops, Nnutr]
    KvecN = args["KN"]  # Half-max degradation rate of nutrient substrates, array of length [Nnutr]

    # Groth due to substrate uptake
    argsN = {'fmax':rmat, 'Khalf':KvecN, 'n':np.ones(len(KvecN))}
    dS_grow = f_Hillf(x,argsN) # Growth    rcx per species
    return( dS_grow )

#### Per-capita death rates, if considered in models
##########################################################################
#### Mortality due to toxic compound
def f_mpc_CR(t, x, args):

    ### Max death rate 'm', half-max constant 'KT', Hill coefficients
    mmat  = args["m"]    # Max death rate of toxic compound, matrix of size [Npops, Ntox]
    KvecT = args["KT"]   # Half-max constant for toxic substrates, array of length [Ntox]
    vecHill = args["nHill"] # Hill coefficient for toxic substrates (subject to change)

    # Death due to toxic compounds
    argsT = {'fmax':mmat, 'Khalf':KvecT, 'n':vecHill}
    dS_mort = f_Hillf(x,argsT) # Growth    rcx per species
    return( dS_mort )

##########################################################################
#### Cost of tradeoff in growth vs degradation
def f_tradeoff(t, x, args):
    ### Per-capita (and per-substrate) cost of growth vs degradation
    ### Max investment 'f', tradeoff constants
    fmat  = args["f"]       # Max investment into degradation [Npops, Ntox]
    u, v  = args["fpowtup"] # Tradeoff constants

    # Tradeoff (1-sum(f^u))^v
    fpowu = np.power(fmat, u)
    fpow_sum = np.sum(fpowu, axis=1)
    fpow_1minus = np.subtract(1.0, fpow_sum)
    f_tradeoff  = np.power(fpow_1minus, v)
    return( f_tradeoff )

def f_percap(t, x_subs, args, funtup):

    ### Collect size of problem and the relevant functions
    Nspec, Nnutr, Ntox = args['size']
    f_grow, f_mort, f_cost = funtup

    ### Collect states
    x_spec = x[:Nspec]
    x_nutr = x[Nspec:(Nspec+Nnutr)]
    x_tox  = x[(Nspec+Nnutr):]

    ### Growth and cost of growth
    dS_gr = f_grow(t, x_nutr, args)
    dS_costly = np.multiply(f_cost(t, x, args), dS_gr)

    ### Death rate (if any)
    dS_mort = f_mort(t, x_tox, args)

    ### Net growth per capita and in total
    dS_percap = np.subtract(dS_costly, dS_mort)
    dS = np.multiply(dS_percap,x_spec)
    return( dS )

##########################################################################
#### BATCH CULTURE MODELS
##########################################################################
#### Consumer-Resource with Hill functions
## Hill functions for growth and death, uptake of nutrient

def f_CRbatchHill_full(t, x, args):
    #### INPUT
    ##  t       the time step of the ODE solver. Not explicitly used.
    ##  x       the ODE state at time t
    ##  args    dictionary with all model parameters as well as community size
    ##  funtup  tuple of function handles (f_growthrate, f_deathrate, f_tradeoff)
    #### OUTPUT
    ##  dX      the change in community state from the present time step to the next

    Npops, Nnutr, Ntox = args['size']

    rmat = args["r"]    # Growth effect of nutrient, matrix of size [Npops, Nnutr]
    mmat = args["m"]    # Mortality of toxins, matrix of size [Npops, Ntox]
    fmat = args["finv"] # Per-species investment [0.0, 1.0] into degradation of each toxin
    fpowu, fpowv = args["fpow"] # Tuple of tradeoff constants (u,v) such that
                        #   the fraction f^uv is invested into degradation of the toxin
                        #   the fraction (1-f^u)^v is invested into growth
    KvecN = args["KN"]  # Half-max degradation rate of nutrient substrates, array of length [Nnutr]
    KvecT = args["KT"]  # Half-max degradation rate of toxic substrates, array of length [Ntox]
    vecHill = args["nHill"] # Hill coefficient for toxic substrates (subject to change)
    Yvec = args["Y"]    # Nutrient Yield [say gram of biomass per mol of nutrient]
    dvec = args["d"]    # Toxin degradation rate [say mol of toxin per gram of biomass]

    ### State variables for abundance and concentrations
    x_spec = x[:Npops]
    c_nutr = x[Npops:(Npops+Nnutr)]
    c_tox  = x[(Npops+Nnutr):]

    ### Parameter combinations
    fmatpow = np.power(fmat,fpowu)
    ftot_spec = np.sum(fmatpow,axis=1)
    fpowuvtr = np.transpose(np.power(fmatpow,fpowv))
    rtr = np.transpose(rmat)

    ## Change in abundance due to growth and mortality
    ## The net growth per capita depends on the model structure, whether ...
    ##    - whether death rates are considered
    ##    - whether there is a tradeoff between investments into growth and degradation
    # dS_percap = f_percap(t, x[Npops:], args, funtup)
    # dS = np.multiply(dS_percap, x_spec)
    ## The above considerations are determined by the choice of growth functions in funtup
    argsN = {'fmax':rmat, 'Khalf':KvecN}
    argsT = {'fmax':mmat, 'Khalf':KvecT, 'n':vecHill}
    dS_grow = np.multiply(f_Monodf(c_nutr, argsN), x_spec) # Growth    rcx per species
    dS_mort = np.multiply(f_Hillf(c_tox, argsT), x_spec) # Mortality sum(mcx) per species
    invest_gr = np.power(np.subtract([1.0]*Npops,ftot_spec), fpowv)
    dS = np.subtract(np.multiply(invest_gr,dS_grow), dS_mort)

    # Change in nutrient concentration
    SY  = np.multiply(x_spec,Yvec)
    rSY    = np.dot(rtr, SY)
    nutrMonod = np.divide(c_nutr,np.add(c_nutr, KvecN))
    dN_nutr = np.negative(np.multiply(rSY, nutrMonod))

    # Change in toxin concentrations
    dT_invest = np.multiply(dvec,dS_grow)
    dT_pertox = np.dot(fpowuvtr,dT_invest)
    dT_tox    = np.negative(np.multiply(c_tox, dT_pertox))

    # Assemble state vector
    dX = np.append(np.append(dS,dN_nutr),dT_tox)
    return( dX )

def f_notroff_death(t, x, args):
    #### INPUT
    ##  t       the time step of the ODE solver. Not explicitly used.
    ##  x       the ODE state at time t
    ##  args    dictionary with all model parameters as well as community size
    ##  funtup  tuple of function handles (f_growthrate, f_deathrate, f_tradeoff)
    #### OUTPUT
    ##  dX      the change in community state from the present time step to the next

    Npops, Nnutr, Ntox = args['size']

    rmat = args["r"]    # Growth effect of nutrient, matrix of size [Npops, Nnutr]
    mmat = args["m"]    # Mortality of toxins, matrix of size [Npops, Ntox]
    fmat = args["finv"] # Per-species investment [0.0, 1.0] into degradation of each toxin
    fpowu, fpowv = args["fpow"] # Tuple of tradeoff constants (u,v) such that
                        #   the fraction f^uv is invested into degradation of the toxin
                        #   the fraction (1-f^u)^v is invested into growth
    KvecN = args["KN"]  # Half-max degradation rate of nutrient substrates, array of length [Nnutr]
    KvecT = args["KT"]  # Half-max degradation rate of toxic substrates, array of length [Ntox]
    vecHill = args["nHill"] # Hill coefficient for toxic substrates (subject to change)
    Yvec = args["Y"]    # Nutrient Yield [say gram of biomass per mol of nutrient]
    dvec = args["d"]    # Toxin degradation rate [say mol of toxin per gram of biomass]

    ### State variables for abundance and concentrations
    x_spec = x[:Npops]
    c_nutr = x[Npops:(Npops+Nnutr)]
    c_tox  = x[(Npops+Nnutr):]

    ### Parameter combinations
    fmatpow = np.power(fmat,fpowu)
    ftot_spec = np.sum(fmatpow,axis=1)
    fpowuvtr = np.transpose(np.power(fmatpow,fpowv))
    rtr = np.transpose(rmat)

    ## Change in abundance due to growth and mortality
    ## The net growth per capita depends on the model structure, whether ...
    ##    - whether death rates are considered
    ##    - whether there is a tradeoff between investments into growth and degradation
    # dS_percap = f_percap(t, x[Npops:], args, funtup)
    # dS = np.multiply(dS_percap, x_spec)
    ## The above considerations are determined by the choice of growth functions in funtup
    argsN = {'fmax':rmat, 'Khalf':KvecN}
    argsT = {'fmax':mmat, 'Khalf':KvecT, 'n':vecHill}
    dS_grow = np.multiply(f_Monodf(c_nutr, argsN), x_spec) # Growth    rcx per species
    dS_mort = np.multiply(f_Hillf(c_tox, argsT), x_spec) # Mortality sum(mcx) per species
    dS = np.subtract(dS_grow, dS_mort)

    # Change in nutrient concentration
    SY  = np.multiply(x_spec,Yvec)
    rSY    = np.dot(rtr, SY)
    nutrMonod = np.divide(c_nutr,np.add(c_nutr, KvecN))
    dN_nutr = np.negative(np.multiply(rSY, nutrMonod))

    # Change in toxin concentrations
    dT_invest = np.multiply(dvec,dS_grow)
    dT_pertox = np.dot(fpowuvtr,dT_invest)
    dT_tox    = np.negative(np.multiply(c_tox, dT_pertox))

    # Assemble state vector
    dX = np.append(np.append(dS,dN_nutr),dT_tox)
    return( dX )

def f_troff_nodeath(t, x, args):
    #### INPUT
    ##  t       the time step of the ODE solver. Not explicitly used.
    ##  x       the ODE state at time t
    ##  args    dictionary with all model parameters as well as community size
    ##  funtup  tuple of function handles (f_growthrate, f_deathrate, f_tradeoff)
    #### OUTPUT
    ##  dX      the change in community state from the present time step to the next

    Npops, Nnutr, Ntox = args['size']

    rmat = args["r"]    # Growth effect of nutrient, matrix of size [Npops, Nnutr]
    fmat = args["finv"] # Per-species investment [0.0, 1.0] into degradation of each toxin
    fpowu, fpowv = args["fpow"] # Tuple of tradeoff constants (u,v) such that
                        #   the fraction f^uv is invested into degradation of the toxin
                        #   the fraction (1-f^u)^v is invested into growth
    KvecN = args["KN"]  # Half-max degradation rate of nutrient substrates, array of length [Nnutr]
    Yvec = args["Y"]    # Nutrient Yield [say gram of biomass per mol of nutrient]
    dvec = args["d"]    # Toxin degradation rate [say mol of toxin per gram of biomass]

    ### State variables for abundance and concentrations
    x_spec = x[:Npops]
    c_nutr = x[Npops:(Npops+Nnutr)]
    c_tox  = x[(Npops+Nnutr):]

    ### Parameter combinations
    fmatpow = np.power(fmat,fpowu)
    ftot_spec = np.sum(fmatpow,axis=1)
    fpowuvtr = np.transpose(np.power(fmatpow,fpowv))
    rtr = np.transpose(rmat)

    ## Change in abundance due to growth and mortality
    ## The net growth per capita depends on the model structure, whether ...
    ##    - whether death rates are considered
    ##    - whether there is a tradeoff between investments into growth and degradation
    # dS_percap = f_percap(t, x[Npops:], args, funtup)
    # dS = np.multiply(dS_percap, x_spec)
    ## The above considerations are determined by the choice of growth functions in funtup
    argsN = {'fmax':rmat, 'Khalf':KvecN}
    dS_grow = np.multiply(f_Monodf(c_nutr, argsN), x_spec) # Growth    rcx per species
    invest_gr = np.subtract([1.0]*Npops,ftot_spec)
    dS = np.multiply(invest_gr,dS_grow)

    # Change in nutrient concentration
    SY  = np.multiply(x_spec,Yvec)
    rSY    = np.dot(rtr, SY)
    nutrMonod = np.divide(c_nutr,np.add(c_nutr, KvecN))
    dN_nutr = np.negative(np.multiply(rSY, nutrMonod))

    # Change in toxin concentrations
    dT_invest = np.multiply(dvec,dS_grow)
    dT_pertox = np.dot(fpowuvtr,dT_invest)
    dT_tox    = np.negative(np.multiply(c_tox, dT_pertox))

    # Assemble state vector
    dX = np.append(np.append(dS,dN_nutr),dT_tox)
    return( dX )

def f_notroff_nodeath(t, x, args):
    #### INPUT
    ##  t       the time step of the ODE solver. Not explicitly used.
    ##  x       the ODE state at time t
    ##  args    dictionary with all model parameters as well as community size
    ##  funtup  tuple of function handles (f_growthrate, f_deathrate, f_tradeoff)
    #### OUTPUT
    ##  dX      the change in community state from the present time step to the next

    Npops, Nnutr, Ntox = args['size']

    rmat = args["r"]    # Growth effect of nutrient, matrix of size [Npops, Nnutr]
    fmat = args["finv"] # Per-species investment [0.0, 1.0] into degradation of each toxin
    fpowu, fpowv = args["fpow"] # Tuple of tradeoff constants (u,v) such that
                        #   the fraction f^uv is invested into degradation of the toxin
                        #   the fraction (1-f^u)^v is invested into growth
    KvecN = args["KN"]  # Half-max degradation rate of nutrient substrates, array of length [Nnutr]
    Yvec = args["Y"]    # Nutrient Yield [say gram of biomass per mol of nutrient]
    dvec = args["d"]    # Toxin degradation rate [say mol of toxin per gram of biomass]

    ### State variables for abundance and concentrations
    x_spec = x[:Npops]
    c_nutr = x[Npops:(Npops+Nnutr)]
    c_tox  = x[(Npops+Nnutr):]

    ### Parameter combinations
    fmatpow = np.power(fmat,fpowu)
    fpowuvtr = np.transpose(np.power(fmatpow,fpowv))
    rtr = np.transpose(rmat)

    ## Change in abundance due to growth and mortality
    ## The net growth per capita depends on the model structure, whether ...
    ##    - whether death rates are considered
    ##    - whether there is a tradeoff between investments into growth and degradation
    # dS_percap = f_percap(t, x[Npops:], args, funtup)
    # dS = np.multiply(dS_percap, x_spec)
    ## The above considerations are determined by the choice of growth functions in funtup
    argsN = {'fmax':rmat, 'Khalf':KvecN}
    dS = np.multiply(f_Monodf(c_nutr, argsN), x_spec) # Growth    rcx per species

    # Change in nutrient concentration
    SY  = np.multiply(x_spec,Yvec)
    rSY    = np.dot(rtr, SY)
    nutrMonod = np.divide(c_nutr,np.add(c_nutr, KvecN))
    dN_nutr = np.negative(np.multiply(rSY, nutrMonod))

    # Change in toxin concentrations
    dT_invest = np.multiply(dvec,dS)
    dT_pertox = np.dot(fpowuvtr,dT_invest)
    dT_tox    = np.negative(np.multiply(c_tox, dT_pertox))

    # Assemble state vector
    dX = np.append(np.append(dS,dN_nutr),dT_tox)
    return( dX )

def f_CRbatchHill(t, x, args):
    Npops, Nnutr, Ntox = args['size']

    rmat = args["r"]    # Growth effect of nutrient, matrix of size [Npops, Nnutr]
    mmat = args["m"]    # Mortality of toxins, matrix of size [Npops, Ntox]
    fmat = args["finv"] # Per-species investment [0.0, 1.0] into degradation of each toxin
    fpow = args["fpow"] # Per-species investment [0.0, 1.0] into degradation of each toxin
    KvecN = args["KN"]  # Half-max degradation rate of nutrient substrates, array of length [Nnutr]
    KvecT = args["KT"]  # Half-max degradation rate of toxic substrates, array of length [Ntox]
    vecHill = args["nHill"] # Hill coefficient for toxic substrates (subject to change)
    Yvec = args["Y"]    # Nutrient Yield [say gram of biomass per mol of nutrient]
    dvec = args["d"]    # Toxin degradation rate [say mol of toxin per gram of biomass]

    ### State variables for abundance and concentrations
    x_spec = x[:Npops]
    c_nutr = x[Npops:(Npops+Nnutr)]
    c_tox  = x[(Npops+Nnutr):]

    ### Parameter combinations
    fmatpow = np.power(fmat,fpow)
    ftot_spec = np.sum(fmatpow,axis=1)
    ftr = np.transpose(fmatpow)
    rtr = np.transpose(rmat)

    # Change in abundance due to growth and mortality
    argsN = {'fmax':rmat, 'Khalf':KvecN, 'n':np.ones(len(KvecN))}
    argsT = {'fmax':mmat, 'Khalf':KvecT, 'n':vecHill}
    dS_grow = np.multiply(f_Hillf(c_nutr,argsN), x_spec) # Growth    rcx per species
    dS_mort = np.multiply(f_Hillf(c_tox, argsT), x_spec) # Mortality sum(mcx) per species
    invest_gr = np.subtract([1.0]*Npops,ftot_spec)
    dS = np.subtract(np.multiply(invest_gr,dS_grow), dS_mort)

    # Change in nutrient concentration
    SY  = np.multiply(x_spec,Yvec)
    rSY    = np.dot(rtr, SY)
    nutrMonod = np.divide(c_nutr,np.add(c_nutr, KvecN))
    dN_nutr = np.negative(np.multiply(rSY, nutrMonod))

    # Change in toxin concentrations
    dT_invest = np.multiply(dvec,dS_grow)
    dT_pertox = np.dot(ftr,dT_invest)
    dT_tox    = np.negative(np.multiply(c_tox, dT_pertox))
    ### MIGHT WANT TO CHANGE ALSO dT

    # Assemble state vector
    dX = np.append(np.append(dS,dN_nutr),dT_tox)
    return( dX )

def jac_CRbatchHill(t, x, args):
    Nspec, Nnutr, Ntox = args['size']
    # N0, T0 = args['nt_inflow']

    rmat = args["r"]    # Growth effect of nutrient, matrix of size [Nspec, Nnutr]
    mmat = args["m"]    # Mortality of toxins, matrix of size [Nspec, Ntox]
    fmat = args["finv"] # Per-species investment [0.0, 1.0] into degradation of each toxin
    fpow = args["fpow"] # Per-species investment [0.0, 1.0] into degradation of each toxin
    KvecN = args["KN"]  # Half-max degradation rate of nutrient substrates, array of length [Nnutr]
    KvecT = args["KT"]  # Half-max degradation rate of toxic substrates, array of length [Ntox]
    vecHill = args["nHill"] # Hill coefficient for toxic substrates (subject to change)
    Yvec = args["Y"]        # Nutrient Yield [say gram of biomass per mol of nutrient]
    dvec = args["d"]        # Toxin degradation rate [say mol of toxin per gram of biomass]

    ### State variables for abundance and concentrations
    x_spec = x[:Nspec]
    c_nutr = x[Nspec:(Nspec+Nnutr)]
    c_tox  = x[(Nspec+Nnutr):]

    #### Prelims
    fmatpow   = np.power(fmat,fpow)
    ftot_spec = np.sum(fmatpow,axis=1)
    ftr = np.transpose(fmatpow)
    rtr = np.transpose(rmat)

    argsN = {'fmax':rmat, 'Khalf':KvecN, 'n':np.ones(len(KvecN))}
    argsT = {'fmax':mmat, 'Khalf':KvecT, 'n':vecHill}
    gr_percap = f_Hillf(c_nutr,argsN) # np.dot(rmat,c_nutr)
    d_percap  = f_Hillf(c_tox,argsT)  # np.dot(mmat,c_tox)
    invest_gr = np.subtract([1.0]*Nspec,ftot_spec)


    #### Compute partial derivatives wrt species abundance
    ## dSdS simple
    dSdS = np.subtract(np.multiply(invest_gr,gr_percap),d_percap)  # -> diag
    ## dSdN a bit more complicated
    XF = np.multiply(invest_gr,x_spec)
    xf_repN = np.transpose(np.reshape([XF]*Nnutr,[Nnutr,Nspec]))
    dHillN_denom = np.power(np.add(c_nutr,KvecN),2.0)
    dHillN = np.divide(KvecN,dHillN_denom)
    dHillN_repS = np.reshape([dHillN]*Nspec,[Nspec,Nnutr])
    rNmut = np.multiply(rmat,dHillN_repS)
    dSdN  = np.multiply(rNmut,xf_repN)
    ## dSdT more complicated
    x_repT  = np.transpose(np.reshape([x_spec]*Ntox,[Ntox,Nspec]))
    KT_powHill = np.power(KvecT,vecHill)
    dHillT_denom = np.add(KT_powHill,np.power(c_tox,vecHill))
    dHillT_div = np.divide(KT_powHill,np.power(dHillT_denom,2.0))
    diffTHill  = np.multiply(vecHill, np.power(c_tox, np.subtract(vecHill,1.0)))
    dTHill_tot = np.multiply(diffTHill,dHillT_div)
    dHillT_repS = np.reshape([dTHill_tot]*Nspec,[Nspec,Ntox])
    mTmut = np.multiply(mmat,dHillT_repS)
    dSdT = np.negative(np.multiply(mTmut,x_repT))

    #### Effect on nutrient dynamics
    ## dNdS
    Y_repN = np.reshape([Yvec]*Nnutr,[Nnutr,Nspec])
    rY = np.multiply(rtr,Y_repN)
    nutrMonod = np.divide(c_nutr,np.add(c_nutr, KvecN))
    CNrepX = np.transpose(np.reshape([nutrMonod]*Nspec,[Nspec,Nnutr]))
    dNdS = np.negative(np.multiply(rY, CNrepX))      # array [Nnutr,Ntox]
    ## dNdN
    XY = np.multiply(Yvec,x_spec)
    rXY = np.dot(rtr,XY)
    diffMonod = np.negative(np.divide(KvecN,np.power(np.add(c_nutr, KvecN), 2.0)))
    dNdN = np.negative(np.multiply(rXY,diffMonod))   # array [Nnutr] -> diag
    ## dNdT
    dNdT = np.zeros((Nnutr,Ntox))                    # array [Nnutr,Ntox]

    ### Effect on toxin
    diffMonod_rep = np.reshape([diffMonod]*Nspec,[Nspec,Nnutr])
    rdiffM = np.multiply(rmat,diffMonod_rep)
    dXvec  = np.multiply(dvec,x_spec) # array [Nspec]
    dX_rep = np.transpose(np.reshape([dXvec]*Nnutr,[Nnutr,Nspec]))  # array [Nspec, Nnutr]
    rdX = np.multiply(dX_rep,rdiffM)  # array [Nspec, Nnutr]
    drN = np.multiply(dvec,gr_percap) # array [Nspec]
    rdXN = np.multiply(drN,x_spec)    # array [Nspec]
    drN_repT = np.reshape([drN]*Ntox,[Ntox,Nspec])                  # array [Ntox, Nspec]

    ctox_repX = np.transpose(np.reshape([c_tox]*Nspec,[Nspec,Ntox]))
    ctox_repN = np.transpose(np.reshape([c_tox]*Nnutr,[Nnutr,Ntox]))
    dTdS = np.negative(np.multiply(ctox_repX, np.multiply(ftr,drN_repT)))  # array [Ntox, Nspec]
    dTdN = np.negative(np.multiply(ctox_repN, np.dot(ftr,rdX)))            # array [Ntox, Nnutr]
    dTdT = np.negative(np.dot(ftr,rdXN))                                   # array [Ntox] -> diag

    # Assemble Jacobian from partial derivatives
    dSdNT = np.hstack((dSdN,dSdT))
    dS    = np.hstack((np.diag(dSdS),dSdNT))
    dNdNT = np.hstack((np.diag(dNdN),dNdT))
    dN    = np.hstack((dNdS,dNdNT))
    dTdNT = np.hstack((dTdN,np.diag(dTdT)))
    dT    = np.hstack((dTdS,dTdNT))
    J_out = np.vstack((dS,dN,dT))
    return J_out

##########################################################################
#### Consumer-Resource model of batch culture with linear growth rates
def f_CRbatch_lingr(t,x, args):
    Npops, Nnutr, Ntox = args['size']

    rmat = args["r"]    # Growth effect of nutrient, matrix of size [Npops, Nnutr]
    mmat = args["m"]    # Mortality of toxins, matrix of size [Npops, Ntox]
    fmat = args["finv"] # Per-species investment [0.0, 1.0] into degradation of each toxin
    Yvec = args["Y"]    # Nutrient Yield [say gram of biomass per mol of nutrient]
    dvec = args["d"]    # Toxin degradation rate [say mol of toxin per gram of biomass]

    ### State variables for abundance and concentrations
    x_spec = x[:Npops]
    c_nutr = x[Npops:(Npops+Nnutr)]
    c_tox  = x[(Npops+Nnutr):]

    ### Parameter combinations
    ftot_spec = np.sum(fmat,axis=1)
    ftr = np.transpose(fmat)
    rtr = np.transpose(rmat)

    # Change in abundance due to growth and mortality
    dS_grow = np.multiply(np.dot(rmat,c_nutr),x_spec) # Growth    rcx per species
    dS_mort = np.multiply(np.dot(mmat,c_tox), x_spec) # Mortality sum(mcx) per species
    invest_gr = np.subtract([1.0]*Npops,ftot_spec)
    dS = np.subtract(np.multiply(invest_gr,dS_grow), dS_mort)

    # Change in nutrient concentration
    SY  = np.multiply(x_spec,Yvec)
    rSY   = np.dot(rtr, SY)
    dN_nutr = np.negative(np.multiply(rSY, c_nutr))

    # Change in toxin concentrations
    dT_invest = np.multiply(dvec,dS_grow)
    dT_pertox = np.dot(ftr,dT_invest)
    dT_tox    = np.negative(np.multiply(c_tox, dT_pertox))

    # Assemble state vector
    dX = np.append(np.append(dS,dN_nutr),dT_tox)
    return dX

def jac_CRbatch_lingr(t, x, args):
    Nspec, Nnutr, Ntox = args['size']
    N0, T0 = args['nt_inflow']

    rmat = args["r"]    # Growth effect of nutrient, matrix of size [Nspec, Nnutr]
    mmat = args["m"]    # Mortality of toxins, matrix of size [Nspec, Ntox]
    fmat = args["finv"] # Per-species investment [0.0, 1.0] into degradation of each toxin
    Yvec = args["Y"]    # Nutrient Yield [say gram of biomass per mol of nutrient]
    dvec = args["d"]    # Toxin degradation rate [say mol of toxin per gram of biomass]

    ### State variables for abundance and concentrations
    x_spec = x[:Nspec]
    c_nutr = x[Nspec:(Nspec+Nnutr)]
    c_tox  = x[(Nspec+Nnutr):]

    #### Prelims
    gr_percap = np.dot(rmat,c_nutr)
    d_percap  = np.dot(mmat,c_tox)
    invest_gr = np.subtract([1.0]*Nspec,np.sum(fmat, axis=1))
    ftr = np.transpose(fmat)
    rtr = np.transpose(rmat)

    #### Compute partial derivatives wrt species abundance
    XF = np.multiply(invest_gr,x_spec)
    xf_repN = np.transpose(np.reshape([XF]*Nnutr,[Nnutr,Nspec]))
    x_repT  = np.transpose(np.reshape([x_spec]*Ntox,[Ntox,Nspec]))
    dSdS = np.subtract(np.multiply(invest_gr,gr_percap),d_percap)  # -> diag
    dSdN = np.multiply(rmat,xf_repN)
    dSdT = np.negative(np.multiply(mmat,x_repT))

    ### Effect on nutrient dynamics
    Y_repN = np.reshape([Yvec]*Nnutr,[Nnutr,Nspec])
    rY = np.multiply(rtr,Y_repN)
    XY = np.multiply(Yvec,x_spec)
    CNrepX = np.transpose(np.reshape([c_nutr]*Nspec,[Nspec,Nnutr]))
    dNdS = np.negative(np.multiply(rY, CNrepX))      # array [Nnutr,Ntox]
    dNdN = np.negative(np.dot(rtr,XY))               # array [Nnutr] -> diag
    dNdT = np.zeros((Nnutr,Ntox))                    # array [Nnutr,Ntox]

    ### Effect on toxin
    dXvec  = np.multiply(dvec,x_spec) # array [Nspec]
    dX_rep = np.transpose(np.reshape([dXvec]*Nnutr,[Nnutr,Nspec]))  # array [Nspec, Nnutr]
    rdX = np.multiply(dX_rep,rmat)    # array [Nspec, Nnutr]
    drN = np.multiply(dvec,gr_percap) # array [Nspec]
    rdXN = np.multiply(drN,x_spec)    # array [Nspec]
    drN_repT = np.reshape([drN]*Ntox,[Ntox,Nspec])                  # array [Ntox, Nspec]

    ctox_repX = np.transpose(np.reshape([c_tox]*Nspec,[Nspec,Ntox]))
    ctox_repN = np.transpose(np.reshape([c_tox]*Nnutr,[Nnutr,Ntox]))
    dTdS = np.negative(np.multiply(ctox_repX, np.multiply(ftr,drN_repT)))  # array [Ntox, Nspec]
    dTdN = np.negative(np.multiply(ctox_repN, np.dot(ftr,rdX)))            # array [Ntox, Nnutr]
    dTdT = np.negative(np.dot(ftr,rdXN))                                   # array [Ntox] -> diag

    # Assemble Jacobian from partial derivatives
    dSdNT = np.hstack((dSdN,dSdT))
    dS    = np.hstack((np.diag(dSdS),dSdNT))
    dNdNT = np.hstack((np.diag(dNdN),dNdT))
    dN    = np.hstack((dNdS,dNdNT))
    dTdNT = np.hstack((dTdN,np.diag(dTdT)))
    dT    = np.hstack((dTdS,dTdNT))
    J_out = np.vstack((dS,dN,dT))
    return J_out

##########################################################################
#### Batch culture with Hill functions for growth and death rates.
##   Function for investment into degradation is also a Hill function:
##   f_ik(T) = f_ik*Tk^h/(Tk^h +K^h)

def f_CRbatchHillInvF(t, x, args):
    #### Co-culture version of late December model candidate
    ##   dS/dt = ((1-f^u)^v*rho(N) -mu(T))S
    ##   dN/dt = -sum( 1/Y*rho(N)*S )
    ##   dT/dt = -sum( d*f^uv*rho(N)*S )
    ##  With Hill functions rho, mu of the form
    ##    rho_i(N) = sum_j( r_ij*Nj^h/(Nj^h +K^h) )
    ##     mu_i(T) = sum_k( m_ik*Tk^h/(Tk^h +K^h) )
    ##  And a no-sum Hill function for degradation investment
    ##        f(T) = f_ik*Tk^h/(Tk^h +K^h)
    Npops, Nnutr, Ntox = args['size']

    rmat, KNr, rHill = args["rtup"]   # Growth effect of nutrient, matrix of size [Npops, Nnutr]
    mmat, KTm, mHill = args["mtup"]   # Mortality of toxins, matrix of size [Npops, Ntox]
    fmat, KTf, fHill = args["ftup"]   # Per-species investment [0.0, 1.0], matrix of size [Npops, Ntox]
    fpowu, fpowv = args["powtup"]     # Investment tradeoff
    Yvec = args["Y"]    # Nutrient Yield [say gram of biomass per mol of nutrient]
    dmat = args["d"]    # Degradation efficiency [say mol of toxic compound per gram of biomass]

    ### State variables for abundance and concentrations
    x_spec = x[:Npops]
    c_nutr = x[Npops:(Npops+Nnutr)]
    c_tox  = x[(Npops+Nnutr):]

    ### Growth rate rho(N)
    argsRho = {'fmax':rmat, 'Khalf':KNr, 'n':rHill}
    rhoN_percap = f_Hillf(c_nutr, argsRho) # Per-capita growth

    ### Death rate mu(T)
    argsMu  = {'fmax':mmat, 'Khalf':KTm, 'n':mHill}
    muT_percap = f_Hillf(c_tox, argsMu)   # Per-capita mortality

    ### Investment function f
    argsInv = {'fmax':fmat, 'Khalf':KTf, 'n':fHill}
    fInv    = f_HillfNotDot(c_tox, argsInv) # Size Nnutr

    ### Parameter combinations
    fInvpow    = np.power(fInv,fpowu)    # f^u, matrix of size [Nspec, Ntox]
    fInvpowpow = np.power(fInvpow,fpowv) # f^uv
    ftot_spec  = np.sum(fInvpow, axis=1)
    ftot_recippow = np.power(np.subtract(1.0,ftot_spec),fpowv) # (1-f^u)^v

    # Change in abundance due to growth and mortality
    dSnet_percap = np.subtract(np.multiply(ftot_recippow,rhoN_percap), muT_percap) # Array length Nspec
    dS = np.multiply(dSnet_percap, x_spec)                                         # Array length Nspec

    # Change in nutrient concentration
    SY = np.multiply(x_spec, Yvec)
    nutrHill = f_HillfNotDot(c_nutr, argsRho)                                      # Array Nspec x Nnutr
    dN_nutr = np.negative(np.dot(SY, nutrHill))                                    # Array Nnutr

    # Change in toxin concentrations
    # argsDel = {'fmax':dmat, 'Khalf':KTd, 'n':dHill}
    # dS_degr = np.multiply(f_HillfNotDot(c_tox, argsDel), x_spec) # Mortality sum(mcx) per species
    dT_pertox = np.multiply(fInvpowpow,dmat)                                       # Array Nspec x Ntox
    dS_grow   = np.multiply(rhoN_percap,x_spec)                                    # Array Nspec
    #dSgrow_repT = np.reshape(np.repeat(dS_grow,Ntox),(Npops,Ntox))                 # Array Nspec x Ntox
    #dT_degr = np.multiply(dSgrow_repT,dT_pertox)
    #dT_tox  = np.negative(np.sum(dT_degr,axis=0))                                  # Array Ntox
    dT_degr = np.dot(dS_grow,dT_pertox)
    dT_tox  = np.negative( dT_degr )                                  # Array Ntox

    # Assemble state vector
    dX = np.zeros((Npops+Nnutr+Ntox))
    dX[:Npops] = dS
    dX[Npops:(Npops+Nnutr)] = dN_nutr
    dX[(Npops+Nnutr):]      = dT_tox
    # dX = np.append(np.append(dS,dN_nutr),dT_tox)
    return( dX )

def jac_CRbatchHillInvF(t, x, args):
    #### Late December model candidate
    ##   dS/dt = ((1-f^u)^v*rho(N) -mu(T))S
    ##   dN/dt = -sum( 1/Y*rho(N)*S )
    ##   dT/dt = -sum( d*f^uv*rho(N)*S )
    ##  With Hill functions rho, mu of the form
    ##    rho_i(N) = sum_j( r_ij*Nj^h/(Nj^h +K^h) )
    ##     mu_i(T) = sum_k( m_ik*Tk^h/(Tk^h +K^h) )
    ##  And a no-sum Hill function for degradation investment
    ##        f(T) = f_ik*Tk^h/(Tk^h +K^h)
    Npops, Nnutr, Ntox = args['size']

    rmat, KNr, rHill = args["rtup"]   # Growth effect of nutrient, matrix of size [Npops, Nnutr]
    mmat, KTm, mHill = args["mtup"]   # Mortality of toxins, matrix of size [Npops, Ntox]
    fmat, KTf, fHill = args["ftup"]   # Per-species investment [0.0, 1.0], matrix of size [Npops, Ntox]
    fpowu, fpowv = args["powtup"]     # Investment tradeoff
    Yvec = args["Y"]    # Nutrient Yield [say gram of biomass per mol of nutrient]
    dmat = args["d"]    # Degradation efficiency [say mol of toxic compound per gram of biomass]

    ### State variables for abundance and concentrations
    x_spec = x[:Npops]
    c_nutr = x[Npops:(Npops+Nnutr)]
    c_tox  = x[(Npops+Nnutr):]

    ### Growth rate rho(N)
    argsRho = {'fmax':rmat, 'Khalf':KNr, 'n':rHill}
    rhoN_percap = f_Hillf(c_nutr, argsRho)   # np.dot(rmat,c_nutr)
    drho_dN = dfHilldx(c_nutr, argsRho)

    ### Death rate mu(T)
    argsMu  = {'fmax':mmat, 'Khalf':KTm, 'n':mHill}
    muT_percap  = f_Hillf(c_tox, argsMu)     # np.dot(mmat,c_tox)
    dmu_dT  = dfHilldx(c_tox, argsMu)

    ### Investment function f(T)
    argsInv = {'fmax':fmat, 'Khalf':KTf, 'n':fHill}
    fT_inv  = f_HillfNotDot(c_tox, argsInv) # Size [Nspec, Ntox]
    dfinv_dT = dfHilldx(c_tox, argsInv)

    #### Powers of f investment
    fmatpow = np.power(fT_inv,fpowu)
    ftot_spec = np.sum(fmatpow, axis=1)     # First idea: f^u degradation, 1-sum(f^u) growth
    fmatpowpow = np.power(fmatpow,fpowv)    # Second idea: f^v degradation, (1-sum(f))^v growth
    ftot_powpow = np.power(ftot_spec,fpowv) # (sum(f^u))^v
    ftot_recippow = np.power(np.subtract(1.0,ftot_spec),fpowv)

    #### Compute partial derivatives wrt species abundance
    ## dSdS simple
    dSdS = np.subtract(np.multiply(ftot_recippow,rhoN_percap),muT_percap) # array [Nspec] -> diag
    ## dSdN needs derivative of Hill function
    XF = np.multiply(ftot_recippow,x_spec)                                          # Array [Nspec]
    XF_repN = np.reshape(np.repeat(XF,Nnutr),(Npops,Nnutr))
    dSdN = np.multiply(drho_dN,XF_repN)                                                  # Array [Nspec x Nnutr]
    ## dSdT needs derivative of Hill function
    # Step 1: Find overall shape of (1-sum(f^u))^v derivative
    dfudT = np.multiply(fpowu, np.power(fT_inv,fpowu-1.0))                          # Array [Nspec x Ntox]
    dfrecipvdT = np.multiply(fpowv, np.power(np.subtract(1.0,ftot_spec),fpowv-1.0)) # Array [Nspec]
    dfvdT_rep  = np.reshape(np.repeat(dfrecipvdT,Ntox),(Npops,Ntox))                #   ->  [Nspec x Ntox]
    dpowuv  = np.multiply(dfudT, dfvdT_rep)                                         # Array [Nspec x Ntox]
    dfpowdT = np.negative(np.multiply(dpowuv, dfinv_dT))                            # Array [Nspec x Ntox]
    # Step 2: Find dmudT rho(N) [length Nspec] to size Nspec x Ntox
    ##### Already computed above as dmu_dT
    # Step 3a: reshape rho(N) and X_spec [both are length Nspec] to size Nspec x Ntox
    x_repT = np.reshape(np.repeat(x_spec,Ntox),(Npops,Ntox))
    rhoN_repT = np.reshape(np.repeat(rhoN_percap,Ntox),(Npops,Ntox))
    # Step 3b: multiply, subtract and multiply again
    dnetgrow = np.subtract(np.multiply(dfpowdT,rhoN_repT),dmu_dT)
    dSdT = np.multiply(dnetgrow,x_repT)

    #### Effect on nutrient dynamics
    ## dNdS
    Y_repN = np.reshape(np.repeat(Yvec,Nnutr),(Npops,Nnutr))     # Array length Nspec -> [Nspec x Nnutr]
    nutrHill = f_HillfNotDot(c_nutr, argsRho)                    # Array [Nspec x Nnutr]
    YtimesHill = np.multiply(Y_repN, nutrHill)                   # Array [Nspec x Nnutr]
    dNdS = np.negative(np.transpose(YtimesHill))                 #   ->  [Nnutr x Nspec]
    ## dNdN
    drho_trans = np.transpose(drho_dN)                           # Array [Nnutr x Nspec]
    XY = np.multiply(Yvec,x_spec)                                # Array [Nspec]
    dNdN = np.negative(np.dot(drho_trans,XY))                    # array [Nnutr] -> diag
    ## dNdT
    dNdT = np.zeros((Nnutr,Ntox))                                # array [Nnutr,Ntox]

    ### Effect on toxin
    degrFpow = np.multiply(dmat, fmatpowpow)                            # array [Nspec x Ntox]
    dTdS_transposed = np.multiply(degrFpow, rhoN_repT)
    dTdS = np.negative(np.transpose(dTdS_transposed))                   # array [Ntox x Nspec]
    degrFpow_trans = np.transpose(degrFpow)                             # array [Ntox x Nspec]
    XrepN = np.reshape(np.repeat(x_spec,Nnutr),(Npops,Nnutr))
    drhoX = np.multiply(drho_dN,XrepN)                                  # array [Nspec x Nnutr]
    dTdN = np.negative(np.dot(degrFpow_trans, drhoX))                   # array [Ntox, Nnutr]

    fpowuv = fpowu*fpowv
    uvfpowuvminus1 = np.multiply(fpowuv,np.power(fT_inv,fpowuv-1.0))    # array [Nspec x Ntox]
    dDeltaF = np.multiply(uvfpowuvminus1,dfinv_dT)
    rhoX = np.multiply(rhoN_percap,x_spec)                              # array [Nspec]
    rX_repT = np.reshape(np.repeat(rhoX,Ntox),(Npops,Ntox))
    dT_grdeg = np.multiply(dmat, rX_repT)                               # array [Nspec x Ntox]
    dTmultsum = np.sum(np.multiply(dDeltaF,dT_grdeg),axis = 0)
    dTdT = np.negative(dTdotprod)                                       # array [Ntox] -> diag

    # Assemble Jacobian from partial derivatives
    dSdNT = np.hstack((dSdN,dSdT))
    dS    = np.hstack((np.diag(dSdS), dSdNT))
    dNdNT = np.hstack((np.diag(dNdN),dNdT))
    dN    = np.hstack((dNdS,dNdNT))
    # dN    = np.transpose(np.vstack((dNdS,np.transpose(dNdNT))))
    dTdNT = np.hstack((dTdN,np.diag(dTdT)))
    dT    = np.hstack((dTdS,dTdNT))
    # dT    = np.transpose(np.vstack((dTdS,np.transpose(dTdNT))))
    J_out = np.vstack((dS,dN,dT))
    return( J_out )

###############################################################
def f_CRbatchHillDegrMat_mono(t, x, args):
    #### Late December model candidate
    ##   dS/dt = ((1-f^u)^v*rho(N) -mu(T))S
    ##   dN/dt = -sum( 1/Y*rho(N)*S )
    ##   dT/dt = -sum( d*f^uv*rho(N)*S )
    ##  With Hill functions rho, mu of the form
    ##    rho_i(N) = sum_j( r_ij*Nj^h/(Nj^h +K^h) )
    ##     mu_i(T) = sum_k( m_ik*Tk^h/(Tk^h +K^h) )
    ##  And a no-sum Hill function for degradation investment
    ##        f(T) = f_ik*Tk^h/(Tk^h +K^h)
    Npops, Nnutr, Ntox = args['size']

    rmat, KNr, rHill = args["rtup"]   # Growth effect of nutrient, matrix of size [Npops, Nnutr]
    mmat, KTm, mHill = args["mtup"]   # Mortality of toxins, matrix of size [Npops, Ntox]
    dmat, KTd, dHill = args["dtup"]   # Per-species investment [0.0, 1.0], matrix of size [Npops, Ntox]
    fpowu, fpowv = args["powtup"]     # Investment tradeoff
    Yvec = args["Y"]    # Nutrient Yield [say gram of biomass per mol of nutrient]
    fmat = args["fmat"]    # Degradation efficiency [say mol of toxic compound per gram of biomass]

    ### State variables for abundance and concentrations
    x_spec = x[0]
    c_nutr = x[1:(1+Nnutr)]
    c_tox  = x[(1+Nnutr):]

    ### Degradation function d
    ### Investment function f(T)
    argsDegr = {'fmax':dmat, 'Khalf':KTd, 'n':dHill}
    degT_inv = f_HillfNotDot(c_tox, argsDegr) # Size [Nspec, Ntox]

    ### Parameter combinations
    fInvpow = np.power(fmat,fpowu)       # f^u, matrix of size [Nspec, Ntox]
    fInvpowpow = np.power(fInvpow,fpowv) # f^uv
    ftot_spec  = np.sum(fInvpow) # np.sum(fInvpow, axis=1)
    ftot_recippow = np.power(np.subtract(1.0,ftot_spec),fpowv) # (1-f^u)^v

    # Change in abundance due to growth and mortality
    argsRho = {'fmax':rmat, 'Khalf':KNr, 'n':rHill}
    argsMu  = {'fmax':mmat, 'Khalf':KTm, 'n':mHill}
    dS_grow = np.multiply(f_Hillf(c_nutr,argsRho), x_spec) # Growth    rcx per species
    dS_mort = np.multiply(f_Hillf(c_tox, argsMu), x_spec) # Mortality sum(mcx) per species
    dS = np.subtract(np.multiply(ftot_recippow,dS_grow), dS_mort)

    # Change in nutrient concentration
    SY = np.multiply(x_spec, Yvec)
    nutrHill = f_HillfNotDot(c_nutr, argsRho) # Size Nnutr
    dN_nutr = np.negative(np.multiply(SY, nutrHill))

    # Change in toxin concentrations
    # argsDel = {'fmax':dmat, 'Khalf':KTd, 'n':dHill}
    # dS_degr = np.multiply(f_HillfNotDot(c_tox, argsDel), x_spec) # Mortality sum(mcx) per species
    dT_pertox = np.multiply(fInvpowpow,degT_inv)
    dT_degr = np.multiply(dS_grow,dT_pertox)
    dT_tox  = np.negative(dT_degr)

    # Assemble state vector
    dX = np.append(np.append(dS,dN_nutr),dT_tox)
    return( dX )

def jac_CRbatchHillDegrMat_mono(t, x, args):
    #### Late December model candidate
    ##   dS/dt = ((1-f^u)^v*rho(N) -mu(T))S
    ##   dN/dt = -sum( 1/Y*rho(N)*S )
    ##   dT/dt = -sum( d*f^uv*rho(N)*S )
    ##  With Hill functions rho, mu of the form
    ##    rho_i(N) = sum_j( r_ij*Nj^h/(Nj^h +K^h) )
    ##     mu_i(T) = sum_k( m_ik*Tk^h/(Tk^h +K^h) )
    ##  And a no-sum Hill function for degradation
    ##        d(T) = d_ik*Tk^h/(Tk^h +K^h)
    Npops, Nnutr, Ntox = args['size']

    rmat, KNr, rHill = args["rtup"]   # Growth effect of nutrient, matrix of size [Npops, Nnutr]
    mmat, KTm, mHill = args["mtup"]   # Mortality of toxins, matrix of size [Npops, Ntox]
    dmat, KTd, dHill = args["dtup"]   # Per-species investment [0.0, 1.0], matrix of size [Npops, Ntox]
    fpowu, fpowv = args["powtup"]     # Investment tradeoff
    Yvec = args["Y"]    # Nutrient Yield [say gram of biomass per mol of nutrient]
    fmat = args["fmat"]    # Degradation efficiency [say mol of toxic compound per gram of biomass]

    ### State variables for abundance and concentrations
    x_spec = x[0]
    c_nutr = x[1:(1+Nnutr)]
    c_tox  = x[(1+Nnutr):]

    ### Growth rate rho(N)
    argsRho = {'fmax':rmat, 'Khalf':KNr, 'n':rHill}
    rhoN_percap = f_Hillf(c_nutr,argsRho) # np.dot(rmat,c_nutr)
    drho_dN = dfHilldx(c_nutr,argsRho)

    ### Death rate mu(T)
    argsMu  = {'fmax':mmat, 'Khalf':KTm, 'n':mHill}
    muT_percap  = f_Hillf(c_tox,argsMu)   # np.dot(mmat,c_tox)
    dmu_dT  = dfHilldx(c_tox,argsMu)

    ### Investment function f(T)
    argsDegr = {'fmax':dmat, 'Khalf':KTd, 'n':dHill}
    degT_inv = f_HillfNotDot(c_tox, argsDegr) # Size [Nspec, Ntox]
    ddeg_dT  = dfHilldx(c_tox, argsDegr)

    #### Powers of f investment
    fmatpow = np.power(fmat,fpowu)        # First idea: f^u degradation, 1-sum(f^u) growth
    fmatpowpow = np.power(fmatpow,fpowv)    # Second idea: f^v degradation, (1-sum(f))^v growth
    ftot_spec  = np.sum(fmatpow)
    ftot_powpow = np.power(ftot_spec,fpowv) # (sum(f^u))^v
    ftot_recippow = np.power(np.subtract(1.0,ftot_spec),fpowv)

    #### Compute partial derivatives wrt species abundance
    dS_growth = np.multiply(ftot_recippow,rhoN_percap)
    ## dSdS simple
    dSdS = np.subtract(dS_growth,muT_percap) # array [Nspec] -> diag
    ## dSdN needs derivative of Hill function
    XF = np.multiply(ftot_recippow,x_spec)
    dSdN  = np.multiply(drho_dN,XF)
    ## dSdT needs derivative of Hill function
    x_repT = np.array([x_spec]*Ntox)
    dSdT = np.negative( np.multiply(dmu_dT, x_repT) )

    #### Effect on nutrient dynamics
    ## dNdS
    Y_repN = np.repeat(Yvec,Nnutr)
    nutrHill = f_HillfNotDot(c_nutr, argsRho) # Single Si
    # nutrHill = np.transpose(f_HillfNotDot(c_nutr, argsRho)) # Transpose for multiple Si
    dNdS = np.negative(np.multiply(Y_repN, nutrHill))
    ## dNdN
    XY = np.multiply(Yvec,x_spec)
    rXY = np.multiply(drho_dN,XY) # Repeat for multiple Si
    dNdN = np.negative(rXY)       # array [Nnutr] -> diag
    ## dNdT
    dNdT = np.zeros((Nnutr,Ntox)) # array [Nnutr,Ntox]

    ### Effect on toxin
    degrFpow = np.multiply(degT_inv, fmatpowpow)
    dTdS = np.negative(np.multiply(degrFpow, rhoN_percap))               # array [Ntox]        # CHECK
    dTdN = np.negative(np.outer(degrFpow, np.multiply(drho_dN,x_spec)))  # array [Ntox, Nnutr] # CHECK

    dDeltaF = np.multiply(ftot_powpow,ddeg_dT)
    rX_repT = np.repeat(np.multiply(rhoN_percap,x_spec),Ntox)  # array [Nnutr]
    dTdT = np.negative( np.multiply(dDeltaF, rX_repT) )                  # array [Ntox] -> diag

    # Assemble Jacobian from partial derivatives
    dSdNT = np.hstack((dSdN,dSdT))
    dS    = np.hstack(([dSdS,], dSdNT))
    dNdNT = np.hstack((np.diag(dNdN),dNdT))
    dN    = np.transpose(np.vstack((dNdS,np.transpose(dNdNT))))
    dTdNT = np.hstack((dTdN,np.diag(dTdT)))
    dT    = np.transpose(np.vstack((dTdS,np.transpose(dTdNT))))
    J_out = np.vstack((dS,dN,dT))
    return J_out

###############################################################
## Batch culture model, ! monoculture !
def f_CRHill_mono(t, x, args):
    Npops, Nnutr, Ntox = args['size']

    rmat = args["r"]    # Growth effect of nutrient, matrix of size [Npops, Nnutr]
    mmat = args["m"]    # Mortality of toxins, matrix of size [Npops, Ntox]
    fmat = args["finv"] # Per-species investment [0.0, 1.0] into degradation of each toxin
    fpow = args["fpow"] # Per-species investment [0.0, 1.0] into degradation of each toxin
    KvecN = args["KN"]  # Half-max degradation rate of nutrient substrates, array of length [Nnutr]
    KvecT = args["KT"]  # Half-max degradation rate of toxic substrates, array of length [Ntox]
    vecHill = args["nHill"] # Hill coefficient for toxic substrates (subject to change)
    Yvec = args["Y"]    # Nutrient Yield [say gram of biomass per mol of nutrient]
    dvec = args["d"]    # Toxin degradation rate [say mol of toxin per gram of biomass]

    ### State variables for abundance and concentrations
    x_spec = x[0]
    c_nutr = x[1:(1+Nnutr)]
    c_tox  = x[(1+Nnutr):]

    ### Parameter combinations
    fmatpow = np.power(fmat,fpow)
    ftot_spec = np.sum(fmatpow,axis=1)
    ftr = np.transpose(fmatpow)
    # rtr = np.transpose(rmat)

    # Change in abundance due to growth and mortality
    argsN = {'fmax':rmat, 'Khalf':KvecN, 'n':np.ones(len(KvecN))}
    argsT = {'fmax':mmat, 'Khalf':KvecT, 'n':vecHill}
    dS_grow = np.multiply(f_Hillf(c_nutr,argsN), x_spec) # Growth    rcx per species
    dS_mort = np.multiply(f_Hillf(c_tox, argsT), x_spec) # Mortality sum(mcx) per species
    invest_gr = np.subtract([1.0]*Npops,ftot_spec)
    dS = np.subtract(np.multiply(invest_gr,dS_grow), dS_mort)

    # Change in nutrient concentration
    SY  = np.multiply(x_spec,Yvec)
    rSY   = np.multiply(rmat, SY)  # Dot product for Nspec>1
    nutrMonod = np.divide(c_nutr,np.add(c_nutr, KvecN))
    dN_nutr = np.negative(np.multiply(rSY, nutrMonod))

    # Change in toxin concentrations
    dT_invest = np.multiply(dvec,dS_grow)
    dT_pertox = np.dot(ftr,dT_invest)
    dT_tox    = np.negative(np.multiply(c_tox, dT_pertox))
    ### MIGHT WANT TO CHANGE ALSO dT

    # Assemble state vector
    dX = np.append(np.append(dS,dN_nutr),dT_tox)
    return( dX )

###############################################################
## Batch culture model, Hill function for degradation investment
##   ! monoculture !

def f_CRbatchHill_mono(t, x, args):
    #### Standard model candidate
    ##   dS/dt = ((1-f^u)^v*rho(N) -mu(T))S
    ##   dN/dt = -sum( 1/Y*rho(N)*S )
    ##   dT/dt = -sum( d*f^uv*rho(N)*S )
    ##  With Hill functions rho, mu of the form
    ##    rho_i(N) = sum_j( r_ij*Nj^h/(Nj^h +K^h) )
    ##     mu_i(T) = sum_k( m_ik*Tk^h/(Tk^h +K^h) )
    Npops, Nnutr, Ntox = args['size']

    rmat = args["r"]    # Growth effect of nutrient, length [Nnutr]
    mmat = args["m"]    # Mortality of toxins, length [Ntox]
    fmat = args["finv"] # Per-species investment [0.0, 1.0], length [Ntox]
    fpowu, fpowv = args["fpow"] # Investment tradeoff
    Yvec = args["Y"]    # Nutrient Yield [say gram of biomass per mol of nutrient]
    dvec = args["d"]    # Degradation efficiency [say mol of toxic compound per gram of biomass]
    KNr  = args["KN"]   # Half-max coefficient for nutrient uptake
    KTm  = args["KT"]   # Half-max coefficient for tox degradation

    ### State variables for abundance and concentrations
    x_spec = x[0]
    c_nutr = x[Npops:(Npops+Nnutr)]
    c_tox  = x[(Npops+Nnutr):]

    ### Parameter combinations
    fpowu = 1.0
    fpowv = 2.0
    fInvpow = np.power(fmat,fpowu)       # f^u, array of length [Ntox]
    fInvpowpow = np.power(fInvpow,fpowv) # f^uv
    ftot_spec  = np.sum(fInvpow)
    ftot_recippow = np.power(np.subtract(1.0,ftot_spec),fpowv) # (1-f^u)^v

    # Change in abundance due to growth and mortality
    argsRho = {'fmax':rmat, 'Khalf':KNr, 'n':1.0}
    argsMu  = {'fmax':mmat, 'Khalf':KTm, 'n':2.0}
    dS_grow = np.multiply(f_Hillf_mono(c_nutr, argsRho), x_spec) # Growth    rcx per species
    dS_mort = np.multiply(f_Hillf_mono(c_tox, argsMu), x_spec)   # Mortality mcx per species
    dS = np.subtract(np.multiply(ftot_recippow, dS_grow), dS_mort)

    # Change in nutrient concentration
    SY = np.multiply(x_spec, Yvec)
    nutrHill = f_HillfNotDot_mono(c_nutr, argsRho) # Len Nnutr
    dN_nutr  = np.negative(np.multiply(SY, nutrHill))

    # Change in toxin concentrations
    dT_pertox = np.multiply(fInvpowpow,dvec) # Len Ntox
    dT_degr = np.multiply(dS_grow,dT_pertox)
    dT_tox  = np.negative(np.multiply(dT_degr, c_tox))

    # Assemble state vector
    dX = np.append(np.append(dS,dN_nutr),dT_tox)
    return( dX )

def jac_CRbatchHill_mono(t, x, args):
    #### Late December model candidate
    ##   dS/dt = ((1-f^u)^v*rho(N) -mu(T))S
    ##   dN/dt = -sum( 1/Y*rho(N)*S )
    ##   dT/dt = -sum( d*f^uv*rho(N)*S )
    ##  With Hill functions rho, mu of the form
    ##    rho_i(N) = sum_j( r_ij*Nj^h/(Nj^h +K^h) )
    ##     mu_i(T) = sum_k( m_ik*Tk^h/(Tk^h +K^h) )
    ##  And a no-sum Hill function for degradation investment
    ##        f(T) = f_ik*Tk^h/(Tk^h +K^h)
    Npops, Nnutr, Ntox = args['size']

    rmat = args["r"]    # Growth effect of nutrient, length [Nnutr]
    mmat = args["m"]    # Mortality of toxins, length [Ntox]
    fmat = args["finv"] # Per-species investment [0.0, 1.0], length [Ntox]
    fpowu, fpowv = args["fpow"] # Investment tradeoff
    Yvec = args["Y"]    # Nutrient Yield [say gram of biomass per mol of nutrient]
    dvec = args["d"]    # Degradation efficiency [say mol of toxic compound per gram of biomass]
    KNr  = args["KN"]   # Half-max coefficient for nutrient uptake
    KTm  = args["KT"]   # Half-max coefficient for tox degradation

    ### State variables for abundance and concentrations
    x_spec = x[0]
    c_nutr = x[1:(1+Nnutr)]
    c_tox  = x[(1+Nnutr):]

    ### Growth rate rho(N)
    argsRho = {'fmax':rmat, 'Khalf':KNr, 'n':1.0}
    rhoN_percap = f_Hillf_mono(c_nutr,argsRho) # np.dot(rmat,c_nutr)
    drho_dN = dfHilldx(c_nutr,argsRho)

    ### Death rate mu(T)
    argsMu  = {'fmax':mmat, 'Khalf':KTm, 'n':2.0}
    muT_percap  = f_Hillf_mono(c_tox,argsMu)   # np.dot(mmat,c_tox)
    dmu_dT  = dfHilldx(c_tox,argsMu)

    # fT_inv = fmat
    # dfinv_dT = 0

    #### Powers of f investment
    fmatpow = np.power(fmat,fpowu)        # First idea: f^u degradation, 1-sum(f^u) growth
    fmatpowpow = np.power(fmatpow,fpowv)    # Second idea: f^v degradation, (1-sum(f))^v growth
    ftot_spec  = np.sum(fmatpow)
    ftot_powpow = np.power(ftot_spec,fpowv) # (sum(f^u))^v
    ftot_recippow = np.power(np.subtract(1.0,ftot_spec),fpowv)

    #### Compute partial derivatives wrt species abundance
    ## dSdS simple
    dSdS = np.subtract(np.multiply(ftot_recippow,rhoN_percap),muT_percap) # array [Nspec] -> diag
    ## dSdN needs derivative of Hill function
    XF = np.multiply(ftot_recippow,x_spec)
    dSdN = np.multiply(drho_dN,XF)
    ## dSdT needs derivative of Hill function
    dSdT = np.negative( np.multiply(x_spec, dmu_dT) )

    #### Effect on nutrient dynamics
    ## dNdS
    nutrHill = f_HillfNotDot_mono(c_nutr, argsRho) # Single Si
    # nutrHill = np.transpose(f_HillfNotDot(c_nutr, argsRho)) # Transpose for multiple Si
    dNdS = np.negative(np.multiply(Yvec, nutrHill))
    ## dNdN
    XY = np.multiply(Yvec,x_spec)
    rXY = np.multiply(drho_dN,XY) # Repeat for multiple Si
    dNdN = np.negative(rXY)       # array [Nnutr] -> diag
    ## dNdT
    dNdT = np.zeros((Nnutr,Ntox))                    # array [Nnutr,Ntox]

    ### Effect on toxin
    degrFpow = np.multiply(dmat, fmatpowpow)                                                # CHECK
    dftimesT = np.multiply(degrFpow,c_tox)
    dTdS = np.negative(np.multiply(degrFpow, rhoN_percap))              # array [Ntox]        # CHECK
    dTdN = np.negative(np.outer(degrFpow, np.multiply(drho_dN,x_spec)))  # array [Ntox, Nnutr] # CHECK

    rX = np.multiply(rhoN_percap,x_spec)  # array [Nnutr]
    dT_grdeg = np.multiply(degrFpow, rX)
    dTdT = np.negative(np.multiply(0.0,dT_grdeg))                  # array [Ntox] -> diag

    # Assemble Jacobian from partial derivatives
    dSdNT = np.hstack((dSdN,dSdT))
    dS    = np.hstack(([dSdS,], dSdNT))
    dNdNT = np.hstack((np.diag(dNdN),dNdT))
    dN    = np.transpose(np.vstack((dNdS,np.transpose(dNdNT))))
    dTdNT = np.hstack((dTdN,np.diag(dTdT)))
    dT    = np.transpose(np.vstack((dTdS,np.transpose(dTdNT))))
    J_out = np.vstack((dS,dN,dT))
    return J_out

def f_CRbatchHillInvF_mono(t, x, args):
    #### Late December model candidate
    ##   dS/dt = ((1-f^u)^v*rho(N) -mu(T))S
    ##   dN/dt = -sum( 1/Y*rho(N)*S )
    ##   dT/dt = -sum( d*f^uv*rho(N)*S )
    ##  With Hill functions rho, mu of the form
    ##    rho_i(N) = sum_j( r_ij*Nj^h/(Nj^h +K^h) )
    ##     mu_i(T) = sum_k( m_ik*Tk^h/(Tk^h +K^h) )
    ##  And a no-sum Hill function for degradation investment
    ##        f(T) = f_ik*Tk^h/(Tk^h +K^h)
    Npops, Nnutr, Ntox = args['size']

    rvec, KNr, rHill = args["rtup"]   # Growth effect of nutrient, length [Nnutr]
    mvec, KTm, mHill = args["mtup"]   # Mortality of toxins, length [Ntox]
    fvec, KTf, fHill = args["ftup"]   # Per-species investment [0.0, 1.0], length [Ntox]
    fpowu, fpowv = args["powtup"]     # Investment tradeoff
    Yvec = args["Y"]    # Nutrient Yield [say gram of biomass per mol of nutrient]
    dvec = args["d"]    # Degradation efficiency [say mol of toxic compound per gram of biomass]

    ### State variables for abundance and concentrations
    x_spec = x[0]
    c_nutr = x[1:(1+Nnutr)]
    c_tox  = x[(1+Nnutr):]

    ### Investment function f
    argsInv = {'fmax':fvec, 'Khalf':KTf, 'n':fHill}
    fInv    = f_HillfNotDot_mono(c_tox, argsInv) # Size Nnutr

    ### Parameter combinations
    fInvpow = np.power(fInv,fpowu)       # f^u, array of length [Ntox]
    fInvpowpow = np.power(fInvpow,fpowv) # f^uv
    ftot_spec  = np.sum(fInvpow)
    ftot_recippow = np.power(np.subtract(1.0,ftot_spec),fpowv) # (1-f^u)^v

    # Change in abundance due to growth and mortality
    argsRho = {'fmax':rvec, 'Khalf':KNr, 'n':rHill}
    argsMu  = {'fmax':mvec, 'Khalf':KTm, 'n':mHill}
    dS_grow = np.multiply(f_Hillf_mono(c_nutr, argsRho), x_spec) # Growth    rcx per species
    dS_mort = np.multiply(f_Hillf_mono(c_tox, argsMu), x_spec)   # Mortality mcx per species
    dS = np.subtract(np.multiply(ftot_recippow, dS_grow), dS_mort)

    # Change in nutrient concentration
    SY = np.multiply(x_spec, Yvec)
    nutrHill = f_HillfNotDot_mono(c_nutr, argsRho) # Len Nnutr
    dN_nutr = np.negative(np.multiply(SY, nutrHill))

    # Change in toxin concentrations
    # argsDel = {'fmax':dmat, 'Khalf':KTd, 'n':dHill}
    # dS_degr = np.multiply(f_HillfNotDot(c_tox, argsDel), x_spec) # Mortality sum(mcx) per species
    dT_pertox = np.multiply(fInvpowpow,dvec) # Len Ntox
    dT_degr = np.multiply(dS_grow,dT_pertox)
    dT_tox  = np.negative(dT_degr)

    # Assemble state vector
    dX = np.append(np.append(dS,dN_nutr),dT_tox)
    return( dX )

def jac_CRbatchHillInvF_mono(t, x, args):
    #### Late December model candidate
    ##   dS/dt = ((1-f^u)^v*rho(N) -mu(T))S
    ##   dN/dt = -sum( 1/Y*rho(N)*S )
    ##   dT/dt = -sum( d*f^uv*rho(N)*S )
    ##  With Hill functions rho, mu of the form
    ##    rho_i(N) = sum_j( r_ij*Nj^h/(Nj^h +K^h) )
    ##     mu_i(T) = sum_k( m_ik*Tk^h/(Tk^h +K^h) )
    ##  And a no-sum Hill function for degradation investment
    ##        f(T) = f_ik*Tk^h/(Tk^h +K^h)
    Npops, Nnutr, Ntox = args['size']

    rmat, KNr, rHill = args["rtup"]   # Growth effect of nutrient, matrix of size [Npops, Nnutr]
    mmat, KTm, mHill = args["mtup"]   # Mortality of toxins, matrix of size [Npops, Ntox]
    fmat, KTf, fHill = args["ftup"]   # Per-species investment [0.0, 1.0], matrix of size [Npops, Ntox]
    fpowu, fpowv = args["powtup"]     # Investment tradeoff
    Yvec = args["Y"]    # Nutrient Yield [say gram of biomass per mol of nutrient]
    dmat = args["d"]    # Degradation efficiency [say mol of toxic compound per gram of biomass]

    ### State variables for abundance and concentrations
    x_spec = x[0]
    c_nutr = x[1:(1+Nnutr)]
    c_tox  = x[(1+Nnutr):]

    ### Growth rate rho(N)
    argsRho = {'fmax':rmat, 'Khalf':KNr, 'n':rHill}
    rhoN_percap = f_Hillf_mono(c_nutr,argsRho) # np.dot(rmat,c_nutr)
    drho_dN = dfHilldx(c_nutr,argsRho)

    ### Death rate mu(T)
    argsMu  = {'fmax':mmat, 'Khalf':KTm, 'n':mHill}
    muT_percap  = f_Hillf_mono(c_tox,argsMu)   # np.dot(mmat,c_tox)
    dmu_dT  = dfHilldx(c_tox,argsMu)

    ### Investment function f(T)
    argsInv = {'fmax':fmat, 'Khalf':KTf, 'n':fHill}
    fT_inv  = f_HillfNotDot_mono(c_tox, argsInv) # Size [Nspec, Ntox]
    dfinv_dT = dfHilldx(c_tox,argsInv)

    #### Powers of f investment
    fmatpow = np.power(fT_inv,fpowu)        # First idea: f^u degradation, 1-sum(f^u) growth
    fmatpowpow = np.power(fmatpow,fpowv)    # Second idea: f^v degradation, (1-sum(f))^v growth
    ftot_spec  = np.sum(fmatpow)
    ftot_powpow = np.power(ftot_spec,fpowv) # (sum(f^u))^v
    ftot_recippow = np.power(np.subtract(1.0,ftot_spec),fpowv)

    #### Compute partial derivatives wrt species abundance
    ## dSdS simple
    dSdS = np.subtract(np.multiply(ftot_recippow,rhoN_percap),muT_percap) # array [Nspec] -> diag
    ## dSdN needs derivative of Hill function
    XF = np.multiply(ftot_recippow,x_spec)
    dSdN = np.multiply(drho_dN,XF)
    ## dSdT needs derivative of Hill function
    # x_repT = np.multiply(x_spec,np.ones(Ntox))
    dfrecipvdT = np.multiply(fpowv, np.power(np.subtract(1.0,ftot_spec),fpowv-1.0))
    dfudT      = np.multiply(fpowu, np.power(fT_inv,fpowu-1.0))
    dpowuv = np.multiply(dfrecipvdT, dfudT)             # Len Ntox
    dfpowdT = np.negative(np.multiply(dpowuv,dfinv_dT))
    dnetgrow = np.subtract(np.multiply(dfpowdT,rhoN_percap),dmu_dT)
    dSdT = np.multiply(x_spec, dnetgrow)

    #### Effect on nutrient dynamics
    ## dNdS
    # Y_repN = np.repeat(Yvec,Nnutr) # Y is scalar
    nutrHill = f_HillfNotDot_mono(c_nutr, argsRho) # Single Si
    # nutrHill = np.transpose(f_HillfNotDot(c_nutr, argsRho)) # Transpose for multiple Si
    dNdS = np.negative(np.multiply(Yvec, nutrHill))
    ## dNdN
    XY = np.multiply(Yvec,x_spec)
    rXY = np.multiply(drho_dN,XY) # Repeat for multiple Si
    dNdN = np.negative(rXY)       # array [Nnutr] -> diag
    ## dNdT
    dNdT = np.zeros((Nnutr,Ntox))                    # array [Nnutr,Ntox]

    ### Effect on toxin
    degrFpow = np.multiply(dmat, fmatpowpow)                                                # CHECK
    dTdS = np.negative(np.multiply(degrFpow, rhoN_percap))              # array [Ntox]        # CHECK
    dTdN = np.negative(np.outer(degrFpow, np.multiply(drho_dN,x_spec)))  # array [Ntox, Nnutr] # CHECK

    fpowuv = fpowu*fpowv
    uvfpowuvminus1 = np.multiply(fpowuv,np.power(fT_inv,fpowuv-1.0))
    dDeltaF = np.multiply(uvfpowuvminus1,dfinv_dT)
    #rX_repT = np.repeat(np.multiply(rhoN_percap,x_spec),Ntox)  # array [Nnutr]
    rX = np.multiply(rhoN_percap,x_spec)  # array [Nnutr]
    dT_grdeg = np.multiply(dmat, rX)
    dTdT = np.negative(np.multiply(dDeltaF,dT_grdeg))                  # array [Ntox] -> diag

    # Assemble Jacobian from partial derivatives
    dSdNT = np.hstack((dSdN,dSdT))
    dS    = np.hstack(([dSdS,], dSdNT))
    dNdNT = np.hstack((np.diag(dNdN),dNdT))
    dN    = np.transpose(np.vstack((dNdS,np.transpose(dNdNT))))
    dTdNT = np.hstack((dTdN,np.diag(dTdT)))
    dT    = np.transpose(np.vstack((dTdS,np.transpose(dTdNT))))
    J_out = np.vstack((dS,dN,dT))
    return J_out


##########################################################################
#### CHEMOSTAT MODELS
##########################################################################
### Standard chemostat model, with linear uptake of nutrients

def f_chemostat(t,x, args):
    Nspec, Nnutr, Ntox = args['size']
    N0, T0 = args['nt_inflow']

    aflow = args["a"]   # Flow rate 'alpha' of chemostat
    rmat = args["r"]    # Growth effect of nutrient, matrix of size [Nspec, Nnutr]
    mmat = args["m"]    # Mortality of toxins, matrix of size [Nspec, Ntox]
    fmat = args["finv"] # Per-species investment [0.0, 1.0] into degradation of each toxin
    Yvec = args["Y"]    # Nutrient Yield [say gram of biomass per mol of nutrient]
    dvec = args["d"]    # Toxin degradation rate [say mol of toxin per gram of biomass]

    ### State variables for abundance and concentrations
    x_spec = x[:Nspec]
    c_nutr = x[Nspec:(Nspec+Nnutr)]
    c_tox  = x[(Nspec+Nnutr):]

    ### Parameter combinations
    ftot_spec = np.sum(fmat,axis=1)
    ftr = np.transpose(fmat)
    rtr = np.transpose(rmat)

    # Change in abundance due to growth and mortality
    dS_grow = np.multiply(np.dot(rmat,c_nutr),x_spec) # Growth    rcx per species
    dS_mort = np.multiply(np.dot(mmat,c_tox)+aflow,x_spec)       # Mortality sum(mcx) per species
    invest_gr = np.subtract([1.0]*Nspec,ftot_spec)
    dS = np.subtract(np.multiply(invest_gr,dS_grow), dS_mort)

    # Change in nutrient concentration
    SY  = np.multiply(x_spec,Yvec)
    rSY = np.dot(rtr, SY)
    dN_flow   = aflow*np.subtract(N0,c_nutr)
    dN_nutr   = np.subtract(dN_flow, np.multiply(rSY, c_nutr))

    # Change in toxin concentrations
    dT_flow   = aflow*np.subtract(T0,c_tox)
    dT_invest = np.multiply(dvec,dS_grow)
    dT_pertox = np.dot(ftr,dT_invest)
    dT_tox    = np.subtract(dT_flow, np.multiply(c_tox, dT_pertox))

    # Assemble state vector
    dX = np.append(np.append(dS,dN_nutr),dT_tox)
    return dX

def jac_chemostat(t, x, args):
    Nspec, Nnutr, Ntox = args['size']
    N0, T0 = args['nt_inflow']

    aflow = args["a"]   # Flow rate 'alpha' of chemostat
    rmat = args["r"]    # Growth effect of nutrient, matrix of size [Nspec, Nnutr]
    mmat = args["m"]    # Mortality of toxins, matrix of size [Nspec, Ntox]
    fmat = args["finv"] # Per-species investment [0.0, 1.0] into degradation of each toxin
    Yvec = args["Y"]    # Nutrient Yield [say gram of biomass per mol of nutrient]
    dvec = args["d"]    # Toxin degradation rate [say mol of toxin per gram of biomass]

    ### State variables for abundance and concentrations
    x_spec = x[:Nspec]
    c_nutr = x[Nspec:(Nspec+Nnutr)]
    c_tox  = x[(Nspec+Nnutr):]

    #### Prelims
    gr_percap = np.dot(rmat,c_nutr)
    d_percap  = np.add(np.dot(mmat,c_tox),aflow)
    invest_gr = np.subtract([1.0]*Nspec,np.sum(fmat, axis=1))
    ftr = np.transpose(fmat)
    rtr = np.transpose(rmat)

    #### Compute partial derivatives wrt species abundance
    XF = np.multiply(invest_gr,x_spec)
    xf_repN = np.transpose(np.reshape([XF]*Nnutr,[Nnutr,Nspec]))
    x_repT  = np.transpose(np.reshape([x_spec]*Ntox,[Ntox,Nspec]))
    dSdS =  np.subtract(np.multiply(invest_gr,gr_percap),d_percap)  # -> diag
    dSdN =  np.multiply(rmat,xf_repN)
    dSdT = -np.multiply(mmat,x_repT)

    ### Effect on nutrient dynamics
    Y_repN = np.reshape([Yvec]*Nnutr,[Nnutr,Nspec])
    rY = np.multiply(rtr,Y_repN)
    XY = np.multiply(Yvec,x_spec)
    CNrepX = np.transpose(np.reshape([c_nutr]*Nspec,[Nspec,Nnutr]))
    dNdS = -1.0*np.multiply(rY, CNrepX)                             # array [Nnutr,Ntox]
    dNdN = -1.0*np.add(aflow, np.dot(rtr,XY))                       # array [Nnutr] -> diag
    dNdT = np.zeros((Nnutr,Ntox))                                   # array [Nnutr,Ntox]

    ### Effect on toxin
    #rd  = np.multiply(rvec,dvec)
    dXvec  = np.multiply(dvec,x_spec) # array [Nspec]
    dX_rep = np.transpose(np.reshape([dXvec]*Nnutr,[Nnutr,Nspec]))  # array [Nspec, Nnutr]
    rdX = np.multiply(dX_rep,rmat)    # array [Nspec, Nnutr]
    drN = np.multiply(dvec,gr_percap) # array [Nspec]
    rdXN = np.multiply(drN,x_spec)    # array [Nspec]
    drN_repT = np.reshape([drN]*Ntox,[Ntox,Nspec])                  # array [Ntox, Nspec]

    ctox_repX = np.transpose(np.reshape([c_tox]*Nspec,[Nspec,Ntox]))
    ctox_repN = np.transpose(np.reshape([c_tox]*Nnutr,[Nnutr,Ntox]))
    dTdS = -1.0*np.multiply(ctox_repX, np.multiply(ftr,drN_repT))   # array [Ntox, Nspec]
    dTdN = -1.0*np.multiply(ctox_repN, np.dot(ftr,rdX))             # array [Ntox, Nnutr]
    dTdT = -1.0*np.add(aflow, np.dot(ftr,rdXN))                     # array [Ntox] -> diag

    # Assemble Jacobian from partial derivatives
    dSdNT = np.hstack((dSdN,dSdT))
    dS    = np.hstack((np.diag(dSdS),dSdNT))
    dNdNT = np.hstack((np.diag(dNdN),dNdT))
    dN    = np.hstack((dNdS,dNdNT))
    dTdNT = np.hstack((dTdN,np.diag(dTdT)))
    dT    = np.hstack((dTdS,dTdNT))
    J_out = np.vstack((dS,dN,dT))
    return J_out

##########################################################################
#### Chemostat model. Single nutrient, single toxin compounds
##########################################################################

def f_chem_1N1T(t,x, args):
    Nspec = args['size']
    N0, T0 = args['nt_inflow']

    aflow = args["a"]   # Flow rate 'alpha' of chemostat
    rvec = args["r"]    # Growth effect of nutrient, matrix of size [Nspec, 1]
    mvec = args["m"]    # Mortality of toxins, matrix of size [Nspec, 1]
    fvec = args["finv"] # Per-species investment [0.0, 1.0] into degradation of each toxin
    Yvec = args["Y"]    # Nutrient Yield [say gram of biomass per mol of nutrient]
    dvec = args["d"]    # Toxin degradation rate [say mol of toxin per gram of biomass]

    ### State variables for abundance and concentrations
    x_spec = x[:Nspec]
    c_nutr = x[Nspec]
    c_tox  = x[Nspec+1]

    # Change in abundance due to growth and mortality
    dS_grow = np.multiply(np.multiply(rvec,c_nutr),x_spec)      # Growth    rcx per species
    dS_mort = np.multiply(np.multiply(mvec,c_tox)+aflow,x_spec) # Mortality sum(mcx) per species
    finvest = np.subtract([1.0]*Nspec,fvec)
    dS = np.subtract(np.multiply(finvest,dS_grow), dS_mort)

    # Change in nutrient concentration
    rY  = np.multiply(rvec,Yvec)
    rYX = np.dot(rY, x_spec)
    dN_inflow = np.multiply(aflow,N0)
    dN_out    = np.multiply(np.add(aflow,rYX), c_nutr)
    dN = np.subtract(dN_inflow,dN_out)

    # Change in toxin concentrations
    dT_inflow = np.multiply(aflow,T0)
    dT_invest = np.multiply(dvec,dS_grow)
    dT_pertox = np.dot(fvec,dT_invest)
    dT_out    = np.multiply(np.add(aflow,dT_pertox), c_tox)
    dT = np.subtract(dT_inflow, dT_out)

    # Assemble state vector
    dX = np.append(np.append(dS,dN),dT)
    return( dX )

def jac_chem_1N1T(t, x, args):
    Nspec = args['size']
    N0, T0 = args['nt_inflow']

    aflow = args["a"]   # Flow rate 'alpha' of chemostat
    rvec = args["r"]    # Growth effect of nutrient, matrix of size [Nspec, Nnutr]
    mvec = args["m"]    # Mortality of toxins, matrix of size [Nspec, Ntox]
    fvec = args["finv"] # Per-species investment [0.0, 1.0] into degradation of each toxin
    Yvec = args["Y"]    # Nutrient Yield [say gram of biomass per mol of nutrient]
    dvec = args["d"]    # Toxin degradation rate [say mol of toxin per gram of biomass]

    ### State variables for abundance and concentrations
    x_spec = x[:Nspec]
    c_nutr = x[Nspec]
    c_tox  = x[(Nspec+1)]

    #### Prelims
    gr_percap = np.multiply(rvec,c_nutr)
    d_percap  = np.add(np.multiply(mvec,c_tox),aflow)
    invest_gr = np.subtract([1.0]*Nspec,fvec)

    #### Compute partial derivatives wrt species abundance
    XF   =  np.multiply(invest_gr,x_spec)
    dSdS =  np.subtract(np.multiply(invest_gr,gr_percap),d_percap) # -> diag
    dSdN =  np.multiply(rvec,XF)
    dSdT = -np.multiply(mvec,x_spec)

    ### Effect on nutrient dynamics
    rY   = np.multiply(Yvec,rvec)
    dNdS = -1.0*np.multiply(rY,c_nutr)
    dNdN = -1.0*np.add(aflow, np.dot(rY,x_spec))
    dNdT = [0.0,]

    ### Effect on toxin
    rd   = np.multiply(rvec,dvec)
    rdX  = np.multiply(rd,x_spec)
    drN  = np.multiply(dvec,gr_percap)
    rdXN = np.multiply(drN,x_spec)

    dTdS = -1.0*np.multiply(c_tox, np.multiply(fvec,drN))
    dTdN = -1.0*np.multiply(c_tox, np.dot(fvec,rdX))
    dTdT = -1.0*np.add(aflow, np.dot(fvec,rdXN))

    # Assemble Jacobian from partial derivatives
    dSdNT = np.vstack((dSdN,dSdT))
    dS    = np.transpose(np.vstack((np.diag(dSdS),dSdNT)))
    dNdNT = np.vstack((np.diag(dNdN),dNdT))
    dN    = np.transpose(np.vstack((dNdS,dNdNT)))
    dTdSN = np.vstack((dTdS,dTdN))
    dT    = np.transpose(np.vstack((dTdSN,np.diag(dTdT))))
    J_out = np.vstack((dS,dN,dT))
    return J_out
