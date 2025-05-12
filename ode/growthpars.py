##### BASIC GROWTH PARAMETERS TO INITIALISE COMMUNITY
## Used by Selsim_ ...
## Written by Bjorn.Vessman@unil.ch
## 2020.02.18
import numpy as np

r0 = (0.01,0.1)
m0 = (1.0e-4,1.0e-3)
# f0_pars = (0.01, 0.4)
Y0 = (np.log(1.0e-3),np.log(5.0))
d0 = (np.log(1.0e-4),np.log(5.0))

### For now, the following are not drawn randomly
KN0 = 25.0
KT0 = 25.0
cHill = 2.0
fpow  = (1.0,1.0)
fracsparse = (0.5, 0.5, 0.5)
sparseflag = (True, True, True)

#### Bimodal distribution of parameters from Nov-20
#Y0 = (np.log(1.0e-6),np.log(1.0e1), np.log(1.0e-2),np.log(1.0e1)) # Mean, dev, mean dev
#d0 = (np.log(1.0e-8),np.log(1.0e1), np.log(1.0e-6),np.log(1.0e1))
#fracY1 = 0.6
#fracd1 = 0.6

def f0(Ntox):
    fout = (0.0,np.divide(0.6,Ntox))
    return( fout )
