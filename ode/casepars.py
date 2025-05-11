##### CASE PARAMETERS FOR SELECTION SIMULATION
## Used by Selsim_ ...
## Written by Bjorn.Vessman@unil.ch
## 2020.04.15
import seltheory_ODEmodels as smodels

## TODO:
# Implement tuples of several functions
#  -> need to modify the repeated-selection
# Implement cheater case / mutations to growth rate and mortality

### Choose model case
casestr = 'Full'
if casestr=='Full':
    fgrowth = smodels.f_CRbatchHill_full
elif casestr=='Nodeath':
    fgrowth = smodels.f_troff_nodeath
elif casestr=='Notroff':
    fgrowth = smodels.f_notroff_death
elif casestr=='Notroff-Nodeath':
    fgrowth = smodels.f_notroff_nodeath

### Compute _all_ possible co-cultures of the species?
runall = True

### Compute stability after the transfers?
runstability = True

### Base directory for output
basedir = './'

#### Differentiate cheating potential here: allow r, m to mutate or not
# traittuple = ('finv','m','r')
# casepars = {'flag_randcomm':False}  # Draw community parameters at random or not?
# flag_mutrm = 0                      # Allow r, m to mutate or not (cases b, d)
