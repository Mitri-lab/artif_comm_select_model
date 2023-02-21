Here the model generating the main results from the IBM.

Code is wrote in python, but designed to be run in unix in order to paralelise jobs for different seeds (random species sets). For instance, for 
the current paper it was run in a workstation as follows:

for seed in {22..26}
do
    nohup nice -n 10 python3 221018_main_model_unix.py $seed > errors$seed.txt &
    sleep 1 
done

To run, whithin the python code one must indicate the propagations and whether they are in random or selection treatment (lines 99 and 102). Propagation should be indicated in
the list "list_prop" (disassembly = "d3", propagule = "p", migrant pool = "m", propagule with invasion = "pip", migrant pool with invasion = "mip", no selection = "n"). The 
treament ( selection = "s", random = "r"), should be indicated within the list "list_treat". For each item in "list_prop" there should be an item in "list_treat", 
so the lenght of both lists must be the same.

There are many coments within the python code, but to know more about the model, the different data structures we use and how results are organised check the "ibm_guide" doc.
