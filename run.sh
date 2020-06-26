#!/bin/bash
## Ford Fishman

# Set pathing
sourcepath=/Users/fordfishman/GitHub/CommunityImmunity
# sourcepath=/geode2/home/u030/ffishman/Carbonate/thindrives/fordfishman/GitHub/CommunityImmunity
python=$sourcepath/python
r=$sourcepath/r
# output=$sourcepath/output/latency_test/m105_pS106_aH_107_f1

sims=30
a=1e7
l=0.9
b=1.5
pS=1e-6

settings=${sims}sims_a${a}_b${b}_l${l}_pS${pS}

output=$sourcepath/output/$settings

# remove output directory if it exists
if [ -d $output ]; then rm -r $output; fi

mkdir $output
time python3 $python/main.py -o $output -S $sims -pS $pS -a $a -b $b -l $l

time Rscript $r/plotGraphsMultiSim.R -f "$output/summary.csv" -o $output
# time Rscript $r/plotGraphs1Sim.R -f $output/*.csvq