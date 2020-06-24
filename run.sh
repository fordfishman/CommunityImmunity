#!/bin/bash
## Ford Fishman

# Set pathing
sourcepath=/Users/fordfishman/GitHub/CommunityImmunity
# sourcepath=/geode2/home/u030/ffishman/Carbonate/thindrives/fordfishman/GitHub/CommunityImmunity
python=$sourcepath/python
r=$sourcepath/r
# output=$sourcepath/output/latency_test/m105_pS106_aH_107_f1

settings=a1e7_b1.5_l0.9_pS1e-6
output=$sourcepath/output/$settings

rm -r $output
mkdir $output

sims=30
a=1e7
l=0.9
b=1.5
pS=1e-6
# Run main python script

# time python3 $python/main.py -o $output/ -S $sims
time python3 $python/main.py -o $output/ -S $sims -pS $pS -a $a -b $b -l $l

time Rscript $r/plotGraphsMultiSim.R -f "$output/${sims}sims_${settings}.csv" -o $output
# time Rscript $r/plotGraphs1Sim.R -f $output/*.csvq