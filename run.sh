#!/bin/bash
## Ford Fishman

# Set pathing
sourcepath=/Users/fordfishman/GitHub/CommunityImmunity
# sourcepath=/geode2/home/u030/ffishman/Carbonate/thindrives/fordfishman/GitHub/CommunityImmunity
python=$sourcepath/python
r=$sourcepath/r
# output=$sourcepath/output/latency_test/m105_pS106_aH_107_f1

settings=m0_pS0
output=$sourcepath/output/$settings

rm -r $output
mkdir $output

sims=1000
# Run main python script

time python3 $python/main.py -o $output/ -S $sims
time Rscript $r/plotGraphsMultiSim.R -f "$output/${sims}sims_${settings}.csv" -o $output
# time Rscript $r/plotGraphs1Sim.R -f $output/*.csvq