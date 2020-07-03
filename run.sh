#!/bin/bash
## Ford Fishman

# Set pathing
sourcepath=/Users/fordfishman/GitHub/CommunityImmunity
# sourcepath=/geode2/home/u030/ffishman/Carbonate/thindrives/fordfishman/GitHub/CommunityImmunity
python=$sourcepath/python
r=$sourcepath/r
# output=$sourcepath/output/latency_test/m105_pS106_aH_107_f1

sims=1
a=1e6
l=0.9
m=1e-7
b=1.2
pS=1e-6
d=0.1
adsp=1e-8
beta=100

settings=${sims}sims_a${a}_adsp${adsp}_b${b}_beta${beta}_d${d}_l${l}_m${m}_pS${pS}

output=$sourcepath/output/$settings

# remove output directory if it exists
if [ -d $output ]; then rm -r $output; fi

mkdir $output
# time python3 $python/main.py -o $output -S $sims -pS $pS -a $a -b $b -l $l
time python3 $python/main.py -s True -o $output -pS $pS -a $a -b $b -l $l -m $m -d $d --adsp $adsp --beta $beta


# time Rscript $r/plotGraphsMultiSim.R -f "$output/summary.csv" -o $output
time Rscript $r/plotGraphs1Sim.R -f $output/