#!/usr/bin/env bash
## Ford Fishman

# Set pathing
# sourcepath=/Users/fordfishman/GitHub/CommunityImmunity
sourcepath=$(pwd)
# sourcepath=/geode2/home/u030/ffishman/Carbonate/thindrives/fordfishman/GitHub/CommunityImmunity
python=$sourcepath/python
r=$sourcepath/r
# output=$sourcepath/output/latency_test/m105_pS106_aH_107_f1

sims=1
t=2000
a=1e6
c=0.01
f=0
l=0.9
m=1e-6
b=1.2
pS=1e-6
d=0.1
adsp=1e-8
beta=100
popinit=1e5
phageinit=1e7

settings=${sims}sims_t${t}_a${a}_adsp${adsp}_b${b}_beta${beta}_c${c}_d${d}_f${f}_l${l}_m${m}_pS${pS}_pop${popinit}_phage${phageinit}

output=$sourcepath/output/$settings

# remove output directory if it exists
if [ -d $output ]; then rm -r $output; fi

mkdir $output
# time python3 $python/main.py -o $output -S $sims -pS $pS -a $a -b $b -l $l
time python3 $python/main.py -s True -t $t -o $output -pS $pS -a $a -b $b -c $c -l $l -m $m -f $f -d $d --adsp $adsp --beta $beta --popinit $popinit --phageinit $phageinit


# time Rscript $r/plotGraphsMultiSim.R -f "$output/summary.csv" -o $output
time Rscript $r/plotGraphs1Sim.R -f $output/