#!/bin/bash
## Ford Fishman

# Set pathing
sourcepath=/Users/fordfishman/GitHub/CommunityImmunity
python=$sourcepath/python
r=$sourcepath/r
output=$sourcepath/output/labmeeting20200428/m105_pS106_aH_107_fP

# Run main python script

time python3 $python/main.py -o $output
time Rscript $r/plotGraphs.R -f $output 