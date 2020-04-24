#!/bin/bash
## Ford Fishman

# Set pathing
sourcepath=/Users/fordfishman/GitHub/CommunityImmunity
python=$sourcepath/python
r=$sourcepath/r
output=$sourcepath/output

# Run main python script

time python3 $python/main.py -o $output/moregraphs1
# time Rscript $r/plotGraphs.R -f $output/test5.csv -o $output/test5.png