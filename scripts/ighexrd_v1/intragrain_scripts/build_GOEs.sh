#!/bin/sh
# -*- coding: utf-8 -*-

cd /nfs/chess/aux/user/ken38/Ti7_project/ti7-05-3percent/
source activate hexrd_0.5.29_working
#for first load use initial 11, for last load use final 68
#counter set to increment up to the number of grains total (see grains.out to decide)

echo starting eta_ome_maps initial
python generate_eta_ome_maps_parallel_initial.py GOE_ti705.yml ti7-05 initial 11
echo Done building maps

echo starting GOE builder initial
counter=0
while [ $counter -le 792 ]; do
echo grain $counter
python GOE_builder_one_load.py $counter GOE_ti705.yml ti7-05 initial 11
((counter++))
done
echo Done

echo starting eta_ome_maps final
python generate_eta_ome_maps_parallel_final.py GOE_ti705.yml ti7-05 final 68
echo Done building maps

echo starting GOE builder final
counter=0
while [ $counter -le 792 ]; do
echo grain $counter
python GOE_builder_one_load.py $counter GOE_ti705.yml ti7-05 final 68
((counter++))
done
echo Done
