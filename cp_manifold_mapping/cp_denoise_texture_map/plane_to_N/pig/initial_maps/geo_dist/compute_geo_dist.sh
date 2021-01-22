#!/bin/bash

echo "Starting run at: `date`" 

# this script simply loops over all points to compute pairwise geodesic distance between N points. The total number of vertices in the triangulation in this example is N = 4430.
for i in $(seq 0 4429); do ./geo_dist_pt_to_Npts.out ../../partial_pig_refined1.txt $i >> ../partial_pig_refined1_geo_dist.txt ; echo $i; done

echo "Job finished at: `date`"