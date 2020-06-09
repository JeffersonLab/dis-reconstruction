#!/usr/bin/bash

for i in {10..24}
do
        echo "-----------------------------------"
        echo "Running DJANGOH Simulation number $i for ep Collider!!!"
        echo "..."
        echo ""

	djangoh < ep_hera.in > logfiles/ep_hera_$i.log

        echo "Completed Simulation number $i!!!"
        echo ""

        echo "Making Output ROOT File..."
        root -l -b -q 'make_tree.C("djangoh_hera_evt.dat")'
        echo "-----------------------------------"

        mv outfiles/djangoh_hera_evt.dat outfiles/djangoh_hera_$i.out
        mv outfiles/djangoh_hera_evt.root outfiles/djangoh_hera_$i.root
	mv outfiles/djangoh_hera_out.dat outfiles/djangoh_hera_out_$i.dat

done

echo "Done!!!"

