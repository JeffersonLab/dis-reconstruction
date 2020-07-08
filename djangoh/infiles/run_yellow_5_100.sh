#!/usr/bin/bash

for i in {5..5}
do
        echo "-----------------------------------"
        echo "Running DJANGOH Simulation number $i for ep Collider!!!"
        echo "..."
        echo ""

	djangoh < ep_yellow_5_100.inp > logfiles/ep_yellow_5_100_$i.log

        echo "Completed Simulation number $i!!!"
        echo ""

        echo "Making Output ROOT File..."
        root -l -b -q 'make_tree.C("djangoh.NC.5x100_evt.dat")'
        echo "-----------------------------------"

        mv outfiles/djangoh.NC.5x100_evt.dat outfiles/djangoh.NC.5x100_$i.out
        mv outfiles/djangoh.NC.5x100_evt.root outfiles/djangoh.NC.5x100_$i.root
	mv outfiles/djangoh.NC.5x100_out.dat outfiles/djangoh.NC.5x100_out_$i.dat

done

echo "Done!!!"

