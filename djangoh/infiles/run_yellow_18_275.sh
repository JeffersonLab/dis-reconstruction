#!/usr/bin/bash

for i in {0..9}
do
        echo "-----------------------------------"
        echo "Running DJANGOH Simulation number $i for ep Collider!!!"
        echo "..."
        echo ""

	djangoh < ep_yellow_18_275.inp > logfiles/ep_yellow_18_275_$i.log

        echo "Completed Simulation number $i!!!"
        echo ""

        echo "Making Output ROOT File..."
        root -l -b -q 'make_tree.C("djangoh.NC.18x275_evt.dat")'
        echo "-----------------------------------"

        mv outfiles/djangoh.NC.18x275_evt.dat outfiles/djangoh.NC.18x275_$i.out
        mv outfiles/djangoh.NC.18x275_evt.root outfiles/djangoh.NC.18x275_$i.root
	mv outfiles/djangoh.NC.18x275_out.dat outfiles/djangoh.NC.18x275_out_$i.dat

done

echo "Done!!!"

