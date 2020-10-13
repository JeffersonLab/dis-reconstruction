#!/usr/bin/bash

for i in {0..0}
do
        echo "-----------------------------------"
        echo "Running DJANGOH Simulation number $i for ep Collider!!!"
        echo "..."
        echo ""

	OUTFILE1=outfiles/NC.5x41_pol_neg_smp.dat
	if test -f "$OUTFILE1"; then
        	rm -f "$OUTFILE1"
	fi

	djangoh < ep_yellow_5_41_pol_neg.inp > logfiles/ep_yellow_5_41_pol_neg$i.log

        echo "Completed Simulation number $i!!!"
        echo ""

        echo "Making Output ROOT File..."
        root -l -b -q 'make_tree.C("NC.5x41_pol_neg_evt.dat")'
        echo "-----------------------------------"

        mv outfiles/NC.5x41_pol_neg_evt.dat outfiles/NC.5x41_pol_neg_$i.out
        mv outfiles/NC.5x41_pol_neg_evt.root outfiles/NC.5x41_pol_neg_$i.root
	mv outfiles/NC.5x41_pol_neg_out.dat outfiles/NC.5x41_out_pol_neg_$i.dat

done

echo "Done!!!"

