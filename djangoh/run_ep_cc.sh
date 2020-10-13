#!/usr/bin/bash

echo "-----------------------------------"
echo "Running DJANGOH Simulation for e+p (positron-proton) Collider!!!"
echo "..."
echo ""


OUTFILE1=outfiles/djangoh.CC.Rad.18x275_evt.dat
if test -f "$OUTFILE1"; then
	rm -f "$OUTFILE1"
fi
OUTFILE2=outfiles/djangoh.CC.Rad.18x275_out.dat
if test -f "$OUTFILE2"; then
	rm -f "$OUTFILE2"
fi
OUTFILE3=outfiles/djangoh.CC.Rad.18x275_smp.dat
if test -f "$OUTFILE3"; then
        rm -f "$OUTFILE3"
fi

djangoh < ep.Rad=1.CC.in > logfiles/ep.Rad=1.CC.log

echo "Completed Simulation!!!"
echo ""

echo "Making Output ROOT File..."
root -l -b -q 'make_tree.C("djangoh.CC.Rad.18x275_evt.dat")'
echo "Done!!!"

echo "-----------------------------------"

#echo "Making First Smeared ROOT File..."
#root -l -b -q 'make_smeared_perfect.C("djangoh.CC.Rad.18x275_evt")'
#echo "Done!!!"

#echo "Making Second Smeared ROOT File..."
#root -l -b -q 'make_smeared_central.C("djangoh.CC.Rad.18x275_evt")'
#echo "Done!!!"

#echo "Making Third Smeared ROOT File..."
#root -l -b -q 'make_smeared_handbook.C("djangoh.CC.Rad.18x275_evt")'
#echo "Done!!!"

echo "-----------------------------------"
