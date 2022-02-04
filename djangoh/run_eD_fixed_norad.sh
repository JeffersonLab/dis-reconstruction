#!/usr/bin/bash

echo "-----------------------------------"
echo "Running DJANGOH Simulation for eD Fixed Target!!!"
echo "..."
echo ""


OUTFILE1=outfiles/djangoh.NC.D.fixed_evt.dat
if test -f "$OUTFILE1"; then
	rm -f "$OUTFILE1"
fi
OUTFILE2=outfiles/djangoh.NC.D.fixed_out.dat
if test -f "$OUTFILE2"; then
	rm -f "$OUTFILE2"
fi
OUTFILE3=outfiles/djangoh.NC.D.fixed_smp.dat
if test -f "$OUTFILE3"; then
        rm -f "$OUTFILE3"
fi

djangoh < ep.Rad=0.NC.D.fixed.in > logfiles/ep.Rad=0.NC.D.fixed.log

echo "Completed Simulation!!!"
echo ""

echo "Making Output ROOT File..."
root -l -b -q 'make_tree.C("djangoh.NC.D.fixed_evt.dat")'
echo "Done!!!"

echo "-----------------------------------"

#echo "Making First Smeared ROOT File..."
#root -l -b -q 'make_smeared_perfect.C("djangoh.NC.20x250_evt")'
#echo "Done!!!"

#echo "Making Second Smeared ROOT File..."
#root -l -b -q 'make_smeared_central.C("djangoh.NC.20x250_evt")'
#echo "Done!!!"

#echo "Making Third Smeared ROOT File..."
#root -l -b -q 'make_smeared_handbook.C("djangoh.NC.20x250_evt")'
#echo "Done!!!"

#echo "-----------------------------------"
