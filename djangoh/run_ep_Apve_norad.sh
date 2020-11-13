#!/usr/bin/bash

echo "-----------------------------------"
echo "Running DJANGOH Simulation for ep Collider!!!"
echo "..."
echo ""


OUTFILE1=outfiles/djangoh.NC.Apve.noRad.20x250_evt.dat
if test -f "$OUTFILE1"; then
	rm -f "$OUTFILE1"
fi
OUTFILE2=outfiles/djangoh.NC.Apve.noRad.20x250_out.dat
if test -f "$OUTFILE2"; then
	rm -f "$OUTFILE2"
fi
OUTFILE3=outfiles/djangoh.NC.Apve.noRad.20x250_smp.dat
if test -f "$OUTFILE3"; then
        rm -f "$OUTFILE3"
fi

#Create file for random number generation
./make_random.py

djangoh < ep.Rad=0.NC.Apve.in > logfiles/ep.Rad=0.NC.Apve.log

echo "Completed Simulation!!!"
echo ""

echo "Making Output ROOT File..."
root -l -b -q 'make_tree.C("djangoh.NC.Apve.noRad.20x250_evt.dat")'
echo "Done!!!"

echo "-----------------------------------"

echo "Making First Smeared ROOT File..."
root -l -b -q 'make_smeared_perfect.C("djangoh.NC.Apve.noRad.20x250_evt")'
echo "Done!!!"

echo "Making Second Smeared ROOT File..."
root -l -b -q 'make_smeared_central.C("djangoh.NC.Apve.noRad.20x250_evt")'
echo "Done!!!"

echo "Making Third Smeared ROOT File..."
root -l -b -q 'make_smeared_handbook.C("djangoh.NC.Apve.noRad.20x250_evt")'
echo "Done!!!"

echo "-----------------------------------"
