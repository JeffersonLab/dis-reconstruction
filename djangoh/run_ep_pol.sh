#!/usr/bin/bash

echo "-----------------------------------"
echo "Running DJANGOH Simulation for ep Collider!!!"
echo "..."
echo ""


OUTFILE1=outfiles/djangoh.NC.pol.20x250_evt.dat
if test -f "$OUTFILE1"; then
	rm -f "$OUTFILE1"
fi
OUTFILE2=outfiles/djangoh.NC.pol.20x250_out.dat
if test -f "$OUTFILE2"; then
	rm -f "$OUTFILE2"
fi

djangoh < ep.Rad=1.NC.pol.in > logfiles/ep.Rad=1.NC.pol.log

echo "Completed Simulation!!!"
echo ""

echo "Making Output ROOT File..."
root -l -b -q 'make_tree.C("djangoh.NC.pol.20x250_evt.dat")'
echo "Done!!!"

echo "-----------------------------------"

echo "Making First Smeared ROOT File..."
root -l -b -q 'make_smeared_perfect.C("djangoh.NC.pol.20x250_evt")'
echo "Done!!!"

echo "Making Second Smeared ROOT File..."
root -l -b -q 'make_smeared_central.C("djangoh.NC.pol.20x250_evt")'
echo "Done!!!"

echo "Making Third Smeared ROOT File..."
root -l -b -q 'make_smeared_handbook.C("djangoh.NC.pol.20x250_evt")'
echo "Done!!!"

echo "-----------------------------------"
