#!/usr/bin/bash 

if [ -z "$1" ]
then
	echo "No job number set."
        echo "Please run as ./smear_5_41.sh jobnumber"
	echo "Exiting..."
	exit 1
fi

echo "-----------------------------------"
echo "Detector Smearing PYTHIA Simulation for ep Collider!!!"
echo "-----------------------------------"
echo "Performing Job $1"
echo "..."
echo ""

echo "Making Output ROOT File..."
root -l -b -q 'make_smeared_perfect.C("outfiles/yellow_report/5_41","ep_5_41_newtune_'$1'")'
echo "Done!!!"
echo ""

