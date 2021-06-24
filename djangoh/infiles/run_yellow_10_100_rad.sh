#!/usr/bin/bash
#-----#-----#----#----#----#----#----#----#----#

if [ -z "$1" ]
then
        echo "No job number set."
        echo "Please run as ./run_yellow_10_100_rad.sh jobnumber"
        echo "Exiting..."
	exit 1
fi

#Go into scratch directory
chmod g+rx ${_CONDOR_SCRATCH_DIR}
cd ${_CONDOR_SCRATCH_DIR}

#Make subdirectory and move there
INPUT=$(( 0 + $1 ))
echo $INPUT
DIR=`printf "%04d" $INPUT`
mkdir $DIR
cd $DIR

#Soft links to necessary files
ln -s /eic/data/baraks/dis-reconstruction/djangoh/make_random.py
ln -s /eic/data/baraks/dis-reconstruction/djangoh/make_tree.C
ln -s /eic/data/baraks/dis-reconstruction/djangoh/infiles/ep_yellow_10_100_rad.inp

#Run simulation
echo "start running in directory $PWD"
echo "-----------------------------------"
echo "Running DJANGOH Simulation number $1 for ep Collider!!!"
echo "..."
echo ""

mkdir outfiles

#Make random number file
./make_random.py
cat fort.8
echo ""

djangoh < ep_yellow_10_100_rad.inp > ep_yellow_10_100.log

echo "Completed Simulation number $1!!!"
echo ""

echo "Making Output ROOT File..."
root -l -b -q 'make_tree.C("djangoh.NC.Rad.10x100_evt.dat")'
echo "-----------------------------------"

echo "Done!!!"
echo ""

#Move output files and cleanup
echo "Cleaning Up..."
mv -v ep_yellow_10_100.log /eic/data/baraks/dis-reconstruction/djangoh/logfiles/yellow/Rad/ep_yellow_10_100_${INPUT}.log
mv -v outfiles/djangoh.NC.Rad.10x100_evt.dat   /eic/data/baraks/dis-reconstruction/djangoh/outfiles/yellow/10_100/Rad/djangoh.NC.10x100_${INPUT}.out
mv -v outfiles/djangoh.NC.Rad.10x100_evt.root  /eic/data/baraks/dis-reconstruction/djangoh/outfiles/yellow/10_100/Rad/djangoh.NC.10x100_${INPUT}.root
mv -v outfiles/djangoh.NC.Rad.10x100_out.dat   /eic/data/baraks/dis-reconstruction/djangoh/outfiles/yellow/10_100/Rad/djangoh.NC.10x100_out_${INPUT}.dat



