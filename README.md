# dis-reconstruction

<br/>

eic-smear package
-----------------

Documentation for the eic-smear package can be found [here](https://wiki.bnl.gov/eic/index.php/Monte_Carlo_and_Smearing) and [here](https://wiki.bnl.gov/eic/index.php/Smearing). The most up-to-date version of the eic-smear code is [here](https://gitlab.com/eic/eic-smear).

The EIC versions of various event generators provide output that can be directly fed into the eic-smear software. See [here](https://wiki.bnl.gov/eic/index.php/Simulations#Event_Generators). In particular, look at the information for [PYTHIA](https://wiki.bnl.gov/eic/index.php/PYTHIA) and [DJANGOH](https://wiki.bnl.gov/eic/index.php/DJANGOH).

<br/>

Working on the RACF (BNL) machines
----------------------------------
If you have either a BNL RACF account, you automatically have access to the software packages mentioned above. 

<BR/>

Follow these instructions if you have an EIC account: [EIC Environment Setup](https://wiki.bnl.gov/eic/index.php/Computing).

The EIC environment setup described on the EIC wiki page above will link to a ROOT5 build. It is better instead to link to ROOT6 by setting your environment as follows:

setenv EIC_LEVEL pro

source /afs/rhic.bnl.gov/eic/restructured/etc/eic_cshrc.csh

(N.B. The above commands should also allow you to setup the EIC environment if you are working on a STAR account.)

<br/>

Follow these steps if you have an sPHENIX account (or a PHENIX account with sPHENIX permissions): [sPHENIX Enviroment Setup](https://wiki.bnl.gov/sPHENIX/index.php/Setup).

If working on with the Sphenix setup, the default LHAPDF directory is

/cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/lhapdf-5.9.1/share/lhapdf/PDFsets

This directory only contains a small number of PDF datasets. It is best to set the LHAPATH variable to

/afs/rhic/eic/share/lhapdf/PDFsets

<br/>

If you have some difficulty setting up either enviroment using the above instructions, please contact me and I can send you my environmental setup scripts for each configuration.

<br/>

Working on the JLAB (ifarm) machines
-----------------------------------
If you have access to the Jefferson Lab farm (ifarm), you should source the 'setup_jlab.csh' script provided in this repository. This will allow you to run the EIC versions of the Pythia6 and the DJANGOH event generators on the ifarm; it will also provide access to the eic-smear library installed on the ifarm.

<br/>


Working with the Singularity container or Docker image
------------------------------------------------------
You can also use the singularity container (https://github.com/EIC-Detector/Singularity) to run the software -- this works best on linux. I was able to get "Option-2" to work on an Ubuntu box. If you have MacOS or Windows, it may be better to use the virtual box as they suggest (https://github.com/EIC-Detector/Singularity/blob/master/VirtualBox.md). I was able to get that running on a Windows10 machine.

Information on obtaining and using the EIC software Docker image can be found [here](https://eic.gitlab.io/documents/quickstart/#ESCalate).

N.B. It seems that not all the simulation packages listed above are provided by default on the Singularity container or Docker image. I'll contact the EIC software group to resolve this ASAP.

<br/>


Files provided on this repository
---------------------------------
Example detectors are provided in the "detectors" folder.

In the "pythia" folder, the script "run_ep.sh" will run a PYTHIA6 simulation and then smear the output for each detector in the "detectors" folder. The simulation output files are saved in the "outfiles" subfolder. Depending on which account the user is working on(i.e. sphenix, eic, jlab), the user will have to link to the correct detector directory in "make_smeared_...". The file "pythia_ep.sub" allows the user to run the simulation on the RACF (BNL) batch farm.

The "analysis" folder gives several examples of analyzing the smeared (and generator-level) output. Simple ROOT macros can be found in the top directory; a more complex analysis code is in the "compiled" subfolder. A different "Makefile" is needed depending on the account (i.e. sphenix, eic, jlab) the user is working on.

<br/>

Contact
--------
Author: [Barak Schmookler](mailto:barak.schmookler@stonybrook.edu)



