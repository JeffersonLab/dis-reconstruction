# dis-reconstruction


eic-smear package
-----------------

Documentation for the eic-smear package can be found [here](https://wiki.bnl.gov/eic/index.php/Computing) and [here](https://wiki.bnl.gov/eic/index.php/Smearing). The most up-to-date version of the eic-smear code is [here](https://gitlab.com/eic/eic-smear).

The EIC versions of various event generators provide output that can be directly fed into the eic-smear software. See [here](https://wiki.bnl.gov/eic/index.php/Simulations#Event_Generators). In particular, look at the information for [PYTHIA](https://wiki.bnl.gov/eic/index.php/PYTHIA) and [DJANGOH](https://wiki.bnl.gov/eic/index.php/DJANGOH).


Working on the RACF (BNL) machines
----------------------------------
If you have either an EIC RACF account or a sPHENIX RACF account, you automatically have access to the software packages mentioned above. Follow these instructions if you have an EIC account: [EIC Environment Setup](https://wiki.bnl.gov/eic/index.php/Computing).

Follow the steps here if you have an sPHENIX account (or a PHENIX account with sPHENIX permissions): [sPHENIX Enviroment Setup](https://wiki.bnl.gov/sPHENIX/index.php/Setup).

If working on with the Sphenix setup, the default LHAPDF directory is

/cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/lhapdf-5.9.1/share/lhapdf/PDFsets

This directory only contains a small number of PDF datasets. It is best to set the LHAPATH variable to

/afs/rhic/eic/share/lhapdf/PDFsets


If you have some difficulty setting up either enviroment using the above instructions, please contact me and I can send you my environmental setup scripts for each configuration.


Working on the JLAB (ifarm) machines
-----------------------------------
If you have access to the Jefferson Lab farm (ifarm), you should source the setup_jlab.sh script provided in this repository. This will allow you to run the eic-smear package. (Need to add instructions for PYTHIA and DJANGOH.)


Working with the Singularity container or Docker image
------------------------------------------------------
You can also use the singularity container (https://github.com/EIC-Detector/Singularity) to run the software -- this works best on linux. I was able to get "Option-2" to work on an Ubuntu box. If you have MacOS or Windows, it may be better to use the virtual box as they suggest (https://github.com/EIC-Detector/Singularity/blob/master/VirtualBox.md). I was able to get that running on a Windows10 machine.

Information on obtaining and using the EIC software Docker image can be found [here](https://eic.gitlab.io/documents/quickstart/#ESCalate).

N.B. It seems that not all the simulation packages listed above are provided by default on the Singularity container or Docker image. I'll contact the EIC software group to resolve this ASAP.


Files provided on this repository
---------------------------------


Contact
--------
Author: [Barak Schmookler](barak.schmookler@stonybrook.edu)



