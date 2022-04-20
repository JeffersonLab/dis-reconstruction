# dis-reconstruction

<br/>

eic-smear package
-----------------

Documentation for the eic-smear package can be found [here](https://eic.github.io/software/eicsmear.html) and [here](https://wiki.bnl.gov/eic/index.php/Smearing). The most up-to-date version of the eic-smear code is [here](https://gitlab.com/eic/eic-smear).

The EIC versions of various event generators provide output that can be directly fed into the eic-smear software. See [here](https://wiki.bnl.gov/eic/index.php/Simulations#Event_Generators). In particular, look at the information for [PYTHIA](https://eic.github.io/software/pythia6.html) and [DJANGOH](https://eic.github.io/software/djangoh.html).

<br/>

Working on the RACF (BNL) machines
----------------------------------
If you have either a BNL RACF account, you automatically have access to the software packages mentioned above. 

<BR/>

Follow these instructions if you have an EIC account: [EIC Environment Setup](https://wiki.bnl.gov/eic/index.php/Computing).

The EIC environmental setup described on the EIC wiki page above will link to the 'pro' setup. To link to the 'dev' version of eic-smear along with ROOT6, set your environment as follows:

> setenv EIC_LEVEL dev 
> source /cvmfs/eic.opensciencegrid.org/x8664_sl7/MCEG/releases/etc/eic_cshrc.csh

(N.B. The above commands should also allow you to setup the EIC environment if you are working on a STAR account.)

<br/>

Follow these steps if you have an sPHENIX account (or a PHENIX account with sPHENIX permissions): [sPHENIX Enviroment Setup](https://wiki.bnl.gov/sPHENIX/index.php/Setup).

If working on with the Sphenix setup, the default LHAPDF directory is

> /cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/lhapdf-5.9.1/share/lhapdf/PDFsets

This directory only contains a small number of PDF datasets. It is best to set the LHAPATH variable to

> /afs/rhic/eic/share/lhapdf/PDFsets

<br/>

If you have some difficulty setting up either enviroment using the above instructions, please contact me and I can send you my environmental setup scripts for each configuration.

<br/>

Working on the JLAB (ifarm) machines
-----------------------------------
To run on the JLAB (ifarm) machines, you also simply need to set the environment as above:

> setenv EIC_LEVEL dev   
> source /cvmfs/eic.opensciencegrid.org/x8664_sl7/MCEG/releases/etc/eic_cshrc.csh

For some versions of the ifarm, you may need to first access the singularity container:

> module load singularity   
> singularity shell -B /cvmfs:/cvmfs /cvmfs/eic.opensciencegrid.org/singularity/rhic_sl7_ext   
> export EIC_LEVEL="dev"   
> source /cvmfs/eic.opensciencegrid.org/x8664_sl7/MCEG/releases/etc/eic_bash.sh  

Some additional information can be found [here](https://eic.github.io/software/escalate_singularity_1.html).

~~If you have access to the Jefferson Lab farm (ifarm), you should source the 'setup_jlab.csh' script provided in this repository. This will allow you to run the EIC versions of the Pythia6 and the DJANGOH event generators on the ifarm; it will also provide access to the eic-smear library installed on the ifarm.~~

<br/>


Working with the Singularity container or Docker image
------------------------------------------------------
You can also use a Singularity container to access the EIC software. General instructions can be found [here](https://eic.github.io/software/eicsmear_generators_singularity.html).

A ready-made image is provided by the sPHENIX and ECCE groups [here](https://github.com/eic/Singularity) to run the software. This works best on linux -- and 'option 1' is preferable if you just want to run the EIC generators. If you have MacOS or Windows, it may be better to use the [virtual box](https://github.com/eic/Singularity/blob/master/VirtualBox.md). I was able to get that running on a Windows10 machine.

Information on accessing the EIC software through the ESCalate package can be found [here](https://eic.github.io/software/escalate.html).

<br/>

Step-by-Step simulation tutorial
--------------------------------
A step-by-step simulation tutorial can be found [here](https://github.com/cipriangal/eicGenTutorials). Follow the [documentation](https://drive.google.com/file/d/1RiiveVGhMEIzmtNQYr4a2U0ghGo0w3Hl/view?usp=sharing) that is linked to in that repository.

<br/>

Files provided in this repository
---------------------------------
Example detectors are provided in the "detectors" folder. Many are taken from the detectors provided [here](https://github.com/eic/eicsmeardetectors).

In the "pythia" folder, the script "run_ep.sh" will run a PYTHIA6 simulation and then smear the output for each detector in the "detectors" folder. The simulation output files are saved in the "outfiles" subfolder. Depending on which account the user is working on (i.e. sphenix, eic, jlab), the user will have to link to the correct detector directory in "make_smeared_...". The file "pythia_ep.sub" allows the user to run the simulation on the RACF (BNL) batch farm. An example script that shows how to smear multiple generator-level files is also provided in the "smear_matrix" subfolder.

The "djangoh" folder has a similar structure as the "pythia" folder.

The "analysis" folder gives several examples of analyzing the smeared (and generator-level) output. Simple ROOT macros can be found in the top directory; a more complex analysis code is in the "compiled" subfolder. A different "Makefile" is needed depending on the account (i.e. sphenix, eic, jlab) the user is working on. The "cross_section" and "kinematic_maps" subfolders provide scripts used for studies included in the EIC yellow report.

<br/>


Minimum bias simulation data located on the RACF
------------------------------------------------
|Data Set| Generator | Beam Energies        | Run Information                                                                 | Number of Events | Int. Luminosity       |
|:-------|:---------:|:--------------------:|:-------------------------------------------------------------------------------:|:----------------:|:---------------------:| 
|1       | Pythia6   | 5x41   GeV e-p       | Q<sup>2</sup> > 0.5 GeV<sup>2</sup>; NC unpolarized; QED Radiation OFF          |  100 million     | 0.14 fb<sup>-1</sup>  |
|2       | Pythia6   | 5x41   GeV e-p       | Q<sup>2</sup> > 3.0 GeV<sup>2</sup>; NC unpolarized; QED Radiation OFF          |  100 million     | 0.96 fb<sup>-1</sup>  |
|3       | Pythia6   | 5x100  GeV e-p       | Q<sup>2</sup> > 0.5 GeV<sup>2</sup>; NC unpolarized; QED Radiation OFF          |  15  million     | 0.016 fb<sup>-1</sup> |
|4       | Pythia6   | 10x100 GeV e-p       | Q<sup>2</sup> > 0.5 GeV<sup>2</sup>; NC unpolarized; QED Radiation OFF          |  11  million     | 9.9e-3 fb<sup>-1</sup>|
|5       | Pythia6   | 10x110 GeV e-p       | Q<sup>2</sup> > 0.5 GeV<sup>2</sup>; NC unpolarized; QED Radiation OFF          |  15  million     | 0.013 fb<sup>-1</sup> |
|6       | Pythia6   | 18x110 GeV e-p       | Q<sup>2</sup> > 0.5 GeV<sup>2</sup>; NC unpolarized; QED Radiation OFF          |  15  million     | 0.011 fb<sup>-1</sup> |
|7       | Pythia6   | 18x275 GeV e-p       | Q<sup>2</sup> > 0.5 GeV<sup>2</sup>; NC unpolarized; QED Radiation OFF          |  15  million     | 9.0e-3 fb<sup>-1</sup>|
|8       | Pythia6   | 27.5x920 GeV e+p     | Q<sup>2</sup> > 1.5 GeV<sup>2</sup>; NC unpolarized; QED Radiation OFF          |  10  million     | 0.011 fb<sup>-1</sup> |
|9       | Djangoh   | 5x41   GeV e-p       | Q<sup>2</sup> > 0.5 GeV<sup>2</sup>; NC unpolarized; QED Radiation OFF          |  ~10 million     | 0.014 fb<sup>-1</sup> |
|10      | Djangoh   | 5x100  GeV e-p       | Q<sup>2</sup> > 0.5 GeV<sup>2</sup>; NC unpolarized; QED Radiation OFF          |  ~10 million     | 0.011 fb<sup>-1</sup> |
|11      | Djangoh   | 10x100 GeV e-p       | Q<sup>2</sup> > 0.5 GeV<sup>2</sup>; NC unpolarized; QED Radiation OFF          |  ~10 million     | 9.1e-3 fb<sup>-1</sup>|
|12      | Djangoh   | 18x275 GeV e-p       | Q<sup>2</sup> > 0.5 GeV<sup>2</sup>; NC unpolarized; QED Radiation OFF          |  ~10 million     | 6.6e-3 fb<sup>-1</sup>|
|13      | Djangoh   | 27.6x920 GeV e+p     | Q<sup>2</sup> > 1.5 GeV<sup>2</sup>; NC unpolarized; QED Radiation OFF          |  ~2.5 million    | 3.5e-3 fb<sup>-1</sup>|
|14      | Pythia6   | 5x41   GeV e-p       | Q<sup>2</sup> down to photo-production limit; NC unpolarized; QED Radiation OFF |  500 million     | 6.3e-3 fb<sup>-1</sup>|
|15      | Pythia6   | 10x100 GeV e-p       | Q<sup>2</sup> down to photo-production limit; NC unpolarized; QED Radiation OFF |  300 million     | 2.3e-3 fb<sup>-1</sup>|
|16      | Pythia6   | 18x275 GeV e-p       | Q<sup>2</sup> down to photo-production limit; NC unpolarized; QED Radiation OFF |  300 million     | 1.7e-3 fb<sup>-1</sup>|
|17      | Djangoh   | 10x100 GeV e-p       | Q<sup>2</sup> > 0.5 GeV<sup>2</sup>; NC unpolarized; QED Radiation ON           |  ~15 million     | 0.013 fb<sup>-1</sup> |
|18      | Pythia8   | 18x275 GeV e-p       | Q<sup>2</sup> > 1.0 GeV<sup>2</sup>; NC unpolarized                             |  ~15 million     | 0.022 fb<sup>-1</sup> |
  
  
Some comments on the generator kinematic limits when QED radiation is OFF. Pythia6 always applies a hard cut of W > 2 GeV. Djangoh requires Q<sup>2</sup> > 0.2 GeV<sup>2</sup>, and it always applies a hard cut of W > ~3.38 GeV. Pythia8 seems to require Q<sup>2</sup> > 1.0 GeV<sup>2</sup>.

<br/>

The paths on the RACF machines to the above data-sets are as follows:

|Data Set| Path                                                              |
|:-------|:------------------------------------------------------------------|
|1       |/gpfs02/eic/baraks/pythia/outfiles/yellow_report/5_41              |
|2       |/gpfs02/eic/baraks/pythia/outfiles/other_studies/5_41/higher_Q2    |
|3       |/gpfs02/eic/baraks/pythia/outfiles/yellow_report/5_100             |
|4       |/gpfs02/eic/baraks/pythia/outfiles/yellow_report/10_100            |
|5       |/gpfs02/eic/baraks/pythia/outfiles/yellow_report/10_110            |
|6       |/gpfs02/eic/baraks/pythia/outfiles/yellow_report/18_110            |
|7       |/gpfs02/eic/baraks/pythia/outfiles/yellow_report/18_275            |
|8       |/gpfs02/eic/baraks/pythia/outfiles/yellow_report/hera              |
|9       |/gpfs02/eic/baraks/djangoh/outfiles/yellow/5_41                    |
|10      |/gpfs02/eic/baraks/djangoh/outfiles/yellow/5_100                   |
|11      |/gpfs02/eic/baraks/djangoh/outfiles/yellow/10_100                  |
|12      |/gpfs02/eic/baraks/djangoh/outfiles/yellow/18_275                  |
|13      |/gpfs02/eic/baraks/djangoh/outfiles/hera                           |
|14      |/gpfs02/eic/baraks/pythia/outfiles/other_studies/5_41              |
|15      |/gpfs02/eic/baraks/pythia/outfiles/other_studies/10_100/local_build|
|16      |/gpfs02/eic/baraks/pythia/outfiles/other_studies/18_275/fullQ2     |
|17      |/gpfs02/eic/baraks/djangoh/outfiles/yellow/10_100/Rad              |
|18      |/gpfs02/eic/baraks/pythia8/basic_DIS/output                        |

<br/>

Contact
--------
Author: [Barak Schmookler](mailto:barak.schmookler@stonybrook.edu)



