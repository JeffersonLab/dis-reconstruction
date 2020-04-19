#!/usr/bin/csh

# For EIC-SMEAR
source /apps/root/PRO/setroot_CUE
setenv ROOT_INCLUDE_PATH /work/halla/gmp12/baraks/eic-smear/include
setenv LD_LIBRARY_PATH .:/work/halla/gmp12/baraks/eic-smear/build:${LD_LIBRARY_PATH}

# For DJANGOH
setenv LD_LIBRARY_PATH /work/JAM/apps/lhapdf5/lib:${LD_LIBRARY_PATH}
setenv PATH /work/JAM/apps/DJANGOH-4.6.10/bin/:${PATH}

# For pythiaeRHIC
setenv LHAPDF5 /work/JAM/apps/lhapdf5/lib
setenv LD_LIBRARY_PATH /work/halla/gmp12/baraks/PYTHIA-RAD-CORR/install/lib:${LD_LIBRARY_PATH}
setenv PATH /work/halla/gmp12/baraks/PYTHIA-RAD-CORR/install/bin:${PATH}
setenv LHAPATH /work/JAM/apps/lhapdf5/share/lhapdf

#This is the sPHENIX version of pythiaeRHIC. It compiles but gives run-time error
#setenv LD_LIBRARY_PATH /work/halla/gmp12/baraks/pythia6/pythiaeRHIC/install/lib:${LD_LIBRARY_PATH}
#setenv PATH /work/halla/gmp12/baraks/pythia6/pythiaeRHIC/install/bin:${PATH}


