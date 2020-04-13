#!/usr/bin/bash

# For EIC-SMEAR
source /apps/root/PRO/setroot_CUE
setenv ROOT_INCLUDE_PATH /work/halla/gmp12/baraks/eic-smear/include
setenv LD_LIBRARY_PATH /work/halla/gmp12/baraks/eic-smear/build:${LD_LIBRARY_PATH}

# For PythiaeRHIC


