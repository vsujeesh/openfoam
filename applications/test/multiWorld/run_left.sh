#!/bin/bash
. /home/mattijs/OpenFOAM/OpenFOAM-plus.feature-localWorld/etc/bashrc 
cd /home/mattijs/OpenFOAM/OpenFOAM-plus.feature-localWorld/applications/test/multiWorld
/home/mattijs/OpenFOAM/OpenFOAM-plus.feature-localWorld/platforms/linux64GccDPInt32Debug/bin/laplacianFoam -case ./left -world LEFT 2>&1 | tee run_left.log
read dummy
