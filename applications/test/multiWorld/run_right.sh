#!/bin/bash
. /home/mattijs/OpenFOAM/OpenFOAM-plus.feature-localWorld/etc/bashrc 
cd /home/mattijs/OpenFOAM/OpenFOAM-plus.feature-localWorld/applications/test/multiWorld
/home/mattijs/OpenFOAM/OpenFOAM-plus.feature-localWorld/platforms/linux64GccDPInt32Debug/bin/laplacianFoam -case ./right -world RIGHT 2>&1 | tee run_right.log
read dummy
