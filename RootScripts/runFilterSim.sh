#!/bin/sh

Pwd=`pwd`

# Comment out all of either (1) or (2)

### (1) execute with arguments
#abspath=`realpath $2`
#root -b -q filterSim.C+
#cd $ARASIM #must be in AraSim directory at execution time to prevent crash
#root -b -q $Pwd/filterSim.C+\(\"$1\",\"$abspath\"\)

### (2) execute without arguments
afile=/data/user/aschultz/filters/Similarity/Simulations/randLocs/output/condensedAraSimOut.root
#afile=/data/user/pfendner/AraSim/A23/A2_final/AraOut.setup_station2_E17.5.run21.root
anOutDir=./
abspath=`realpath $anOutDir`
root -b -q filterSim.C+
cd $ARASIM #must be in AraSim directory at execution time to prevent crash
root -b -q $Pwd/filterSim.C+\(\"$afile\",\"$abspath\"\)
