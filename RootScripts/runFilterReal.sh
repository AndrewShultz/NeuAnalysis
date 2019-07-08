#!/bin/sh

afile=/data/exp/ARA/2014/unblinded/L1/ARA02/0622/run003799/event003799.root
anOutDir=./

#root -b -q filterReal.C+\(\"$1\",\"$2\"\)
root -b -q filterReal.C+\(\"$afile\",\"$anOutDir\"\)
