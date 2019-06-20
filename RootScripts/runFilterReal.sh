#!/bin/sh

afile=/data/exp/ARA/2014/unblinded/L1/ARA02/0622/run003799/event003799.root

root -b -q filter.C+\(\"$1\"\)
