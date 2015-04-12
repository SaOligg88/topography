#!/bin/bash

# Usage: ./x.1_prepSurf.sh [input surface path] [output prefix]

# SurfPatch

mris_convert ${1} ./tmp.asc
tail -n +2 ./tmp.asc > ./tmp2.asc
sed -e 's/\( 0\)*$//g' ./tmp2.asc > ${2}.asc
rm -f ./tmp.asc ./tmp2.asc

