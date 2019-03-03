# !/usr/bin/bash
mkdir Data
mkdir Data/ScatRates
mkdir Data/ScatRates/gamma
mkdir Data/ScatRates/L
mkdir Data/ScatRates/X
mkdir _Data

make

./MonteCarlo

mv Data _Data/$1_$2