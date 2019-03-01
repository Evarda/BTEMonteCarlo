# !/usr/bin/bash
mkdir Data
mkdir Data/gamma
mkdir Data/L
mkdir Data/X
mkdir _Data

make

./MonteCarlo

mv Data _Data/$1_$2