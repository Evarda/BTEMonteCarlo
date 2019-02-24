# !/usr/bin/bash
mkdir Data
mkdir _Data

make

./MonteCarlo

mv Data _Data/$1_$2