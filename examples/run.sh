#!/bin/bash

bin=/home/tienhung/Softwares/CG_MD

$bin/runmd -p 1l2x.pdb -b hbond.dat -k stack.dat -m $bin/src/maxi_explicit -u $bin/uvv/pmf_Mg_P_1264 -T 25 -M 0.001 -K 0.1 -s 5000000000 -v "700 700 700" -a 20000 -e 20000 -C 30 -S $RANDOM -z 29 -d 10
