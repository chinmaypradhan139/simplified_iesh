#!/bin/bash
mkdir 1
cp a.out 1
cp job.sh 1
cp fort.23 1
cp raw_x.txt 1
cp raw_w.txt 1
mv input.txt 1
cd 1
qsub job.sh	
cd .. 
