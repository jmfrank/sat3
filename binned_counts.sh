#!/bin/bash

# Input 1 is chromosome sizes. 
# Input 2 is bed file. Might need to change default column number? 


bedtools makewindows -g $1 -w 2000 | bedtools map -a - -b $2 -o sum
