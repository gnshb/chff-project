#!/bin/bash

g++ -lm -fopenmp main.cpp;
./a.out | tee ad.txt;
bash new.sh ad.txt;
