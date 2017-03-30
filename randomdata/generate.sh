#!/bin/bash
djenrandom -b -s -k 10 -m correlated --correlation=-0.1 > corrm0p1.bin
djenrandom -b -s -k 10 -m correlated --correlation=-0.2 > corrm0p2.bin
djenrandom -b -s -k 10 -m correlated --correlation=-0.3 > corrm0p3.bin
djenrandom -b -s -k 10 -m correlated --correlation=-0.4 > corrm0p4.bin
djenrandom -b -s -k 10 -m correlated --correlation=-0.5 > corrm0p5.bin
djenrandom -b -s -k 10 -m correlated --correlation=-0.6 > corrm0p6.bin
djenrandom -b -s -k 10 -m correlated --correlation=-0.7 > corrm0p7.bin
djenrandom -b -s -k 10 -m correlated --correlation=-0.8 > corrm0p8.bin
djenrandom -b -s -k 10 -m correlated --correlation=-0.9 > corrm0p9.bin

djenrandom -b -s -k 10 -m correlated --correlation=0.1 > corr0p1.bin
djenrandom -b -s -k 10 -m correlated --correlation=0.2 > corr0p2.bin
djenrandom -b -s -k 10 -m correlated --correlation=0.3 > corr0p3.bin
djenrandom -b -s -k 10 -m correlated --correlation=0.4 > corr0p4.bin
djenrandom -b -s -k 10 -m correlated --correlation=0.5 > corr0p5.bin
djenrandom -b -s -k 10 -m correlated --correlation=0.6 > corr0p6.bin
djenrandom -b -s -k 10 -m correlated --correlation=0.7 > corr0p7.bin
djenrandom -b -s -k 10 -m correlated --correlation=0.8 > corr0p8.bin
djenrandom -b -s -k 10 -m correlated --correlation=0.9 > corr0p9.bin

djenrandom -b -s -k 10 -m biased --bias=0.1 > bias0p1.bin
djenrandom -b -s -k 10 -m biased --bias=0.2 > bias0p2.bin
djenrandom -b -s -k 10 -m biased --bias=0.3 > bias0p3.bin
djenrandom -b -s -k 10 -m biased --bias=0.4 > bias0p4.bin
djenrandom -b -s -k 10 -m biased --bias=0.5 > bias0p5.bin
djenrandom -b -s -k 10 -m biased --bias=0.6 > bias0p6.bin
djenrandom -b -s -k 10 -m biased --bias=0.7 > bias0p7.bin
djenrandom -b -s -k 10 -m biased --bias=0.8 > bias0p8.bin
djenrandom -b -s -k 10 -m biased --bias=0.9 > bias0p9.bin

djenrandom -b -s -k 10 > uniform.bin

