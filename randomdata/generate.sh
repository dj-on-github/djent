#!/bin/bash
djrandom -b -s -k 10 -m correlated --correlation=0.1 > corr0p1.bin
djrandom -b -s -k 10 -m correlated --correlation=0.2 > corr0p2.bin
djrandom -b -s -k 10 -m correlated --correlation=0.3 > corr0p3.bin
djrandom -b -s -k 10 -m correlated --correlation=0.4 > corr0p4.bin
djrandom -b -s -k 10 -m correlated --correlation=0.5 > corr0p5.bin
djrandom -b -s -k 10 -m correlated --correlation=0.6 > corr0p6.bin
djrandom -b -s -k 10 -m correlated --correlation=0.7 > corr0p7.bin
djrandom -b -s -k 10 -m correlated --correlation=0.8 > corr0p8.bin
djrandom -b -s -k 10 -m correlated --correlation=0.9 > corr0p9.bin

djrandom -b -s -k 10 -m biased --bias=0.1 > bias0p1.bin
djrandom -b -s -k 10 -m biased --bias=0.2 > bias0p2.bin
djrandom -b -s -k 10 -m biased --bias=0.3 > bias0p3.bin
djrandom -b -s -k 10 -m biased --bias=0.4 > bias0p4.bin
djrandom -b -s -k 10 -m biased --bias=0.5 > bias0p5.bin
djrandom -b -s -k 10 -m biased --bias=0.6 > bias0p6.bin
djrandom -b -s -k 10 -m biased --bias=0.7 > bias0p7.bin
djrandom -b -s -k 10 -m biased --bias=0.8 > bias0p8.bin
djrandom -b -s -k 10 -m biased --bias=0.9 > bias0p9.bin

djrandom -b -s -k 10 > uniform.bin

