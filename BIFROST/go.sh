#!/bin/sh
# V-scan, 181 points over 90 degrees
mcrun --mpi=16 Full_simple.instr INCOH=1 -n1e7  Omega=0,90 -N181 -dV_scan
# Magnon-scan, 181 points over 90 degrees
mcrun --mpi=16 Full_simple.instr INCOH=0 -n1e8  Omega=0,90 -N181 -dMagnon_scan
