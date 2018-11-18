#!/bin/sh
# V-run
mcrun --mpi=16 Full_simple.instr INCOH=1 -n1e8  Omega=0 -dVrun
mcrun --mpi=16 Full_simple.instr INCOH=1 INCOHspread=1 -n1e8 -dVrun_inel
# Magnon-scan, 46 points over 90 degrees
mcrun --mpi=16 Full_simple.instr INCOH=0 -n2e7  Omega=0,90 -N46 -dMagnon_rerun
