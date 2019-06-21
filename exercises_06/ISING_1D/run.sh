#!/usr/bin/env bash

for t in $(seq 0.5 0.01 2.0)
 do
  ./Monte_Carlo_ISING_1D.exe $t
done
