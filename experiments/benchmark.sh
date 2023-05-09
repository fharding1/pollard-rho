#!/bin/bash

for r in {3..13}
do
  (while IFS=, read -r p n g h x
  do
      out=$(./rho $p $n $g $h 20)
      bc -l <<< "${out}"
  done < "dlps/zp-${r}-digit-prime-subgroups-10000.csv" ) | datamash mean 1
done
