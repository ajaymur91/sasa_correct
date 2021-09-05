#!/bin/bash
GRO=$1
n=$(awk 'NR==2 {print $1}' $GRO)
N_ions=$(echo $n/2 | bc)

sed -n '1,2p' $GRO > rev1$GRO
sed '1,2d; $d' $GRO | tac >> rev1$GRO
sed -n '$p' $GRO >> rev1$GRO
echo system | gmx trjconv -f rev1$GRO -o rev$GRO -s rev1$GRO
rm -rf rev1$GRO
