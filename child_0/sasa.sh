#!/bin/bash
cat << 'EOF2' > vdwradii.dat 
???  Na1   0.35
???  Cl1   0.35
???  Na    0.35
???  Cl    0.35
EOF2

GRO=$1
Xx=$2
echo "atomname $Xx" | gmx select -f $GRO -s $GRO -on index.ndx &> /dev/null
echo "0" | gmx sasa -f $GRO -s $GRO -probe 0 -ndots 1024 -tv V_sasa.xvg -n index.ndx &> /dev/null
tail -n 1 V_sasa.xvg | awk '{print $2}' #| tee volume_"$GRO".txt
rm -rf area.xvg vdwradii.dat V_sasa.xvg
