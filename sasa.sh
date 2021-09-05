#!/bin/bash
cat << 'EOF2' > vdwradii.dat 
???  Na1   0.35
???  Cl1   0.35
???  Na    0.35
???  Cl    0.35
EOF2

GRO=$1
Xx=$2
echo "$Xx" | gmx sasa -f $GRO -s $GRO -probe 0 -ndots 1024 -tv V_sasa.xvg &> /dev/null
tail -n 1 V_sasa.xvg | awk '{print $2}' #| tee volume_"$GRO".txt
rm -rf area.xvg vdwradii.dat V_sasa.xvg
