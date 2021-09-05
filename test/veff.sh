#!/bin/bash
#PROC=1
#TRIALS=${2:-10}
WDIR=$(pwd)
GRO=$1
shape=${3:-triclinic}
Si=${2:-10000}
export PLUMED_MAXBACKUP=-1
export GMX_MAXBACKUP=-1
#echo "$shape box"


mkdir -p tempdir
cd tempdir

######################################################################

cp $WDIR/$GRO $WDIR/tempdir/$GRO

# make index file
echo q | gmx make_ndx -f $GRO &> /dev/null

# pbc Wrap cluster 
echo "$(cat << EOF 
Ion: GROUP NDX_FILE=index.ndx NDX_GROUP=Ion
WRAPAROUND ATOMS=Ion AROUND=1
DUMPATOMS FILE=dump.gro ATOMS=Ion
EOF
)" | plumed driver --igro $GRO --plumed /dev/stdin &> /dev/null 

# Center cluster in box
gmx trjconv -f dump.gro -s $GRO -center -o box.gro &> /dev/null << EOF
Ion
Ion
EOF

# Pad box (for inserting a pair of ions)
echo "System" | gmx editconf -princ -f box.gro -bt $shape -d 0.7 -o pad.gro &> /dev/null
sed '1,2d; $d' pad.gro | awk '{print $4"\t"}' | tr -d '\n' | awk '{print $0 "\n"}' > $WDIR/X
sed '1,2d; $d' pad.gro | awk '{print $5"\t"}' | tr -d '\n' | awk '{print $0 "\n"}' > $WDIR/Y
sed '1,2d; $d' pad.gro | awk '{print $6"\t"}' | tr -d '\n' | awk '{print $0 "\n"}' > $WDIR/Z
n=$(awk 'NR==2 {print $1}' pad.gro)
Lx=$(tail -n 1 pad.gro | awk '{print $1}')
Ly=$(tail -n 1 pad.gro | awk '{print $2}')
Lz=$(tail -n 1 pad.gro | awk '{print $3}')

echo $n > $WDIR/Np
tail -n 1 pad.gro | awk '{print $0}' > $WDIR/L
############################## SASA NA
cd $WDIR
Rscript --vanilla veff.R $Si $((n+2)) $Si | tee volume_$GRO.txt
########################################

# make child gro
n=$(awk 'NR==2 {print $1}' tempdir/pad.gro)

sed -n '1p' tempdir/pad.gro > child_$GRO
echo $((n+2)) >> child_$GRO
sed '1,2d; $d' tempdir/pad.gro >> child_$GRO
  
Na=$(head -n 1 -q {X..Z}C | awk '{print $1"\t"}' | tr -d '\n' | awk '{print $0"\n"}')
Cl=$(head -n 1 -q {X..Z}C | awk '{print $2"\t"}' | tr -d '\n' | awk '{print $0"\n"}')

printf "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" $((n+1)) Na+ Na  $((n+1)) $Na >> child_$GRO
printf "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" $((n+2)) Cl- Cl  $((n+2)) $Cl >> child_$GRO

sed -n '$p' tempdir/pad.gro >> child_$GRO

#rm -rf tempdir {X..Z} {X..Z}C Np L 
###########################################
