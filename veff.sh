#!/bin/bash
#PROC=1
#TRIALS=${2:-10}
WDIR=$(pwd)
GRO=$1
rev=$3
shape=${4:-triclinic}
Si=${2:-1}
export PLUMED_MAXBACKUP=-1
export GMX_MAXBACKUP=-1
#echo "$shape box"



#mkdir -p tempdir
tempdir=$(mktemp -d -t -p $WDIR tempdir-XXXXXXXXXX)
echo $tempdir
cd $tempdir

######################################################################

# Create Na.gro
cat << 'EOF' > rev_gro.sh
#!/bin/bash
GRO=$1
n=$(awk 'NR==2 {print $1}' $GRO)
N_ions=$(echo $n/2 | bc)

sed -n '1,2p' $GRO > rev1$GRO
sed '1,2d; $d' $GRO | tac >> rev1$GRO
sed -n '$p' $GRO >> rev1$GRO
echo system | gmx trjconv -f rev1$GRO -o rev$GRO -s rev1$GRO
rm -rf rev1$GRO
EOF

cat << 'EOF' > sasa.sh 
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
EOF
######################################################################

cp $WDIR/$GRO $tempdir/$GRO
cp $WDIR/vquick.R $tempdir/vquick.R

#rev=$(Rscript -e 'cat(sample(c(0,1),1,replace=T))')
#rev=1
#if [ $rev -eq 1 ]
#then
#echo reversed
#bash rev_gro.sh $GRO &> /dev/null
#cp $GRO orig$GRO
#cp rev$GRO $GRO
#fi

# make index file
echo q | gmx make_ndx -f $GRO &> /dev/null

# pbc Wrap cluster 
echo "$(cat << EOF 
Ion: GROUP NDX_FILE=index.ndx NDX_GROUP=Ion
com: CENTER ATOMS=Ion
WRAPAROUND ATOMS=Ion AROUND=com
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

if [ $rev -eq 1 ]
then
echo reversed
bash rev_gro.sh pad.gro &> /dev/null
cp revpad.gro pad.gro
fi

sed '1,2d; $d' pad.gro | awk '{print $4"\t"}' | tr -d '\n' | awk '{print $0 "\n"}' > X
sed '1,2d; $d' pad.gro | awk '{print $5"\t"}' | tr -d '\n' | awk '{print $0 "\n"}' > Y
sed '1,2d; $d' pad.gro | awk '{print $6"\t"}' | tr -d '\n' | awk '{print $0 "\n"}' > Z
n=$(awk 'NR==2 {print $1}' pad.gro)
Lx=$(tail -n 1 pad.gro | awk '{print $1}')
Ly=$(tail -n 1 pad.gro | awk '{print $2}')
Lz=$(tail -n 1 pad.gro | awk '{print $3}')

echo $n > Np
tail -n 1 pad.gro | awk '{print $0}' > L



############################## SASA NA
#cd $WDIR
Rscript --vanilla vquick.R $Si $((n+1)) $Si &> /dev/null

########################################

# make child gro
n=$(awk 'NR==2 {print $1}' pad.gro)

sed -n '1p' pad.gro > child_$GRO
echo $((n+1)) >> child_$GRO
sed '1,2d; $d' pad.gro >> child_$GRO

if [ $((n%2)) -eq 0 ]
then
Xx=$(sed -n '3p' pad.gro | awk '{print $2}')
Xy=$(sed -n '4p' pad.gro | awk '{print $2}')
Xxc=$(sed -n '3p' pad.gro | awk '{print $1}' | cut -c 2-)
else
Xx=$(sed -n '4p' pad.gro | awk '{print $2}')
Xy=$(sed -n '3p' pad.gro | awk '{print $2}')
Xxc=$(sed -n '4p' pad.gro | awk '{print $1}' | cut -c 2-)
fi

mkdir -p child_$rev
rm -rf child_$rev/*
cp sasa.sh child_$rev/
L=$(wc -l < XC)

for i in `seq 1 $L`
do
cp child_$GRO child_$rev/$i.gro
crd=$(echo {X..Z}C | xargs -n 1 sed -n "$i p" | awk '{print $1"\t"}' | tr -d '\n' | awk '{print $0"\n"}')
printf "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" $((n+1)) $Xxc $Xx  $((n+1)) $crd >> child_$rev/$i.gro
sed -n '$p' pad.gro | awk '{print $1+0.35, $2+0.35, $3+0.35}' >> child_$rev/$i.gro
done

cd child_$rev
for i in `seq 1 $L`
do
bash sasa.sh $i.gro $Xx >> volume.txt
done
cd ..
bash sasa.sh $GRO $Xy | tee volume_$Xy.txt
#cp M* child_$rev
cp T* child_$rev
cp Cl.txt child_$rev
cp volume_$Xy.txt child_$rev/
rm -rf $WDIR/child_$rev
mv child_$rev $WDIR
cd $WDIR
#cp volume_$GRO.txt volume.txt
#rm -rf $tempdir 
###########################################
