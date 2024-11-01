#!/bin/bash

##root -l -b -q 'get_npe.C("ch1_point2_U53V_T15")'

ch=1
U=53
T=15
#T=45
#list=( 1 2 3 4 5 6 7 8 9 10 11 12 )
#list=( 2 6 )
#list=( 2 4 9 )
#list=( 2 7 8 9 )
#list=( 2 9 )
list=( 2 )

for i in ${list[@]}
#for (( i = 1; i <= 16; i++ ))
do 
   echo "U = $U, T = $T, i = $i" 
   qsub -shell n -b y -V -N npe_res_ch${ch}_point${i}_${U}V_T${T} -cwd root -l -b -q 'get_npe.C("ch'${ch}'_point'${i}'_U'${U}'V_T'${T}'")'
done
