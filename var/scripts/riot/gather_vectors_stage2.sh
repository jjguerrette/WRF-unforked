#!/bin/bash

#----------------------------------------------
#Gather relevant vectors into single directory
#----------------------------------------------
cd $CWD

if [ $rand_type -eq 1 ]; then
   vectors=("yhat_obs.e*.p")
   nvec=("$NENS")
fi
if [ $rand_type -eq 2 ]; then
   vectors=("ahat.e*.p" "ghat.p")
   nvec=("$NENS" "1")
fi
if [ $rand_type -eq 6 ]; then
   vectors=("omega.e*.p" "yhat.e*.p" "ghat.p")
   nvec=("$NENS" "$NENS" "1")
fi
if [ $rand_type -eq 10 ]; then
   vectors=("omega.e*.p" "yhat.e*.p")
   nvec=("$NENS" "$NENS" "1")
fi

vcount=0
for var in ${vectors[@]}
do
   echo "Working on $var files"
   #Test for the presence of each vector type
   ls ../run.*/$var"0000"
   dummy=`ls ../run.*/$var"0000" | wc -l`
   echo "$dummy present of ${nvec[$vcount]} $var files"
   if [ $dummy -ne ${nvec[$vcount]} ]; then 
      echo "ERROR: Missing or extra $var""0000"
      echo $((rand_stage*100+vcount*10+4)); exit $((rand_stage*100+vcount*10+4));
   fi
   mv -v ../run.*/$var* ../vectors_$it0
   ln -sf ../vectors_$it0/$var* ./
   vcount=$((vcount+1))
done
