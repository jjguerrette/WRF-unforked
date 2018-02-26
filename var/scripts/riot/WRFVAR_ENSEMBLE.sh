#!/bin/bash
SECONDS=0
#========================================================
#WRFVAR_ENSEMBLE.sh is used to run the ensemble for:
#stage 1 (rand_type==[3,6]) of WRFDA_RIOT.sh
#stage 3 (rand_type==[3]) of WRFDA_RIOT.sh
#========================================================

#Probably should add a check for each of this I/O variables to make sure they are not empty
it="$1"; iSAMP0="$2"; NSAMP="$3"; NJOBS="$4"; rand_stage="$5"; innerit="$6"

if [ $NJOBS -gt $((NSAMP+1)) ]; then echo "ERROR: NJOBS should never be more than NSAMP+1"; echo $((rand_stage*100+4)); exit $((rand_stage*100+4)); fi


set -m          # allow for job control
EXIT_CODE=0     # exit code of overall script
function job_exit_codes() {
    for job in `jobs -p`; do
        echo "PID => ${job}"
         CODE=0;
         wait ${job} || CODE=$?
         if [[ "${CODE}" != "0" ]]; then
            echo "At least one WRFVAR_MEMBER failed with exit code => ${CODE}" ;
            EXIT_CODE=1;
         fi
    done
}
trap 'job_exit_codes' CHLD


proc_i=0
for (( iSAMP = $iSAMP0 ; iSAMP <= $NJOBS ; iSAMP++))
#for (( iSAMP = $NJOBS ; iSAMP >= $iSAMP0 ; iSAMP--))
do
   ii=$iSAMP
   if [ $iSAMP -lt 10 ]; then ii="0"$ii; fi
   if [ $iSAMP -lt 100 ]; then ii="0"$ii; fi
   if [ $iSAMP -lt 1000 ]; then ii="0"$ii; fi

   cd ../run.$ii

   echo "Starting job number $iSAMP for $NSAMP samples"

   ## Assign the processes to the hostlist for the current ensemble member
   npiens=$nproc_local
   if [ $iSAMP -eq $((NSAMP+1)) ]; then
      npiens=$nproc_local_grad
   fi

   # Forward consecutive process placement
   ## tail -$((NUMPROC-proc_i)) $PBSNODE0 | head -$npiens > hostlist

   # Reverse consecutive process placement
   head -$((NUMPROC-proc_i)) $PBSNODE0 | tail -$npiens > hostlist

   proc_i=$((proc_i+npiens))

   ./WRFVAR_MEMBER.sh "$iSAMP" "$npiens" "$it" "$rand_stage" "$innerit" &
#   wait_pids+=($!)
   wait_pids2[$iSAMP]=$!
done

## Are there limitations to running large numbers of background processes?
# Could also use GNU "parallel", where npiens must be determined internally from values of $NSAMP and $iSAMP
# (on small scale tests, this does not seem to give an advantage; need to try larger-scale tests)
# parallel -j $NJOBS ./WRFVAR_MEMBER.sh ::: `seq 1 $NJOBS` ::: "$NSAMP" ::: "$it" ::: "$rand_stage" ::: "$innerit"

##WAIT for all ensembles to finish
#echo "wait ${wait_pids[@]}"
#wait "${wait_pids[@]}"

EXIT_CODE=0
for (( iSAMP = $iSAMP0 ; iSAMP <= $NJOBS ; iSAMP++))
do
   wait "${wait_pids2[$iSAMP]}"
   errcode=$?
   if [ $errcode -ne 0 ]; then 
      echo "WRFVAR_MEMBER $iSAMP with PID ${wait_pids2[$iSAMP]} failed with exit code $errcode"
      EXIT_CODE=1
   fi
done

#Reset PBS_NODEFILE
export PBS_NODEFILE=$PBSNODE0

hr0=$(($SECONDS / 3600))
if [ $hr0 -lt 10 ]; then hr0="0"$hr0; fi
min0=$((($SECONDS / 60) % 60))
if [ $min0 -lt 10 ]; then min0="0"$min0; fi
sec0=$(($SECONDS % 60))
if [ $sec0 -lt 10 ]; then sec0="0"$sec0; fi
echo "WRFVAR_ENSEMBLE time: $hr0:$min0:$sec0"

echo "WRFVAR_ENSEMBLE EXIT_CODE => $EXIT_CODE"
exit "$EXIT_CODE"
