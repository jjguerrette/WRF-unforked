#!/bin/bash
SECONDS=0
#========================================================
#WRFVAR_MEMBER.sh is used to run a single ensemble member
#========================================================
iSAMP="$1"; it="$2"; rand_stage="$3"; innerit="$4"

ii=$iSAMP
if [ $iSAMP -lt 10 ]; then ii="0"$ii; fi
if [ $iSAMP -lt 100 ]; then ii="0"$ii; fi
if [ $iSAMP -lt 1000 ]; then ii="0"$ii; fi

cd $MEMBERPREFIX$ii

it0=$it
if [ $it -lt 10 ]; then it0="0"$it0; fi

it0_last=$((it-1))
if [ $it0_last -lt 10 ]; then it0_last="0"$it0_last; fi

innerit0=$innerit
if [ $innerit -lt 10 ]; then innerit0="0"$innerit0; fi
if [ $innerit -lt 100 ]; then innerit0="0"$innerit0; fi
if [ $innerit -lt 1000 ]; then innerit0="0"$innerit0; fi

innerit0_last=$((innerit-1))
if [ $((innerit-1)) -lt 10 ]; then innerit0_last="0"$innerit0_last; fi
if [ $((innerit-1)) -lt 100 ]; then innerit0_last="0"$innerit0_last; fi
if [ $((innerit-1)) -lt 1000 ]; then innerit0_last="0"$innerit0_last; fi

# Archive old log files
TEMPTIME=$SECONDS
if [ $innerit -eq 1 ]; then
   if [ $rand_stage -eq 1 ]; then
      if [ $it -gt 1 ] && \
         ([ $RIOT_RESTART -eq 0 ] || [ $RIOT_RESTART -eq 2 ]); then
         mkdir oldrsl_$it0_last
         mv rsl.* oldrsl_$it0_last
      fi
   else
      mkdir oldrsl_$it0".iter"$innerit0".stage1"
      mv rsl.* oldrsl_$it0".iter"$innerit0".stage1"
   fi
else
   mkdir oldrsl_$it0".iter"$innerit0_last".stage3"
   mv rsl.* oldrsl_$it0".iter"$innerit0_last".stage3"
fi
TEMPTIME=$((SECONDS-TEMPTIME))
hr_temp1=$(($TEMPTIME / 3600))
if [ $hr_temp1 -lt 10 ]; then hr_temp1="0"$hr_temp1; fi
min_temp1=$((($TEMPTIME / 60) % 60))
if [ $min_temp1 -lt 10 ]; then min_temp1="0"$min_temp1; fi
sec_temp1=$(($TEMPTIME % 60))
if [ $sec_temp1 -lt 10 ]; then sec_temp1="0"$sec_temp1; fi


#Handle external checkpoint and obs associated files that change each outer iteration
TEMPTIME=$SECONDS
if [ $innerit -eq 1 ]; then
if [ $rand_stage -eq 1 ] && [ $rand_type -eq 3 ]; then
   if [ $CPDT -gt 0 ]; then
      rm wrf_checkpoint_d01_*
      if [ $WRF_MET -gt 0 ]; then
         rm xtraj_for_obs_d01_*
      fi
      if [ $WRF_CHEM -gt 0 ]; then
         rm SURFACE_Hx_y*
         rm AIRCRAFT_Hx_y*
      fi
   fi
else if ([ $rand_stage -eq 3 ] && [ $rand_type -eq 3 ]) || \
        ([ $rand_stage -eq 1 ] && [ $rand_type -eq 6 ]); then
   if [ $CPDT -gt 0 ]; then
      ln -sf $CWD_rel/wrf_checkpoint_d01_* ./
      if [ $WRF_MET -gt 0 ]; then
        ln -sf $CWD_rel/xtraj_for_obs_d01_* ./
      fi
   fi
   if [ $WRF_CHEM -gt 0 ]; then
      ln -sf $CWD_rel/SURFACE_Hx_y* ./
      ln -sf $CWD_rel/AIRCRAFT_Hx_y* ./
   fi
fi
fi
fi

#Link cvt file for current outer iteration
if [ $it -gt 1 ] && [ $rand_stage -eq 1 ] && ([ $RIOT_RESTART -eq 0 ] || [ $RIOT_RESTART -eq 2 ]); then
   #Test for the presence of cvt
   if [ $(ls "$CWD"/cvt.it"$it0_last".* | wc -l) -eq 0 ]; then echo "ERROR: Missing cvt.it$it0_last"; echo $((rand_stage*100+3)); exit $((rand_stage*100+3)); fi
   ln -sf $CWD_rel/cvt.it$it0_last* ./
#   ln -sf $CWD_rel/cvt.* ./
fi

TEMPTIME=$((SECONDS-TEMPTIME))
hr_temp2=$(($TEMPTIME / 3600))
if [ $hr_temp2 -lt 10 ]; then hr_temp2="0"$hr_temp2; fi
min_temp2=$((($TEMPTIME / 60) % 60))
if [ $min_temp2 -lt 10 ]; then min_temp2="0"$min_temp2; fi
sec_temp2=$(($TEMPTIME % 60))
if [ $sec_temp2 -lt 10 ]; then sec_temp2="0"$sec_temp2; fi


TEMPTIME=$SECONDS
#Large source of load imbalance as written.  Need to rewrite.
if [ $innerit -eq 1 ] && [ $rand_stage -eq 3 ] && [ $rand_type -eq 3 ]; then
   #Check for presence of qhat vectors for this ensemble member (cv space)
   vectors=("qhat.e")
   vectors2=(".iter$innerit0")

   vcount=0
   for var in ${vectors[@]}
   do
      suffix=${vectors2[$vcount]}
#      echo "Checking for $var$ii.*$suffix files"

      #Test for the presence of qhat for this ensemble member
      if [ $(ls ../vectors_$it0/$var$ii.*$suffix | wc -l) -eq 0 ]; then echo "ERROR: Missing $var$ii.*$suffix"; echo $((rand_stage*100+2)); exit $((rand_stage*100+2)); fi

      vcount=$((vcount+1))
   done
fi

TEMPTIME=$((SECONDS-TEMPTIME))
hr_temp3=$(($TEMPTIME / 3600))
if [ $hr_temp3 -lt 10 ]; then hr_temp3="0"$hr_temp3; fi
min_temp3=$((($TEMPTIME / 60) % 60))
if [ $min_temp3 -lt 10 ]; then min_temp3="0"$min_temp3; fi
sec_temp3=$(($TEMPTIME % 60))
if [ $sec_temp3 -lt 10 ]; then sec_temp3="0"$sec_temp3; fi

hr0=$(($SECONDS / 3600))
if [ $hr0 -lt 10 ]; then hr0="0"$hr0; fi
min0=$((($SECONDS / 60) % 60))
if [ $min0 -lt 10 ]; then min0="0"$min0; fi
sec0=$(($SECONDS % 60))
if [ $sec0 -lt 10 ]; then sec0="0"$sec0; fi
echo "WRFVAR_MEMBER.$ii time, TOTAL-$hr0:$min0:$sec0; "$'\n'"ARCHIVE-$hr_temp1:$min_temp1:$sec_temp1; "$'\n'"LINK-$hr_temp2:$min_temp2:$sec_temp2; "$'\n'"CHECK-$hr_temp3:$min_temp3:$sec_temp3; "

exit 0

