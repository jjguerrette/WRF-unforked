#!/bin/bash


SECONDS=0
iteration_sec=0
cleanup_sec=0
echo ""
echo "********************************************************"
echo ""
echo "         WRFDA-RIOT script with randomised SVD"
echo ""
echo "********************************************************"
echo ""

#=======================================================================================
# INSTRUCTIONS:
#=======================================================================================
# Place this script and RIOT_settings.sh in the same directory where you would normally 
# execute da_wrfvar.exe. Ensure that you are able to successfully run da_wrfvar.exe for 
# a 4D-Var inversion before using RIOT, including generation of cost_fn and 
# wrfvar_output files.
# The same namelist.input file used for standard 4D-Var will be used for RIOT 4D-Var,
#   with additional modifications automated by this script. You must have access to the
#   parent directory within which this working directory is located.  This script will
#   generate directories such as ../run.0001, ../run.0002, etc. that are essential.

export MEMBERPREFIX="../run."

# Fill in the environment variables in RIOT_settings.sh and the User Options section. 
#   - nout_RIOT takes the place of max_ext_its in namelist.input as the # outer iterations
#   - ntmax(1) in namelist.input is replaced by NBLOCK
#   - ntmax(2:nout_RIOT) should be set in namelist.input
#   - Follow the instructions for individual options below
# 
# Disclaimer: This script has been tested on a very small number of platforms. Each
#  platform will have a unique combination of MPI/MPT implementations and job queue 
#  management tools (e.g., PBS).  Specific solutions may be required for your 
#  environment. For questions/comments, please contact:
#
#  JJ Guerrette 
#  jonathan.guerrette@noaa.gov
#  jonathan.guerrette@colorado.edu
#
#=======================================================================================
if [ $(ls da_wrfvar.exe | wc -l) -eq 0 ]; then echo "da_wrfvar.exe required for WRFDA-RIOT"; echo 101; exit 101; fi
if [ $(ls namelist.input | wc -l) -eq 0 ]; then echo "namelist.input required for WRFDA-RIOT"; echo 102; exit 102; fi

if [ -z "$MPICALL" ]; then #Could be defined externally
   ##Check the manual for your compute system in order to choose one of the following:
   #MPICALL=mpiexec #(e.g., NASA Pleiades)
   MPICALL=mpirun #(e.g., NOAA Theia)

   # Note: the script performance has been verified on the platforms in parentheses;
   #       please add more as they are confirmed
fi
echo "MPI Calling Wrapper = "
echo "==> "$MPICALL

#=======================================================================================
# Begin User Options
#=======================================================================================
#####################################################################################
#Modify RIOT_settings.sh to control RIOT behavior

#Read in baseline RIOT settings
if [ -z "$RIOT_SETTINGS_CALLED" ]; then
   echo "Running default settings file: RIOT_settings.sh"
   . ./RIOT_settings.sh
fi

if [ -z "$WRF_CHEM" ]; then #Could be defined externally
   # Set to 0 (default) to not conduct CHEM DA
   # 1 to conduct chemical emission inversion (requires WRFDA-Chem)
   export WRF_CHEM=0
fi
if [ $WRF_CHEM -eq 1 ]; then
   if [ -z $FORCE_SURFACE ]; then #Could be defined externally
      #Set to 0 (default) or 1 to conduct chemical emission inversion (requires WRFDA-Chem)
      export FORCE_SURFACE=1
   fi
   if [ -z $FORCE_AIRCRAFT ]; then #Could be defined externally
      #Set to 0 (default) or 1 to conduct chemical emission inversion (requires WRFDA-Chem)
      export FORCE_AIRCRAFT=1
   fi
fi
if [ -z "$WRF_MET" ]; then #Could be defined externally
   # Set to 0 (default) to not conduct MET DA
   # 1 to conduct MET DA
   export WRF_MET=0
fi
if [ $WRF_CHEM -eq 0 ] && [ $WRF_MET -eq 0 ]; then
   echo "ERROR: Either WRF_CHEM or WRF_MET must be set to 1"
   echo 1; exit 1
fi

if [ -z "$nin_RIOT" ]; then #Could be defined externally
   export nin_RIOT=1
fi

#Manually turn on/off gathering of subprocedure times
export SUBTIMING=1 #0 or 1
#####################################################################################

#Additional local settings
OUTDIR="ALL_RIOT_OUTPUT"

#MPI options
DEBUGSTR=
#DEBUGSTR="-show"
#DEBUGSTR="-verbose"

export MPI_VERBOSE=1
RIOT_EXECUTABLE="./da_wrfvar.exe"

#------------------------------------
CLEANUP=2 #Takes extra time or space:
# 0 - leave temp files as-is
# 1 - store intermediate and end temp files in large tars
# 2 - store intermediate temp files
# >2 - remove temp files
#------------------------------------

#All of these must be set to 1 to perform RIOT - toggle for debugging
STAGE0=1
STAGE1=1 #Set to 1 to perform parallel A * OMEGA calculations
STAGE2=1 #Set to 1 to perform final SVD and increment CVs

#=======================================================================================
# End User Options
#=======================================================================================

echo "=================================================="
echo " Critical RIOT Options"
echo "=================================================="
echo "(1) rand_type=$rand_type"
echo " * valid rand_type options: "
echo "   + 6-RSVD5.6"
echo "   + 3-Block Lanczos"
echo " * rand_type=[3,6] are functional in WRFDA for CHEM and non-CHEM variables"
if [ $rand_type -ne 6 ] && [ $rand_type -ne 2 ] && [ $rand_type -ne 3 ]; then echo "ERROR: unknown rand_type=$rand_type"; echo 2; exit 2; fi
echo ""
echo "(2) RIOT_RESTART=$RIOT_RESTART"
echo ""
echo "(3) ALT_START=$ALT_START"
echo ""
echo "(4) MAX processors per job: NPpJMAX=$NPpJMAX"
echo ""
#Could use these two lines instead:
#NOUTER=`grep max_ext_its namelist.input`
#NOUTER=${NOUTER//[!0-9]/}
NOUTER=$nout_RIOT # overwrites NOUTER with environment variable
NINNER=$nin_RIOT # overwrites NINNER with environment variable
if [ $rand_type -ne 3 ]; then NINNER=1; fi

echo "(5) # of outer iterations: nout_RIOT=$NOUTER"
ex -c :"/max_ext_its" +:":s/=.*/=1,/" +:wq namelist.input

#Retrieve ntmax array from namelist - set NSAMP for it>1
ntmax_all=`grep ntmax namelist.input`
ntmax_all=${ntmax_all#*=}
IFS=',' read -ra ntmax_array <<< "$ntmax_all"

ntmax_array[0]=$NBLOCK
for it in $(seq 1 $NOUTER)
do
   echo "# samples in outer iteration $it: ntmax($it)=${ntmax_array[$((it-1))]}"

   #Hard constraint that NSAMP must be <= NBLOCK (necessary?)
   if [ ${ntmax_array[$((it-1))]} -gt $NBLOCK ]; then
      echo "WARNING ntmax($it)=$NSAMP > NBLOCK=$NBLOCK, setting ntmax($it)=$NBLOCK"
      ${ntmax_array[$((it-1))]}=$NBLOCK
   fi
   ntmax=$ntmax${ntmax_array[$((it-1))]}","
done
ntmax=$ntmax${ntmax_array[$((NOUTER-1))]}","
#ex -c :"/ntmax" +:":s/=.*/=$ntmax/" +:wq namelist.input
ex -c :":%s/ntmax.*/ntmax=$ntmax/" +:wq namelist.input

NSAMP=$NBLOCK
# Set number of parallel jobs
NJOBS=$((NSAMP+1)) #Extra parallel ensemble member for gradient member
if [ $rand_type -eq 3 ]; then  NJOBS=$NSAMP; fi #Not for Block Lanczos

echo ""
echo "(7) LRA-LRU Adaptation: ADAPT_SVD=$ADAPT_SVD"
echo ""
if [ -z "$GLOBAL_OPT" ]; then
   GLOBAL_OPT="false" #default choice
fi
GLOBAL_OPT="false" #default choice

echo "(9) GLOBAL_OPT=$GLOBAL_OPT"
if [ "$GLOBAL_OPT" == "true" ]; then
   echo "For very large size(cv), global cv I/O requires"
   echo "1 - additional wall-time during WRFDA gather step"
   echo "2 - larger minimum memory on head node, for same step"
   echo "Consider these factors if da_wrfvar.exe hangs or quits unexpectedly during cv I/O"
   echo "Turn off global I/O with, e.g., export GLOBAL_OPT=\"false\""
   echo ""
else if [ "$GLOBAL_OPT" == "false" ]; then
   echo "Note: Global I/O is necessary to change number of processors per ensemble member between outer iterations"
else
   echo "ERROR: GLOBAL_OPT must either be true or false"; echo 5; exit 5
fi
fi

echo "================================================="
echo "Characerize the nodes and processors being used"
echo "================================================="
PPN=  #Use this to override the number of processors per node
if [ -z "$PBS_NODEFILE" ]; then
   echo "ERROR: Using SVD requires PBS_NODEFILE to be set to your host file."
   echo "ERROR: This tool requires a large number of processors (> 100)"
   echo 6; exit 6
fi
PBSNODE0=$PBS_NODEFILE
NUMNODES=`cat $PBSNODE0 | uniq | wc -l`
if [ $NUMNODES -lt $NJOBS ]; then
   echo "ERROR: NUMNODES must be set >= NJOBS (i.e., $NJOBS)"
   echo "For RSVD (rand_type=6), NJOBS=NBLOCK+1"
   echo "For Block Lanczos (rand_type=3), NJOBS=NBLOCK"
   echo 7; exit 7
fi

if [ -n "$PPN" ]; then
   echo "Using user defined cores per node: PPN=$PPN"
   NUMPROC=$((NUMNODES*PPN))
else
   NUMPROC=`cat $PBSNODE0 | wc -l`
   PPN=$((NUMPROC/NUMNODES))
   echo "Automatically determining cores per node: PPN=$PPN"
fi

#Manually set the maximum number of processes per job 
# - limited by WRF patch overlap
# - depends on domain (nx x ny) and PPN
# - critical for speeding up single-job (serial) portions of RIOT 
#    (e.g., STAGE 2 of RSVD5.6 and Block Lanczos, and STAGE 4 of Block Lanczos)
# - A good rule of thumb is (nx/10 * ny/10) <= NPpJMAX << (nx/5 * ny/5)
# - Requirement: NPpJMAX <= NUMPROC
if [ -z "$NPpJMAX" ]; then
#   export NPpJMAX=64

   nx=`grep e_we namelist.input`
   nx=${nx#*=}
   nx=${nx%,*}

   ny=`grep e_sn namelist.input`
   ny=${ny#*=}
   ny=${ny%,*}

   NPFAC=8
#   NPFAC=10 #Will be required in WRF V4.0
   export NPpJMAX=$(($((nx*ny))/$((NPFAC*NPFAC))))

   #Use integer number of nodes (optional)
   NODEMAX=$((NPpJMAX/PPN))
   export NPpJMAX=$((NODEMAX*PPN))
fi

if [ $NPpJMAX -gt $NUMPROC ]; then
   NPpJMAX=$NUMPROC
   export NPpJMAX=$NPpJMAX

   #Use integer number of nodes (optional)
   NODEMAX=$((NPpJMAX/PPN))
   export NPpJMAX=$((NODEMAX*PPN))
fi
#Use as many of the cores as possible for ensemble members
#NODES_all_ens=$((NUMNODES/NSAMP))
#NPpJ=$((NODES_all_ens*PPN))

NPpJ=$(($((NUMNODES/NSAMP))*PPN)) 
if [ $NPpJ -gt $NPpJMAX ]; then NPpJ=$(($((NPpJMAX/PPN))*PPN)); fi
NODES_all_ens=$((NPpJ/PPN*NSAMP))

#Use remaining cores for gradient member (if needed)
NODES_grad=$((NUMNODES-NODES_all_ens))
NPpJ_grad=$((NODES_grad*PPN))
if [ $NJOBS -gt $NSAMP ]; then
   #Fail-safe for when mod(NUMNODES,NSAMP)==0
   if [ $NODES_grad -lt 1 ]; then
      #NODES_all_ens=$(($((NUMNODES-1))/NSAMP))
      #NPpJ=$((NODES_all_ens*PPN))
      NPpJ=$(($(($((NUMNODES-1))/NSAMP))*PPN)) 
      NODES_all_ens=$((NPpJ/PPN*NSAMP))
      NODES_grad=$((NUMNODES-NODES_all_ens))
   fi
   NPpJ_grad=$((NODES_grad*PPN))
   if [ $NPpJ_grad -gt $NPpJMAX ]; then NPpJ_grad=$(($((NPpJMAX/PPN))*PPN)); fi
   NODES_grad=$((NPpJ_grad/PPN))

   if [ "$GLOBAL_OPT" == "false" ] && [ $NPpJ_grad -ne $NPpJ ]; then
       echo "ERROR: nproc_local and nproc_local_grad must be equal when GLOBAL_OPT==false"
       echo "NPpJ=$NPpJ"
       echo "NPpJ_grad=$NPpJ_grad"
       echo 9; exit 9
   fi
fi

echo "Cores per ensemble member job: "$NPpJ
echo "Cores per non-global gradient job: "$NPpJ_grad
if [ $((NPpJ % PPN)) -ne 0 ] || [ $((NPpJ_grad % PPN)) -ne 0 ] || [ $NPpJ -eq 0 ] || ([ $NPpJ_grad -eq 0 ] && [ $NJOBS -gt $NSAMP ]); then
    echo "ERROR: cores per ensemble member must be equal to PPN times a positive integer multiplier"
    echo "NPpJ=$NPpJ"
    echo "NPpJ_grad=$NPpJ_grad"
    echo "PPN=$PPN"
    echo 10; exit 10
fi

# Set processor counts for: 
#(1) stage 0
nproc_max=$NUMPROC
if [ $nproc_max -gt $NPpJMAX ]; then nproc_max=$NPpJMAX; fi

#(2) ensemble members
nproc_local=$NPpJ

#(3) gradient member
nproc_local_grad=$NPpJ_grad

#(4) global gathering stages
if [ "$GLOBAL_OPT" == "false" ]; then
   #FOR NON-GLOBAL SIM
   nproc_global=$nproc_local
else if [ "$GLOBAL_OPT" == "true" ]; then
   #FOR GLOBAL ENABLED SIM
   nproc_global=$nproc_max
fi
fi

echo "LOCAL PROCS = $nproc_local"
echo "GRAD PROCS = $nproc_local_grad"
echo "GLOBAL PROCS = $nproc_global"

echo "=================================================================="
echo " Automatically initialize 4D-Var namelist.input settings for RIOT"
echo "=================================================================="
#------------------------------------------------------------------
# Set checkpoint interval and OSSE from environment variables
# - only meaningful for CHEM so far
# - set to 0 otherwise
# - if non-zero, namelist values should already be set
#------------------------------------------------------------------
CPDT=$checkpoint_interval 
#quick_rand="true"
if [ -z "$CPDT" ]; then
   if [ $WRF_CHEM -eq 0 ]; then
      CPDT=0 #Or set to some other default for WRF_CHEM==0
   else
      #Thus far WRF_CHEM==1 requires CPDT>0
      CPDT=180 # (or manually, e.g., 5, 6, 10, 12, 15, 20, 30, 60, 180)
   fi
fi
#if [ $CPDT -le 0 ]; then quick_rand="false"; fi
echo "checkpoint interval = "
echo "==> "$CPDT

#CHEM OSSE

OSSE=$OSSE_CHEM # Set OSSE (from env var)
if [ -z "$OSSE" ] || [ $WRF_CHEM -eq 0 ]; then
   OSSE=0 #Default if env var missing
fi
echo "OSSE = "
echo "==> "$OSSE

# Turn off lanczos
if (grep -q use_lanczos namelist.input); then
   ex -c :"/use_lanczos" +:":s/use_lanczos.*/use_lanczos=false,/" +:wq namelist.input
fi

# Add or modify RIOT namelist options
SVD_VARS=("rand_outer" "rand_minimise" "ensmember" "rand_stage" "use_randomblock" "rand_type" "adapt_svd" "svd_p" "read_omega" "use_global_cv_io" "rand_inner_it" "max_rand_inner")
SVD_VALS=("1" "true" "0" "0" "true" "$rand_type" "$ADAPT_SVD" "$svd_p" "$GLOBAL_OMEGA" "$GLOBAL_OPT" "1" "$NINNER")
ivar=0
for var in ${SVD_VARS[@]}
do
   if (grep -q $var namelist.input); then # check if NL option is present already
      ex -c :"/$var" +:":s/=.*/=${SVD_VALS[$ivar]},/" +:wq namelist.input
   else
      ex -c :"/wrfvar6" +:":s/\(wrfvar6\)\n/\1\r$var=${SVD_VALS[$ivar]},\r/" +:wq namelist.input
   fi
   ivar=$((ivar+1))
done

echo "=================================================="
echo " Setup all working directories"
echo "=================================================="
CWD=$(pwd)
echo "WORKING DIRECTORY: "$CWD
cd ../
PARENTDIR=$(pwd)
CWD_rel0=${CWD#$PARENTDIR"/"}
CWD_rel="../"$CWD_rel0 #Use for linking, more stable
cd $CWD

#Cleanup Parent
dummy=$(ls ../cost_fn*)
if [ $? -eq 0 ]; then  rm -r ../cost_fn*; fi
dummy=$(ls ../grad_fn*)
if [ $? -eq 0 ]; then  rm -r ../grad_fn*; fi

#Cleanup CWD
if [ $(ls -d RIOT_stash | wc -l) -gt 0 ]; then rm -rv RIOT_stash; fi
if [ $(ls -d *rsl* | wc -l) -gt 0 ] || [ $(ls -d *_fn | wc -l) -gt 0 ]; then 
   mkdir RIOT_stash
   mv -v rsl.* RIOT_stash
   mv -v oldrsl* RIOT_stash
   mv -v cost_fn RIOT_stash
   mv -v grad_fn RIOT_stash
fi

if [ $RIOT_RESTART ] && [ $RIOT_RESTART -eq 1 ]; then
#----------------------------------------------------------------
# Option to perform posterior error calculation w/o minimization
#----------------------------------------------------------------
# Either setup links here or externally
   ln -sfv $ALT_START ./fg

   itstart=2
   NOUTER=$((NOUTER+1))
   ex -c :"/rand_minimise" +:":s/=.*/=false,/" +:wq namelist.input
else if [ $RIOT_RESTART ] && [ $RIOT_RESTART -eq 2 ]; then
#------------------------------------------------------------------
# Option to restart RIOT minimization from previous outer iteration
#------------------------------------------------------------------
# Either setup links here or externally
   ln -sfv $ALT_START ./fg

   itstart=$ALT_it1
   if [ -z "$itstart" ]; then
      itstart=2
   fi
   ii=$((itstart-1))
   if [ $ii -lt 10 ]; then ii="0"$ii; fi

   cvt_current="./cvt.it"$ii".p0000"
   xhat_current="./xhat.it"$ii".p0000"
   wrfvar_current="./wrfvar_output_$ii"

   if [ "$cvt_current" -ne "$ALT_CVT" ]; then ln -sfv $ALT_CVT ./cvt.it"$ii".p0000; fi
   if [ "$xhat_current" -ne "$ALT_XHAT" ]; then ln -sfv $ALT_XHAT ./xhat.it"$ii".p0000; fi
   if [ "$wrfvar_current" -ne "$ALT_START" ]; then cp -v $ALT_START ./wrfvar_output_$ii; fi

#   if [ $GRAD_PRECON -gt 0 ] || [ $SPECTRAL_PRECON -gt 0 ]; then
#      cd ../
#      ln -sfv $ALT_hess_dir/hessian_eigenpairs.it* ./
##      ln -sfv $ALT_hess ./hessian_eigenpairs.it"$ii".0000
#      cd $CWD
#   fi
   else
      itstart=1 #default

      if [ $(ls ./fg_orig | wc -l) -eq 0 ]; then cp -v fg fg_orig; fi
      if [ $(ls ./fg_orig | wc -l) -eq 0 ]; then echo 17; exit 17;fi
      ln -sfv ./fg_orig ./fg
   fi
fi
ln -sfv ./fg ./wrfinput_d01


#Currently only available for WRF_CHEM==1
if [ $OSSE -gt 0 ]; then
   echo "Generating/linking OSSE files..."
   if [ $OSSE -eq 1 ]; then
      ex -c :"/osse_chem" +:":s/=.*,/=false,/g" +:wq namelist.input
      ex -c :"/init_osse_chem" +:":s/=.*,/=true,/g" +:wq namelist.input
      mpistring="$MPICALL $DEBUGSTR -np $nproc_max $RIOT_EXECUTABLE"
      #COULD MAKE THIS FASTER (MORE MEMORY) by distributing across more nodes (currently chooses first $npiens processors in $PBS_NODEFILE)
      
      echo "$mpistring"
      eval "$mpistring"

      mpireturn=$?
      echo "WRFDA return value: $mpireturn"

      mkdir rsl_init_osse
      cp -v rsl.* rsl_init_osse
   fi    
   if [ $OSSE -eq 2 ]; then
      if [ -z "$DIRpert" ]; then
         DIRpert="BCLARGE_osse_init"
         DIRpert="BCLARGE_osse_init_xb=0"
      fi    
      fg_osse="$CWD/../../"$DIRpert"/$CWD_rel0/fg_osse"
      osse_obs="$CWD/../../"$DIRpert"/$CWD_rel0/AIRCRAFT_DATA_*.txt"
      osse_var="$CWD/../../"$DIRpert"/$CWD_rel0/AIRCRAFT_MODEL_VARIANCE_*.txt"

      echo "fg_osse="$fg_osse
      echo "osse_obs="$osse_obs
      echo "osse_var="$osse_var

      cp -v $fg_osse ./fg_osse
      cp -v $osse_obs ./
      cp -v $osse_var ./
   fi    
   if [ $RIOT_RESTART -ne 2 ] || [ $itstart -eq 1 ]; then
      mv fg fg_0

      cp -v ./fg_osse ./fg
   fi

   ex -c :"/osse_chem" +:":s/=.*,/=true,/g" +:wq namelist.input
   ex -c :"/init_osse_chem" +:":s/=.*,/=false,/g" +:wq namelist.input
fi 

#Populate ensemble member directories
if [ $(ls $MEMBERPREFIX* | wc -l) -gt 0 ]; then  rm -r $MEMBERPREFIX[0-9]*; fi
if [ $WRF_CHEM -gt 0 ]; then rm *Hx*; fi

for (( iSAMP = 1 ; iSAMP <= $NJOBS ; iSAMP++))
do
   ii=$iSAMP
   if [ $iSAMP -lt 10 ]; then ii="0"$ii; fi
   if [ $iSAMP -lt 100 ]; then ii="0"$ii; fi
   if [ $iSAMP -lt 1000 ]; then ii="0"$ii; fi

   mkdir $MEMBERPREFIX$ii
   cd $MEMBERPREFIX$ii

   # link all files necessary to run da_wrfvar.exe
   ln -s $CWD_rel/* ./
   rm namelist.input
   rm rsl.out.*
   rm rsl.error.*
done
cd $CWD

hr0=$(($SECONDS / 3600))
if [ $hr0 -lt 10 ]; then hr0="0"$hr0; fi
min0=$((($SECONDS / 60) % 60))
if [ $min0 -lt 10 ]; then min0="0"$min0; fi
sec0=$(($SECONDS % 60))
if [ $sec0 -lt 10 ]; then sec0="0"$sec0; fi
echo "RIOT initialization time: $hr0:$min0:$sec0"

#-----------------------------------------
#Export variables used in WRFVAR_ENSEMBLE
#-----------------------------------------
export CWD
export CWD_rel
export rand_type
export GLOBAL_OMEGA
export GLOBAL_OPT
export CPDT
export WRF_CHEM
export WRF_MET
export RIOT_RESTART
export nproc_local
export nproc_local_grad
export nproc_global
export NUMPROC
export PBSNODE0
export MPICALL
export DEBUGSTR
export RIOT_EXECUTABLE
export CLEANUP

## BEGIN Outer Loop
for (( it = $itstart ; it <= $NOUTER ; it++))
do
   SECONDS=0
   echo ""
   echo ""
   echo "===================================================="
   echo "===================================================="
   echo "==                                                =="
   echo "==  Starting outer loop iteration $it of $NOUTER  =="
   echo "==                                                =="
   echo "===================================================="
   echo "===================================================="
   echo ""
   echo ""
   ex -c :"/rand_outer" +:":s/=.*/=$it,/" +:wq namelist.input
   ex -c :"/rand_inner_it" +:":s/=.*/=1,/" +:wq namelist.input
   it0=$it
   if [ $it -lt 10 ]; then it0="0"$it0; fi

   cd $CWD

   if [ $(ls -d "../vectors_"$it0* | wc -l) -gt 0 ]; then rm -rv ../vectors_$it0*; fi
   mkdir -v ../vectors_$it0
   rm omega.e*.p*
   rm yhat.e*.p*

   NSAMP=${ntmax_array[$((it-1))]}
#NOTE: The following constraint has been moved to start of script.
#   #Hard constraint that NSAMP must be <= NBLOCK (necessary?)
#   if [ $NSAMP -gt $NBLOCK ]; then
#      echo "WARNING NBLOCK=$NBLOCK, ntmax(it)=$NSAMP, setting NSAMP=$NBLOCK"
#      NSAMP=$NBLOCK
#      ntmax_array[$((it-1))]=$NSAMP
#      ntmax=
#      for iit in $(seq 1 $NOUTER)
#      do
#         echo "# samples in outer iteration $iit: ntmax($iit)=${ntmax_array[$((iit-1))]}"
#         ntmax=$ntmax${ntmax_array[$((iit-1))]}","
#      done
#      ntmax=$ntmax${ntmax_array[$((NOUTER-1))]}","
##      ex -c :"/ntmax" +:":s/=.*/=$ntmax/" +:wq namelist.input
#      ex -c :":%s/ntmax.*/ntmax=$ntmax/" +:wq namelist.input
#   fi
   echo ""
   echo "Using ensemble size $NSAMP in iteration $it, according to ntmax in namelist.input."
   echo ""
   NJOBS=$((NSAMP+1))
   if [ $rand_type -eq 3 ]; then  NJOBS=$NSAMP; fi # no. of ensemble members = NSAMP for Block Lanczos
   ii_grad=$NJOBS
   if [ $NJOBS -lt 10 ]; then ii_grad="0"$ii_grad; fi
   if [ $NJOBS -lt 100 ]; then ii_grad="0"$ii_grad; fi
   if [ $NJOBS -lt 1000 ]; then ii_grad="0"$ii_grad; fi

   RSLTIME=0
   FILETIME=0
   if [ $STAGE0 -gt 0 ]; then
#====================================================================================
#====================================================================================
      echo "============="
      echo "SVD STAGE 0"
      echo "============="
      # This stage only necessary to generate checkpoint files or global OMEGA vectors
      ex -c :"/rand_stage" +:":s/=.*/=0,/" +:wq namelist.input
      if [ $CPDT -gt 0 ] || [ "$GLOBAL_OMEGA" == "true"  ] ; then
         echo "Generating checkpoint files and/or global omega..."

#         if [ "$GLOBAL_OMEGA" == "true" ] && [ "$GLOBAL_OPT" == "false" ]; then
         if [ "$GLOBAL_OPT" == "false" ]; then #(For cvt reading and GLOBAL_OMEGA writing)
            npiens=$nproc_local
         else
            npiens=$nproc_global
            #Could make this faster (more memory) by distributing across more nodes 
            # - currently chooses first $npiens processors in $PBS_NODEFILE
         fi

         mpistring="$MPICALL $DEBUGSTR -np $npiens $RIOT_EXECUTABLE"
#         mpistring="$MPICALL $DEBUGSTR -np $nproc_max $RIOT_EXECUTABLE"
         echo "$mpistring"
         eval "$mpistring"

         mpireturn=$?
         echo "WRFDA return value: $mpireturn"

         if [ $SUBTIMING -eq 1 ]; then grep da_end_timing rsl.out.0000 > ../bench_time_stage0.it$it0; fi
         TEMPTIME=$SECONDS
         mkdir oldrsl_$it0".stage0"
         mv rsl.* oldrsl_$it0".stage0"
         RSLTIME=$((RSLTIME+$((SECONDS-TEMPTIME))))
      fi
   fi

   if [ $STAGE1 -le 0 ]; then
      echo "EXIT: STAGE1 > 0 required for multiple outer iterations"; echo 11; exit 11;
   fi

   if [ $rand_type -ne 3 ]; then
      #-------------------------------------------------------------------
      # Check for presence of checkpoint files
      #-------------------------------------------------------------------
      if [ $CPDT -gt 0 ]; then
         if [ $(ls "$CWD"/wrf_checkpoint_d01* | wc -l) -eq 0 ]; then echo "ERROR: Missing checkpoint files"; echo 13; exit 13; fi
         if [ $WRF_MET -gt 0 ]; then
            if [ $(ls "$CWD"/xtraj_for_obs_d01* | wc -l) -eq 0 ]; then echo "ERROR: Missing xtraj_for_obs files"; echo 14; exit 14; fi
         fi
      fi

      #-------------------------------------------------------------------
      # Check for presence of obs output files (CHEM only)
      #-------------------------------------------------------------------
      if [ $WRF_CHEM -gt 0 ]; then
         if [ $FORCE_SURFACE -eq 1 ] && [ $(ls "$CWD"/SURFACE_Hx_y* | wc -l) -eq 0 ]; then echo "ERROR: Missing SURFACE_Hx_y files"; echo 15; exit 15; fi
         if [ $FORCE_AIRCRAFT -eq 1 ] && [ $(ls "$CWD"/AIRCRAFT_Hx_y* | wc -l) -eq 0 ]; then echo "ERROR: Missing AIRCRAFT_Hx_y files"; echo 16; exit 16; fi
      fi
   fi


#====================================================================================
#====================================================================================
   echo "========================"
   echo "Prepare for SVD STAGE 1"
   echo "========================"
   cd $CWD
   ex -c :"/rand_stage" +:":s/=.*/=1,/" +:wq namelist.input

   if [ $it -gt $itstart ] || ([ $it -gt 1 ] && [ $RIOT_RESTART -eq 2 ]); then
      dependent_grad=0
      if [ $rand_type -eq 3 ]; then dependent_grad=1; fi

      #----------------------------------------------------
      # Redistribute cores among samples as NSAMP changes
      #----------------------------------------------------
      # When max(ntmax_array(iteration>=it)) < max(ntmax_array(iteration<it)), 
      # or ntmax_array(it) ~= ntmax_array(it-1), there are several
      # options for improved efficiency or performance:

      ##--------------------------------------------------------------------------
      ## 1 - Release nodes that would go unused (narrow the job width on the fly)
      ##  --> Possible on NAS Pleiades with pbs_release_nodes...still need to script
      ##  --> Available on any other systems?
      ##--------------------------------------------------------------------------

      ##--------------------------------------------------------------------------
      ##2 - Adjust cores per job to fit NUMPROC as NSAMP changes
      #      (possible only with global CV I/O in WRFDA)
      if [ "$GLOBAL_OPT" == "true" ]; then # && [ $NSAMP -ne $NSAMP_last ]; then
         if [ $dependent_grad -eq 1 ]; then
            NODES_all_ens=$NUMNODES
         else
            NODES_all_ens=$((NUMNODES-NODES_grad))
         fi
         NPpJ=$(($((NODES_all_ens/NSAMP))*PPN)) #Use remaining cores for gradient member
         if [ $((NPpJ/PPN)) -eq 0 ]; then
            echo "ERROR: cores per ensemble member less than cores per node"
            echo "$NUMPROC, $NODES_all_ens, $PPN, $NPpJ"
            echo 18; exit 18
         fi
         if [ $NPpJ -gt $NPpJMAX ]; then NPpJ=$(($((NPpJMAX/PPN))*PPN)); fi

         echo "Cores per ensemble member job: "$NPpJ
         export nproc_local=$NPpJ

         dummy=$(($((nproc_local*NSAMP))+$(($((NJOBS-NSAMP))*nproc_local_grad))))
         if [ $dummy -gt $NUMPROC ]; then
            echo "ERROR: Too many processors requested in it=$it, NUMPROC_request=$dummy, NUMPROC_avail=$NUMPROC"; echo 19; exit 19
         fi
      fi
      ##--------------------------------------------------------------------------
   fi

#====================================================================================
#====================================================================================
   if [ $rand_type -eq 6 ] || [ $rand_type -eq 2 ]; then
      echo "============================================================================="
      echo "SVD STAGE 1: Multiply Hessian by $NSAMP OMEGA samples and calculate GRAD(J)"
      echo "============================================================================="
   fi
   if [ $rand_type -eq 3 ]; then
      echo "======================================================="
      echo "SVD STAGE 1: Generate $NSAMP Gradient Realizations"
      echo "======================================================="
   fi

   ./WRFVAR_ENSEMBLE.sh "$it" "1" "$NSAMP" "$NJOBS" "1" "1"
   err=$?; if [ $err -ne 0 ]; then echo $err; exit $err; fi

   if [ $STAGE2 -le 0 ]; then
      echo "EXIT: STAGE2 > 0 required for multiple outer iterations"; echo 20; exit 20;
   fi

   #Store cost_fn and grad_fn values
   cd $CWD
   if [ $it -eq $itstart ]; then
      if [ $rand_type -eq 3 ]; then
         for (( iSAMP = 1 ; iSAMP <= $NSAMP ; iSAMP++))
         do
            ii=$iSAMP
            if [ $iSAMP -lt 10 ]; then ii="0"$ii; fi
            if [ $iSAMP -lt 100 ]; then ii="0"$ii; fi
            if [ $iSAMP -lt 1000 ]; then ii="0"$ii; fi
            cat $MEMBERPREFIX$ii/cost_fn > ../cost_fn.$ii
            cat $MEMBERPREFIX$ii/grad_fn > ../grad_fn.$ii
         done
      else
         cat $MEMBERPREFIX$ii_grad/cost_fn > ../cost_fn
         cat $MEMBERPREFIX$ii_grad/grad_fn > ../grad_fn
      fi
   else
      if [ $rand_type -eq 3 ]; then
         for (( iSAMP = 1 ; iSAMP <= $NSAMP ; iSAMP++))
         do
            ii=$iSAMP
            if [ $iSAMP -lt 10 ]; then ii="0"$ii; fi
            if [ $iSAMP -lt 100 ]; then ii="0"$ii; fi
            if [ $iSAMP -lt 1000 ]; then ii="0"$ii; fi
            grep [0-9] $MEMBERPREFIX$ii/cost_fn >> ../cost_fn.$ii
            grep [0-9] $MEMBERPREFIX$ii/grad_fn >> ../grad_fn.$ii
         done
      else
         # Use these lines if cost_fn file is overwritten each outer iteration
         grep [0-9] $MEMBERPREFIX$ii_grad/cost_fn >> ../cost_fn
         grep [0-9] $MEMBERPREFIX$ii_grad/grad_fn >> ../grad_fn
      fi
   fi


#===================================================================================
#===================================================================================

   if [ $rand_type -eq 6 ] || [ $rand_type -eq 2 ]; then
      echo "================================================="
      echo "SVD STAGE 2: Perform Eigen Decomp + Increment CVs"
      echo "================================================="
   fi
   if [ $rand_type -eq 3 ]; then
      echo "=============================================================="
      echo "SVD STAGE 2: Orthogonalize Gradient Realizations, generate Q"
      echo "=============================================================="

      #--------------------------------------------------------------------
      # Establish identical checkpointed trajectories for stage > 1
      #--------------------------------------------------------------------
      if [ $CPDT -gt 0 ]; then
         if [ $(ls "$MEMBERPREFIX"0001/wrf_checkpoint_d01* | wc -l) -eq 0 ]; then echo "ERROR: Missing checkpoint files"; echo 13; exit 13; fi
         rm wrf_checkpoint_d01*
         mv "$MEMBERPREFIX"0001/wrf_checkpoint_d01* ./
         rm $MEMBERPREFIX*/wrf_checkpoint_d01*

         if [ $WRF_MET -gt 0 ]; then
            if [ $(ls "$MEMBERPREFIX"0001/xtraj_for_obs_d01* | wc -l) -eq 0 ]; then echo "ERROR: Missing xtraj_for_obs files"; echo 14; exit 14; fi
            rm xtraj_for_obs_d01*
            mv "$MEMBERPREFIX"0001/xtraj_for_obs_d01* ./
            rm $MEMBERPREFIX*/xtraj_for_obs_d01*
         fi
      fi

      #--------------------------------------------------------------------
      # Establish identical checkpointed CHEM obs for stage > 1
      #--------------------------------------------------------------------
      if [ $WRF_CHEM -gt 0 ]; then
         if [ $FORCE_SURFACE -eq 1 ] && [ $(ls "$MEMBERPREFIX"0001/SURFACE_Hx_y* | wc -l) -eq 0 ]; then echo "ERROR: Missing SURFACE_Hx_y files"; echo 15; exit 15; fi
         rm SURFACE_Hx_y*
         mv "$MEMBERPREFIX"0001/SURFACE_Hx_y* ./
         rm $MEMBERPREFIX*/SURFACE_Hx_y*
         if [ $FORCE_AIRCRAFT -eq 1 ] && [ $(ls "$MEMBERPREFIX"0001/AIRCRAFT_Hx_y* | wc -l) -eq 0 ]; then echo "ERROR: Missing AIRCRAFT_Hx_y files"; echo 16; exit 16; fi
         rm AIRCRAFT_Hx_y*
         mv "$MEMBERPREFIX"0001/AIRCRAFT_Hx_y* ./
         rm $MEMBERPREFIX*/AIRCRAFT_Hx_y*
      fi
   fi

   ./check_vectors.sh "$it" "$NSAMP" "2" "1"
   err=$?; if [ $err -ne 0 ]; then echo $err; exit $err; fi

   cd $CWD
   ex -c :"/rand_stage" +:":s/=.*/=2,/" +:wq namelist.input

   if [ "$GLOBAL_OPT" == "false" ]; then #(For cvt reading/writing, ghat+yhat+omega reading, and qhat writing)
      npiens=$nproc_local
   else
      npiens=$nproc_global
      #Could make this faster (more memory) by distributing across more nodes 
      # - currently chooses first $npiens processors in $PBS_NODEFILE
   fi

   mpistring="$MPICALL $DEBUGSTR -np $npiens $RIOT_EXECUTABLE"
   #Could make this faster (more memory) by distributing across more nodes 
   # - currently chooses first $npiens processors in $PBS_NODEFILE

   echo "$mpistring"
   eval "$mpistring"

   mpireturn=$?
   echo "WRFDA return value: $mpireturn"

   if ([ $rand_type -eq 6 ] || [ $rand_type -eq 3 ]) && [ $SUBTIMING -eq 1 ]; then
      grep da_end_timing $MEMBERPREFIX*/rsl.out.0000 > ../bench_time_hess-vec.it$it0
      grep da_end_timing rsl.out.0000 > ../bench_time_finalize-riot.it$it0
   fi

#===================================================================================
#===================================================================================
   if [ $rand_type -eq 3 ]; then
      cd $CWD
      TEMPTIME=$SECONDS
      mkdir oldrsl_$it0".iter0001.stage2"
      mv rsl.* oldrsl_$it0".iter0001.stage2"
      RSLTIME=$((RSLTIME+$((SECONDS-TEMPTIME))))

      hr0=$(($SECONDS / 3600))
      if [ $hr0 -lt 10 ]; then hr0="0"$hr0; fi
      min0=$((($SECONDS / 60) % 60))
      if [ $min0 -lt 10 ]; then min0="0"$min0; fi
      sec0=$(($SECONDS % 60))
      if [ $sec0 -lt 10 ]; then sec0="0"$sec0; fi
      echo "Iteration $it gradient realization time: $hr0:$min0:$sec0"

      innerSECONDS=$SECONDS
      innerSECONDS0=$SECONDS
      for (( innerit = 1 ; innerit <= $NINNER ; innerit++))
      do
         echo ""
         echo "========================================================="
         echo "==                                                     =="
         echo "==  Starting inner loop iteration $innerit of $NINNER  =="
         echo "==                                                     =="
         echo "========================================================="
         echo ""
         innerit0=$innerit
         if [ $innerit -lt 10 ]; then innerit0="0"$innerit0; fi
         if [ $innerit -lt 100 ]; then innerit0="0"$innerit0; fi
         if [ $innerit -lt 1000 ]; then innerit0="0"$innerit0; fi
         ex -c :"/rand_inner_it" +:":s/=.*/=$innerit,/" +:wq namelist.input

         echo "====================================================="
         echo "SVD STAGE 3: Multiply Hessian by Q, inner = $innerit"
         echo "====================================================="
         ex -c :"/rand_stage" +:":s/=.*/=3,/" +:wq namelist.input

         ./WRFVAR_ENSEMBLE.sh "$it" "1" "$NSAMP" "$NJOBS" "3" "$innerit"
         err=$?; if [ $err -ne 0 ]; then echo $err; exit $err; fi         

         echo "==========================================================="
         echo "SVD STAGE 4: Finish Block Lanczos inner iteration $innerit"
         echo "==========================================================="

         ./check_vectors.sh "$it" "$NSAMP" "4" "$innerit"
         err=$?; if [ $err -ne 0 ]; then echo $err; exit $err; fi

         cd $CWD
         ex -c :"/rand_stage" +:":s/=.*/=4,/" +:wq namelist.input

         if [ "$GLOBAL_OPT" == "false" ]; then #(For cvt reading, ghat+yhat_i+qhat_i reading, and qhat_i+1 writing)
            npiens=$nproc_local
         else
            npiens=$nproc_global
            #Could make this faster (more memory) by distributing across more nodes 
            # - currently chooses first $npiens processors in $PBS_NODEFILE
         fi

         mpistring="$MPICALL $DEBUGSTR -np $npiens $RIOT_EXECUTABLE"
         #Could make this faster (more memory) by distributing across more nodes 
         # - currently chooses first $npiens processors in $PBS_NODEFILE

         echo "$mpistring"
         eval "$mpistring"

         mpireturn=$?
         echo "WRFDA return value: $mpireturn"

         if [ $SUBTIMING -eq 1 ]; then
            grep da_end_timing $MEMBERPREFIX*/rsl.out.0000 > ../bench_time_hess-vec.it$it0"."$innerit0
            grep da_end_timing rsl.out.0000 > ../bench_time_finalize-riot.it$it0"."$innerit0
         fi

         if [ $innerit -lt $NINNER ]; then
            TEMPTIME=$SECONDS
            mkdir oldrsl_$it0".iter"$innerit0".stage4"
            mv rsl.* oldrsl_$it0".iter"$innerit0".stage4"
            RSLTIME=$((RSLTIME+$((SECONDS-TEMPTIME))))
         fi

         hr0=$(($((SECONDS-innerSECONDS)) / 3600))
         if [ $hr0 -lt 10 ]; then hr0="0"$hr0; fi
         min0=$((($((SECONDS-innerSECONDS)) / 60) % 60))
         if [ $min0 -lt 10 ]; then min0="0"$min0; fi
         sec0=$(($((SECONDS-innerSECONDS)) % 60))
         if [ $sec0 -lt 10 ]; then sec0="0"$sec0; fi
         echo "Inner iteration $innerit compute time: $hr0:$min0:$sec0"

         innerSECONDS=$SECONDS

         hr0=$(($((SECONDS-innerSECONDS0)) / 3600))
         if [ $hr0 -lt 10 ]; then hr0="0"$hr0; fi
         min0=$((($((SECONDS-innerSECONDS0)) / 60) % 60))
         if [ $min0 -lt 10 ]; then min0="0"$min0; fi
         sec0=$(($((SECONDS-innerSECONDS0)) % 60))
         if [ $sec0 -lt 10 ]; then sec0="0"$sec0; fi
         echo "Total accum. Inner iteration compute time: $hr0:$min0:$sec0"
      done

      hr0=$(($innerSECONDS0 / 3600))
      if [ $hr0 -lt 10 ]; then hr0="0"$hr0; fi
      min0=$((($innerSECONDS0 / 60) % 60))
      if [ $min0 -lt 10 ]; then min0="0"$min0; fi
      sec0=$(($innerSECONDS0 % 60))
      if [ $sec0 -lt 10 ]; then sec0="0"$sec0; fi
      echo "Iteration $it gradient realization time: $hr0:$min0:$sec0"
   fi


   if [ $it -lt $NOUTER ]; then
#==============================================================================
      echo ""
      echo "Extract analysis for current outer loop, and initialize the next one..."
      echo ""
#==============================================================================

      TEMPTIME=$SECONDS
      if [ $it -eq $itstart ]; then cp -v fg fg_it0; fi
      mv -v wrfvar_output wrfvar_output_$it0
      ln -sfv ./wrfvar_output_$it0 ./fg 

      if [ $WRF_MET -gt 0 ]; then
         it1=$it
         if [ $it -lt 10 ]; then it1="0"$it1; fi
         mv -v wrfbdy_d01 wrfbdy_d01_$it1
         mv -v wrfvar_bdyout wrfvar_bdyout_$it0
         ln -sfv ./wrfvar_bdyout_$it0 ./wrfbdy_d01
      fi
      FILETIME=$((FILETIME+$((SECONDS-TEMPTIME))))

      TEMPTIME=$SECONDS
      mkdir oldrsl_$it0
      mv rsl.* oldrsl_$it0
      RSLTIME=$((RSLTIME+$((SECONDS-TEMPTIME))))
   fi
   if [ $WRF_CHEM -gt 0 ]; then
      if [ $it -eq $itstart ]; then mkdir oldhxy; fi
      mv *Hx_y* oldhxy/
   fi

   hr0=$(($SECONDS / 3600))
   if [ $hr0 -lt 10 ]; then hr0="0"$hr0; fi
   min0=$((($SECONDS / 60) % 60))
   if [ $min0 -lt 10 ]; then min0="0"$min0; fi
   sec0=$(($SECONDS % 60))
   if [ $sec0 -lt 10 ]; then sec0="0"$sec0; fi
   echo "Iteration $it compute time: $hr0:$min0:$sec0"

   echo "Iteration $it rsl archive time: $RSLTIME sec."
   echo "Iteration $it wrfinput,wrfbdy update time: $FILETIME sec."

   iteration_sec=$((iteration_sec+SECONDS))
   hr0=$(($iteration_sec / 3600))
   if [ $hr0 -lt 10 ]; then hr0="0"$hr0; fi
   min0=$((($iteration_sec / 60) % 60))
   if [ $min0 -lt 10 ]; then min0="0"$min0; fi
   sec0=$(($iteration_sec % 60))
   if [ $sec0 -lt 10 ]; then sec0="0"$sec0; fi
   echo "Total accum. compute time: $hr0:$min0:$sec0"

   if [ $CLEANUP -gt 0 ]; then
      cd ../
      SECONDS=0

      if [ $CLEANUP -eq 1 ] || [ $CLEANUP -eq 2 ]; then
         tar -czf vectors_$it0.tar.gz vectors_$it0
         echo "Completed tar of vectors_$it0.tar.gz"
      fi 
      rm -r vectors_$it0

      hr0=$(($SECONDS / 3600))
      if [ $hr0 -lt 10 ]; then hr0="0"$hr0; fi
      min0=$((($SECONDS / 60) % 60))
      if [ $min0 -lt 10 ]; then min0="0"$min0; fi
      sec0=$(($SECONDS % 60))
      if [ $sec0 -lt 10 ]; then sec0="0"$sec0; fi
      echo "Iteration $it cleanup time: $hr0:$min0:$sec0"

      cleanup_sec=$((cleanup_sec+SECONDS))
      hr0=$(($cleanup_sec / 3600))
      if [ $hr0 -lt 10 ]; then hr0="0"$hr0; fi
      min0=$((($cleanup_sec / 60) % 60))
      if [ $min0 -lt 10 ]; then min0="0"$min0; fi
      sec0=$(($cleanup_sec % 60))
      if [ $sec0 -lt 10 ]; then sec0="0"$sec0; fi
      echo "Total accum. cleanup time: $hr0:$min0:$sec0"

   fi
   cd $CWD
done
## END Outer Loop

SECONDS=0

# The remainder of the script should be modified to meet the needs of your application.

#Extract important output files
OUTLOC="../../$OUTDIR/$PBS_JOBNAME"
mkdir ../../$OUTDIR
mkdir $OUTLOC
cp -v wrfvar_out* $OUTLOC/
if [ $WRF_CHEM -gt 0 ]; then cp -v AminusB* $OUTLOC/; fi
cp -v ../cost_fn $OUTLOC/
cp -v ../grad_fn $OUTLOC/
cp -r oldrsl_* $OUTLOC/
cp -v rsl.* $OUTLOC/
if [ $WRF_CHEM -gt 0 ]; then cp -v oldhxy/* $OUTLOC/; fi

if [ $CLEANUP -gt 0 ]; then
   echo "CLEANING UP SVD $(date)"
   cd $CWD/../
   if [ $CLEANUP -eq 1 ]; then 
      mkdir oldvectors
      mv vectors_0* oldvectors/
      tar -czf oldvectors.tar.gz oldvectors
      rm -r oldvectors/

      mkdir oldrundirs
      mv run.0* oldrundirs/
      tar -czf oldrundirs.tar.gz oldrundirs
      rm -r oldrundirs/
   fi
   echo "Completed Final Cleanup"
fi

hr0=$(($SECONDS / 3600))
if [ $hr0 -lt 10 ]; then hr0="0"$hr0; fi
min0=$((($SECONDS / 60) % 60))
if [ $min0 -lt 10 ]; then min0="0"$min0; fi
sec0=$(($SECONDS % 60))
if [ $sec0 -lt 10 ]; then sec0="0"$sec0; fi
echo "Current cleanup time: $hr0:$min0:$sec0"

cleanup_sec=$((cleanup_sec+SECONDS))
hr0=$(($cleanup_sec / 3600))
if [ $hr0 -lt 10 ]; then hr0="0"$hr0; fi
min0=$((($cleanup_sec / 60) % 60))
if [ $min0 -lt 10 ]; then min0="0"$min0; fi
sec0=$(($cleanup_sec % 60))
if [ $sec0 -lt 10 ]; then sec0="0"$sec0; fi
echo "Total accum. cleanup time: $hr0:$min0:$sec0"

echo "FINISHED RIOT $(date)"

exit 0
