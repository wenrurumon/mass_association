#!/bin/csh
#         <------ Setup Parameters ------>
#
#SBATCH -J mfpca2expr1         # Job name
#SBATCH -N 11                 # Total number of nodes (16 cores/node)
#SBATCH -n 167                  # Total number of tasks
#SBATCH -p normal              # Queue name
#SBATCH -o mfpca2expr1.o%j     # Name of stdout output file (%j expands to jobid)
#SBATCH -t 47:59:59            # Run time (hh:mm:ss)
#SBATCH -A Epistasis-influncing  #Specify the allocation used
#------------------------------------------------------

module load launcher
module load Rstats
module load intel/15.0.2  mvapich2/2.1
setenv EXECUTABLE     $TACC_LAUNCHER_DIR/init_launcher 
setenv CONTROL_FILE   paramlist_mfpca2expr1
setenv WORKDIR        /work/03963/zhu2/rushdata_qtl/code
#--------------------------------------------------------
#--------------------------------------------------------

#----------------
# Error Checking
#----------------

if ( ! -e $WORKDIR ) then
        echo " "
	echo "Error: unable to change to working directory."
	echo "       $WORKDIR"
	echo " "
	echo "Job not submitted."
	exit
endif

if ( ! -f $EXECUTABLE ) then
	echo " "
	echo "Error: unable to find launcher executable $EXECUTABLE."
	echo " "
	echo "Job not submitted."
	exit
endif

if ( ! -f $WORKDIR/$CONTROL_FILE ) then
	echo " "
	echo "Error: unable to find input control file $CONTROL_FILE."
	echo " "
	echo "Job not submitted."
	exit
endif


#----------------
# Job Submission
#----------------

cd $WORKDIR/
echo " WORKING DIR:   $WORKDIR/"

$TACC_LAUNCHER_DIR/paramrun $EXECUTABLE $CONTROL_FILE

echo " "
echo " Parameteric Job Complete"
echo " "
