# get current working directory, i.e. directory this script is called from
WD=$(pwd)
# get this directory, i.e. directory this scripts is placed in, i.e. the SANS directory
PTH=$(dirname $0)

# run SANS with check_N option
$PTH/SANS "$@" -M "$PTH/makefile"

# if Ns okay, we're done. Re-compile and re-run if necessary, i.e. if exit code is 3
if [ $? == "3" ]
then
	# go to SANS directory and re-compile
	echo "Re-compile with makefile_autoN in $PTH ..."
	cd $PTH
	make -f makefile_autoN
	# return to working directory and re-run
	cd $WD
	echo "Re-run in working directory $WD ..."
	$PTH/SANS "$@"
fi
