# qsub -q high.q -pe smp 12 getTn.sh

START=$(date +%s) 
R CMD BATCH --no-save getTn.R getTn.Rout
END=$(date +%s)
DIFF=$(( ${END} - ${START} ))
echo "getTn.R took ${DIFF} seconds"