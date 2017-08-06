#!/usr/bin/zsh

PM_DIR="~/Projects/3rd/panicmage-master"
PGS_DIR="$PWD/.."

cd PGS_DIR

integer i
ITERATIONS=100
for (( i = 0; i < $ITERATIONS; i++ ));
do
	$PGS_DIR/pangenomesim -o $PGS_DIR/safeguard --param gene-length=100 --param num-genomes=10,core-size=200 -p theta=200,rho=0.2 -p seed=$i 2> /dev/null
	$PM_DIR/panicmage $PGS_DIR/safeguard/coalescent.newick $PGS_DIR/safeguard/genefrequency.gfs 10 2> /dev/null
done |
	grep 'p-value' |
	tee safeguard.txt |
	awk '/probab/{c++}END{print c/NR}'
