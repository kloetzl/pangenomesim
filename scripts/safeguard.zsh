#!/usr/bin/zsh

integer i
ITERATIONS=100
for (( i = 0; i < $ITERATIONS; i++ )); do
	./pangenomesim -o safeguard --model IMG --param loci_length=100 --param num_genomes=10,img_core_size=200 -p img_theta=200,img_rho=0.2 -p seed=$i 2> /dev/null
	cd ~/Projects/3rd/panicmage-master
	./panicmage ../../pangenomesim/safeguard/coalescent.newick ../../pangenomesim/safeguard/genefrequency.gfs 10 2> /dev/null
	cd -
done | grep 'p-value' | tee safeguard.txt | awk '/probab/{c++}END{print c/NR}'
