#!/usr/bin/zsh

# cd bar

genomes=(genome*)
shift genomes


cat genome001.fasta > synthetic.fasta

for file in $genomes; do
	echo $file

	makeblastdb -dbtype nucl -in synthetic.fasta &> /dev/null

	blastn -db synthetic.fasta -outfmt 7 -query $file |
		grep '^# 0 hits' -B 2 |
		awk '/Query/{print $NF}' > non_matches.txt

	xargs -a non_matches.txt -d'\n' echo |
		tr ' ' '|' |
		xargs -I '{}' getSeq -s '{}' $file >> synthetic.fasta
done
