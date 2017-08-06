# PanGenomeSim

Simulate bacterial pangenomes according to the 'infinitely many genes' model.

## Installation

    autoreconf -fi
    ./configure
    make
    make install  # optional, may require sudo

Verify that the installation was successful by running the following.

    pangenomesim --help

This will give you a list of available options and parameters.

## Usage

    pangenomesim --out-dir sim

A new directory `sim` will be created and filled with various files. If no directory parameter was given, the current directory will be populated.

- `ref.fasta`: A technical reference genome for simulation; Union of core and accessory genes.
- `genome*.fasta`: The individual simulated bacterial genomes.
- `accessory.fasta`: The dispensable part of the reference genome.
- `core.fasta`: Genes shared by all genomes.
- `alignment.maf`: Sequence alignment of all genes.
- `coalescent.newick`: Evolutionary history of the genomes.
- `distance.mat`: Matrix of evolutionary distances.
- `genefrequency.gfs`: Distribution of gene frequencies.
- `reproducible.seed`: Parameters to reproduce this simulation.

Different parameters can be used to change the simulation. Use `-p KEY=VALUE` or `--param KEY=VALUE` to set a new value.

- `core-size`: Simulate additional core genes.
- `gene-length`: Set the length of a gene.
- `mut-rate`: Change the nucleotide substitution rate.
- `num-genomes`: Set the total number of genomes.
- `rho`: Rate of gene loss; See IMG model.
- `seed`: Seed for RNG.
- `theta`: Rate of gene gain; See IMG model.

The following two calls are equivalent.

    pangenomesim --out-dir sim -p core-size=10 --param gene-length=2000
    pangenomesim -o sim -p core-size=10,gene-length=2000

## License

BSD 2-Clause License

Copyright (c) 2017, Fabian Klötzl <fabian-pangenomesim@kloetzl.info>

See the License file for full detail.
