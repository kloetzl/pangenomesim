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

    pangenomesim -o dir

A new directory `dir` will be created and filled with various files.

- `ref.fasta`: A technical reference genome for simulation; Union of core and accessory genes.
- `genome*.fasta`: The individual simulated bacterial genomes.
- `accessory.fasta`: The dispensable part of the reference genome.
- `core.fasta`: Genes shared by all genomes.
- `alignment.maf`: Sequence alignment of all genes.
- `coalescent.newick`: Evolutionary history of the genomes.
- `distance.mat`: Matrix of evolutionary distances.
- `genefrequency.gfs`: Distribution of gene frequencies.
- `reproducible.seed`: Parameters to reproduce this simulation.

## License

BSD 2-Clause License

Copyright (c) 2017, Fabian Klötzl <fabian-pangenomesim@kloetzl.info>

See [License](License) for full detail.
