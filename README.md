# Genomic UV Susceptibility Model

Implementation of the genomic UV susceptibility model from the Ultraviolet
Germicidal Irradiation Handbook (Kowalski, 20009), chapter 4: UV Rate Constants.

Corrections from scattering and hyperchromaticity are not included.

## Usage

Simply pass the sequence file, or the name of the directory where the sequence files 
are stored. Note that only one value of `num_strands` is possible for now, so if you
wish to calculate the dimerization probabilities of a mix of single- and double-stranded
viruses, you should separate them into separate directories.

`python dimerization_probability.py --input_path my_viral_genome.fasta --num_strands 1`

You may also optionally indicate an `--output_path` argument

## TODO:
1. Proper logging
2. Proper exception handling
3. Hyperchromaticity correction
4. Cleaner file organization/packaging
5. Stricter filetype handling
6. Scattering correction
