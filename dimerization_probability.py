import click
from pathlib import Path
from Bio import SeqIO

"""
Implementation of the genomic UV susceptibility model from the Ultraviolet
Germicidal Irradiation Handbook (Kowalski, 20009), chapter 4: UV Rate Constants

Corrections from scattering are not included. Hyperchromaticity corrections
are also not included. 

TODO:
1. Proper logging
2. Proper exception handling
3. Hyperchromaticity correction
4. Scattering correction
5. Cleaner file organization/packaging
6. Stricter filetype handling
"""


PAIRS = {"A": "T", "T": "A", "G": "C", "C": "G"}

FLANKED_PURINES = [
    "TTA",
    "ATT",
    "TTG",
    "GTT",
    "CTA",
    "ACT",
    "CTG",
    "GCT",
    "TCA",
    "ATC",
    "TCG",
    "GTC",
    "CCA",
    "ACC",
    "CCG",
    "GCC",
]


class DimerCounter:

    """
    Implementation of the genomic UV susceptibility model from the
    Ultraviolet Germicidal Irradiation Handbook (Kowalski, 20009),
    Chapter 4: UV Rate Constants

    Default values for Fa, Fb, and Fc were chosen empirically, and represent
    the relative prevalnce of TT, CT, and CC dimers.

    """

    def __init__(
        self,
        genome: str,
        num_strands: int,
        Fa: float = 0.1,
        Fb: float = 6,
        Fc: float = 4,
    ):

        # type and value checks. Should eventually be replaced by proper exception handling.
        assert (
            type(genome) == str
        ), "Genome must be a string"  # verify that this is a string
        assert type(num_strands) == int, "num_strands must be an int"
        assert num_strands in [1, 2], "num_strands must be either 1 or 2"
        assert set(genome.upper()) == set(
            "ACGT"
        ), "Genome must contain only A, C, G, and T"

        # set constants
        self.genome = genome.upper()
        self.num_strands = num_strands
        self.Fa = Fa
        self.Fb = Fb
        self.Fc = Fc

    # transcription tools
    def transcribe(self, genome: str):
        return "".join(list(map(self._lookup, self.genome)))

    def _lookup(self, char: str):
        return PAIRS[char]

    # counting basics
    def tt(self, genome: str) -> int:
        return genome.count("TT")

    def ct(self, genome: str) -> int:
        return genome.count("CT") + genome.count("TC")

    def cc(self, genome: str) -> int:
        return genome.count("CC")

    def purines(self, genome: str) -> int:
        return sum(genome.count(dimer) for dimer in FLANKED_PURINES)

    def _count_dimers(self, genome: str) -> float:
        """count pyrimidine dimers, weighted by their relative proportions"""
        return (
            self.tt(genome)
            + self.Fa * self.ct(genome)
            + self.Fb * self.cc(genome)
            + self.Fc * self.purines(genome)
        )

    def calculate(self) -> float:

        sum_dimers = self._count_dimers(self.genome)

        if self.num_strands == 2:
            transcription = self.transcribe(self.genome)
            transcribed_dimers = self._count_dimers(transcription)
            sum_dimers = sum_dimers + transcribed_dimers

        numerator = sum_dimers**0.5
        denominator = len(self.genome) ** (2 / 3)
        return numerator / denominator


def get_genome(file) -> str:
    """read the genome from a SeqIO-compatible file."""
    filetype = file.suffix[1:]
    with open(file) as input_handle:
        genome = str(list(SeqIO.parse(input_handle, filetype))[0].seq)
    return genome


def calculate_probability(file, num_strands: int, silent: bool = False) -> float:
    """Load the genome file and calculate its dimerization probability"""
    genome = get_genome(file)
    counter = DimerCounter(genome=genome, num_strands=num_strands)
    probability = counter.calculate()
    if not silent:
        print(file, round(probability,3), sep="\t")
    return probability


def write_file(file, probs: list) -> None:
    """write the list of dimerization probabilities to an output file"""
    with open(file, "w") as fp:
        for item in probs:
            # write each item on a new line
            fp.write("%s\n" % item)


@click.command()
@click.option(
    "-i",
    "--input_path",
    "input_path",
    type=click.Path(),
    required=True,
    help="""Either a single fasta file, or a directory. Note that only one 
    value for num_strands is possible, so if you have both single-stranded 
    and double-stranded viruses, you should separate them into separate
    directories.""",
)
@click.option(
    "-o", "--output_path", "output_path", type=click.Path(), default="./probability.txt"
)
@click.option(
    "-n",
    "--num_strands",
    "num_strands",
    type=click.Choice([1, 2]),
    default=1,
    help="""Whether the virus is single- or double-stranded. Default is 1, 
    which will produce conservative estimates if virus is actually 
    double-stranded.""",
)
@click.option(
    "-s",
    "--silent",
    "silent",
    type=bool,
    default=False,
)
def main(input_path, output_path, num_strands: int, silent: bool):

    path = Path(input_path)
    if path.is_dir():
    # if directory, iterate through directories
        probabilities = []
        for file in path.iterdir():
            probabilities.append(
                calculate_probability(file=file, num_strands=num_strands, silent=silent)
            )
            write_file(file=output_path, probs=probabilities)
    elif path.is_file():
        probability = calculate_probability(
            file=path, num_strands=num_strands, silent=silent
        )
        write_file(file=output_path, probs=[probability])
    else:
        raise TypeError("Could not parse input_path")


if __name__ == "__main__":

    main()
