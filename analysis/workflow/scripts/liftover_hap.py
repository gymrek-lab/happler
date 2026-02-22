#!/usr/bin/env python

from pathlib import Path

import click
from liftover import ChainFile


@click.command()
@click.argument("hap", type=click.Path(exists=True, path_type=Path))
@click.argument("chain", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--update-varIDs",
    is_flag=True,
    default=False,
    show_default=True,
    help="Also update the variant IDs to use CHROM:POS format",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("/dev/stdout"),
    show_default="stdout",
    help="A haplotypes file with its coordinates lifted over",
)
def main(
    hap: Path,
    chain: Path,
    update_varids: bool = False,
    output: Path = Path("/dev/stdout"),
):
    """
    Lifts over the coordinates in a .hap file using a .chain.gz file.
    """

    converter = ChainFile(chain, one_based=True)

    with open(hap, "r") as hap_file:
        with open(output, "w") as out_file:
            chrom = None
            # Iterate over each line in the input file
            for line in hap_file:
                first_char = line[0]
                line = line.strip("\n")
                if not first_char == "#" and (
                    first_char == "H" or first_char == "V" or first_char == "R"
                ):
                    line_parts = line.split("\t")
                    if first_char == "H":
                        chrom = line_parts[1]
                    line_parts[2] = str(converter[chrom][int(line_parts[2])][0][1])
                    line_parts[3] = str(converter[chrom][int(line_parts[3])][0][1])
                    if update_varids and (first_char == "V" or first_char == "R"):
                        line_parts[4] = f"{chrom}:{line_parts[2]}"
                    line = "\t".join(line_parts)
                out_file.write(line + "\n")


if __name__ == "__main__":
    main()
