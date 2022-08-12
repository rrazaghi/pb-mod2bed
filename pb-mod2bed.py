#!/usr/bin/env python

import click
import pysam
import sys


def binerize_mod_call(call, min_prob=0.5, max_prob=0.5):
    prob = call / 255
    if prob < min_prob:
        return "unmethylated"
    elif prob > max_prob:
        return "methylated"
    else:
        return "unknown"


def get_read_names(file):
    names = set()
    file = open(file)
    for line in file:
        l = line.strip()
        names.add(l)
    # file.close()

    return names


@click.command()
@click.argument("bam", nargs=1, type=click.Path(exists=True), required=True)
@click.option(
    "-n",
    "--read_names",
    is_flag=False,
    default=None,
    required=False,
    type=click.Path(exists=True),
    help="filter analysis based on file containing read names per line",
)
@click.option(
    "-u",
    "--can_prob",
    is_flag=False,
    default=0.5,
    type=float,
    help="probability threshold for canonical bases",
)
@click.option(
    "-m",
    "--mod_prob",
    is_flag=False,
    default=0.5,
    type=float,
    help="probability threshold for modified bases",
)
@click.option(
    "-o",
    "--out",
    required=False,
    type=click.File("w"),
    default=sys.stdout,
    help="output path",
)
def pbmod2bed(bam, read_names, can_prob, mod_prob, out):
    """
    This script converts pacbio modified bam files to expanded bed files in the following format:

    read_name start end methylation_status pass_tag

    """
    if read_names:
        click.echo("\nParsing file containing read names...")
        names = get_read_names(read_names)
        num_reads = len(names)
        click.echo(f"{num_reads} reads were found!\n")
    samfile = pysam.AlignmentFile(bam, "rb", check_sq=False)
    click.echo("\nProcessing BAM...\n")
    # out_file = open(out, "w")
    i = 0
    match = 0
    for read in samfile.fetch(until_eof=True):
        i += 1
        if i % 100000 == 0:
            click.echo(f"{i} reads processed...")

        if read_names:
            if read.qname not in names:
                continue
            else:
                match += 1
        mod_dict = read.modified_bases
        pass_tag = read.get_tag("np")
        for mod_key, mod_list in mod_dict.items():
            for mod in mod_list:
                pos = mod[0]
                status = binerize_mod_call(mod[1], min_prob=can_prob, max_prob=mod_prob)
                print(
                    "\t".join(
                        [read.qname, str(pos), str(pos + 1), status, str(pass_tag)]
                    ),
                    end="\n",
                    file=out,
                )

    samfile.close()
    # out_file.close()
    click.echo(f"Total of {i} reads were processed.\n")
    if read_names:
        click.echo(f"{match} reads from read_name file were found in BAM.\n")
    click.echo("Done!\n")


if __name__ == "__main__":
    pbmod2bed()
