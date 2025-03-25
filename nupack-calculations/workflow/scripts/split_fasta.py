from Bio import SeqIO
from pathlib import Path
from collections import defaultdict
import sys


def split_monomers(input_fasta, output_directory: Path):
    output_directory.mkdir(exist_ok=True)

    # split into 10 sequences
    to_write = []
    max_len = 0
    too_large = []
    counter = 1
    for fasta in SeqIO.parse(open(input_fasta), 'fasta'):
        fasta_len = len(fasta.seq)
        if fasta_len > 40_000:
            too_large.append(fasta)
            continue
        to_write.append(fasta)
        max_len = max(max_len, fasta_len)
        if len(to_write) >= 10:
            SeqIO.write(to_write, output_directory / f"{counter}_{max_len}.fa", 'fasta')
            to_write.clear()
            counter += 1
            max_len = 0

    # last group
    if len(to_write) > 0:
        SeqIO.write(to_write, output_directory / f"{counter}_{max_len}.fa", 'fasta')
    SeqIO.write(too_large, output_directory / f"too_large.fa.skipped", 'fasta')


def split_from_csv(input_fasta, dimer_csv, output_directory: Path):
    output_directory.mkdir(exist_ok=True)

    # read in fasta to pull dimer pairs from
    fastas = {
        fasta.name: fasta
        for fasta in SeqIO.parse(open(input_fasta), 'fasta')
    }

    # split into 10 sequences
    to_write = []
    max_len = 0
    too_large = []
    counter = 1
    with open(dimer_csv) as pairs:
        header = pairs.readline()
        num_seqs = 2
        if header.startswith('genes'):
            num_seqs = 1
        for line in pairs:

            tokens = line.strip().split(',')
            lengths = [len(fastas[seq]) for seq in tokens[:num_seqs]]
            if sum(lengths) > 40_000:
                for seq in tokens[:num_seqs]:
                    too_large.append(fastas[seq])
                continue

            for seq in tokens[:num_seqs]:
                to_write.append(fastas[seq])

            # max length is the total length of the dimers
            max_len = max(max_len, sum(lengths))
            if len(to_write) >= 10:
                SeqIO.write(to_write, output_directory / f"{counter}_{max_len}.fa", 'fasta')
                to_write.clear()
                counter += 1
                max_len = 0

        # last group
        if len(to_write) > 0:
            SeqIO.write(to_write, output_directory / f"{counter}_{max_len}.fa", 'fasta')
        SeqIO.write(too_large, output_directory / f"too_large.fa.skipped", 'fasta')


if __name__ == '__main__':
    if len(sys.argv) == 3:
        split_monomers(sys.argv[1], Path(sys.argv[2]))
    elif len(sys.argv) == 4:
        split_from_csv(sys.argv[1], sys.argv[2], Path(sys.argv[3]))
    else:
        raise ValueError('Need 2 or 3 arguments')
