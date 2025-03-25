import nupack
from Bio import SeqIO
import sys

def main(input_fasta, output_file):
    with open(output_file, 'w') as output:
        for fasta1, fasta2 in grouper(SeqIO.parse(open(input_fasta), 'fasta'), 2):
            name1, sequence1 = fasta1.id, str(fasta1.seq).replace('N', '')
            name2, sequence2 = fasta2.id, str(fasta2.seq).replace('N', '')
            print(len(sequence1))
            print(len(sequence2))
            deltaG = get_dimer_landscape_FE(sequence1, sequence2)
            output.write(f"{name1}\t{name2}\t{deltaG}\n")


def get_dimer_landscape_FE(seq1, seq2):

    nupack_model = nupack.Model(material='rna',
                                kelvin=310.15,
                                )
    nupack_seqs = [
        nupack.Strand(seq, name=f'seq{e}')
        for e, seq in enumerate((seq1, seq2), start=1)
    ]

    to_exclude = [
        nupack.Complex([nupack_seqs[0]], name='(seq1)'),
        nupack.Complex([nupack_seqs[1]], name='(seq2)')
        ]

    to_exclude += [nupack.Complex([nupack_seqs[1], nupack_seqs[1]], name='(seq2+seq2)')]
    if seq1 == seq2:  # Then we care about the homodimer
        to_exclude += [nupack.Complex([nupack_seqs[0], nupack_seqs[1]], name='(seq1+seq2)')]
    else:  # Then we care about the heterodimer
        to_exclude += [nupack.Complex([nupack_seqs[0], nupack_seqs[0]], name='(seq1+seq1)')]
       
    cset1 = nupack.ComplexSet(
        strands=nupack_seqs,
        complexes=nupack.SetSpec(max_size=2, exclude=to_exclude),
    )
    complex_analysis_result = nupack.complex_analysis(
        cset1, model=nupack_model, compute=['pfunc'])

    if seq1 == seq2:  # Then we are actually considering homodimers
        return complex_analysis_result['(seq1+seq1)'][2]
    else:
        heterodimer_name = [i.name for i in cset1 if
                            'seq1' in i.name and 'seq2' in i.name][0]  # can be either '(seq1+seq2)' or '(seq2+seq1)'
        return complex_analysis_result[heterodimer_name][2]


def grouper(iterable, n):
    "Collect data into non-overlapping fixed-length chunks or blocks."
    # grouper('ABCDEFG', 3, incomplete='strict') → ABC DEF ValueError
    iterators = [iter(iterable)] * n
    return zip(*iterators)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
