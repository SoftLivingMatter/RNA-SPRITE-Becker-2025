import nupack
from Bio import SeqIO
import sys


def get_monomer_landscape_FE(seq: str):
    cset1 = nupack.ComplexSet(
        strands=[nupack.Strand(seq, name='seq1')],
        complexes=nupack.SetSpec(max_size=1),
    )
    complex_analysis_result = nupack.complex_analysis(
        cset1,
        model=nupack.Model(material='rna', kelvin=310.15),
        compute=['pfunc'],
    )

    monomer_FE = complex_analysis_result[cset1.complexes[0]].free_energy

    return monomer_FE

def main(input_fasta, output_file):
    with open(output_file, 'w') as output:
        for fasta in SeqIO.parse(open(input_fasta), 'fasta'):
            name, sequence = fasta.id, str(fasta.seq).replace('N', '')
            print(len(sequence))
            deltaG = get_monomer_landscape_FE(sequence)
            output.write(f"{name}\t{deltaG}\n")

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
