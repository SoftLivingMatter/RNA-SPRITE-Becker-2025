import pysam
import re
import argparse


def main():
    opts = parse_args()

    bids = generate_bid_map(opts.bid)
    add_tag(opts.input_bam, opts.out_bam, bids)


def parse_args():
    parser = argparse.ArgumentParser(description='Add tags to bam file')
    parser.add_argument('-b', '--bid', type=str, required=True,
                        help='bID file for translating barcode to sequences')
    parser.add_argument('-i', '--input-bam', type=str, required=True,
                        help='BAM path to which to add XT tag from chrom field')
    parser.add_argument('-o', '--out-bam', type=str, required=True,
                        help='Out BAM path')

    args = parser.parse_args()
    return args

def generate_bid_map(bid_file):
    with open(bid_file, 'r') as bid:
        bid.readline()  # header
        bid.readline()  # blank line
        result = {}
        for line in bid:
            tokens = line.split()
            result[tokens[1]] = tokens[2]

    return result

def add_tag(input_bam, output_bam, bids):
    tag_regex = re.compile(r'\]\[|\[|\]')
    with pysam.AlignmentFile(input_bam, "rb") as in_bam, \
        pysam.AlignmentFile(output_bam, "wb", template=in_bam) as out_bam:

        for read in in_bam.fetch(until_eof = True):
            tags = tag_regex.split(read.query_name)[1:-1]
            seq = ''.join(bids[tag] for tag in tags)
            read.tags += [('BC', seq)]
            out_bam.write(read)

if __name__ == "__main__":
    main()
