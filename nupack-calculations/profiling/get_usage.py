import sys


def main(infile):
    with open(infile) as input:
        header = input.readline().strip()
        print(f"{header}\tlengths")
        for line in input:
            line = line.strip()
            slurm_file = line.split()[0]
            print(line + '\t' + ','.join(str(s) for s in read_slurm_stats(slurm_file)))


def read_slurm_stats(slurm_file):
    with open(slurm_file) as slurm:
        lengths = []
        recording = False
        for line in slurm:
            if recording and line.startswith('['):
                return lengths
                return len(lengths), min(lengths), max(lengths), sum(lengths)
            if recording:
                lengths.append(int(line))
            if line.startswith('Activating conda'):
                recording = True

if __name__ == '__main__':
    main(sys.argv[1])
