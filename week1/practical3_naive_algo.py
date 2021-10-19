import random
from file_reading import *


def naive(pattern, text):
    occurrences = []
    for i in range(len(text) - len(pattern) + 1):
        match = True
        for j in range(len(pattern)):
            if text[i + j] != pattern[j]:
                match = False
                break
        if match:
            occurrences.append(i)
    return occurrences


def generate_reads(genome, num_reads, read_len):
    ''' Generate reads from random positions in the given genome. '''
    reads = []
    for _ in range(num_reads):
        start = random.randint(0, len(genome) - read_len) - 1
        reads.append(genome[start: start + read_len])
    return reads


def reverse_complement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t


if __name__ == "__main__":
    reference_genome = read_genome('..\\phix.fa')
    phix_reads, phix_quals = read_seq_and_qual('..\\ERR266411_1.first1000.fastq')

    t = 'AGCTTAGATAGC'
    p = 'AG'
    pattern_occurences = naive(p, t)
    print(f'Pattern "AG" appears in sequence "AGCTTAGATAGC" {len(pattern_occurences)} times at '
          f'positions {pattern_occurences}')

    # Generate 100 reads of length 100
    random_reads_from_genome = generate_reads(reference_genome, 100, 100)

    # Count how many randomly generated reads match the genome exactly (should be 100%)
    num_matched = 0
    for r in random_reads_from_genome:
        matches = naive(r, reference_genome)
        if len(matches) > 0:
            num_matched += 1
    print('%d / %d reads randomly generated from genom matched the genome exactly!' % (num_matched,
                                                                                       len(random_reads_from_genome)))

    num_matched = 0
    n = 0
    for seq in phix_reads:
        matches = naive(seq, reference_genome)
        n += 1
        if len(matches) > 0:
            num_matched += 1
    print('%d / %d sequences matched the reference genome exactly!' % (num_matched, n))

    # Now let's try matching just the first 30 bases of each read
    num_matched = 0
    n = 0
    for seq in phix_reads:
        seq = seq[:30]  # just taking the first 30 bases
        matches = naive(seq, reference_genome)
        n += 1
        if len(matches) > 0:
            num_matched += 1
    print('%d / %d sequences of first 30 bases from each sequence matched the reference genome '
          'exactly!' % (num_matched, n))

    # Now let's try matching just the first 30 bases of each read and their reverse complement
    num_matched = 0
    n = 0
    for seq in phix_reads:
        seq = seq[:30]  # just taking the first 30 bases
        matches = naive(seq, reference_genome)
        matches.extend(naive(reverse_complement(seq), reference_genome))
        n += 1
        if len(matches) > 0:
            num_matched += 1
    print('%d / %d sequences or reverse compelent of first 30 bases from each sequence matched the '
          'reference genome exactly!' % (num_matched, n))