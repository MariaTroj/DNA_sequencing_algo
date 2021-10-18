import random
from practical1 import readFastq


def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


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


def generateReads(genome, numReads, readLen):
    ''' Generate reads from random positions in the given genome. '''
    reads = []
    for _ in range(numReads):
        start = random.randint(0, len(genome) - readLen) - 1
        reads.append(genome[start: start + readLen])
    return reads


def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t


if __name__ == "__main__":
    reference_genome = readGenome('..\\phix.fa')
    phix_reads, phix_quals = readFastq('..\\ERR266411_1.first1000.fastq')

    t = 'AGCTTAGATAGC'
    p = 'AG'
    pattern_occurences = naive(p, t)
    print(f'Pattern "AG" appears in sequence "AGCTTAGATAGC" {len(pattern_occurences)} times at '
          f'positions {pattern_occurences}')

    # Generate 100 reads of length 100
    random_reads_from_genome = generateReads(reference_genome, 100, 100)

    # Count how many randomly generated reads match the genome exactly (should be 100%)
    numMatched = 0
    for r in random_reads_from_genome:
        matches = naive(r, reference_genome)
        if len(matches) > 0:
            numMatched += 1
    print('%d / %d reads randomly generated from genom matched the genome exactly!' % (numMatched,
                                                          len(random_reads_from_genome)))

    numMatched = 0
    n = 0
    for seq in phix_reads:
        matches = naive(seq, reference_genome)
        n += 1
        if len(matches) > 0:
            numMatched += 1
    print('%d / %d sequences matched the reference genome exactly!' % (numMatched, n))

    # Now let's try matching just the first 30 bases of each read
    numMatched = 0
    n = 0
    for seq in phix_reads:
        seq = seq[:30]  # just taking the first 30 bases
        matches = naive(seq, reference_genome)
        n += 1
        if len(matches) > 0:
            numMatched += 1
    print('%d / %d sequences of first 30 bases from each sequence matched the reference genome '
          'exactly!' % (numMatched, n))

    # Now let's try matching just the first 30 bases of each read and their reverse complement
    numMatched = 0
    n = 0
    for seq in phix_reads:
        seq = seq[:30]  # just taking the first 30 bases
        matches = naive(seq, reference_genome)
        matches.extend(naive(reverseComplement(seq), reference_genome))
        n += 1
        if len(matches) > 0:
            numMatched += 1
    print('%d / %d sequences or reverse compelent of first 30 bases from each sequence matched the '
          'reference genome exactly!' % (numMatched, n))