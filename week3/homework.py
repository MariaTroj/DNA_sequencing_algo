from file_reading import *


def local_alignment(pattern, text):
    '''
    Dynamic programming to find approximate occurrences of a pattern in a text
    The minimal value in the bottom row is the edit distance of the closest match between P and T
    '''
    # Create distance matrix
    D = []
    for i in range(len(pattern) + 1):
        D.append([0] * (len(text) + 1))
    # Initialize first row and column of matrix
    for i in range(len(pattern) + 1):
        D[i][0] = i
    for i in range(len(text) + 1):
        D[0][i] = 0
    # Fill in the rest of the matrix
    for i in range(1, len(pattern) + 1):
        for j in range(1, len(text) + 1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if pattern[i - 1] == text[j - 1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the value in the bottom right corner of the matrix
    return min(D[-1])


def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match


def overlap_all_pairs(reads, min_overlap_length):
    kmers = {}
    overlap_pairs = []
    for read in reads:
        for kmer_start in range(0, len(read) - min_overlap_length + 1):
            kmer = read[kmer_start:kmer_start + min_overlap_length]
            if kmer not in kmers.keys():
                kmers[kmer] = set()
            kmers[kmer].add(read)
    for read in reads:
        suffix = read[-min_overlap_length:]
        for possible_overlaping_read in kmers[suffix]:
            if possible_overlaping_read is not read:
                if overlap(read, possible_overlaping_read, min_overlap_length) != 0:
                    overlap_pairs.append((read, possible_overlaping_read))
    return overlap_pairs


if __name__ == "__main__":
    reads, qualities = read_seq_and_qual("..\\ERR266411_1.for_asm.fastq")
    overlap_pairs = overlap_all_pairs(reads, 30)
    print(len(overlap_pairs))