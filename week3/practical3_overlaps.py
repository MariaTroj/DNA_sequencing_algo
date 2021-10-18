from itertools import permutations


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


def naive_overlap_map(reads, k):
    olaps = {}
    for a, b in permutations(reads, 2):
        olen = overlap(a, b, min_length=k)
        if olen > 0:
            olaps[(a, b)] = olen
    return olaps


if __name__ == "__main__":
    print("Overlap length of strings TTACGT and CGTGTGC is", overlap('TTACGT', 'CGTGTGC'))
    print("Overlap length of strings TTACGT and GTGTCG is", overlap('TTACGT', 'GTGTGC'))


    reads = ['ACGGATC', 'GATCAAGT', 'TTCACGGA']
    print("", naive_overlap_map(reads, 3))