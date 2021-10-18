import itertools
import timeit

import week4.homework


def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1  # move just past previous match


def scs(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = None
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss) - 1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i + 1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i + 1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
    return shortest_sup  # return shortest


def scs_list(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup_list = []
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss) - 1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i + 1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i + 1][olen:]
        shortest_sup_list.append(sup)  # found shorter superstring
    shortest_sup_list = sorted(shortest_sup_list, key=lambda x: len(x))
    shortest_sup_list = [sup for sup in shortest_sup_list if len(sup) == len(shortest_sup_list[0])]
    return shortest_sup_list  # return shortest


def read_sequences(filename):
    sequences = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            fh.readline()  # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
    return sequences


# class GreedySCS:
#
#     def __init__(self, reads, min_overlap_length):
#         self.kmers: dict = {}
#         self.overlap_pairs: list = []
#         self.reads = reads
#         self.min_overlap_length = min_overlap_length
#
#     def prepare_kmers(self):
#         for read in self.reads:
#             for kmer_start in range(0, len(read) - self.min_overlap_length + 1):
#                 kmer = read[kmer_start:kmer_start + self.min_overlap_length]
#                 self.kmers.setdefault(kmer, set()).add(read)
#
#     def update_kmers(self, read_a, read_b, new_read):
#         for read in (read_a, read_b):
#             if read not in self.reads:
#                 for kmer_start in range(0, len(read) - self.min_overlap_length + 1):
#                     kmer = read[kmer_start:kmer_start + self.min_overlap_length]
#                     try:
#                         self.kmers[kmer].remove(read)
#                     except:
#                         continue
#
#         for kmer_start in range(0, len(new_read) - self.min_overlap_length + 1):
#             kmer = new_read[kmer_start:kmer_start + self.min_overlap_length]
#             self.kmers.setdefault(kmer, set()).add(new_read)
#
#     def overlap(self, a, b):
#         """ Return length of longest suffix of 'a' matching
#             a prefix of 'b' that is at least 'min_length'
#             characters long.  If no such overlap exists,
#             return 0. """
#         start = 0  # start all the way at the left
#         while True:
#             start = a.find(b[:self.min_overlap_length], start)  # look for b's suffx in a
#             if start == -1:  # no more occurrences to right
#                 return 0
#             # found occurrence; check for full suffix/prefix match
#             if b.startswith(a[start:]):
#                 return len(a) - start
#             start += 1  # move just past previous match
#
#     def overlap_all_pairs(self):
#         for read in self.reads:
#             suffix = read[-self.min_overlap_length:]
#             for possible_overlaping_read in self.kmers[suffix]:
#                 if possible_overlaping_read != read:
#                     overlap_len = self.overlap(read, possible_overlaping_read)
#                     if overlap_len != 0:
#                         self.overlap_pairs.append((read, possible_overlaping_read, overlap_len))
#         self.overlap_pairs = sorted(self.overlap_pairs, key=lambda x: x[2])
#
#     def update_overlap_pairs(self, read_a, read_b, new_read):
#         for read in (read_a, read_b):
#             if read not in self.reads:
#                 pairs_to_remove = []
#                 for pair in self.overlap_pairs:
#                     if pair[0] is read or pair[1] is read:
#                         pairs_to_remove.append(pair)
#                 for pair in pairs_to_remove:
#                     self.overlap_pairs.remove(pair)
#             else:
#                 pairs_to_remove = set()
#                 for pair in self.overlap_pairs:
#                     if pair[0] == read or pair[1] == read:
#                         pairs_to_remove.add(pair)
#                 for pair in pairs_to_remove:
#                     self.overlap_pairs.remove(pair)
#
#         suffix = new_read[-self.min_overlap_length:]
#         for possible_overlaping_read in self.kmers[suffix]:
#             if possible_overlaping_read != new_read:
#                 overlap_len = self.overlap(new_read, possible_overlaping_read)
#                 if overlap_len != 0:
#                     self.overlap_pairs.append((new_read, possible_overlaping_read, overlap_len))
#
#         self.overlap_pairs = sorted(self.overlap_pairs, key=lambda x: x[2])
#
#     def pick_maximal_overlap(self):
#         """ Return a pair of reads from the list with a
#             maximal suffix/prefix overlap >= k.  Returns
#             overlap length 0 if there are no such overlaps."""
#         try:
#             return self.overlap_pairs[-1]
#         except:
#             return ('', '', 0)
#
#     def greedy_scs(self):
#         """ Greedy shortest-common-superstring merge.
#             Repeat until no edges (overlaps of length >= k)
#             remain. """
#         self.prepare_kmers()
#         self.overlap_all_pairs()
#         read_a, read_b, olen = self.pick_maximal_overlap()
#         while olen > 0:
#             try:
#                 self.reads.remove(read_a)
#                 self.reads.remove(read_b)
#             except:
#                 print(read_a, read_b)
#             self.reads.append(read_a + read_b[olen:])
#             self.update_kmers(read_a, read_b, read_a + read_b[olen:])
#             self.update_overlap_pairs(read_a, read_b, read_a + read_b[olen:])
#             read_a, read_b, olen = self.pick_maximal_overlap()
#         return ''.join(self.reads)

class GreedySCS:
    def __init__(self, reads, k):
        self.reads = reads
        self.k = k
        self.kmers: dict[str: set] = {}
        self.overlap_pairs: dict[tuple[str, str]: int] = {}
        self.overlap_all_pairs()

    def overlap(self, a, b):
        """ Return length of longest suffix of 'a' matching
            a prefix of 'b' that is at least 'min_length'
            characters long.  If no such overlap exists,
            return 0. """
        start = 0  # start all the way at the left
        while True:
            start = a.find(b[:self.k], start)  # look for b's prefix in a
            if start == -1:  # no more occurrences to right
                return 0
            # found occurrence; check for full suffix/prefix match
            if b.startswith(a[start:]):
                return len(a) - start
            start += 1  # move just past previous match

    def overlap_all_pairs(self):
        self.kmers = {}
        self.overlap_pairs = {}
        for read in self.reads:
            for kmer_start in range(0, len(read) - self.k + 1):
                kmer = read[kmer_start:kmer_start + self.k]
                if kmer not in self.kmers.keys():
                    self.kmers[kmer] = set()
                self.kmers[kmer].add(read)
        for read in self.reads:
            suffix = read[-self.k:]
            for possible_overlaping_read in self.kmers[suffix]:
                if possible_overlaping_read != read:
                    overlap_len = self.overlap(read, possible_overlaping_read)
                    if overlap_len != 0:
                        self.overlap_pairs[(read, possible_overlaping_read)] = overlap_len

    def update_kmers_and_overlap_pairs(self, read_a, read_b, new_read):
        """ Every time after you merge read a and b, just update this dictionary by pop() any
        entry with either a or b, and add entries of (merge_of_ab, read):olen and (read,
        merge_of_ab) of each read. """
        for read in read_a, read_b:
            for kmer_start in range(0, len(read) - self.k + 1):
                kmer = read[kmer_start:kmer_start + self.k]
                try:
                    self.kmers[kmer].remove(read)
                except:
                    continue
        for kmer_start in range(0, len(new_read) - self.k + 1):
            kmer = new_read[kmer_start:kmer_start + self.k]
            self.kmers[kmer].add(new_read)

        self.overlap_pairs = {key: value for key, value in self.overlap_pairs.items()
                              if ((read_a not in key) and (read_b not in key))}
        suffix = new_read[-self.k:]
        for possible_overlaping_read in self.kmers[suffix]:
            if possible_overlaping_read is not new_read:
                overlap_len = self.overlap(new_read, possible_overlaping_read)
                if overlap_len != 0:
                    self.overlap_pairs[(new_read, possible_overlaping_read)] = overlap_len
        prefix = new_read[:self.k]
        for possible_overlaping_read in self.kmers[prefix]:
            if possible_overlaping_read is not new_read:
                overlap_len = self.overlap(possible_overlaping_read, new_read)
                if overlap_len != 0:
                    self.overlap_pairs[(possible_overlaping_read, new_read)] = overlap_len

    def pick_maximal_overlap_new(self):
        """ Return a pair of reads from the list with a
            maximal suffix/prefix overlap >= k.  Returns
            overlap length 0 if there are no such overlaps."""
        try:
            read_a, read_b = max(self.overlap_pairs, key=self.overlap_pairs.get)
            return read_a, read_b, self.overlap_pairs[(read_a, read_b)]
        except:
            return None, None, 0

    def greedy_scs_new(self):
        """ Greedy shortest-common-superstring merge.
            Repeat until no edges (overlaps of length >= k)
            remain. """
        self.overlap_all_pairs()
        read_a, read_b, olen = self.pick_maximal_overlap_new()
        while olen > 0:
            self.reads.remove(read_a)
            self.reads.remove(read_b)
            self.reads.append(read_a + read_b[olen:])
            self.update_kmers_and_overlap_pairs(read_a, read_b, read_a + read_b[olen:])
            read_a, read_b, olen = self.pick_maximal_overlap_new()
        return ''.join(self.reads)

    def pick_maximal_overlap(self):
        """ Return a pair of reads from the list with a
            maximal suffix/prefix overlap >= k.  Returns
            overlap length 0 if there are no such overlaps."""
        reada, readb = None, None
        best_olen = 0
        for a, b in itertools.permutations(self.reads, 2):
            olen = overlap(a, b, min_length=self.k)
            if olen > best_olen:
                reada, readb = a, b
                best_olen = olen
        return reada, readb, best_olen

    def greedy_scs(self):
        """ Greedy shortest-common-superstring merge.
            Repeat until no edges (overlaps of length >= k)
            remain. """
        read_a, read_b, olen = self.pick_maximal_overlap()
        while olen > 0:
            self.reads.remove(read_a)
            self.reads.remove(read_b)
            self.reads.append(read_a + read_b[olen:])
            read_a, read_b, olen = self.pick_maximal_overlap()
        return ''.join(self.reads)


if __name__ == "__main__":
    # print(len(reads), len(set(reads)))
    # print(len(scs(['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT'])))
    # print(len(scs_list(['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT'])))
    setup = '''text = 'Czekolada z solonymi cząstkami karmelu, z solonymi kawałkami karmelu i orzechami laskowymi'
    reads = [text[i:i+7] for i in range(len(text) - 6)]
    gscs = GreedySCS(reads, 6)'''
    code_to_execute = '''gscs.greedy_scs()'''
    gscs = timeit.repeat(stmt=code_to_execute.replace('    ', ''), setup=setup.replace('    ', ''),
                         globals=globals(), number=1000, repeat=100)
    print('GSCS average: ', sum(gscs) / len(gscs))

    code_to_execute = '''gscs.greedy_scs_new()'''
    gscs = timeit.repeat(stmt=code_to_execute.replace('    ', ''), setup=setup.replace('    ', ''),
                         globals=globals(), number=1000, repeat=100)
    print('GSCS new average: ', sum(gscs) / len(gscs))

    setup = '''text = 'Czekolada z solonymi cząstkami karmelu, z solonymi kawałkami karmelu i orzechami laskowymi'
        reads = [text[i:i+7] for i in range(len(text) - 6)]
        gscs = GreedySCS(reads, 3)'''
    code_to_execute = '''gscs.greedy_scs_new()'''
    gscs = timeit.repeat(stmt=code_to_execute.replace('    ', ''), setup=setup.replace('    ', ''),
                         globals=globals(), number=1000, repeat=100)
    print('GSCS new average: ', sum(gscs) / len(gscs))



    # setup = 'reads = read_sequences("ads1_week4_reads.fq")'
    # code_to_execute = '''
    # gscs = GreedySCS(reads, 30)
    # assembled_genome = gscs.greedy_scs()
    # print(len(assembled_genome))
    # print(assembled_genome.count('A'))
    # print(assembled_genome.count('T'))
    # '''
    # k_30 = timeit.repeat(stmt=code_to_execute.replace('    ', ''), setup=setup.replace('    ',''),
    #                      globals=globals(), repeat=1)
    # print('k = 30: ', k_30, ' Average: ', sum(k_30) / len(k_30))
    #
    # code_to_execute = '''
    #     gscs = GreedySCS(reads, 70)
    #     assembled_genome = gscs.greedy_scs()
    #     print(len(assembled_genome))
    #     print(assembled_genome.count('A'))
    #     print(assembled_genome.count('T'))
    #     '''
    # k_70 = timeit.repeat(stmt=code_to_execute.replace('    ', ''), setup=setup.replace('    ', ''),
    #                      globals=globals(), repeat=1)
    # print('k = 70: ', k_70, ' Average: ', sum(k_70) / len(k_70))

    # reads = read_sequences("ads1_week4_reads.fq")
    # gscs = GreedySCS(reads, 80)
    # assembled_genome = gscs.greedy_scs()
    # print(len(assembled_genome))
    # print(assembled_genome.count('A'))
    # print(assembled_genome.count('T'))
