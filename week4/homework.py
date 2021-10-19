import itertools
import timeit
from week4.practical1_shortest_common_superstring import overlap, scs
from file_reading import read_seq_and_qual


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


class GreedySCS:
    def __init__(self, reads: list[str], k: int):
        self.reads = reads
        self.k = k
        self.kmers: dict[str: set] = {}
        self.overlap_pairs: dict[tuple[str, str]: int] = {}
        self.overlap_all_pairs()

    def overlap(self, a: str, b: str) -> int:
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

    def pick_maximal_overlap(self) -> (str, str, int):
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

    def greedy_scs(self) -> str:
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

    def overlap_all_pairs(self):
        """
        kmers dictionary is filled according to the rule for each read in reads:
         - key: k-mer from read,
         - value: set of reads containing this k-mer
        overlap_pairs dictionary is filled according to the rule for each read in reads and for
        each next_read from kmers[suffix of read]:
         - key: (read, next_read),
         - value: overlap length. Entry key: value is added only if overlap length is not 0
        """
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

    def update_kmers(self, read_a: str, read_b: str, new_read: str):
        """
        Every time after merge read_a and read_b, kmers dict is updated:
         - read_a and read_b are deleted from kmers for k-mers from read_a and read_b
         - new_read is added to kmers for k-mer from new_read
        """
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

    def update_overlap_pairs(self, read_a: str, read_b: str, new_read: str):
        """
        Every time after merge read_a and read_b, kmers dict is updated:
         - every entry with either read_a or read_b in the key is deleted
         - entries: (new_read, read_with_overlaping_prefix): olen and
         (read_with_overlaping_suffix, new_read): olen are added.
         """
        self.overlap_pairs = {key: value for key, value in self.overlap_pairs.items()
                              if ((read_a not in key) and (read_b not in key))}
        suffix = new_read[-self.k:]
        for possible_overlaping_read in self.kmers[suffix]:
            if possible_overlaping_read != new_read:
                overlap_len = self.overlap(new_read, possible_overlaping_read)
                if overlap_len != 0:
                    self.overlap_pairs[(new_read, possible_overlaping_read)] = overlap_len
        prefix = new_read[:self.k]
        for possible_overlaping_read in self.kmers[prefix]:
            if possible_overlaping_read != new_read:
                overlap_len = self.overlap(possible_overlaping_read, new_read)
                if overlap_len != 0:
                    self.overlap_pairs[(possible_overlaping_read, new_read)] = overlap_len

    def pick_maximal_overlap_new(self) -> (str, str, int):
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
            self.update_kmers(read_a, read_b, read_a + read_b[olen:])
            self.update_overlap_pairs(read_a, read_b, read_a + read_b[olen:])
            read_a, read_b, olen = self.pick_maximal_overlap_new()
        return ''.join(self.reads)


if __name__ == "__main__":
    print("Text: 'Czekolada z solonymi cząstkami karmelu, z solonymi kawałkami karmelu i "
          "orzechami laskowymi'\nReads length: 7")
    setup = '''text = 'Czekolada z solonymi cząstkami karmelu, z solonymi kawałkami karmelu i orzechami laskowymi'
    reads = [text[i:i+7] for i in range(len(text) - 6)]
    gscs = GreedySCS(reads, 6)'''
    code_to_execute = '''gscs.greedy_scs()'''
    gscs = timeit.repeat(stmt=code_to_execute.replace('    ', ''), setup=setup.replace('    ', ''),
                         globals=globals(), number=1000, repeat=100)
    print('GSCS average (minimum overlap length = 6): ', sum(gscs) / len(gscs))

    code_to_execute = '''gscs.greedy_scs_new()'''
    gscs = timeit.repeat(stmt=code_to_execute.replace('    ', ''), setup=setup.replace('    ', ''),
                         globals=globals(), number=1000, repeat=100)
    print('GSCS new average (minimum overlap length = 6): ', sum(gscs) / len(gscs))

    setup = '''text = 'Czekolada z solonymi cząstkami karmelu, z solonymi kawałkami karmelu i orzechami laskowymi'
        reads = [text[i:i+7] for i in range(len(text) - 6)]
        gscs = GreedySCS(reads, 3)'''
    code_to_execute = '''gscs.greedy_scs_new()'''
    gscs = timeit.repeat(stmt=code_to_execute.replace('    ', ''), setup=setup.replace('    ', ''),
                         globals=globals(), number=1000, repeat=100)
    print('GSCS new average (minimum overlap length = 3): ', sum(gscs) / len(gscs))

    # setup = 'reads, quals = read_seq_and_qual("ads1_week4_reads.fq")'
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

    # reads, quals = read_seq_and_quals("ads1_week4_reads.fq")
    # gscs = GreedySCS(reads, 80)
    # assembled_genome = gscs.greedy_scs()
    # print(len(assembled_genome))
    # print(assembled_genome.count('A'))
    # print(assembled_genome.count('T'))
