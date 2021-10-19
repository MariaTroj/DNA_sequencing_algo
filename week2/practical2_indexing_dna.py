import bisect

class Index(object):
    def __init__(self, text: str, k: int):
        ''' Create index from all substrings of size 'length' '''
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(text) - k + 1):  # for each k-mer
            self.index.append((text[i:i + k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer

    def query(self, pattern):
        ''' Return index hits for first k-mer of P '''
        kmer = pattern[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits


def query_index(pattern: str, text: str, index: Index):
    k = index.k
    offsets = []
    for i in index.query(pattern):
        if pattern[k:] == text[i + k:i + len(pattern)]:  # verify that rest of P matches
            offsets.append(i)
    return offsets

if __name__ == "__main__":
    t = 'ACTTGGAGATCTTTGAGGCTAGGTATTCGGGATCGAAGCTCATTTCGGGGATCGATTACGATATGGTGGGTATTCGGGA'
    p = 'GGTATTCGGGA'
    index = Index(t, 4)
    print(query_index(p, t, index))