import collections
import random
import matplotlib.pyplot as plt
from file_reading import read_seq_and_qual

# import wget
# url = 'http://d28rh4a8wq0iu5.cloudfront.net/ads1/data/SRR835775_1.first1000.fastq'
# url = 'http://d28rh4a8wq0iu5.cloudfront.net/ads1/data/phix.fa'
# filename = wget.download(url)

genome = [random.choice('ACGT') for i in range(0, 1000)]
base_counts = collections.Counter(genome)


def phred_33_to_q(qual):
    return ord(qual) - 33

def create_hist(quality_strings):
    # Create a histogram of quality scores
    hist = [0]*50
    for read in quality_strings:
        for phred in read:
            q = phred_33_to_q(phred)
            hist[q] += 1
    return hist


if __name__ == '__main__':

    seqs, quals = read_seq_and_qual('..\\SRR835775_1.first1000.fastq')

    count = collections.Counter()
    for seq in seqs:
        count.update(seq)
    print(count)

    h = create_hist(quals)
    # Plot the histogram
    plot1 = plt.figure(1)
    plt.title("Quality Scores")
    plt.plot(range(len(h)), h)

    plt.show()

