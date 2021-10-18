import collections
import random
import matplotlib.pyplot as plt


# import wget
# url = 'http://d28rh4a8wq0iu5.cloudfront.net/ads1/data/SRR835775_1.first1000.fastq'
# url = 'http://d28rh4a8wq0iu5.cloudfront.net/ads1/data/phix.fa'
# filename = wget.download(url)

genome = [random.choice('ACGT') for i in range(0, 1000)]
base_counts = collections.Counter(genome)

def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline() # skip name line
            seq = fh.readline().rstrip() # read base sequence
            fh.readline() # skip placeholder line
            qual = fh.readline().rstrip() #base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

def phred33ToQ(qual):
    return ord(qual) - 33

def createHist(qualityStrings):
    # Create a histogram of quality scores
    hist = [0]*50
    for read in qualityStrings:
        for phred in read:
            q = phred33ToQ(phred)
            hist[q] += 1
    return hist


if __name__ == 'main':

    seqs, quals = readFastq('..\\SRR835775_1.first1000.fastq')

    count = collections.Counter()
    for seq in seqs:
        count.update(seq)
    print(count)

    h = createHist(quals)
    # Plot the histogram
    plot1 = plt.figure(1)
    plt.title("Quality Scores")
    plt.plot(range(len(h)), h)

    plt.show()

