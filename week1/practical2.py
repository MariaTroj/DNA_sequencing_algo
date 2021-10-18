from practical1 import readFastq
import matplotlib.pyplot as plt


def findGCByPos(reads):
    ''' Find the GC ratio at each position in the read '''
    # Keep track of the number of G/C bases and the total number of bases at each position
    gc = [0] * 100
    totals = [0] * 100
    for read in reads:
        for i in range(len(read)):
            if read[i] == 'C' or read[i] == 'G':
                gc[i] += 1
            totals[i] += 1
    # Divide G/C counts by total counts to get the average at each position
    for i in range(len(gc)):
        if totals[i] > 0:
            gc[i] /= float(totals[i])
    return gc


if __name__ == "__main__":

    seqs, quals = readFastq('..\\SRR835775_1.first1000.fastq')

    gc = findGCByPos(seqs)
    # Plot Gs and Cs content
    plot2 = plt.figure(2)
    plt.plot(range(len(gc)), gc)
    plt.title("GC content")

    plt.show()