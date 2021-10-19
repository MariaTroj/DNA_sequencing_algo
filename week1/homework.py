from practical1 import phred_33_to_q
from practical3_naive_algo import reverse_complement
from file_reading import *
import matplotlib.pyplot as plt


def naive_with_reverse_complement(pattern, text):
    occurrences = []
    reverse_pattern = reverse_complement(pattern)
    for i in range(len(text) - len(pattern) + 1):
        match = True
        reverse_match = not (pattern == reverse_pattern)
        for j in range(len(pattern)):
            if text[i + j] != pattern[j]:
                match = False
            if text[i + j] != reverse_pattern[j]:
                reverse_match = False
            if not (match or reverse_match):
                break
        if match:
            occurrences.append(i)
        if reverse_match:
            occurrences.append(i)
    return occurrences


def naive_with_mismatches(pattern, text, n_of_mismatches = 2):
    occurrences = []
    for i in range(len(text) - len(pattern) + 1):
        mismaches = 0
        for j in range(len(pattern)):
            if text[i + j] != pattern[j]:
                mismaches += 1
                if mismaches > n_of_mismatches:
                    break
        if mismaches <= n_of_mismatches:
            occurrences.append(i)
    return occurrences


if __name__ == "__main__":
    sequence = 'AGGAGGTT'
    reference_genome = read_genome('..\\lambda_virus.fa')

    occurrences = naive_with_mismatches(sequence, reference_genome)
    print(occurrences)

    seqs, quals = read_seq_and_qual("..\\ERR037900_1.first1000.fastq")

    average_quality = [0]*len(quals[0])

    for qual in quals:
        for i, phred in enumerate(qual):
            average_quality[i] += phred_33_to_q(phred)

    plot1 = plt.figure(1)
    plt.title("Quality Scores")
    plt.plot(range(len(average_quality)), average_quality)
    plt.show()