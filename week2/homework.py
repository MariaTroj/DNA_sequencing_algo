from bm_preproc import BoyerMoore
from kmer_index import *

def boyer_moore_with_counts(pattern: str, p_bm: BoyerMoore, text: str) -> tuple[list, int, int]:
    """ Do Boyer-Moore matching
    count and return (a) the number of character comparisons performed and (b) the number of alignments tried"""
    i = 0
    occurrences: list = []
    num_alignments: int = 0
    num_character_comparisons: int = 0
    while i < len(text) - len(pattern) + 1:
        shift = 1
        mismatched = False
        num_alignments += 1
        for j in range(len(pattern)-1, -1, -1):
            num_character_comparisons += 1
            if pattern[j] != text[i+j]:
                skip_bc = p_bm.bad_character_rule(j, text[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences, num_alignments, num_character_comparisons


def naive_with_counts(pattern, text) -> tuple[list, int, int]:
    occurrences: list = []
    num_alignments: int = 0
    num_character_comparisons: int = 0
    for i in range(len(text) - len(pattern) + 1):
        match = True
        num_alignments += 1
        for j in range(len(pattern)):
            num_character_comparisons += 1
            if text[i + j] != pattern[j]:
                match = False
                break
        if match:
            occurrences.append(i)
    return occurrences, num_alignments, num_character_comparisons


def index_assisted_approximate_matching(pattern: str, k_mers_index: Index, text: str,
                                        max_n_mismatches:
int):
    """
    Function that, given a length-24 pattern P and given an Index object
    built on 8-mers, finds all approximate occurrences of P within T with up to 2 mismatches.
    Insertions and deletions are not allowed. Don't consider any reverse complements """
    segment_length = int(round(len(pattern) / (max_n_mismatches+1)))
    index_hits = 0
    all_matches = set()
    for i in range(max_n_mismatches+1):
        start = i*segment_length
        end = min((i+1)*segment_length, len(text))
        matches = k_mers_index.query(pattern[start:end])
        index_hits += len(matches)
        # Extend matching segments to see if whole p matches
        for m in matches:
            if m < start or m-start+len(pattern) > len(text):
                continue
            mismatches = 0
            for j in range(0, start):
                if not pattern[j] == text[m-start+j]:
                    mismatches += 1
                    if mismatches > max_n_mismatches:
                        break
            for j in range(end, len(pattern)):
                if not pattern[j] == text[m-start+j]:
                    mismatches += 1
                    if mismatches > max_n_mismatches:
                        break
            if mismatches <= max_n_mismatches:
                all_matches.add(m - start)
    return list(all_matches), index_hits


def query_subseq(pattern: str, text: str, subseq_ind: SubseqIndex, max_n_mismatches: int):
    ival = subseq_ind.ival
    all_matches = set()
    index_hits = 0
    for i in range(ival):
        matches = subseq_ind.query(pattern[i:])
        index_hits += len(matches)
        # Extend matching segments to see if whole p matches
        for m in matches:
            if m < i or m - i + len(pattern) > len(text):
                continue
            mismatches = 0
            for j in range(ival):
                if j == i: continue
                for k in range(j, len(pattern), ival):
                    if not pattern[k] == text[m - i + k]:
                        mismatches += 1
                        if mismatches > max_n_mismatches:
                            break
            if mismatches <= max_n_mismatches:
                all_matches.add(m - i)
    return list(all_matches), index_hits


if __name__ == "__main__":
    pattern = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
    text = ''
    alphabet = 'ACGT '
    p_bm = BoyerMoore(pattern, alphabet)

    with open("..\\chr1.GRCh38.excerpt.fasta", 'r') as file:
        for line in file:
            if line.startswith('>'):
                continue
            text += line.strip()

    # Question 1, 2, 3
    occurrences_bm, num_alignments_bm, num_character_comparisons_bm = boyer_moore_with_counts(
        pattern, p_bm, text)
    occurrences_naive, num_alignments_naive, num_character_comparisons_naive = naive_with_counts(
        pattern, text)
    print(len(occurrences_naive), num_alignments_naive, num_character_comparisons_naive)
    print(len(occurrences_bm), num_alignments_bm, num_character_comparisons_bm)

    # Question 4, 5
    pattern = 'GGCGCGGTGGCTCACGCCTGTAAT'
    k_mers_index = Index(text, 8)
    approximate_occurences, index_hits = index_assisted_approximate_matching(pattern,
                                                                              k_mers_index, text, 2)
    print(index_hits, len(approximate_occurences), approximate_occurences)

    total_num_alignments_bm = 0
    for i in range(0, 3):
        p_bm = BoyerMoore(pattern[i*8:(i+1)*8], alphabet)
        occurrences_bm, _, _ = boyer_moore_with_counts(pattern[i*8:(i+1)*8], p_bm, text)
        total_num_alignments_bm += len(occurrences_bm)
    print(total_num_alignments_bm)

    # Question 6
    pattern = 'GGCGCGGTGGCTCACGCCTGTAAT'
    subseq_ind = SubseqIndex(text, 8, 3)
    occurrences, num_index_hits = query_subseq(pattern, text, subseq_ind, 2)
    print(occurrences, num_index_hits)