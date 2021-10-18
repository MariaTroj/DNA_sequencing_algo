from homework import *


def test_naive_with_counts():
    pattern = 'word'
    text = 'there would have been a time for such a word'
    occurrences, num_alignments, num_character_comparisons = naive_with_counts(pattern, text)
    print(occurrences, num_alignments, num_character_comparisons)
    assert (occurrences, num_alignments, num_character_comparisons) == ([40], 41, 46)

    pattern = 'needle'
    text = 'needle need noodle needle'
    occurrences, num_alignments, num_character_comparisons = naive_with_counts(pattern, text)
    assert (occurrences, num_alignments, num_character_comparisons) == ([0, 19], 20, 35)


def test_boyer_moore_with_counts():
    pattern = 'word'
    text = 'there would have been a time for such a word'
    lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz '
    p_bm = BoyerMoore(pattern, lowercase_alphabet)
    occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(pattern,
                                                                                     p_bm, text)
    assert (occurrences, num_alignments, num_character_comparisons) == ([40], 12, 15)

    pattern = 'needle'
    text = 'needle need noodle needle'
    p_bm = BoyerMoore(pattern, lowercase_alphabet)
    occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(pattern,
                                                                                     p_bm, text)
    assert (occurrences, num_alignments, num_character_comparisons) == ([0, 19], 5, 18)


def test_query_subseq():
    t = 'to-morrow and to-morrow and to-morrow creeps in this petty pace'
    p = 'to-morrow and to-morrow '
    subseq_ind = SubseqIndex(t, 8, 3)
    occurrences, num_index_hits = query_subseq(p, t, subseq_ind, 2)
    assert occurrences == [0, 14]
    assert num_index_hits == 6
