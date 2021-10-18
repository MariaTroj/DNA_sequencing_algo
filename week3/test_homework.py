from week3.homework import local_alignment, overlap_all_pairs


def test_local_alignment():
    pattern = 'GCGTATGC'
    text = 'TATTGGCTATACGGTT'
    assert local_alignment(pattern, text) == 2


def test_overlap_all_pairs():
    reads = ['ABCDEFG', 'EFGHIJ', 'HIJABC']
    assert overlap_all_pairs(reads, 3) == [('ABCDEFG', 'EFGHIJ'), ('EFGHIJ', 'HIJABC'), ('HIJABC',
                                                                                      'ABCDEFG')]
    assert overlap_all_pairs(reads, 4) == []

    reads = ['CGTACG', 'TACGTA', 'GTACGT', 'ACGTAC', 'GTACGA', 'TACGAT']
    pairs = overlap_all_pairs(reads, 4)
    for pair in [('CGTACG', 'TACGTA'),
                 ('CGTACG', 'GTACGT'),
                 ('CGTACG', 'GTACGA'),
                 ('CGTACG', 'TACGAT'),
                 ('TACGTA', 'ACGTAC'),
                 ('TACGTA', 'CGTACG'),
                 ('GTACGT', 'TACGTA'),
                 ('GTACGT', 'ACGTAC'),
                 ('ACGTAC', 'GTACGA'),
                 ('ACGTAC', 'GTACGT'),
                 ('ACGTAC', 'CGTACG'),
                 ('GTACGA', 'TACGAT')]:
        assert pair in pairs
        pairs.remove(pair)
    assert len(pairs) == 0
    pairs = overlap_all_pairs(reads, 5)
    for pair in [('CGTACG', 'GTACGT'),
                 ('CGTACG', 'GTACGA'),
                 ('TACGTA', 'ACGTAC'),
                 ('GTACGT', 'TACGTA'),
                 ('ACGTAC', 'CGTACG'),
                 ('GTACGA', 'TACGAT')]:
        assert pair in pairs
        pairs.remove(pair)
    assert len(pairs) == 0

