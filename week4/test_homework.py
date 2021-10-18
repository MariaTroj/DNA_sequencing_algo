from week4.homework import scs, scs_list


def test_scs_list():
    strings = ['ABC', 'BCA', 'CAB']
    # Returns just one shortest superstring
    assert scs(strings) == 'ABCAB'
    # Returns list of all superstrings that are tied for shorest
    assert sorted(scs_list(strings)) == sorted(['ABCAB', 'BCABC', 'CABCA'])

    strings = ['GAT', 'TAG', 'TCG', 'TGC', 'AAT', 'ATA']
    # Returns just one shortest superstring
    assert scs(strings) == 'TCGATGCAATAG'
    # Returns list of all superstrings that are tied for shorest
    assert sorted(scs_list(strings)) == ['AATAGATCGTGC', 'AATAGATGCTCG', 'AATAGTCGATGC',
                                         'AATCGATAGTGC', 'AATGCTCGATAG', 'TCGAATAGATGC',
                                         'TCGATAGAATGC', 'TCGATGCAATAG', 'TGCAATAGATCG',
                                         'TGCAATCGATAG']
