'''
Functions:
 - z_array
 - n_array
 - big_l_prime_array
 - smal_l_prime-array
 - good_suffix_table
 - dense_bad_char_tab
are preprocessing pattern and preparing table used by teh bad character rule and good suffix rule
'''
def z_array(s):
    """ Use Z algorithm (Gusfield theorem 1.4.1) to preprocess s """
    assert len(s) > 1
    z = [len(s)] + [0] * (len(s) - 1)
    # Initial comparison of s[1:] with prefix
    for i in range(1, len(s)):
        if s[i] == s[i - 1]:
            z[1] += 1
        else:
            break
    r, l = 0, 0
    if z[1] > 0:
        r, l = z[1], 1
    for k in range(2, len(s)):
        assert z[k] == 0
        if k > r:
            # Case 1
            for i in range(k, len(s)):
                if s[i] == s[i - k]:
                    z[k] += 1
                else:
                    break
            r, l = k + z[k] - 1, k
        else:
            # Case 2
            # Calculate length of beta
            nbeta = r - k + 1
            zkp = z[k - l]
            if nbeta > zkp:
                # Case 2a: Zkp wins
                z[k] = zkp
            else:
                # Case 2b: Compare characters just past r
                nmatch = 0
                for i in range(r + 1, len(s)):
                    if s[i] == s[i - k]:
                        nmatch += 1
                    else:
                        break
                l, r = k, r + nmatch
                z[k] = r - k + 1
    return z


def n_array(s):
    """ Compile the N array (Gusfield theorem 2.2.2) from the Z array """
    return z_array(s[::-1])[::-1]


def big_l_prime_array(p, n):
    """ Compile L' array (Gusfield theorem 2.2.2) using p and N array.
        L'[i] = largest index j less than n such that N[j] = |P[i:]| """
    lp = [0] * len(p)
    for j in range(len(p) - 1):
        i = len(p) - n[j]
        if i < len(p):
            lp[i] = j + 1
    return lp


def big_l_array(p, lp):
    """ Compile L array (Gusfield theorem 2.2.2) using p and L' array.
        L[i] = largest index j less than n such that N[j] >= |P[i:]| """
    l = [0] * len(p)
    l[1] = lp[1]
    for i in range(2, len(p)):
        l[i] = max(l[i - 1], lp[i])
    return l


def small_l_prime_array(n):
    """ Compile lp' array (Gusfield theorem 2.2.4) using N array. """
    small_lp = [0] * len(n)
    for i in range(len(n)):
        if n[i] == i + 1:  # prefix matching a suffix
            small_lp[len(n) - i - 1] = i + 1
    for i in range(len(n) - 2, -1, -1):  # "smear" them out to the left
        if small_lp[i] == 0:
            small_lp[i] = small_lp[i + 1]
    return small_lp


def good_suffix_table(p):
    """ Return tables needed to apply good suffix rule. """
    n = n_array(p)
    lp = big_l_prime_array(p, n)
    return lp, big_l_array(p, lp), small_l_prime_array(n)


def dense_bad_char_tab(p, amap):
    """ Given pattern string and list with ordered alphabet characters, create
        and return a dense bad character table.  Table is indexed by offset
        then by character. """
    tab = []
    nxt = [0] * len(amap)
    for i in range(0, len(p)):
        c = p[i]
        assert c in amap
        tab.append(nxt[:])
        nxt[amap[c]] = i + 1
    return tab


class BoyerMoore(object):
    """ Encapsulates pattern and associated Boyer-Moore preprocessing. """

    def __init__(self, pattern: str, alphabet='ACGT'):
        self.pattern = pattern
        self.alphabet = alphabet
        # Create map from alphabet characters to integers
        self.amap = dict((letter, i) for i, letter in enumerate(alphabet))
        # Make bad character rule table
        self.bad_char = dense_bad_char_tab(pattern, self.amap)
        # Create good suffix rule table
        _, self.big_l, self.small_l_prime = good_suffix_table(pattern)

    def bad_character_rule(self, offset: int, character: str) -> int:
        """ Return # skips given by bad character rule at offset """
        assert character in self.amap
        char_index = self.amap[character]
        assert offset > (self.bad_char[offset][char_index] - 1)
        return offset - (self.bad_char[offset][char_index] - 1)

    def good_suffix_rule(self, offset: int) -> int:
        """ Given a mismatch at offset, return amount to shift
            as determined by (weak) good suffix rule. """
        length = len(self.big_l)
        assert offset < length
        if offset == length - 1:
            return 0
        offset += 1  # i points to leftmost matching position of P
        if self.big_l[offset] > 0:
            return length - self.big_l[offset]
        return length - self.small_l_prime[offset]

    def match_skip(self) -> int:
        """ Return amount to shift in case where P matches T """
        return len(self.small_l_prime) - self.small_l_prime[1]


def boyer_moore(pattern: str, p_bm: BoyerMoore, text: str) -> list:
    """ Do Boyer-Moore matching """
    i = 0
    occurrences = []
    while i < len(text) - len(pattern) + 1:
        shift = 1
        mismatched = False
        for j in range(len(pattern)-1, -1, -1):
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
    return occurrences

if __name__ == "__main__":
    t = 'GCTAGCTCTACGAGTCTA'
    p = 'TCTA'

    p_bm = BoyerMoore(p, alphabet='ACGT')
    print(p_bm.bad_char)
    print(p_bm.bad_character_rule(1, 'T'))
    print(boyer_moore(p, p_bm, t))