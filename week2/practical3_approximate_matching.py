from practical1_boyer_moore import *


def approximate_match(pattern: str, text: str, n: int):
    segment_length = int(round(len(pattern) / (n + 1)))
    all_matches = set()
    for i in range(n+1):
        start = i*segment_length
        end = min((i+1) * segment_length, len(pattern))
        p_bm = BoyerMoore(pattern[start:end], alphabet='ACGT')
        matches = boyer_moore(pattern[start:end], p_bm, text)
        # Extend matching segments to see if whole p matches
        for m in matches:
            if m < start or m-start+len(pattern) > len(text):
                continue
            mismatches = 0
            for j in range(0, start):
                if not pattern[j] == text[m - start + j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            for j in range(end, len(pattern)):
                if not pattern[j] == text[m - start + j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            if mismatches <= n:
                all_matches.add(m - start)
    return list(all_matches)


if __name__ == '__main__':
    pattern = 'AACTTG'
    text = 'CACTTAATTTG'
    print(approximate_match(pattern, text, 2))
