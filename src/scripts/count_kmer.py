from collections import Counter


def count_kmers(seq, k=6):
    cnt = Counter()
    if len(seq) < k:
        return cnt

    kmer = seq[:k]
    cnt[kmer] += 1
    for i in range(k, len(seq)):
        kmer = kmer[1:] + seq[i]
        cnt[kmer] += 1

    return cnt


def test_no_kmers():
    expected = Counter()
    assert count_kmers("AAA") == expected


def test_small_seq():
    expected = Counter({"AAACGT": 1, "AACGTG": 1})
    assert count_kmers("AAACGTG") == expected


if __name__ == "__main__":
    test_no_kmers()
    test_small_seq()
