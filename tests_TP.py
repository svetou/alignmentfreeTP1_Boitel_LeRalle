from TP.kmers import kmer2str, encode_kmer_forward_sense, encode_kmer_reverse_complement, encode_kmer
from TP.kmers import stream_kmers
from TP.__main__ import jaccard # in a big project we would move tested functions out of __main__
# to avoid importing from it, but in our case this seemed simpler than changing the TP file structure.

# Single kmer encoding
assert kmer2str(encode_kmer_forward_sense("ATATCGCAA", 7),7)      == "ATATCGC"
assert kmer2str(encode_kmer_reverse_complement("ATATCGCAA", 7),7) == "GCGATAT"
assert kmer2str(encode_kmer("ATATCGCAA", 7),7)                    == "ATATCGC"
assert kmer2str(encode_kmer("GCGATATAA", 7),7)                    == "ATATCGC"


# Kmer stream encoding
seq = "ATATCGCAAAAAGAGAAATTCAGAGAGACT"
for k in range(1, len(seq) + 1):
    encoded_stream = list(stream_kmers(seq, k))
    encoded_one_by_one = [encode_kmer(seq[i:i+k], k) for i in range(len(seq) + 1 - k)]
    assert encoded_stream == encoded_one_by_one # Python list equality => all elements equal
    
assert(list(stream_kmers("ATC", 42)) == [])

# Jaccard distance testing
## Sanity checks
seq1 = "ATATCGCAAAAAGAGAAATTCAGAGAGACTATATCGCAAAAAGAGAAATTCAGAGAGACT"
assert jaccard([seq1], [seq1], 5) == 1.0
assert jaccard([seq1], [""], 5) == 0.0
seq2 = "GGGAGGGTGGCCGAGAAAAAAAAGTGTGAGC"
assert jaccard([seq1], [seq2], 5) == jaccard([seq2], [seq1], 5)
assert jaccard([seq1, seq1], [seq2, seq2], 5) == jaccard([seq1], [seq2], 5)
assert jaccard([seq1], [seq2], len(seq2)) == 0.0

## Simple manual tests
assert jaccard(["AAAT"], ["AAGGG"], 2) == 1 / 6
assert jaccard(["TTTA"], ["AAGGG"], 2) == 1 / 6

print("All tests passed.")
