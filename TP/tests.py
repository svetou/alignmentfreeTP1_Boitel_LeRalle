from kmers import kmer2str, encode_kmer_forward_sense, encode_kmer_reverse_complement, encode_kmer
from kmers import stream_kmers

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


print("All tests passed.")
