from kmers import kmer2str, encode_kmer_forward_sense, encode_kmer_reverse_complement, encode_kmer

# Single kmer encoding
assert kmer2str(encode_kmer_forward_sense("ATATCGCAA", 7),7)      == "ATATCGC"
assert kmer2str(encode_kmer_reverse_complement("ATATCGCAA", 7),7) == "GCGATAT"
assert kmer2str(encode_kmer("ATATCGCAA", 7),7)                    == "ATATCGC"
assert kmer2str(encode_kmer("GCGATATAA", 7),7)                    == "ATATCGC"

print("All tests passed.")
