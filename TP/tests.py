from kmers import kmer2str, encode_kmer_positive_sense, encode_kmer_reverse_complement, encode_kmer

assert kmer2str(encode_kmer_positive_sense("ATATCGC", 7),7)     == "ATATCGC"
assert kmer2str(encode_kmer_reverse_complement("ATATCGC", 7),7) == "GCGATAT"
assert kmer2str(encode_kmer("ATATCGC", 7),7)                    == "ATATCGC"
assert kmer2str(encode_kmer("GCGATAT", 7),7)                    == "ATATCGC"

print("All tests passed.")
