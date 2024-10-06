
def kmer2str(val, k):
    """ Transform a kmer integer into a its string representation
    :param int val: An integer representation of a kmer
    :param int k: The number of nucleotides involved into the kmer.
    :return str: The kmer string formatted
    """
    letters = ['A', 'C', 'T', 'G']
    str_val = []
    for _ in range(k):
        str_val.append(letters[val & 0b11])
        val >>= 2

    str_val.reverse()
    return "".join(str_val)

encode_nuc = { 'A': 0, 'C': 1, 'T': 2, 'G': 3 }

def encode_kmer_forward_sense(seq, k):
    kmer = 0
    for nuc in seq[0:k]:
        kmer <<= 2
        kmer |= encode_nuc[nuc]
    return kmer

complement_nuc = { 'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C' }

def encode_kmer_reverse_complement(seq, k):
    kmer = 0
    for nuc in seq[k-1::-1]: # first k nucleotides of seq in reverse order
        kmer <<= 2
        kmer |= encode_nuc[complement_nuc[nuc]]
    return kmer

def encode_kmer(seq, k):
    return min(encode_kmer_forward_sense(seq, k), encode_kmer_reverse_complement(seq, k))

def stream_kmers(text, k):
    mask = ( 1 << 2 * (k - 1)) - 1

    pass
