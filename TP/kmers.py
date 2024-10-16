import heapq
import hashlib

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

def stream_kmers(seq, k):
    if len(seq) < k:
        return # don't yield anything for sequences shorter than k
    mask_forward = ( 1 << (2 * (k - 1))) - 1
    # There can be inversions within the sequence so we cannot just compute the strand sense once,
    # but we need to keep track of forward and reverse kmers to yield the correct canonical ones
    kmer_forward = encode_kmer_forward_sense(seq, k)
    kmer_reverse = encode_kmer_reverse_complement(seq, k)
    for nuc in seq[k:]:
        yield min(kmer_forward, kmer_reverse)
        kmer_forward &= mask_forward
        kmer_forward <<= 2
        kmer_forward |= encode_nuc[nuc]

        # no need to mask kmer_reverse as extra bits on the right will be discarded
        kmer_reverse >>= 2
        kmer_reverse |= (encode_nuc[complement_nuc[nuc]] << (2 * (k - 1)))
    yield min(kmer_forward, kmer_reverse)

def hash_kmers(kmer_list):
 
    hashed_kmers = []
    for kmer in kmer_list:
        kmer_str = kmer2str(kmer, len(kmer_list[0]))  # Convert kmer to string for hashing
        hashed_value = int(hashlib.md5(kmer_str.encode()).hexdigest(), 16)
        hashed_kmers.append(hashed_value)
    return hashed_kmers

def filter_smallests(seq, k, s):

    # Generate all k-mers with our function stram_kmers
    kmers = list(stream_kmers(seq, k))
    hashed_kmers = hash_kmers(kmers)
    heap = []
    
    for hashed_kmer in hashed_kmers:
        if len(heap) < s:
            heapq.heappush(heap, -hashed_kmer)  # Negative for max-heap behavior
        else:
            max_heap = -heap[0]
            if hashed_kmer < max_heap:
                heapq.heapreplace(heap, -hashed_kmer)
    
    # Return the sorted sketch
    return sorted(-x for x in heap)
