from TP.loading import load_directory
from TP.kmers import stream_kmers, kmer2str, encode_kmer
from sys import stderr


def jaccard(fileA, fileB, k):

    union_count = 0
    intersection_count = 0
    kmer_counts_A = {}

    for seqA in fileA:
        for kmer_A in stream_kmers(seqA, k):
            if kmer_A not in kmer_counts_A:
                kmer_counts_A[kmer_A] = 0
            kmer_counts_A[kmer_A] += 1
            union_count += 1

    for seqB in fileB:
        for kmer_B in stream_kmers(seqB, k):
            if kmer_B in kmer_counts_A and kmer_counts_A[kmer_B] > 0:
                kmer_counts_A[kmer_B] -= 1
                intersection_count += 1
            else:
                union_count += 1

    if union_count == 0: # checking to avoid an error with empty/short sequences
        print("Warning: no k-mers to compare.", file=stderr)
        return 0

    return intersection_count / union_count


if __name__ == "__main__":

    print("Computation of Jaccard similarity between files")

    # Load all the files in a dictionary
    files = load_directory("data")
    k = 21
    
    print("Computing Jaccard similarity for all pairs of samples")
    filenames = list(files.keys())
    for i in range(len(files)):
        for j in range(i+1, len(files)):
            
            # --- Complete here ---
            # TODO what is there to complete here ?

            distance = jaccard(files[filenames[i]], files[filenames[j]], k)
            print(filenames[i], filenames[j], distance)
