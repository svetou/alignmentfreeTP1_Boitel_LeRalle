from TP.loading import load_directory
from TP.kmers import stream_kmers, kmer2str, encode_kmer, hash_kmers, filter_smallests
from sys import stderr
import heapq
import hashlib
import itertools


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
    
def jaccard2(sketch1, sketch2):
    i, j = 0, 0
    count_union = 0
    count_intersection = 0

    while i < len(sketch1) and j < len(sketch2):
        if sketch1[i] == sketch2[j]:
            count_intersection += 1
            count_union += 1
            i += 1
            j += 1
        elif sketch1[i] < sketch2[j]:
            count_union += 1
            i += 1
        else:
            count_union += 1
            j += 1

    # Add remaining elements from either sketch
    count_union += (len(sketch1) - i) + (len(sketch2) - j)

    return count_intersection / count_union


if __name__ == "__main__":

    print("Computation of Jaccard similarity between files")

    # Load all the files in a dictionary
    files = load_directory("data")
    k = 21
    s = 10000 # number of elements in a sketch

    filenames = list(files.keys())

    print("Precomputing sketches.") # To save time as we use the same file[i] in all i,j pairs for comparison
    sketches = []
    for i in range(len(files)):
        kmers_it = itertools.chain.from_iterable(stream_kmers(seq, k) for seq in files[filenames[i]])
        sketches.append(filter_smallests(kmers_it, s))
        print(f"{filenames[i]} sketch computed.")

    print("Computing Jaccard similarity for all pairs of samples")
    distances = { i: { i: 1.0 } for i in range(len(filenames)) }

    for i in range(len(files)):
        for j in range(i+1, len(files)):            
            distances[i][j] = jaccard2(sketches[i], sketches[j])
            print(filenames[i], filenames[j], distances[i][j])

    # Generate the Markdown similarity matrix for the report
    print("Similarity matrix:")
    filenames_no_ext = [ '.'.join(fn.split('.')[:-1]) for fn in filenames ]
    max_col_length = max(len(fn) for fn in filenames_no_ext)
    print("|  " + " " * max_col_length + "  |  " + "  |  ".join(filenames_no_ext) + "  |  ") # column names
    print("|  " + ("-" * max_col_length + "  |  ") * (len(files) + 1)) # markdown ---- line
    for i in range(len(files)):
        print("|  " + filenames_no_ext[i], end="  | ") # row names
        for j in range(len(files)):
            print(f"{distances[min(i,j)][max(i,j)]: {max_col_length-4}.{max_col_length-2}f}", end="  | ")
        print("")

    
