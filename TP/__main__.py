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
    distances = { i: { i: 1.0 } for i in range(len(filenames)) }

    for i in range(len(files)):
        for j in range(i+1, len(files)):            
            distances[i][j] = jaccard(files[filenames[i]], files[filenames[j]], k)
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
