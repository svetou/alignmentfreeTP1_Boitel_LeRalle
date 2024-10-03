from TP.loading import load_directory
from TP.kmers import stream_kmers, kmer2str, encode_kmer
#from pysam import FastaFile


def jaccard(fileA, fileB, k):
    #open file : 

    print(fileA)
    

    j = 0
    
    return j



if __name__ == "__main__":
    print(kmer2str(encode_kmer("ATATCGCG", 7),7))
    # print("Computation of Jaccard similarity between files")

    # # Load all the files in a dictionary
    # files = load_directory("data")
    # k = 21
    
    # print("Computing Jaccard similarity for all pairs of samples")
    # filenames = list(files.keys())
    # for i in range(len(files)):
    #     for j in range(i+1, len(files)):
            
    #         # --- Complete here ---

    #         j = jaccard(files[filenames[i]], files[filenames[j]], k)
    #         print(filenames[i], filenames[j], j)
