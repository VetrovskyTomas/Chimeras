from multiprocessing import Pool
from itertools import chain
from collections import deque, Counter
import re
import time
import pickle
import gzip
from Bio import SeqIO
import os
from statistics import median
import sys
import numpy as np

def extract_kmers(seq, k_size): # this is somehow faster than this: [seq[r:r+k_size] for r in range(len(seq)-k_size+1)]
    """
    Extracts k-mers (substrings of length k) from a sequence with their starting indices.
    Args:
        seq (str): Input sequence.
        k_size (int): Length of each k-mer.
    Returns:
        list of tuple: List of (k-mer, index) or empty list if seq is shorter than k_size.
    Example:
        >>> extract_kmers("ACGTACGT", 3)
        [('ACG', 0), ('CGT', 1), ('GTA', 2), ('TAC', 3), ('ACG', 4), ('CGT', 5)]
    """
    if len(seq) < k_size:
        return []  # Skip short sequences
    kmers = []
    window = deque(seq[:k_size], maxlen=k_size)  # Initialize rolling window
    
    kmers.append("".join(window))
    for i in range(k_size, len(seq)):  # Slide the window
        window.append(seq[i])  # Fast O(1) update
        kmers.append("".join(window))
    #print("kmers","len:",len(kmers))
    return [(k,i) for k,i in zip(kmers,list(range(len(kmers))))]


def process_chunk(para_input):
    k, v, k_size, out_dir,compress_or_return = para_input
    #print(len(v))
    #if len(v) < 2: return
    extended_list_of_all_kmers_with_index = []
    for seq in v: extended_list_of_all_kmers_with_index.extend(extract_kmers(seq,k_size))
    #print(extended_list_of_all_kmers_with_index)
    #print(list_of_kmers)
    all_kmers_set = set([k[0] for k in extended_list_of_all_kmers_with_index])
    all_kmers_index_dic = {k:[] for k in all_kmers_set}
    for kmer, i in extended_list_of_all_kmers_with_index:
        all_kmers_index_dic[kmer].append(i)
    #all_kmers_index_dic = {kmer: (median(indices),len(indices)) for kmer, indices in all_kmers_index_dic.items()}
    #print({k:all_kmers_index_dic})
    # Convert all_kmers_index_dic to a memory-efficient structure
    kmer_array = np.array(list(all_kmers_index_dic.keys()), dtype=f'U{k_size}')  # Array of k-mers
    index_array = np.array([np.array(indices, dtype=np.uint16) for indices in all_kmers_index_dic.values()], dtype=object)  # Array of lists of indices

    # Save the numpy arrays
    if compress_or_return == "c":
        np.savez_compressed(f'{out_dir}{k}.npz', kmer_array=kmer_array, index_array=index_array)
    if compress_or_return == "r":
        return kmer_array, index_array
    #with gzip.open(f'{out_dir}{k}.pickle.gzip', 'wb') as f:
    #    pickle.dump({k:all_kmers_index_dic}, f, protocol=pickle.HIGHEST_PROTOCOL)
    
def get_kmers_para_pickle(ref_seqs,k_size=50,threads = 1,name_patter=None,out_dir=None,compress_or_return = "c"): # dont need par_of_fasta_id
    """
    Takes list of biopython reads from UNITE database and saves kmers present in a taxonomic group [genus, family, class, order] in .../all{k_mer_size}.pickle.gzip
    compress_or_return = "c" for compressing and saving each file separetly, "r" for returning the kmers as a dictionary and saving them at once
    """
    start = time.time()
    ref_seqs = [i for i in SeqIO.parse(ref_seqs,"fasta")]
    pat = None
    if not name_patter: print("No name pattern!") ; return
    elif "genus" in name_patter.lower(): pat = "g"
    elif "family" in name_patter.lower(): pat = "f"
    elif "class" in name_patter.lower(): pat = "c"
    elif "order" in name_patter.lower(): pat = "o"
    if pat: name_patter = f"({pat}__).+?(;)"
    
    if not out_dir : print("No output directory!") ; return
    if out_dir[-1] != "/": out_dir = out_dir+"/"
    out_dir = f"{out_dir}all{k_size}mers/"
    os.makedirs(out_dir, exist_ok=True)
    
    ref_seqs = [(re.search(name_patter,i.id).group().strip(";"),i.seq) for i in ref_seqs] # makes list of tuple(family-name , sequence) from list of SeqRecords
    names = set([i[0] for i in ref_seqs])
    print(f"{len(ref_seqs)} of sequences, {len(names)} of unique taxa, K size: {k_size}")
    dic_seqs = {k:[] for k in names}
    for i,s in ref_seqs: dic_seqs[i].append(s)
    del(ref_seqs)
    for k in dic_seqs.keys():
        if "unidentified" in k.lower():
            dic_seqs.pop(k)
            break
    lenka = len(names)
    if compress_or_return == "c":
        with Pool(threads) as pool:
            para_input = [(k,v,k_size,out_dir,"c") for k,v in dic_seqs.items()]
            for i, result in enumerate(pool.imap_unordered(process_chunk, para_input)):
                if i % 10 == 0: print(f"{i}/{lenka}      ", end="\r")
                #if i ==3:break
    if compress_or_return == "r":
            with Pool(threads) as pool:
                para_input = [(k,v,k_size,out_dir,"r") for k,v in dic_seqs.items()]
                results = []
                for i, result in enumerate(pool.imap_unordered(process_chunk, para_input)):
                    if i % 500 == 0: print(f"{i}/{lenka}      ", end="\r")
                    results.append(result)
                all_kmers_index_dic = {k: v for k, v in zip(names, results)}
            with open(f"{out_dir}all_kmers_index_dic.pickle", 'wb') as f:
                pickle.dump(all_kmers_index_dic, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"Execution Time: {time.time() - start}")

