# sorted inputs ?? will it help?? # nestarea se o "unidentified"
from itertools import chain
from collections import deque, Counter
import re
import time

def get_kmers(ref_seqs,k_size=50,kmer_count_threshold=10,name_patter = "(f__).+?(;)"): # dont need par_of_fasta_id
    """
    Takes list of biopython reads from UNITE database and returns dictionary of kmers present in each family
    """
    ref_seqs = [(re.search(name_patter,i.id).group().strip(";"),i.seq) for i in ref_seqs] # makes list of tuple(family-name , sequence) from list of SeqRecords
    names = set([i[0] for i in ref_seqs])
    print(f"{len(ref_seqs)} of sequences, {len(names)} of unique taxa")
    dic_seqs = {k:[] for k in names}
    for i,s in ref_seqs: dic_seqs[i].append(s)
    del(ref_seqs)
    for k in dic_seqs.keys():
        if "unidentified" in k.lower():
            dic_seqs.pop(k)
            break
    dic_kmers = {k:[] for k in names}

    def extract_kmers(seq, k_size): # this is somehow faster than this: [seq[r:r+k_size] for r in range(len(seq)-k_size+1)]
        if len(seq) < k_size:
            return []  # Skip short sequences
        kmers = []
        window = deque(seq[:k_size], maxlen=k_size)  # Initialize rolling window
        kmers.append("".join(window))
        for i in range(k_size, len(seq)):  # Slide the window
            window.append(seq[i])  # Fast O(1) update
            kmers.append("".join(window))
        return kmers

    for k,v in dic_seqs.items():
        #kmers = list(chain.from_iterable([[seq[r:r+k_size] for r in range(len(seq)-k_size+1)] for seq in v]))
        
        kmers = list(chain.from_iterable([extract_kmers(seq,k_size) for seq in v]))
        #print(kmers)
        if kmer_count_threshold > 1: com_all_kmers = set([kmer for kmer, count in Counter(kmers).items() if count > kmer_count_threshold])
        elif kmer_count_threshold == 1: com_all_kmers = set(kmers)
        dic_kmers[k] = com_all_kmers
    return dic_kmers