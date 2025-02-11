from multiprocessing import Pool
from itertools import chain
from collections import deque, Counter
import re
import time

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


def process_chunk(para_input):
    chunk_dict, k_size, kmer_count_threshold = para_input
    res_dic = {k:[] for k in chunk_dict.keys()}
    for k,v in chunk_dict.items():
        kmers = list(chain.from_iterable([extract_kmers(seq,k_size) for seq in v]))
        if kmer_count_threshold > 1: com_all_kmers = set([kmer for kmer, count in Counter(kmers).items() if count > kmer_count_threshold])
        elif kmer_count_threshold == 1: com_all_kmers = set(kmers)
        res_dic[k] = com_all_kmers
        
    return res_dic

def get_kmers_para2(ref_seqs,k_size=50,kmer_count_threshold=10,threads = 1,name_patter = "(f__).+?(;)"): # dont need par_of_fasta_id
    """
    Takes list of biopython reads from UNITE database and returns dictionary of kmers present in each family
    """
    start = time.time()
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
            
    keys = list(dic_seqs.keys())
    chunk_size = len(keys) // threads  
    key_chunks = [keys[i:i + chunk_size] for i in range(0, len(keys), chunk_size)]

    if len(keys) > chunk_size:
        key_chunks[-1].extend(keys[chunk_size*threads:])

    para_input = [({k:dic_seqs[k] for k in k_c},k_size,kmer_count_threshold) for k_c in key_chunks]
    with Pool(threads) as pool:
        results = pool.map(process_chunk, para_input)

    final_result = {}
    for r in results:
        final_result.update(r)
    print(f"Execution Time: {time.time() - start}")
    return final_result