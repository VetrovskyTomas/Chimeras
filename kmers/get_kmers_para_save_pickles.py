from multiprocessing import Pool
from itertools import chain
from collections import deque, Counter
import re
import time
import pickle
import gzip
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
    k, v, k_size, out_dir = para_input
    kmers = list(chain.from_iterable([extract_kmers(seq,k_size) for seq in v]))
    with gzip.open(f'{out_dir}{k}.pickle.gzip', 'wb') as f:
        pickle.dump({k:kmers}, f, protocol=pickle.HIGHEST_PROTOCOL)

def get_kmers_para_pickle(ref_seqs,k_size=50,threads = 1,name_patter=None,out_dir=None): # dont need par_of_fasta_id
    """
    Takes list of biopython reads from UNITE database and saves kmers present in a taxonomic group [genus, family, class, order] in .../all{k_mer_size}.pickle.gzip
    """
    start = time.time()
    ref_seqs = [i for i in SeqIO.parse(ref_seqs,"fasta")]

    if not name_patter: print("No name pattern!") ; return
    elif "g" in name_patter.lower(): pat = "g"
    elif "f" in name_patter.lower(): pat = "f"
    elif "c" in name_patter.lower(): pat = "c"
    elif "o" in name_patter.lower(): pat = "o"
    name_patter = f"({pat}__).+?(;)"
    
    if not out_dir : print("No output directory!") ; return
    if out_dir[-1] != "/": out_dir = out_dir+"/"
    out_dir = f"{out_dir}all{k_size}mers/"
    os.makedirs(out_dir, exist_ok=True)
    
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
    lenka = len(names)
    with Pool(threads) as pool:
        para_input = [(k,v,k_size,out_dir) for k,v in dic_seqs.items()]
        for i, result in enumerate(pool.imap_unordered(process_chunk, para_input)):
            if i % 10 == 0: print(f"{i}/{lenka}      ", end="\r")
                
    print(f"Execution Time: {time.time() - start}")