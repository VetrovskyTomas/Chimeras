# detekce kmeru 

import time
from collections import Counter
from multiprocessing import Pool

def process_chunk(para_input):
    df_col,seq = para_input
    hits = []
    #print(df_col);return
    for n,row in enumerate(df_col):
        for kmer in row:
            if kmer in seq.seq:
                hits.append(n)
    if len(set(hits)) > 1:
        fung_counts = [(fung,count) for fung,count in Counter(hits).items() if count > 0] # tahle 0 jde ladit
        if len(fung_counts) > 1:
            #print(seq.id,fung_counts)
            return (seq.id,fung_counts) # make dict insted of tuple? faster?
def detect_chim(my_seqs = None, kmer_df = None, threads = 4, columns = None):
    """
    Takes sequences you want to detect chimeras in and dataframe of unique Kmers per fungal taxa
    If no columns are provided, takes all columns in the dataframe
    Returns dictionary with { column : suspected chimeras = (seq.id , kmers from different taxa found in this sequence) }
    The dataframe has fungal taxonomic groups as rows (no names needed) and columns with various Kmer sizes - each cell in this DF will contain set of kmers unique for this taxa or an empty set
    """
    start = time.time()

    if not columns: columns = df.columns 
    chim_per_col = {k:[] for k in columns}
    for col in columns:
        chims = []
        with Pool(threads) as pool:
            para_input = [(df[col],seq) for seq in my_seqs]
            for i, result in enumerate(pool.imap_unordered(process_chunk, para_input)):
                if i % 500 == 0: print(f"{i}       ", end="\r")
                if result:chims.append(result)
        chim_per_col[col] = chims
    print(time.time() - start)
    return chim_per_col