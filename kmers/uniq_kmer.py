# read the pickles 
import os
import sys
import gzip
import pickle
from collections import Counter
import pandas as pd
from multiprocessing import Pool
import time
import re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

def process_file(args):
    k, pathh, thres, seqs_per_this_tax = args
    with gzip.open(pathh+k+".pickle.gzip", 'rb') as f:
        a = pickle.load(f)
    kmers = list(a.values())[0]
    thres_counts = list(map(lambda x: round(x*seqs_per_this_tax),thres))
    all_thres = {tr:[] for tr in thres}
    all_counts = [(kmer,count) for kmer, count in Counter(kmers).items()]
    for tr,tr_c in zip(thres,thres_counts):
        com_all_kmers = set([kmer for kmer, count in all_counts if count > tr_c])
        all_thres[tr] = com_all_kmers
    
    return k, all_thres

def parallel_process(filenames, pathh, thres, seqs_per_tax, n_cores=6):
    with Pool(n_cores) as pool:
        # Create args list for each file
        args = [(k, pathh, thres, seqs_per_tax[k]) for k in filenames]
        results = []
        for i, result in enumerate(pool.imap_unordered(process_file, args)):
            if i % 10 == 0: print(f"{i}/{len(filenames)}      ", end="\r")
            results.append(result)
            #sys.exit()
    return dict(results)
    
def uniq_kmers(ref_seqs, out_dir, name_patter, k_vals = [], thres = [],n_cores=6):
    start = time.time()
    
    if len(k_vals) == 0: print("No list of kmer sizes provided!") ; return
    if not out_dir: print("No output directory provided!") ; return
    if out_dir[-1] != "/": out_dir = out_dir+"/" 
    
    if name_patter == None: print("No name pattern!") ; return
    elif "g" in name_patter.lower(): pat = "g"
    elif "f" in name_patter.lower(): pat = "f"
    elif "c" in name_patter.lower(): pat = "c"
    elif "o" in name_patter.lower(): pat = "o"
    name_patter = f"({pat}__).+?(;)"
    id_list = [re.search(name_patter,i.id).group().strip(";") for i in SeqIO.parse(ref_seqs,"fasta")]
    #print(*id_list[:100],sep="\n")
    seqs_per_tax = {name:count for name, count in Counter(id_list).items()}
    #print(seqs_per_tax)
    dfs = []
    for k_val in k_vals:
        pathh = f"{out_dir}all{k_val}mers/"
        #os.makedirs(pathh, exist_ok=True)
        filenames = os.listdir(pathh)
        filenames = list(map(lambda x:x.split(".pickle")[0],filenames))
        print("files:",len(filenames))

        print(f"KMER: {k_val}, THRESHOLD: {thres}")
        dicc = {k:[] for k in filenames}

        dicc = parallel_process(filenames, pathh, thres, seqs_per_tax, n_cores=n_cores)
        
        df = pd.DataFrame(dicc).T
        for col in list(df.columns):
            df.rename({col:str(k_val)+"mer_"+str(col)+"thres"},inplace=True,axis="columns")
        dfs.append(df)

    
    fdf = pd.concat(dfs,axis="columns")
    ofdf = pd.DataFrame()
    for c in fdf.columns:
        all_kmers = []
        for j in [list(i) for i in list(fdf[c])]: all_kmers.extend(j)
        rare_kmers = set([kmer for kmer, count in Counter(all_kmers).items() if count > 1])
        ofdf[c] = fdf[c].apply(lambda x: x & rare_kmers)
        
    print("Execution Time:",time.time()-start)

    return {"raw":fdf,"uniq":ofdf}