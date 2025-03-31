# read the pickles !!!
import os
import sys
import gzip
import pickle
from collections import Counter
import pandas as pd
from multiprocessing import Pool
import time
pathh = "/home/mg/ChimeraProject/pickles/all50mers/"
filenames = os.listdir(pathh)
filenames = list(map(lambda x:x.split(".pickle")[0],filenames))
print("files:",len(filenames))

def process_file(args):
    k, pathh, thres = args
    with gzip.open(pathh+k+".pickle.gzip", 'rb') as f:
        a = pickle.load(f)
        kmers = list(a.values())[0]
        if thres > 1: com_all_kmers = set([kmer for kmer, count in Counter(kmers).items() if count > thres])
        elif thres == 1: com_all_kmers = set(kmers)
        return k, com_all_kmers

def parallel_process(filenames, pathh, thres, n_cores=6):
    with Pool(n_cores) as pool:
        # Create args list for each file
        args = [(k, pathh, thres) for k in filenames]
        results = []
        for i, result in enumerate(pool.imap_unordered(process_file, args)):
            if i % 10 == 0: print(f"{i}/{len(filenames)}      ", end="\r")
            results.append(result)
    return dict(results)


for k_val,thres in zip([50]*11 ,[5,10,20,30,40,50,70,100,150,200,300,500]):
    print(f"KMER: {k_val}, THRESHOLD: {thres}\n")
    dicc = {k:[] for k in filenames}
    start = time.time()
    """
    for n,k in enumerate(filenames):
        if n % 10 == 0: print(f"{n}/{len(filenames)}      ",end="\r")
        with gzip.open(pathh+k+".pickle.gzip", 'rb') as f:
            a = pickle.load(f)
            kmers = list(a.values())[0]
            if thres > 1: com_all_kmers = set([kmer for kmer, count in Counter(kmers).items() if count > thres])
            elif thres == 1: com_all_kmers = set(kmers)
            dicc[k] = com_all_kmers
            #print(len(com_all_kmers),len(kmers))  
    print(time.time()-start)
    """
    
    dicc = parallel_process(filenames, pathh, thres)
    print(time.time()-start)
    #kkmm_par_unite = get_kmers_para2(uchime_ref_data,k_val,thres,threads=7)
    with gzip.open(f"/home/mg/ChimeraProject/pickles/kmers_unite_k{k_val}_t{thres}.pickle.gzip", 'wb') as f:
        pickle.dump(dicc, f, protocol=pickle.HIGHEST_PROTOCOL)
    
    log = [f"{k_val}-MERS with threshold: {thres} occurences per family"]
    new_kkmm = {}
    for k,v in dicc.items():
        if len(v) == 0: continue
        new_kkmm[k] = v
    log.append( f"number of all families:\t{len(dicc.keys())}")
    log.append(f"number of families with Kmers above threashold:\t{len(new_kkmm.keys())}")
    
    avg_val_len = []
    for v in new_kkmm.values():
        avg_val_len.append(len(v))
    log.append(f"average count of Kmers per family:\t{sum(avg_val_len)/len(avg_val_len)}")
    log.append(f"max and min count of Kmers per family:\t{max(avg_val_len),min(avg_val_len)}")
    all_seqs_list = []
    for k,v in new_kkmm.items():
        for seq in v:
            all_seqs_list.append(seq)
    log.append(f"number of all Kmers (kmers per family are in a set):\t{len(all_seqs_list)}")
    log.append(f"set of all those Kmers (all family kmer sets joined into one set):\t{len(set(all_seqs_list))}")
    counts = [(kmer,count) for kmer, count in Counter(all_seqs_list).items()]
    uniq_counts = [kmer for kmer, count in Counter(all_seqs_list).items() if count == 1] ######!!!!!!!!!!!!!
    log.append(f"number of all Unique Kmers for a family (these Kmers are Unique to only one family):\t{len(uniq_counts)}")
    with open(f"/home/mg/ChimeraProject/pickles/kmers_unite_k{k_val}_t{thres}.info","w") as f:
        f.write("\n".join(log))
