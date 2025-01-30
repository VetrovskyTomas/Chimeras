# new function to make chimeras
# TODO: prefix_range/suffix_range variable does nothing, does not take fastq, 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random 

def make_chimeras(fasta_file,prefix_range = None, suffix_range = None,file_type="fasta",split_fn="random",normal_split_pars=[0.5,0.5],seed=None,repeat=1):
    if fasta_file is list: 
        orig_seqs = fasta_file # takes biopython records
        orig_len = len(orig_seqs)
        print(orig_len)

    else: 
        orig_len = 0
        if file_type == "fasta":
            with open(fasta_file,"r") as f:
                for i in f:
                    if i.startswith(">"): orig_len+=1
        print(orig_len)
        orig_seqs = SeqIO.parse(fasta_file,file_type)
    
    if seed != None: random.seed(seed) # sets a seed if provided otherwise random

    #two ways of getting prefixes and suffixes to make chimera from the input:
    def random_split(seq,lenka): # random range
        return (seq[:int(lenka*random.random())] , seq[int(lenka*random.random()):])
    
    def gauss_split(seq,lenka): # splits have normal distribution
        sigma = ((lenka / 2) * normal_split_pars[1]) # if normal_split_pars[1] is 0.1, sigma will be 10% of mu
        split_point = int(random.gauss(normal_split_pars[0]*lenka,sigma)) # mu and stddev
        #print(split_point,lenka)
        split_point = max(1, min(lenka - 1, split_point))  # index must be in range of the sequence
        return (seq[:split_point], seq[split_point:])
            
    if split_fn.lower() == "random":  split_fn = random_split
    elif split_fn.lower() == "gauss": split_fn =  gauss_split
            
    idx = 0
    prefixes, suffixes = [None] * (orig_len * repeat), [None] * (orig_len * repeat) # preallocating list
    for n,f in enumerate(orig_seqs):
        #if n % 100000 == 0:print(n,"/",orig_len,end="\r")
        lenka = len(f.seq)
        for r in range(repeat): # if repeat > 1 -> it makes multiple preffixes and suffixes of the sequence
            seq_p, seq_s = split_fn(str(f.seq),lenka)
            prefixes[idx] = (seq_p, f"{f.id}_ARTIF_CHIM-0-{len(seq_p)}")
            suffixes[idx] = (seq_s, f"{f.id}_ARTIF_CHIM-{lenka-len(seq_s)}-{lenka}")
            idx += 1

            if seed != None and repeat > 1: random.seed(seed+1+r)
        #if n > 1000:break
    random.shuffle(prefixes) # randomizes the order
    random.shuffle(suffixes)
    chimeras_recs = [SeqRecord(Seq("".join([p[0],s[0]])),id="___".join([p[1],s[1]]),description="") for p,s in zip(prefixes,suffixes)]

    return chimeras_recs
