# new function to make chimeras
# TODO: prefix_range/suffix_range variable does nothing - make random_split take it as input
def make_chimeras(fasta_file,out_file,prefix_range = None, suffix_range = None,file_type="fasta",split_fn="random",normal_split_pars=(0.5,0.2),seed=None,repeat=1):
    if fasta_file is list: orig_seqs = fasta_file # takes biopython records
    else: orig_seqs = [i for i in SeqIO.parse(fasta_file,file_type)] # or reads the file
    
    if seed != None: random.seed(seed) # sets a seed if provided otherwise random

    #two ways of getting prefixes and suffixes to make chimera from the input:
    def random_split(seq,lenka,is_prefix): # random range
        if is_prefix: return seq[:int(lenka*random.random())]
        else: return seq[int(lenka*random.random()):]
    
    def gauss_split(seq,lenka,is_prefix): # splits have normal distribution
        split_point = int(random.gauss(normal_split_pars[0]*lenka,normal_split_pars[1])) # mu and stderror
        split_point = max(1, min(lenka - 1, split_point))  # index must be in range of the sequence
        if is_prefix: return seq[:split_point]
        else: return seq[split_point:]
            
    if split_fn.lower() == "random":  split_fn = random_split
    elif split_fn.lower() == "gauss": split_fn =  gauss_split

    prefixes, suffixes = [],[]    
    for f in orig_seqs:
        lenka = len(f.seq)
        for r in range(repeat): # if repeat > 1 -> it makes multiple preffixes and suffixes of the sequence
            seq = split_fn(str(f.seq),lenka,True)
            prefixes.append((seq,"".join([f.id,f"_ARTIF_CHIM-0-{len(seq)}"])))
            seq = split_fn(str(f.seq),lenka,False)
            suffixes.append((seq,"".join([f.id,f"_ARTIF_CHIM-{lenka-len(seq)}-{lenka}"])))
            
            if seed != None and repeat > 1: random.seed(seed+1+r)
                
    random.shuffle(prefixes) # randomizes the order
    random.shuffle(suffixes)
    chimeras_recs = [SeqRecord(Seq("".join([p[0],s[0]])),id="___".join([p[1],s[1]])) for p,s in zip(prefixes,suffixes)]

    return chimeras_recs


in_file = "/home/mg/ChimeraProject/my_data/kratky_100000_ITS2_TEST.fa"
out_file = in_file.replace(".fa","_artif_chims.fa")
artif_chims = make_chimeras(in_file,out_file,split_fn="gauss")