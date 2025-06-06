
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import random 



def cut_fasta_in_two_and_shuff(region_fasta,orig_fasta=None,out="region",shift_list=[-100,-50,0,50,100],file_type="fasta",rand=0,seed=666):
    """
    Takes path to records of extracted region (by ITSx) from it (ITS1,ITS2,5.8S etc..)
    and optionally takes path to original file 
    Input file of extracted regions must be subset of records from the original file and must have record names generated by ITSx convention
    IF out == "orig" or 0 ---> 
        cuts the original sequences in the middle of the target region shifted by 0, 50, -100 etc. bases (depending on your shift_list)
        and output dict of lists of randomly joined 1st and 2nd parts of these cut sequences
    IF out == "region" or 1 ->
        reeds only file with the extracted region (from ITSx), cuts it in middle shifted by 0, 50, -100 etc. bases (depending on your shift_list)
        and output dict of lists of randomly joined 1st and 2nd parts of these cut sequences
    """
    random.seed(seed)
    shifted_halves = {str(i):{"first":[],"second":[]} for i in shift_list}
    for shift_from_half in shift_list:
        frst,scnd=[],[]

        if out in ["orig",0]:
            f_iter = SeqIO.parse(region_fasta, file_type)
            f=next(f_iter)
            f_id = f.id[:re.search("\|.*\|",f.id).start()]
            for t in SeqIO.parse(orig_fasta, file_type): 
                if t.id != f_id: continue
                else:
                    span = re.search(str(f.seq),str(t.seq)).span()
                    sum_span = sum(span)
                    #if sum_span % 2 != 0:print(f"WARNING: this 5.8S (index:{n}) length is.. odd, span:{span}, length:{len(str(f.seq))}")
                    if rand:cut_pos = int(sum_span/2)+shift_from_half+random.randint(-rand,rand)
                    else:cut_pos = int(sum_span/2)+shift_from_half
                    if cut_pos < 0: cut_pos=0
                    elif cut_pos > len(str(t.seq)): cut_pos=len(str(t.seq))
                    shifted_halves[str(shift_from_half)]["first"].append((f_id,str(t.seq)[:cut_pos]))
                    shifted_halves[str(shift_from_half)]["second"].append((f_id,str(t.seq)[cut_pos:]))   
                    try:
                        f=next(f_iter)
                        f_id = f.id[:re.search("\|.*\|",f.id).start()]
                    except StopIteration: break
                
        elif out in ["region",1]:
            for f in SeqIO.parse(region_fasta, file_type):
                try:f_id = f.id[:re.search("\|.*\|",f.id).start()]
                except:f_id = f.id
                lenka = len(f.seq)
                if rand:cut_pos = lenka//2 + shift_from_half+random.randint(-rand,rand)
                else:cut_pos = lenka//2 + shift_from_half
                if cut_pos < 0: cut_pos=0
                elif cut_pos > lenka: cut_pos = lenka
                shifted_halves[str(shift_from_half)]["first"].append((f_id,str(f.seq)[:cut_pos]))
                shifted_halves[str(shift_from_half)]["second"].append((f_id,str(f.seq)[cut_pos:]))  
                
    out_recs = {}
    for k, v in shifted_halves.items():
        random.shuffle(v["first"])
        random.shuffle(v["second"])
        new = ["".join([f[1], s[1]]) for f, s in zip(v["first"], v["second"])]
        ids = [f"{f[0]}-FIRST_PART___" + f"{s[0]}-SECOND_PART" for f, s in zip(v["first"], v["second"])]
        out_rec = [SeqRecord(Seq(seq), id=idd) for idd, seq in zip(ids, new)]
        #print(ids);break
        #yield {k:out_recs}
        out_recs[k] = out_rec
    return out_recs
    
    
chimeras = cut_fasta_in_two_and_shuff("../ITSx_1.1.3/outputs/five8_test_j/five8.ITS1.fasta","../ITSx_1.1.3/test_seq/test_j.fasta",rand=6)

print(chimeras.keys(),chimeras["0"][0].id,str(chimeras["0"][0].seq),len(chimeras["0"]))
    
