# make dataset seq+abund for SIMERA 
# old?? does it even work???
from collections import Counter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import time
start = time.time()
unite = [i for i in SeqIO.parse("/home/mg/ChimeraProject/datasets/unite_latest/UNITE_public_21.04.2024.fasta","fasta")]
unite_seq = [i.seq for i in unite]
unite_counts = [(seq,count) for seq,count in Counter(unite_seq).items()]
unite_orig_dic = {i.seq:i.id for i in unite}

simera_abund = []
for u_c in unite_counts:
    seq,count = u_c
    record = SeqRecord(
        seq,
        id=unite_orig_dic[seq]+f"_{str(count)}",
        name="",
        description="",
    )
    simera_abund.append(record)

with open("/home/mg/ChimeraProject/simera/unite_with_abund_simera.fasta","w") as f:
    SeqIO.write(simera_abund,f,"fasta")
print(time.time() - start)