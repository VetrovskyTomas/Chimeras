import argparse
from time import sleep
import sys
import random
parser = argparse.ArgumentParser(
                    prog='put_abund_lable',
                    description=' Takes fasta from Simera and converts it into fasta with abundances for vschearch uchime -i input fasta file from simera, output is stdout',
                    epilog='''nic''')


parser.add_argument('-i', '--input')   
#parser.add_argument('-l', '--label')      
args = parser.parse_args()
input_file = args.input


with open(input_file,"r") as f:
    a=f.read()

all_seqs = list(filter(lambda x: len(x)>1,a.split(">")))
all_seqs2=list(map(lambda x: ">"+"_".join(x.split("\n")[0].split("_")[:-1])+";size="+x.split("\n")[0].split("_")[-1]+";\n"+"".join(x.split("\n")[1:]),all_seqs))
all_seqs2 = list(filter(lambda x: "size=0" not in x,all_seqs2))
print("\n".join(all_seqs2))

