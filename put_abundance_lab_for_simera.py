# puts "1" (or something else) to read name so it can be accepted by Simera
import argparse
from time import sleep
import sys
import random
parser = argparse.ArgumentParser(
                    prog='put_abund_lable',
                    description='-i input fasta file with only one line sequences, -l puts "_{num}" to the end of each read name, ',
                    epilog='''nic''')


parser.add_argument('-i', '--input')   
parser.add_argument('-l', '--label')      
args = parser.parse_args()
input_file = args.input
label = ""
if args.label: label = str(args.label)
#print("label:",label,"input:",input_file,sep="\n")

with open(input_file,"r") as f:
    a=f.read()

all_seqs = list(filter(lambda x: len(x)>1,a.split(">")))
all_seqs2=list(map(lambda x: ">"+x.split("\n")[0]+label+"\n"+"".join(x.split("\n")[1:]),all_seqs))
print("\n".join(all_seqs2))

