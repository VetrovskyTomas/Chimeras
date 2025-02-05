### DATASETS

**ITS**

ITS full: ITS_ ITS9mun_ITS4ngsUn_full_149350seq.7z

raw ITS2 (gITS7/ITS4)
https://www.biomed.cas.cz/mbu/lbwrf/seed/archive/ITS_EXAMPLE_DATA.zip

ITS2: ITS2_gITS7_ITS4_joined_68000seq.zip


**16S**

raw 16S short (515F/806R)
https://www.biomed.cas.cz/mbu/lbwrf/seed/archive/BAC_EXAMPLE_DATA.zip

16S full: 16S_27F_1429R_full_90000seq.7z

### Chimeras
Reference dataset for UCHIME: 10.1264/jsme2.ME14121 - just fungi from unite

### REFERENCE BASED

ITS2_PERMANENT_CLUSTERS_SEEDs.fa https://www.biomed.cas.cz/mbu/lbwrf/seed/archive/ITS2_PERMANENT_CLUSTERS_SEEDs.7z

ITS2_TEST.fa https://www.biomed.cas.cz/mbu/lbwrf/seed/archive/ITS2_TEST.7z

ITS2_TEST_chimeric.zip [ chimeric = chimclean.fasta(output) - ITS2_TEST.fa(input) ]

`vsearch --uchime_ref ITS2_TEST.fa --db ITS2_PERMANENT_CLUSTERS_SEEDs.fa --nonchimeras chimclean.fasta`

vsearch v2.21.2_linux_x86_64, 5948.9GB RAM, 1152 cores

https://github.com/torognes/vsearch


Reading file REFERENCE.fa 100%

110417743 nt in 633435 seqs, min 40, max 1415, avg 174

Masking 100%

Counting k-mers 100%

Creating k-mer index 100%

Detecting chimeras 100%

Found 347395 (6.5%) chimeras, 4935564 (92.6%) non-chimeras, and 45776 (0.9%) borderline sequences in 5328735 unique sequences.

Taking abundance information into account, this corresponds to 347395 (6.5%) chimeras, 4935564 (92.6%) non-chimeras, and 45776 (0.9%) borderline sequences in 5328735 total sequences.

Simera folder uploaded here was taken from https://github.com/bnichols1979/Simera?tab=readme-ov-file and Simera.h file was modified to work with modern GCC:
  Tested on WSL1 Debian with gcc (Debian 12.2.0-14) 12.2.0 , GNU Make 4.3 - Built for x86_64-pc-linux-gnu , and GSL (libgsl-dev) (2.7.1+dfsg-5+deb12u1)
  
./Simera -i tutorial/example_input.fa -o michal_examples/example_output -n 25 -s 10000 -l 0.00005 -c 10000 -f GTGNCAGCMGCCGCGGTAA -r GGACTACHVGGGTWTCTAAT -x 

========================================
  Simera - Simulates Chimera Formation
  Use option -h for help
========================================
Generated 10,000 potential chimeras...
PCR Round 25...
Done.
