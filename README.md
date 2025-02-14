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

## Simera ( doi.org/10.1101/072447 )

Simera folder uploaded here was taken from https://github.com/bnichols1979/Simera?tab=readme-ov-file and Simera.h file was modified to work with modern GCC:
  Tested on WSL1 Debian with gcc (Debian 12.2.0-14) 12.2.0 , GNU Make 4.3 - Built for x86_64-pc-linux-gnu , and GSL (libgsl-dev) (2.7.1+dfsg-5+deb12u1)
  
./Simera -i tutorial/example_input.fa -o michal_examples/example_output -n 25 -s 10000 -l 0.00005 -c 10000 -f GTGNCAGCMGCCGCGGTAA -r GGACTACHVGGGTWTCTAAT -x 

\========================================

  Simera - Simulates Chimera Formation
  
  Use option -h for help
  
\========================================

Generated 10,000 potential chimeras...

PCR Round 25...

Done.

--------------------

cat summary.txt

\# Sequences = 7031

\# Chimeras = 6050

Max length = 491

Total abundance = 38444725

Chimera abundance = 4136442

--------------------

cat samp_summary.txt

\# Sequences = 1271

\# Chimeras = 490

Max length = 301

Total abundance = 10000

Chimera abundance = 1064


## Detection Tools

Pintail (WigeoN) (Ashelford et al., 2005) (WigeoN)
	Detects all 16S anomalies, most of them are chimeras, work only on 16S
	Complicated 
Runs Global Allignment against some database and detect general anomalies
	– reimplemented as WigeoN by Haas 2011
 
Bellerophon (Huber, Faulkner and Hugenholtz, 2004)           (DeSantis et al., 2006)
	Reimplemented by Haas 2011
Checks divergence between the sequence and potential perents and if it meets minimal threshold it is labeled as chimera

KmerGenus (Haas et al., 2011)
computed a catalog of all overlaping 50-mers unique to each genus withina reference 16S sequence set. 
Make set of taxon specific Kmers.
Those matching multiple taxa are chimeras

Chimera slayer (Haas et al., 2011)
	“is sensitive to chimeras between closely related 16S genes”
	Tested on 454 and Sanger
Procedure:
30% of the length from each end) were searched against database of reference chimera-free 16S sequences to identify potential parents of a chimera 
candidate parents of a chimera were selected in NAST 
Scoring higher in alignment with potential parents than any other real reference sequence – it is a chimera. It uses function similar to CHECK_CHIMERA in (Komatsoulis and Waterman, 1997) – also complicated
Three-way alignment 
	CS recognized >87% of chimeras with a minimum of 4% chimera-pair divergence
Creation of Artificial Chimeras: 
Join two sequences with a random break point, each being at least 50bp long
Parent sequence divergence must be at least 10%

Chimera-Checker (Nilsson et al., 2010)
	Only ITS, reference based
	Useful only on whole ITS region (ITS1 + ITS2 ± 100bp)
Extracts ITS1 and ITS2 using HMMER
Blasts against INSD
If ITS1 belongs to one taxa and ITS2 to other -> chimera

Perseus (Quince et al., 2011)
	De novo
	Designed for 454, but might work on other sequences (says Edgar)
	Also a Noise-removing tool
Assumes lower frequency of chimeras compared to their parents due to less PCR cycles and each parent will be present with a frequency at least of equal to chimeras
Pairwise alignment of all sequences and sequences with equal or greater abundance

UCHIME (Edgar et al., 2011)
Reference based:
Sequence is divided into 4 chunks and each chunk is searched in a reference database
Two (or more?) candidate parents are identified
Global-X alignment
If segments extracted from parents have higher identity to the sequence, than to them selves -> chimera
De novo:
	Database is build on the fly and is sorted by decreasing abundance 
	Candidate Parents must have abundance of at least 2x abundant than the chimera

DECIPHER 
makes 30-mers from a sequence starting every 5 nucleotides and then matches a query sequence to a set of these Kmers.
There are sets for fungal groups ...

CATCh
Combines DECIPHER, ChimeraSlayer and UCHIME

