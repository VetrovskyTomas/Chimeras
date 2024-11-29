### DATASETS

ITS

ITS long:

raw ITS2 (gITS7/ITS4)
https://www.biomed.cas.cz/mbu/lbwrf/seed/archive/ITS_EXAMPLE_DATA.zip

ITS2: ITS2_gITS7_ITS4_joined_68000seq.zip


16S

raw 16S short (515F/806R)
https://www.biomed.cas.cz/mbu/lbwrf/seed/archive/BAC_EXAMPLE_DATA.zip


### Chimeras
Reference dataset for UCHIME: 10.1264/jsme2.ME14121 - just fungi from unite

### REFERENCE BASED

ITS2_PERMANENT_CLUSTERS_SEEDs.fa https://www.biomed.cas.cz/mbu/lbwrf/seed/archive/ITS2_PERMANENT_CLUSTERS_SEEDs.fa.gz

`vsearch --uchime_ref TEST.fa --db ITS2_PERMANENT_CLUSTERS_SEEDs.fa --nonchimeras chimclean.fasta`

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

