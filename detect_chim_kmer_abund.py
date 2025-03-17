import gzip
import pickle
from Bio import SeqIO
import numpy as np

def detect_kmers(fasta_file, pickle_gzip_files):
    # Read the FASTA file
    fasta_reads = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_reads[record.id] = str(record.seq)
    
    # Load all pickle.gzip files into memory
    fungi_data = {}
    for file in pickle_gzip_files:
        with gzip.open(file, 'rb') as f:
            data = pickle.load(f)
            fungi_data.update(data)
    
    # Process each read in the FASTA file
    kmer_size = 10
    results = {}
    for read_id, sequence in fasta_reads.items():
        read_results = {}
        
        # Generate kmers for the sequence
        for i in range(len(sequence) - kmer_size + 1):
            kmer = sequence[i:i + kmer_size]
            
            # Check for kmer matches in fungi data
            for fungi, kmer_dict in fungi_data.items():
                if kmer in kmer_dict:
                    median_value, count_of_kmers = kmer_dict[kmer]
                    if fungi not in read_results:
                        read_results[fungi] = []
                    read_results[fungi].append((i, median_value, count_of_kmers))
        
        # Calculate mean and median for each fungi
        fungi_summary = {}
        for fungi, matches in read_results.items():
            positions = [match[0] for match in matches]
            medians = [match[1] for match in matches]
            counts = [match[2] for match in matches]
            
            fungi_summary[fungi] = {
                "positions": positions,
                "mean_median_value": np.mean(medians),
                "median_median_value": np.median(medians),
                "mean_count_of_kmers": np.mean(counts),
                "median_count_of_kmers": np.median(counts)
            }
        
        results[read_id] = fungi_summary
    
    return results
