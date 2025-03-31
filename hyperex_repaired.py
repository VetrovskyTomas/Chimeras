def hyperex_repaired(fasta_seq,forward_primer,reverse_primer):
    from time import time
    import os
    import subprocess
    import re
    from collections import Counter
    from Bio import SeqIO
    start = time()
    # Initialize variables
    results = []  # List to store the content of hyperex_out.fa files
    warnings = []  # List to store the warnings
    errors = []  # List to store the errors
    hyperex_out_fa = "hyperex_temp_out.fa"
    hyperex_log = "hyperex.log"
    hyperex_gff = "hyperex_temp_out.gff"
    hyperex_infile = "infile.fa"
    # Clean up any existing temporary files
    if os.path.exists(hyperex_log): os.remove(hyperex_log)
    if os.path.exists(hyperex_gff): os.remove(hyperex_gff)
    if os.path.exists(hyperex_out_fa): os.remove(hyperex_out_fa)
    if os.path.exists(hyperex_infile): os.remove(hyperex_infile)
    if os.path.exists("hyperex_temp_out_2.fa"): os.remove("hyperex_temp_out_2.fa")
    if os.path.exists("hyperex_temp_out_2.gff"): os.remove("hyperex_temp_out_2.gff")
    if os.path.exists("temp_fasta_input_file_2.fasta"): os.remove("temp_fasta_input_file_2.fasta")
    
    # Command template for hyperex

    # Function to count occurrences of the substring in the log file
    def count_log_occurrences(log_file, substring):
        count = 0
        with open(log_file, "r") as log:
            for line in log:
                if substring in line:
                    count += 1
        return count

    # Process the fasta sequences
    lenka = len(fasta_seq)
    start_index = 0
    forward_primer = forward_primer.strip()
    reverse_primer = reverse_primer.strip()
    command = ["hyperex", "-f", forward_primer, "-r", reverse_primer, "-p", "hyperex_temp_out"]
    #command = ["hyperex", "-f", "GTGARTCATCGARTCTTTG", "-r", "TCCTCCGCTTATTGATATGC", "-p", "hyperex_temp_out"]

    #this loop searches for bad reads
    bad_indices = []
    while start_index < lenka:
        print(f"{start_index}/{lenka}", end="\r")
        fasta_data = "".join([f">{record.id}\n{str(record.seq)}\n" for record in fasta_seq[start_index:]])
        process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate(input=fasta_data)
        warnings.extend(list(map(lambda x: x.replace("\x1b[0m]", ""), re.findall("WARN(.+)", stdout))))
        errors.append(stderr)
        #print(stderr)
        # Read results if the output file exists
        #if os.path.exists(hyperex_out_fa):
        #    results.extend(list(SeqIO.parse(hyperex_out_fa,"fasta")))
        
        # Check for errors and adjust the start index
        if stderr:
            if os.path.exists(hyperex_log):
                log_count = count_log_occurrences(hyperex_log, "[hyperex::utils][INFO] Sequence type is DNA")
                start_index += log_count
                bad_indices.append(start_index-1) # -1
                #print("STDERRR!!!",stderr)
                #print("start_index!!!",start_index,"!!!")
            if len(stderr) == 0:
                print("\nno more errors")
                break
        else:break
        if os.path.exists(hyperex_out_fa): os.remove(hyperex_out_fa)
        if os.path.exists(hyperex_log): os.remove(hyperex_log)
        if os.path.exists(hyperex_gff): os.remove(hyperex_gff)
        if os.path.exists(hyperex_infile): os.remove(hyperex_infile)

    if os.path.exists(hyperex_out_fa): os.remove(hyperex_out_fa)
    if os.path.exists(hyperex_log): os.remove(hyperex_log)
    if os.path.exists(hyperex_gff): os.remove(hyperex_gff)
    if os.path.exists(hyperex_infile): os.remove(hyperex_infile)



    results = []  # List to store the content of hyperex_out.fa files
    #warnings = []  # List to store the warnings
    #errors = []  # List to store the errors
    
        #this runs hyperex on fasta without bad reads
        # Remove bad indices from the fasta sequence
    fasta_seq = [record for i, record in enumerate(fasta_seq) if i not in bad_indices]
    # Save the filtered fasta_seq to a temporary file
    command = ["hyperex", "-f", forward_primer, "-r", reverse_primer, "-p", "hyperex_temp_out_2","temp_fasta_input_file_2.fasta"]
    with open("temp_fasta_input_file_2.fasta", "w") as temp_file:
        SeqIO.write(fasta_seq, temp_file, "fasta-2line")
    subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)

    results = list(SeqIO.parse("hyperex_temp_out_2.fa", "fasta"))

    
    # Clean up any remaining temporary files
    if os.path.exists(hyperex_out_fa): os.remove(hyperex_out_fa)
    if os.path.exists("hyperex_temp_out_2.fa"): os.remove("hyperex_temp_out_2.fa")
    if os.path.exists("hyperex_temp_out_2.gff"): os.remove("hyperex_temp_out_2.gff")
    if os.path.exists("temp_fasta_input_file_2.fasta"): os.remove("temp_fasta_input_file_2.fasta")
    if os.path.exists(hyperex_log): os.remove(hyperex_log)
    if os.path.exists(hyperex_gff): os.remove(hyperex_gff)
    if os.path.exists(hyperex_infile): os.remove(hyperex_infile)
    

    # Print results
    real_results = [i for i in results if len(i) > 0]
    
    warning_counts = Counter(warnings)
    #error_counts = Counter(errors)
    print("WARNINGS:")
    for k,v in warning_counts.items():
        print(v,"-times ",k,sep="")
    print("number of real results:", len(real_results)," "*10, end=";  ")
    print("ERROR COUNT:",len(errors), sep="  ",end=";  ")
    print("SKIPPED READS:",lenka-len(real_results), sep="  ",end=";  ")
    print("Execution time:", time()-start)
    return real_results