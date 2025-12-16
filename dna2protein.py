#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script name: 
    dna2protein.py
Input:
    fasta file "DNA_seq.fasta"

Output:
    text file "translated_seq.txt"
    
Description: 
    The task is to translate a DNA nucleotide sequence into a protein amino acid sequence using the standard genetic code.
    First, I need to create a "fasta" or similar text file with different sequences. (Of course one can also use already available files).
    The file can contain not only "A,G,T,C" but also "N". Unknown bases are translated into "X" amino acid, and stop codons should
    be displayed as "*" in the protein sequence.  
        
Procedure:
    1. Preparation of the script
    2. Handling the fasta file
    3. Translating the sequence
    4. Catching errors, including conditions (if/else)
    5. Write to output file
 
User defined functions: 
    read_fasta, translate_sequence
    
Usage: 
    python3 dna2protein.py DNA_seq.fasta translated_seq.txt

Version: 1.0
Date 2025-10-18
Author: Ariane Neumann
"""
#------------------------------------------------------------------------------
import sys
import os

# first setting the arguments for my code
if len(sys.argv) != 3: # script is 0 + input + output = 3
    sys.argv = ["dna2protein.py", "DNA_seq.fasta", "translated_seq.txt"]  # I modified a fasta file from my old project, so i kept the ".fasta"
# for sharing the script with someone, I might need to add a "sys.exit" check here. But for whatever reason it always crashes my code. So for now, I will leave it out.
DNAseq = sys.argv[1]
seq_translated = sys.argv[2]

# Validate my input and output files, to be sure they open and exist
# Input
if not os.path.exists(DNAseq): 
    print(f"Error: '{DNAseq}' does not exist.")
    sys.exit(1)
if not DNAseq.lower().endswith(('.fna', '.fa', '.fasta')):
    print(f"Warning: '{DNAseq}' is NOT a FASTA file.")
# Output
output_dir = os.path.dirname(seq_translated)
if output_dir and not os.path.exists(output_dir):
    print(f"Error: '{output_dir}' does not exist.")
    sys.exit(1)
     
#------------------------------------------------------------------------------
# Creating a dictionary for the standard genetic code. Would be great ot have this within the "read_fasta" function, but did not manage
genetic_code = {                    
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

#------------------------------------------------------------------------------

# Defining a function for reading the fasta file
def read_fasta(file_path): 
    sequences = {}  # Dictionary for sequences with headers
    
    # Error checking if DNA strand is codig or non-coding
    strand_type = input("Is the input the coding strand? (yes/no): ").strip().lower()

    if strand_type == "yes": # checking if this is coding strand
        print("--> Super, let's go with the coding strand")
    elif strand_type == "no":
        print("--> Ok, first convert T to U.") # if user says non-coding strand, then T needs to be converted to U first. 
        for seq_id in sequences:
            sequences[seq_id] = sequences[seq_id].replace('T', 'U')
        print("Conversion complete. Proceeding with translation.")
    else:
        print("Invalid input. Please answer 'yes' or 'no'.")
        sys.exit(1)

    with open(file_path, 'r') as a:
        current_id = ''         # current sequence ID
        found_sequence = False  # Inbuilt check with boolean to make sure that sequence is found

        for line in a:
            line = line.strip()  # Remove whitespaces in beginning and end of sequence. 

            if not line:
                continue  # Skip empty lines. Might not be needed here, since I generated the file myself. But can happen in manually generated files

            if line.startswith('>'): # The sequences always start with a header line ">", this prompt tells the code where to start reading
                current_id = line[1:]
                sequences[current_id] = '' # empty string for this sequence. Also making sure that multi-line sequences are part of the current sequence and are appended
            else:
                if current_id == '':
                    # Sequence line appears before any header — malformed FASTA
                    print("Error: FASTA file appears malformed. Sequence found before any header.")
                    sys.exit(1) # exit silently if error occurs

                # Add sequence data to the current sequence in case of multi-line sequences
                sequences[current_id] += line.upper() # ensures consistent formatting. 
                found_sequence = True

    # Final checks after reading the file to make sure they exist
    if not sequences:
        print("Error: No sequences found in the FASTA file.")
        sys.exit(1)
    if not found_sequence:
        print("Error: Headers found but no sequence data.")
        sys.exit(1)
    return sequences  # Return the dictionary of sequences

    
#------------------------------------------------------------------------------
 
# Defining the function for the output file
def translate_sequence(dna, warned):
    protein = '' #empty string for translated protein
    i = 0 # Start at index 0 of DNA sequence
    
    if len(dna) % 3 != 0 and not warned: # Check if the sequence length is divisible by 3
        print("Warning: DNA sequence cannot be divided by 3.")
    
    # Loop through the DNA sequence in steps of 3 codons
    while i + 3 <= len(dna): 
        codon = dna[i:i+3] # 3 nucleotides = 1 codon
        
        # Convert T to U, in case of non-coding strand
        codon = codon.replace('U', 'T')

        if 'N' in codon: # if DNA sequence contains "N", replace in protein sequence with "X". This can be changed to any letter,e.g. "R" in the DNA seq
            protein += 'X'
        else:
            protein += genetic_code.get(codon, 'X') # also translate into "X" if the codon is not correct
        i += 3
    return protein, i // 3, warned # integer division operator, as discussed in lecture 1. Dividing the sequence by 3, but returns only whole numbers (meaning it does not takes something that cannot be divided by 3)

#------------------------------------------------------------------------------

# This is the actual translation step, calling the functions defined above
sequences = read_fasta(DNAseq) # this calls the function defined above "read_fasta"
warned = False
with open(seq_translated, 'w') as out: # writing into the seq_translated output file with "with open", safe way.
    for seq_id, dna_seq in sequences.items(): # starting a FOR loop for the sequences and the IDs
        total_codons = 0 # counts the total codons translated
        protein_seq, codon_count, warned = translate_sequence(dna_seq, warned) # this calls the "translate_sequence" function defined above.
        total_codons += codon_count
        out.write(f">{seq_id}\n{protein_seq}\n")
        warned = True
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# These statements will be printed only if the script runs without crashing
print(f"Translation complete. Output written to '{seq_translated}'")
print(f"Total sequences translated: {len(sequences)}")
print(f"Total codons processed: {total_codons}")

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# Optional: view the translated sequences, only runs if script did not crash
# wanted to avoid printing again everything to the screen as in assignment 1
try: # Another error check to see if the whole code worked.
    show_sequences = input("Do you want to view the translated sequences? (yes/no): ").strip().lower()
    if show_sequences == "yes":
        with open(seq_translated, 'r') as result:
            print("\nTranslated Sequences:\n")
            print(result.read())

except Exception:
    print("This did not work. Try again")









