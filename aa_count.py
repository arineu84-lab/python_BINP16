#!/usr/bin/env python3
# -*- coding: utf-8 -*- 
"""
Script name: 
    aa_count.py

Description: 
    The task is to count the absolute abundance for the amino acids in a file. 
    As sub-tasks, I ask for the count of aa in each sequence in the file, in addition to 
    the total abundance count of aa in all sequences together.
    In order to count something in Python, a counted needs to be imported from the collections
    before starting the code. 
        
User defined functions: 
    read_txt, count_amino_acids, get_count, sort_counts

Procedure:
    1. Preparation of the script
    2. Opening and handling the text file containing amino acids sequences
    3. Count the amino acid abundance
    4. Write to output file
 
Input:
    text file "amino.faa"

Output:
    text file "counted_aa.txt"
    
Usage: 
    python3 aa_count.py amino.faa counted_aa.txt

Version: 1.0
Date 2025-10-18
Author: Ariane Neumann
"""
#------------------------------------------------------------------------------

import sys
import os
from collections import Counter # to use a counter, it first needs to be imported

# Setting the arguments for my code. I use the output file from the previous section (dna2protein) as input file here
sys.argv = ["aa_count.py", "amino.faa", "counted_aa.txt"]

try:
    if len(sys.argv) >= 3:
        aa_input = sys.argv[1]
        counted_aa = sys.argv[2]
    else:
        print("Please enter file path")

# Validate input file
    if not os.path.exists(aa_input):
        print(f"Error: '{aa_input}' does not exist.")
        sys.exit(1)
    if not os.path.isfile(aa_input):
        print(f"Error: '{aa_input}' is not a file.")
        sys.exit(1)
        if not aa_input.lower().endswith('.txt'):
            print(f"Warning: '{aa_input}' is NOT a text file.")

# Validate output path
    output_dir = os.path.dirname(counted_aa)
    if output_dir and not os.path.exists(output_dir):
        print(f"Error: '{output_dir}' does not exist.")
        sys.exit(1)
except Exception as e:
    print(f"Error for input:{e}")
    sys.exit(1)
    
# Now I create a set for the amino acids 
amino_acids = set("QRIKLMNACDYEPVWSTFGH")

#------------------------------------------------------------------------------

# Defining the function to read the txt input file
def read_txt(file_path):
    try:
        sequences = [] # creating an empty list for the sequences I will look at 
        headers = []
        
        with open(file_path, 'r') as a: # open in "read" mode
            current_seq = ''# giving an empty string to the current sequence as placeholder
            current_header = ''
            for line in a: # looking into the file
                line = line.strip() # remove whitespace in beginning and end of sequence
                
                if line.startswith(">"): # to make sure that I only consider a sequence when it starts with ">"
                    if current_seq:
                        sequences.append(current_seq) # While still in the FOR loop, I will only append to the current sequence. Once I jump out of the loop, and iterate to the next sequence, I will stop added this that sequence and move to the next one
                        current_seq = ''
                    current_header = line[1:] # removes ">" from the header line
                    
                else:
                    current_seq += line.upper() # check to have all letters in capital
            if current_seq:
                sequences.append(current_seq) # after FOR loop, add last sequence
                headers.append(current_header)

        return sequences, headers # returns the former empty list "sequences" without printing to screen
    
    except Exception as e:
        print(f"Error: cannot reader fasta file: {e}")
        sys.exit(1)
#------------------------------------------------------------------------------

# Defining function to count the amino acids
def count_amino_acids(sequences):
    try: 
        per_sequence_counts = [] # creating an empty list for counting aa per sequence storing individual counter
        total_counts = Counter() # while it would be more intuitiv to start with this before "per_seq", it might mix up the counts from the individual sequences.
    
        for seq in sequences: # starting the FOR loop
            seq_counter = Counter() # counting within one loop iteration
        
            for aa in seq: # one more level in, counting within the sequence
                if aa == '*' :
                    continue
                
                if aa in amino_acids: # if it is a standard aa, add 1
                    seq_counter[aa] += 1
                    total_counts[aa] += 1
                else:
                    seq_counter['X'] += 1 # even if it is an unknown aa, count 1
                    total_counts['X'] += 1
                    per_sequence_counts.append(seq_counter) # append or add all counted aa for the current sequence
    
        return per_sequence_counts, total_counts 
    except Exception as e:
        print(f" Warning: Cannot count aa: {e}")
#------------------------------------------------------------------------------

# Defining function to count each amino_acid/count pair
def get_count(pair): # this needs to be defined before I can use "get_count" in the next function
    return pair[1]

# Defining the sorting function
def sort_counts(counts):
    return sorted(counts.items(), key=get_count, reverse=True)

#------------------------------------------------------------------------------

# Calling the functions and counting the amino acids
sequences,headers = read_txt(aa_input) # reading the input file
per_seq_counts, total_counts = count_amino_acids(sequences) # counting the amino acids
sorted_total = sort_counts(total_counts) # sorting the amino acids

# This will be in the output file "count_out", starting with the total aa sorted and then adding the per sequence sorted.
output_lines = ["# Total amino acid counts (sorted by abundance):"]
for aa, count in sorted_total:
    output_lines.append(f"{aa}\t{count}")

# Open with in "write" mode, saving the output into the text file
with open(counted_aa, 'w') as out:
    out.write('\n'.join(output_lines) + '\n')

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# Optional: In case I want to see the output printed on screen
show_output = input("Do you want to print the contents of the output file to screen? (yes/no): ").strip().lower()
if show_output == 'yes':
    with open(counted_aa, 'r') as out:  # Use the correct variable name here
        print("\n--- Output File Content ---")
        print(out.read())

