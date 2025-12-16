#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 11:29:10 2025

@author: Ariane Neumann
"""
'''
Script name: 
malaria.py

Description: 
The task is to combine to files in a certain way. A fasta file and a 
blastx file. From the blastx file, I want to take the protein description, which 
is missing in the fasta file. Based on the gene_id, I want to add the protein description 
from the blastx file then to the matching gene_id in the fastafile. 
        
User defined functions: 
fasta_malaria, blastx_malaria, output_txt, blast_pos_hits, column, gene_id, first_query_pos,
protein_description, fasta_header, fasta_sequences, seq_fragments, current_header, sequence, 
matched, sequence, seq_fragments, full_seq

Procedure:
1. Preparation of the script
2. Handling the blastx file
3. Handling the fasta file
4. Write the output file
 
Input:
    fasta file "malaria.fna"
    blastx file "malaria.blastx.tab"

Output:
    text file "output.txt"
    
Usage: 
    python3 malaria.py malaria.fna malaria.blastx.tab output.txt

Version: 1.1
Date 2025-10-13
Author: Ariane Neumann

#remarks: This is the updated version, as i thought in-code documentation is enough.
So here is now formatted documentation section.
'''
#%% Preparation of the script

import sys

if len(sys.argv) != 4:
    sys.argv = ["malaria.py", "malaria.fna", "malaria.blastx.tab", "output.txt"] 

fasta_malaria = sys.argv[1]
blastx_malaria = sys.argv[2]
output_txt = sys.argv[3]
 
print("Creating a dictionary") # good to write a header for my own thought process and to divide into sections
blast_pos_hits ={} # creates an empty dictionary
print("="*80) # 50 was suggested in the lecture, but 80 looks better for me.

print("="*80)

#%% Blastx files handling

print("Handling the blastx file")
try: 
    with open(blastx_malaria, "r") as r:  # this is softcoding
        header = r.readline() # best idea for bigger files
        counter = 1
       
        for line in r: #starting the FOR loop
            print(f"{counter}, {line.strip()}") #debugging and remove empty whitespace
            
            column = line.strip().split("\t")  
            gene_id = column[0] # In column 0 is the gene_id
            first_query_pos = column[2] #first_query_pos is in column 2
            protein_description = column[9] #protein description is in column 9
           
            if first_query_pos.lower() != "null": #if the first query position has not "null" as an entry, the hit is valid. "null" indicates that no hit is found. If i put parenthesis here, it will think it is a function and therefore cannot call it. So i removed the () and then the code runs.
                print(f"I found a valid hit for: {gene_id}")
                blast_pos_hits[gene_id] = protein_description # this connects the dictionary with the key word "gene_id" to the value "protein descriptio"
            else:
                print(f"I found no valid hit for: {gene_id}")    
                
            counter = counter + 1  #here the loop ends
        
except FileNotFoundError: # be careful here with the indent. "except" needs to align with "try", otherwise it gives a syntax error. Try and except is again a debugging technique.
    print(f" Error: The blastx file '{blastx_malaria}' was not found")
    sys.exit(1)
            
print("\nI finished reading blastx_malaria.") 
print(f"My collected {len(blast_pos_hits)} protein descriptions.\n")
        
print("="*80)

#%% Handling the fasta file

fasta_header = [] # create a list to store headers
fasta_sequences = [] # creates a list to store the sequences
seq_fragments = [] # creates a temporary list for multi line sequence assembly
current_header = ''
counter = 1

with open(fasta_malaria, "r") as r: # with open safely closes files even if errors occur
    for line in r: # FOR loop
        print(f"{counter}, {line.strip()}") #debugging print statement
        
        if line.startswith('>'): #fasta files start with this sign
            print("Header line") #debugging statement, confirming header is found
            
            if seq_fragments:
                print("Please combine previous sequence")
                full_seq = ''.join(seq_fragments)
                fasta_sequences.append(full_seq)
                seq_fragments = []
                
            current_header = line.strip()
            fasta_header.append(current_header)
        else:
            print("A sequence line detected")
            seq_fragments.append(line.strip())
            
        counter += 1
        
    if seq_fragments:
        print("Please combine last sequence")
        full_seq = ''.join(seq_fragments)
        fasta_sequences.append(full_seq)
        
print("\nI finished reading malaria.fna")
print(f"Total sequence read is: {len(fasta_header)}\n")

print("="*80)

#%% Write the output file

with open(output_txt, "w") as w:
    matched = 0 # initialising the counter, setting it to zero
    for i in range(len(fasta_header)):
        header = fasta_header[i]
        sequence = fasta_sequences[i]
        
        gene_id = header.split()[0][1:]
        
        if gene_id in blast_pos_hits:
            updated_header = f"{header}\tprotein={blast_pos_hits[gene_id]}"
            w.write(updated_header + '\n')
            w.write(sequence + '\n')
            print(f"Matching gene {gene_id} with protein: {blast_pos_hits[gene_id]}")
            matched += 1
        else:
            print(f"Skip gene with no blastx hit: {gene_id}")

print(f"\nI finished writing to my {output_txt}")
print(f"The code matched {matched} out of {len(fasta_header)} sequences.")

print("=" * 80)

print("Super cool, your code runs without breaking down")

































