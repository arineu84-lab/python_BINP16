#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script name: 
    CalculateHapmap.py

Description: 
    The task is to produce a haplotype map for mtDNA and Y mutations. For this,
    SNPs need to be identified from the different input files with different 
    individuals providing their genetic material. Next, major and minor alleles 
    need to be determined and the minor allele frequency needs to be calculated. 
    The final results will be written in an output file.

User defined functions: clean_seq, read_seq, analyse_SNPs, align_haplotype_map

Procedure: 
    1. Read sequences 
    2. Clean up input file, headers and sequences
    3. Identify SNPs and determine alleles 
    4. Align to make output look clean
    5. Write output file

Input: 
    GeneticData - 5.txt 

Output: 
    mtDNA_hapmap.txt Y_hapmap.txt

Usage: 
    python3 CalculateHapmap.py input_file chromosome_name output_file

Version: 1.0 Date 2025-10-23 Author: Ariane Neumann


"""
#==============================================================================
# Setting up the script and cleaning up the input file
#==============================================================================

# Setting up the parameters and importing packages that will be required
from collections import Counter
import os

# Due to inconsistent formatting in the input file, certain special characters need to be removed from the headers first
# Create a dictionary with keys and values for these characters
# Defining the characters to remove within the dictionary and what to replace them with (empty space)
remove_chars = {'>': '', "'": '', '´': '', 'í': '', '\x92': '', '\u00a0': ' ' }  # \x92 removes right single quote, while \u00a0 replaces non-breaking space with regular space

# Create variable for the input file
try:
    input_file = ("GeneticData - 5.txt")
    # Check if the file exists before trying to open it
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file '{input_file}' not found.")

    # Open input file using ISO-8859-1 encoding, as it struggles to otherwise read these characters
    with open(input_file, "r", encoding="ISO-8859-1") as a:
        content = a.read()

    # Replacing each special character in the content with empty spaces as defined above
    for char, replacement in remove_chars.items():
        content = content.replace(char, replacement)

    # Saving the genetic data file without special symbols as "clean version"
    with open("Input_clean.txt", "w", encoding="utf-8") as a:
        a.write(content)

    print("The file'Input_clean.txt' is saved using UTF-8 encoding.")

except Exception as e:
    print(f"Error while cleaning input file: {e}")


# Creating a list of all valid header (which is the first line of each individual). These header are the cleaned header from the "Input_clean.txt" file
valid_headers = [ "Princess Irene", "Prince Fred", "Nicolas II Romanov", "Alexandra Romanov", "Olga Romanov",
    "Tatiana Romanov", "Maria Romanov", "Alexei Romanov", "Suspected body of Anastasia Romanov",
    "Anastasia1", "Anastasia2", "Anastasia3", "Anastasia4", "Anastasia4 son", "Anastasia5",
    "Farmers daughter", "Farmers grandson", "Grigori Rasputin"]

# Open the input file in read mode.
try:
    with open("Input_clean.txt", "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip()]

    current_header = None #this will be first empty 
    mtDNA_profile = [] # creating an empty list for mtDNA
    Y_chrom_profile = [] # creating an empty list for Y chromosome

    # Starting a FOR loop going through the individiduals line by line
    for i in range(len(lines) - 1):
        line = lines[i]
        next_line = lines[i + 1]
    # including conditions to be certain only to use lines with valid header as assigned above
        if line in valid_headers:
            current_header = line
        elif line == "mtDNA" and current_header: # this will combine the header with mtDNA
            mtDNA_profile.append(f"{current_header}\n{next_line}\n")
        elif line == "Y chromosome" and current_header: # combining header line with Y chromosome
            Y_chrom_profile.append(f"{current_header}\n{next_line}\n")
        else: 
            # Raise error if DNA type appears without a valid header
            if current_header is None:
                raise NameError(f"Missing header before line {i}: '{line}'")

    # Save the newly created sets into 2 different txt files. Be sure that they have now good formatting.
    # Save mtDNA sequences
    with open("mtDNA.txt", "w", encoding="utf-8") as f:
        f.writelines(mtDNA_profile)

    # Save Y chromosome sequences
    with open("Ychrom.txt", "w", encoding="utf-8") as f:
        f.writelines(Y_chrom_profile)

    print(" The files 'mtDNA' and 'Ychrom' have been created.")

except NameError as ne:
    print(f"Header error: {ne}")
except Exception as e:
    print(f"Error while splitting profiles: {e}")
    

# Cleaning up the DNA sequences
# Defining the function to clean sequences by replacing "?" with "N"
def clean_seq(input_file, output_file):
    try: 
        with open(input_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()

        cleaned_lines = []
        for line in lines:
            # Replace '?' with 'N' only in sequence lines
            if line.strip() and not line.startswith(" "):  
                cleaned_lines.append(line.replace('?', 'N'))
            else:
                cleaned_lines.append(line)

        with open(output_file, 'w', encoding='utf-8') as f:
            f.writelines(cleaned_lines)

    except FileNotFoundError:
        print(f"File '{input_file}' not found.")
    except ValueError as ve:
        print(f"Value error while cleaning '{input_file}': {ve}")
    except Exception as e:
        print(f"Unexpected error while cleaning '{input_file}': {e}")

# Apply cleaning to both files
clean_seq("mtDNA.txt", "mtDNA_seq.txt")
clean_seq("Ychrom.txt", "Ychrom_seq.txt")

print("Sequences have been cleaned: '?' replaced with 'N'.")


#==============================================================================
# Actual processing of the data set
#==============================================================================

# Function to read sequences from the cleaned files
def read_seq(filename):
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            lines = f.read().splitlines()
        sequences = []
        for i in range(0, len(lines), 2):
            if i + 1 < len(lines):
                seq = lines[i + 1].strip()
                if seq:
                    sequences.append(seq)
        return sequences
    
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return []
    except Exception as e:
        print(f"Unexpected error while reading '{filename}': {e}")
        return []

# Function to find SNPs and calculate MAF
def analyse_SNPs(sequences, chromosome_label):
    try:
        if not sequences:
            print(f"No sequences found for {chromosome_label}.")
            return []

        snp_results = []
        sequence_length = max(len(seq) for seq in sequences)

        for pos in range(sequence_length):
            alleles = []
            for seq in sequences:
                if pos < len(seq):
                    base = seq[pos]
                    if base in 'ACGT':
                        alleles.append(base)
            if len(set(alleles)) > 1: # SNP detection
                counts = Counter(alleles)
                major = counts.most_common(1)[0][0]
                minor = counts.most_common()[-1][0]
                maf = round(counts[minor] / sum(counts.values()), 2)
                allele_str = "/".join(sorted(set(alleles)))
                snp_results.append([ chromosome_label, str(pos + 1), allele_str, major, minor, f"{maf:.2f}"])
        return snp_results

    except ValueError as ve:
        print(f"ValueError: {ve}")
        return []
    except Exception as e:
        print(f"Unexpected error during SNP analysis: {e}")
        return []

# Function to align columns using fixed-width formatting
def align_haplotype_map(data, output_file):
    try:
        header = ["Chromosome", "Position", "Alleles", "MajorAllele", "MinorAllele", "MinorFreq"]
        lines = [header] + data
        col_widths = [max(len(row[i]) for row in lines) for i in range(len(header))]

        formatted_lines = []
        for row in lines:
            formatted_line = "  ".join(row[i].ljust(col_widths[i]) for i in range(len(row)))
            formatted_lines.append(formatted_line)

        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("\n".join(formatted_lines))

    except Exception as e:
        print(f"Error writing haplotype map to '{output_file}': {e}")

# Read sequences from cleaned files
mtDNA_sequences = read_seq("mtDNA_seq.txt")
y_sequences = read_seq("Ychrom_seq.txt")

# Analyse SNPs
mtDNA_snps = analyse_SNPs(mtDNA_sequences, "mtDNA")
y_snps = analyse_SNPs(y_sequences, "Y")

# Align and save output
align_haplotype_map(mtDNA_snps, "mtDNA_hapmap.txt")
align_haplotype_map(y_snps, "Y_hapmap.txt")

print("Aligned haplotype map files have been created.")

#==============================================================================
# !Optional!
#==============================================================================

# Does the user want to print the haplotype maps to the screen
user_input = input("Do you want to print the haplotype maps to the screen? (yes/no): ").strip().lower()

if user_input == 'yes':
    with open("mtDNA_hapmap.txt", "r", encoding="utf-8") as f:
        print("\nContents of mtDNA_hapmap.txt:")
        print(f.read())

    with open("Y_hapmap.txt", "r", encoding="utf-8") as f:
        print("\nContents of Y_hapmap.txt:")
        print(f.read())