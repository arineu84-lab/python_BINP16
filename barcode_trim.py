#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script name: 
    barcode_trim.py

Description: 
    The task is to trim specific barcode sequences provided from DNA sequences. 
    The DNA sequences used here were provided in a fastq file format and contained 
    4 lines per sequence. The difference of 3 barcodes is to be detected as pattern, removed/trimmed
    and sequences that previously contained these different barcodes should be stored 
    in separate output files. In case no barcode pattern is found in a sequence,
    this sequence should be stored separately as well. The outcome should be 4 
    output files from 1 input file.
        
User defined functions: 
    trim_barcode, process_fastq

Procedure:
    1. Preparation of the script
    2. Handling the fastq file creating different functions
    3. Removing the barcodes
    4. Write to output file
 
Input:
    text file "barcode.fastq"
 
Output:
    text file "trimmed_DNA.txt"
    
Usage: 
    python3 barcode_trim.py barcode.fastq trimmed_DNA.txt

Version: 1.0
Date 2025-10-18
Author: Ariane Neumann
"""
#------------------------------------------------------------------------------

import sys
import os

# Setting the arguments for my code. I use the output file from the previous section (dna2protein) as input file here
sys.argv = ["barcode_trim.py", "barcode.fastq", "trimmed_DNA.txt"]

barcode_in = sys.argv[1]
trimmed = sys.argv[2]

# Validate input file
if not os.path.exists(barcode_in):
    print(f"Error: '{barcode_in}' does not exist.")
    sys.exit(1)
if not os.path.isfile(barcode_in):
    print(f"Error: '{barcode_in}' is not a file.")
    sys.exit(1)
if not barcode_in.lower().endswith('.fastq'):
    print(f"Warning: '{barcode_in}' is NOT a fastq file.")

# Validate output path
output_dir = os.path.dirname(trimmed)
if output_dir and not os.path.exists(output_dir):
    print(f"Error: '{output_dir}' does not exist.")
    sys.exit(1)
   
# Making a dictionary for the barcodes according to the pdf file
barcodes = {
    "TATCCTCT": "sample1.fastq",
    "GTAAGGAG": "sample2.fastq",
    "TCTCTCCG": "sample3.fastq"
}
undetermined_file = "undetermined.fastq" # when no barcode is found, this should be saved in a separate file

#------------------------------------------------------------------------------

# Defining the function for the barcode trimming 
def trim_barcode(seq, qual, barcode): 
    if seq.startswith(barcode): # check if sequence starts with barcode or ....
        return seq[len(barcode):], qual[len(barcode):]
    elif seq.endswith(barcode): # check if sequence ends with barcode. Better to check in both places even though it should only be in the end
        return seq[:-len(barcode)], qual[:-len(barcode)]
    else:
        return None, None  # if no barcode was found and the sequence is undetermined 
#------------------------------------------------------------------------------

# Defining function how to handle the fastq file
def process_fastq(input_file, output_prefix=None):
    output = {}  # Create empty dictionary 
    try:
        for barcode, filename in barcodes.items():
            name = f"{output_prefix}_{filename}" if output_prefix else filename
            output[barcode] = open(name, 'w')
    
            undetermined_name = f"{output_prefix}_{undetermined_file}" if output_prefix else undetermined_file
            undetermined = open(undetermined_name, 'w')
    except IOError as e:
            print(f"File error: {e}")
            sys.exit(1)
    try:
        with open(input_file, 'r') as a:  # Opens fastq input file in reading mode
            while True:  # Starting a WHILE loop, which indicates a condition with a boolean. "while" a condition is "true", the loop continues
                header = a.readline().strip()  # Reading header and removing whitespace 
                seq = a.readline().strip() 
                plus = a.readline().strip()  
                qual = a.readline().strip()  

                if not header:  # If no header is found in the line, this is the end of the file
                    break  # Exits WHILE loop, when boolean "true" no longer holds

                matched = False  # Initializes a flag to track if a barcode match is found

                for barcode in barcodes:  # FOR loop inside WHILE loop. While line starts with header, it will look for matching barcodes
                    try: 
                        trimmed_seq, trimmed_qual = trim_barcode(seq, qual, barcode)  # Calls the trimming function to remove barcode if matched
                    except Exception as e:
                        print(f"Error trimming barcode: {e}")
                        continue

                    if trimmed_seq is not None:  # If a match was found and trimming was successful
                        output[barcode].write(f"{header}\n{trimmed_seq}\n{plus}\n{trimmed_qual}\n")  # Writes trimmed sequence to the corresponding output file
                        matched = True  # Sets the match flag to True
                        break  # Exits FOR loop when match was found

        if not matched:  # If no barcode matched the sequence
            undetermined.write(f"{header}\n{seq}\n{plus}\n{qual}\n")  # Writes the original sequence to the undetermined file
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error while processing FASTQ: {e}")
        sys.exit(1)

    for trimmed_out in output.values():  # Starting another FOR loop to check for trimmed output in the output values
        trimmed_out.close()  # Closes each output file when loop ran

    undetermined.close()  # Closes the undetermined output file

#------------------------------------------------------------------------------

# Run the barcode trimming
process_fastq(barcode_in, trimmed)
print("Barcode trimming complete. Output files created.")

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# Optional: In case I want to see the output printed on screen
show_output = input("Do you want to print the contents of the output files to screen? (yes/no): ").strip().lower()
if show_output == 'yes':
    print("\n--- Content of output files ---")
    for filename in barcodes.values():
        full_name = f"{trimmed}_{filename}"
        if os.path.exists(full_name):
            print(f"\nContents of {full_name}:")
            with open(full_name, 'r') as f:
                print(f.read())
    undetermined_name = f"{trimmed}_{undetermined_file}"
    if os.path.exists(undetermined_name):
        print(f"\nContents of {undetermined_name}:")
        with open(undetermined_name, 'r') as f:
            print(f.read())









