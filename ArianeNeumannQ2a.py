#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script name: 
    ArianeNeumannQ2a.py

Description: 
    The task is to write a function called "gc_content(sequence, window_size)", which 
    calculates the GC content within a sliding window (sliding 1 position each time) 
    across a DNA sequence from an input fasta file provided in "examples" folder.
    The function should return a list of GC content percentages for each window.

User defined functions: 
    gc_content(sequence: str, window_size: int) -> list

Procedure:  
    1. open and run manage_examples.py in Spyder6 to populate directory
    2. Inspect the folder and familiarise yourself with the content
    3. Define function and set parameters
    4. Write output file

Input: 
    input_file 

Output: 
    2a_output_ArianeNeumann.txt

Usage: 
    python3 ArianeNeumannQ2a.py input_file output_file

Version: 1.0 Date 2025-10-30 Author: Ariane Neumann
"""
#========================== Defining functions =================================

# Defining function setting parameters for the GC content
def gc_content(sequence: str, window_size: int) -> list:  # this is a suggestion by Anna to better remember the exact type of the variables
    gc_percent = [] # creating an empty list to enter the percentages
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        gc_count = window.count('G') + window.count('C')
        gc_percentage = (gc_count / window_size) * 100
        gc_percent.append(round(gc_percentage, 1))
    return gc_percent

# Defining function to read the input file
def read_sequence(input_file: str) -> str: # returns function as a string
    with open(input_file, 'r') as a:
        lines = a.readlines() # reading line by line
        sequence = ''.join(line.strip() 
                           for line in lines 
                           if not line.startswith('>')) # removing header lines starting with ">"
    return sequence.upper() # making sure that sequence letters all capital case

# Defining function to write the output file
def write_output(output_file):
    output_file = "2a_output_ArianeNeumann.txt"
    with open(output_file, 'w') as a:
        for value in gc_values:
            a.write(f"{value}\n") # writing values followed by new line
    print(f"GC content values written to {output_file}")
    
#============================= Calling functions ==============================

# User input required entering file path
input_file = input("Enter file path + file name here: ")
sequence = read_sequence(input_file)

# User can decide actively on window size, suggested is a size of 5
window_size = int(input("Enter window size: "))

# Based on user input, calculate the GC content and write output to file
gc_values = gc_content(sequence, window_size)
write_output(gc_values)