#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script name: 
    ArianeNeumannQ2b.py

Description: 
    The task is to write a function find_motif(sequence, motif) that takes 
    a DNA sequence and a short motif as input after user request and returns the starting 
    positions (1-based indexing) of all occurrences of the motif within the sequence. 
    Additionally the positions of the motif within the sequence should be displayed in a plot.

User defined functions: 
    read_sequence, find_motif, plot_positions

Procedure:  
    1. open and run manage_examples.py in Spyder6 to populate directory
    2. Inspect the folder and familiarise yourself with the content
    3. Define function and set parameters
    4. Request user input for sequence file and motif
    4. Display position of motif within sequence in plot

Input: 
    input_file (for test purpose, "gene.fna" from BRCA1_dataset was used, motif was "TCTT")

Output: 
    2b_output_ArianeNeumann.png 

Usage: 
    python3 ArianeNeumannQ2b.py input_file output_file

Version: 1.0 Date 2025-10-30 Author: Ariane Neumann
"""
import matplotlib.pyplot as plt

# ================= Defining the functions =================

# Defining function to read sequence
def read_sequence(filename: str) -> str: # function will return a string
    with open(filename, "r") as a:
        lines = a.readlines() # reading all lines
    seq = "".join(line.strip() # joining all lines, removing whitespace
                  for line in lines 
                  if not line.startswith(">")) # ignoring lines starting with fasta header ">"
    return seq.upper() # returns a sequence in upper case

# Defining function to find motif 1-based. This will start at python index 0
def find_motif(sequence: str, motif: str) -> list: # return function as a list
    motif = motif.upper() # making sure that all letters in motif are written with upper case
    positions = [] # empty list for storing the positions detected within sequence
    for i in range(len(sequence) - len(motif) + 1): 
        if sequence[i:i+len(motif)] == motif: # searching within sequence for motif with specific length
            positions.append(i + 1) # if motif within sequence, it will be appended
    return positions

# Defining function for plot
def plot_positions(positions: list, motif: str):
    plt.figure(figsize=(8, 4))
    plt.bar(positions, [1]*len(positions), color="teal")
    plt.title(f'Occurrences of Motif "{motif}"') # takes user input for motif into title
    plt.xlabel("Position in Sequence")
    plt.ylabel("Occurrence")
    plt.savefig("2b_output_ArianeNeumann.png") # saving plot to user output folder
    plt.show()

# ================= Calling the functions =================

# User input required for sequence
input_file = input("Enter file path for sequence here: ")
sequence = read_sequence(input_file)

# User input required for motif
motif = input("What motif should we look for?: ")
positions = find_motif(sequence, motif)

# Will print all found positions of motif within sequence to screen (no user input required this time)
print(f"Motif positions: {positions}")
if positions:
    plot_positions(positions, motif)
else:
    print("No occurrences of the motif found.")
plot_positions(positions, motif)
