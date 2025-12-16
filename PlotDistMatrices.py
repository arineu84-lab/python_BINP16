#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script name: 
    PlotDistMatrices.py

Description: 
    Based on identity and alignment scores scripted by another student, this 
    code should create 4 dendrograms as well as 2 heatmaps displayed the differences or similarities 
    of the "Romanov" family and other individuals. For training purposes, the 
    script was built with a mock comparison file. In the final runs, this should 
    work with the real input file, which should have the same formatting style.
        
User defined functions: 
    
Procedure:
    1. Import numpy, pandas, matplotlib and scipy
    2. Assess the format of the input file
    3. Identify the two different clusters in the file
    4. Define plot style
    5. Create the actual dendrograms 
 
Input:
    input_file 

Output:
    dendrogram.png   
Usage: 
    python3 PlotDistMatrices.py input_file output_pic1, output_pic2, output_pic3, output_pic4, output_pic5, output_pic6

Version: 1.0
Date 2025-10-23
Author: Ariane Neumann
""" 

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
import seaborn as sb


# Prompt user for input file
input_filename = input("Enter the input filename (e.g., mtDNA_all_vs_all_comparison.tsv): ")

try:
    # Load the input file
    input_data_df = pd.read_csv(input_filename, sep='\t')

    # Rename column for easier access
    input_data_df.rename(columns={'OR Score': 'ORScore'}, inplace=True)

    # Convert IdentityScore to float (strip % if present)
    input_data_df['IdentityScore'] = input_data_df['IdentityScore'].astype(str).str.replace('%', '').astype(float)

    # Extract unique individuals
    individuals = pd.unique(input_data_df[['SampleA', 'SampleB']].values.ravel())
    individuals.sort()

    # Define Romanov identifiers
    romanov_names = [
        'Princess Irene', 'Prince Fred', 'Nicolas II Romanov', 'Alexandra Romanov',
        'Olga Romanov', 'Tatiana Romanov', 'Maria Romanov', 'Alexei Romanov',
        'Suspected body of Anastasia Romanov'
    ]

    # Function to create distance matrix from similarity scores
    def create_distance_matrix(df, score_col):
        matrix = pd.DataFrame(np.ones((len(individuals), len(individuals))) * 100,
                              index=individuals, columns=individuals)
        for _, row in df.iterrows():
            a, b = row['SampleA'], row['SampleB']
            score = row[score_col]
            matrix.loc[a, b] = 100 - score
            matrix.loc[b, a] = 100 - score
        np.fill_diagonal(matrix.values, 0)
        return matrix

    # Function to plot and save dendrogram
    def plot_dendrogram(dist_matrix, title, filename, show_plot):
        condensed = squareform(dist_matrix.values)
        linkage_matrix = linkage(condensed, method='average')

        def label_colors(label):
            return 'teal' if label in romanov_names else 'purple'

        plt.figure(figsize=(10, 10))
        dendro = dendrogram(
            linkage_matrix,
            labels=dist_matrix.index.tolist(),
            leaf_font_size=10,
            leaf_rotation=0,
            orientation='right',
            link_color_func=lambda k: 'black',
            color_threshold=0)

        ax = plt.gca()
        ylbls = ax.get_yticklabels()
        for lbl in ylbls:
            lbl.set_color(label_colors(lbl.get_text()))
            fig = plt.figure(figsize=(12, 10), facecolor='white')  # Set figure background

        plt.savefig(filename, facecolor='white')  # Ensure saved image has white background
        plt.title(title, fontsize=14)
        plt.xlabel("Distance")
        plt.tight_layout()
        plt.savefig(filename)
        if show_plot:
            plt.show()
        plt.close()

    # User option to show plots
    show_plots = False

    # Generate and save dendrograms
    plot_dendrogram(create_distance_matrix(input_data_df, 'ORScore'),
                    "Alignment Score", "dendrogram_alignment.png", show_plots)

    plot_dendrogram(create_distance_matrix(input_data_df, 'IdentityScore'),
                    "Identity Score", "dendrogram_identity.png", show_plots)

    print("Dendrograms saved as PNG files.")

except FileNotFoundError:
    print(f"Error: The file '{input_filename}' was not found.")
except KeyError as e:
    print(f"Error: Missing expected column {e} in the input file.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")
    
# Load the genetic distance matrix
genetic_matrix = pd.read_csv('genetic_distance_matrix.csv', sep='\t', index_col=0)

# Define Romanov identifiers
romanov_cluster = [ 'Princess Irene', 'Prince Fred', 'Nicolas II Romanov', 'Alexandra Romanov',
    'Olga Romanov', 'Tatiana Romanov', 'Maria Romanov', 'Alexei Romanov',
    'Suspected body of Anastasia Romanov']

# Assign labels
labels = ['Romanov' if name in romanov_cluster else 'non-Romanov' for name in genetic_matrix.index]

# Map labels to colours: purple for Romanov, teal for non-Romanov
color_code = {'Romanov': 'purple', 'non-Romanov': '#d0f0ff'}
row_colors = pd.Series(labels, index=genetic_matrix.index).map(color_code)

# Create clustered heatmap with row colors
sb.set(font_scale=1.1)
sb.clustermap(genetic_matrix, cmap='viridis', linewidths=0.5, figsize=(12, 10), row_colors=row_colors)

# Show the plot
plt.show()