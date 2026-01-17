# -*- coding: utf-8 -*-
"""
Created on Sun Jul 27 16:10:42 2025

@author: pablo
"""
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import jaccard
import numpy as np

# --- CONFIGURATION ---
PATHS_BY_CONDITION = {
    '10_per_sp': 'C:/Users/pablo/Desktop/data_mitri/211005_all_communities',
    '150_total': 'C:/Users/pablo/Desktop/jul2025'
}


# Simulation parameters
conditions = ['10_per_sp', '150_total']
seeds = range(22, 27) # Original seed numbers used to find folders
total_species = 15 # The total number of possible species (0-14)
output_filename = 'community_comparison_plot.pdf' # Name for the output PDF file
# --------------------


# --- DATA PROCESSING ---
print("Searching for the best community for each condition...")
best_communities = {}

for condition in conditions:
    for seed in seeds:
        current_base_path = PATHS_BY_CONDITION[condition]
        file_path = os.path.join(current_base_path, str(seed), 'df.csv')
        
        try:
            temp_df = pd.read_csv(file_path)
            best_community_row = temp_df.sort_values(by='deg_score', ascending=False).iloc[0]
            composition_str = best_community_row['com']
            composition_tuple = eval(composition_str)
            best_communities[(condition, seed)] = set(composition_tuple)
            print(f"  -> Best community for '{condition}' (Seed {seed}): {set(composition_tuple)}")
        except FileNotFoundError:
            print(f"  -> WARNING: File not found for {condition}, Seed {seed}")
            best_communities[(condition, seed)] = set()

print("\nProcessing complete.")
# --------------------


# --- COMPARISON AND VISUALIZATION ---
print("\n--- Generating final plot ---")

# Create a figure with 5 subplots (one for each species set) arranged vertically
fig, axes = plt.subplots(nrows=len(seeds), ncols=1, figsize=(16, 12))
fig.suptitle('Comparison of Best Communities Across Species Sets', fontsize=22)

for i, seed in enumerate(seeds):
    species_set_num = i + 1 # Display seeds 22-26 as 1-5
    ax = axes[i] # Select the subplot for the current species set
    
    print(f"\n--- Analysis for Species Set: {species_set_num} (Original Seed: {seed}) ---")
    
    # Get the two sets of species to compare
    composition_10_sp = best_communities.get(('10_per_sp', seed), set())
    composition_150_total = best_communities.get(('150_total', seed), set())
    
    # Create binary vectors (presence/absence)
    vector_10_sp = [1 if sp in composition_10_sp else 0 for sp in range(total_species)]
    vector_150_total = [1 if sp in composition_150_total else 0 for sp in range(total_species)]
    
    # Calculate Jaccard Distance
    distance = jaccard(vector_10_sp, vector_150_total)
    print(f"Composition '10_per_sp':    {composition_10_sp}")
    print(f"Composition '150_total':  {composition_150_total}")
    print(f"Jaccard Distance: {distance:.4f}")
    
    # Prepare data for the heatmap
    heatmap_data = pd.DataFrame(
        [vector_10_sp, vector_150_total],
        columns=[j + 1 for j in range(total_species)], # Use 1-15 for species labels
        index=['10_per_sp', '150_total']
    )
    
    # Generate the heatmap on the specific subplot axis
    sns.heatmap(
        heatmap_data, 
        annot=True,
        cmap=['#E0E0E0', '#616161'],
        linewidths=.5,
        linecolor='black',
        cbar=False,
        ax=ax,
        annot_kws={"size": 14}
    )
    
    # Set the title for the subplot
    ax.set_title(f'Species Set: {species_set_num}   |   Jaccard Distance: {distance:.4f}', fontsize=18)
    
    # Set labels with larger fonts
    ax.set_xlabel('Species', fontsize=14)
    ax.set_ylabel('Seeding Condition', fontsize=14)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12, rotation=0)

# Adjust layout to prevent titles/labels from overlapping
plt.tight_layout(rect=[0, 0, 1, 0.96]) 

# --- SAVE THE PLOT ---
# Construct the full path to save the file
save_path = os.path.join(PATHS_BY_CONDITION['150_total'], output_filename)
plt.savefig(save_path, bbox_inches='tight') # Save the figure to a PDF
print(f"\nPlot saved to: {save_path}")
# ---------------------

# Show the plot on the screen
plt.show()

print("\nAnalysis finished.")
# --------------------