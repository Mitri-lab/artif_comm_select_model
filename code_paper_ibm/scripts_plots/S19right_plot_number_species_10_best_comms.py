# -*- coding: utf-8 -*-
"""
Created on Sun Jul 27 16:39:26 2025

@author: pablo
"""
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import numpy as np

# --- CONFIGURATION ---
PATHS_BY_CONDITION = {
    '10_per_sp': 'C:/Users/pablo/Desktop/data_mitri/211005_all_communities',
    '150_total': 'C:/Users/pablo/Desktop/jul2025'
}


# Simulation parameters
conditions = ['10_per_sp', '150_total']
seeds = range(22, 27) # Original seed numbers used to find folders
output_filename = 'species_richness_comparison.pdf' # Name for the output PDF file
# --------------------


# --- DATA PROCESSING ---
print("Selecting top 10 communities from each run...")
top_communities_data = []

for condition in conditions:
    for seed in seeds:
        current_base_path = PATHS_BY_CONDITION[condition]
        file_path = os.path.join(current_base_path, str(seed), 'df.csv')
        species_set_num = seed - 21 # Convert seed 22-26 to set 1-5
        
        try:
            temp_df = pd.read_csv(file_path)
            top_10_communities = temp_df.sort_values(by='deg_score', ascending=False).head(10)
            
            for index, row in top_10_communities.iterrows():
                composition_str = row['com']
                num_species = len(eval(composition_str))
                
                top_communities_data.append({
                    'condition': condition,
                    'species_set': species_set_num,
                    'num_species': num_species
                })

        except FileNotFoundError:
            print(f"  -> WARNING: File not found for {condition}, Seed {seed}")

df_plot = pd.DataFrame(top_communities_data)
print("Processing complete.")
# --------------------


# --- STATISTICAL ANALYSIS ---
print("\n--- Statistical Test (Mann-Whitney U) ---")

species_counts_10_sp = df_plot[df_plot['condition'] == '10_per_sp']['num_species']
species_counts_150_total = df_plot[df_plot['condition'] == '150_total']['num_species']

stat, p_value = mannwhitneyu(species_counts_10_sp, species_counts_150_total, alternative='two-sided')

print(f"P-value: {p_value:.6f}")
if p_value < 0.05:
    print("The difference in the number of species is statistically significant.")
else:
    print("There is no statistically significant difference in the number of species.")
# --------------------


# --- VISUALIZATION ---
print("\nGenerating plot...")
plt.figure(figsize=(10, 8))

# Use a boxplot to show the distribution summary
sns.boxplot(
    data=df_plot, x='condition', y='num_species',
    color="#cccccc",
    width=0.4,
    showfliers=False,
    boxprops=dict(alpha=.9)
)

# Use a stripplot to show all individual data points
sns.stripplot(
    data=df_plot, 
    x='condition', 
    y='num_species', 
    hue='species_set', 
    palette='colorblind',
    jitter=True,
    dodge=False,
    alpha=0.7,
    s=10 
)

# Add titles and labels
plt.title('Number of Species in Top 10 Degrading Communities', fontsize=18)
plt.ylabel('Number of Species per Community', fontsize=16)
plt.xlabel('Seeding Condition', fontsize=16) 

# --- AESTHETIC MODIFICATIONS ---
# Set axis ticks and labels with larger fonts
plt.xticks(fontsize=14) 
plt.yticks(range(1, 16), fontsize=14)
plt.ylim(bottom=0, top=17)

# Add the p-value text to the plot
plt.text(0.5, 16.1, f"Mann-Whitney U test\np-value = {p_value:.4f}",
         ha='center', va='top', fontsize=12, style='italic',
         bbox=dict(facecolor='white', alpha=0.5, edgecolor='none', pad=2))
# ---------------------

plt.legend(title='Species Set', bbox_to_anchor=(1.02, 1), loc='upper left')
plt.tight_layout(rect=[0, 0, 0.9, 1]) 

# --- SAVE THE PLOT ---
save_path = os.path.join(PATHS_BY_CONDITION['150_total'], output_filename)
plt.savefig(save_path, bbox_inches='tight')
print(f"\nPlot saved to: {save_path}")
# ---------------------

# Show the plot on the screen
plt.show()

print("\nAnalysis finished.")
# --------------------