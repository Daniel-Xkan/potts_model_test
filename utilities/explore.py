import re
import numpy as np
import pandas as pd
import seaborn as sns
import utilities.functions as functions

def calculate_residue_e(pos, seq, J_dict, min_position, max_position):
    energy = 0
    aa1 = seq[pos - min_position]  # Access the amino acid at the given position
    for pos2 in range(min_position, max_position + 1):
        if pos2 != pos:  # Ensure we are not calculating self-interaction
            aa2 = seq[pos2 - min_position]  # Access the amino acid at pos2
            energy += J_dict.get((pos, pos2, aa1, aa2), 0)
    return energy

def calculate_residue_de12(pos, pair1,pair2,seq, J_dict, min_position, max_position):
    aa = seq[pos - min_position]  # Access the amino acid at the given position
    wt1,pos1,mt1 = functions.split_pair(pair1)
    wt2,pos2,mt2 = functions.split_pair(pair2)
    energy_wt = J_dict.get((pos, pos1, aa, wt1), 0) + J_dict.get((pos, pos2, aa, wt2), 0)
    energy_mt = J_dict.get((pos, pos1, aa, mt1), 0) + J_dict.get((pos, pos2, aa, mt2), 0)

    energy_diff = energy_wt - energy_mt
    return energy_diff
    

def calculate_residue_relative_e(consensus_seq,pos, seq, J_dict, min_position, max_position):
    e_consensus = calculate_residue_e(pos, consensus_seq, J_dict, min_position, max_position)
    e_seq = calculate_residue_e(pos, seq, J_dict, min_position, max_position)
    relative_e = e_seq - e_consensus
    return relative_e

#E residue with double mutaiton -E residue without double mutation
def calculate_de_residue(pos, sequence, double_mutation, J_dict, redux):
    single_mutation1, single_mutation2  = functions.split_pairs(double_mutation)
    sm1_reduced = functions.unreduced_to_reduced(redux,single_mutation1)
    sm2_reduced = functions.unreduced_to_reduced(redux,single_mutation2)
    wt1,pos1,mt1 = functions.split_pair(sm1_reduced)
    wt2,pos2,mt2 = functions.split_pair(sm2_reduced)
    seq_wt = list(sequence)
    seq_mt = list(sequence)

    seq_wt[pos1 - 1] = wt1
    seq_wt[pos2 - 1] = wt2

    seq_mt[pos1 - 1] = mt1
    seq_mt[pos2 - 1] = mt2

    seq_wt = "".join(seq_wt)
    seq_mt = "".join(seq_mt)
    # print("seq_wt:", seq_wt)
    # print("seq_mt:", seq_mt)
    # if seq_wt != seq_mt:
    #     print("Differences between seq_wt and seq_mt:")
    #     for i, (a, b) in enumerate(zip(seq_wt, seq_mt), start=1):
    #         if a != b:
    #             print(f"Position {i}: {a} -> {b}")
    # else:
    #     print("No differences between seq_wt and seq_mt.")
    e_wt = calculate_residue_e(pos, seq_wt, J_dict, 1, len(seq_wt))
    # print(e_wt)
    e_mt = calculate_residue_e(pos, seq_mt, J_dict, 1, len(seq_mt))
    # print(e_mt)

    energy_diff = e_wt - e_mt
    return energy_diff


import matplotlib.pyplot as plt
def draw_heatmap(consensus_seq, sequences, J_dict, min_position, max_position):
    heatmap_data = []
    for seq in sequences:
        residue_energies = [
            calculate_residue_relative_e(consensus_seq, pos, seq, J_dict, min_position, max_position)
            for pos in range(min_position, max_position + 1)
        ]
        heatmap_data.append(residue_energies)

    heatmap_data = np.array(heatmap_data)
    plt.figure(figsize=(10, 8))
    ax = sns.heatmap(
        heatmap_data,
        cmap="bwr",
        center=0,
        cbar_kws={'label': 'Residue Energy'},
        xticklabels=False,
        yticklabels=False
    )

    # X-axis: mark positions 10, 20, 30, ... (sequence positions)
    ax.set_xticks([0] + list(np.arange(9, (max_position - min_position + 1), 10)))
    ax.set_xticklabels([1] + list(range(10, max_position + 1, 10)), fontsize=6)
    # return
    # tick_positions = np.arange(9, (max_position - min_position + 1), 10)
    # tick_labels = list(range(10, max_position + 1, 10))
    # ax.set_xticks(tick_positions)
    # ax.set_xticklabels(tick_labels, fontsize=6)

    plt.xlabel("Position")
    plt.ylabel("Sequence Index")
    plt.title("Residue Energy Heatmap")
    plt.show()

def draw_heatmap2(vector, min_position, max_position):
    arr = np.array(vector)
    plt.figure(figsize=(10, 8))
    ax = sns.heatmap(
        arr,
        cmap="bwr",
        center=0,
        vmin=-7,
        vmax=5,
        cbar_kws={'label': 'Value'},
        xticklabels=False,
        yticklabels=False
    )
    ax.set_xticks([0] + list(np.arange(9, (max_position - min_position + 1), 10)))
    ax.set_xticklabels([1] + list(range(10, max_position + 1, 10)), fontsize=6)
    plt.xlabel("Column")
    plt.ylabel("Row")
    plt.title("Heatmap")
    plt.show()
