import re
import numpy as np
import pandas as pd
import math

def load_J_dict(j_file, min_position,max_position):
    J_dict = {}
    J = np.load(j_file)

    row = 0

    for pos1 in range(min_position, max_position + 1):
        for pos2 in range(pos1 + 1, max_position + 1):
            for i, aa1 in enumerate(['A', 'B', 'C', 'D']):
                for j, aa2 in enumerate(['A', 'B', 'C', 'D']):
                    col = i * 4 + j
                    J_dict[(pos1, pos2, aa1, aa2)] = J[row, col]
                    J_dict[(pos2, pos1, aa2, aa1)] = J[row, col]
            row += 1

    return J_dict

def split_pair(pair):
    wildtype = pair[0]
    position = int(pair[1:-1])  # int
    mutate = pair[-1]
    return wildtype, position, mutate


def calculate_delta_e_double(mutation_pair1, mutation_pair2, seq, J_dict,min_position,max_position):
    old_amino_acid1, pos1, new_amino_acid1 = split_pair(mutation_pair1)
    old_amino_acid2, pos2 , new_amino_acid2 = split_pair(mutation_pair2)

    # Check if positions are already mutated
    current_aa1 = seq[pos1 - min_position]
    current_aa2 = seq[pos2 - min_position]

    #if the positions are mutated, swap back to the original amino acid for calculation only
    if current_aa1 != old_amino_acid1 or current_aa2 != old_amino_acid2:
        current_aa1 = old_amino_acid1
        current_aa2 = old_amino_acid2

    energy_old = 0
    energy_new = 0

    # old energy
    for other_pos in range(min_position, max_position+ 1):
        other_aa = seq[other_pos - min_position]
        if other_pos == pos1:
            continue
        if other_pos == pos2:
            energy_old += J_dict.get((pos1, pos2, current_aa1, current_aa2), 0)
        else:
            energy_old += J_dict.get((pos1, other_pos, old_amino_acid1, other_aa), 0)

    for other_pos in range(min_position, max_position + 1):
        other_aa = seq[other_pos - min_position]
        if other_pos == pos2:
            continue
        if other_pos == pos1:
            continue
        else:
            energy_old += J_dict.get((pos2, other_pos, old_amino_acid2, other_aa), 0)

    # new energy
    for other_pos in range(min_position, max_position + 1):
        other_aa = seq[other_pos - min_position]
        if other_pos == pos1:
            continue
        if other_pos == pos2:
            energy_new += J_dict.get((pos1, pos2, new_amino_acid1, new_amino_acid2), 0)
        else:
            energy_new += J_dict.get((pos1, other_pos, new_amino_acid1, other_aa), 0)

    for other_pos in range(min_position, max_position + 1):
        other_aa = seq[other_pos - min_position]
        if other_pos == pos2:
            continue
        if other_pos == pos1:
            continue
        else:
            energy_new += J_dict.get((pos2, other_pos, new_amino_acid2, other_aa), 0)

    de12 = energy_old - energy_new
    return de12


def calculate_double_mutant_probablity(seq, mutation_pair1, mutation_pair2, J_dict, min_position, max_position):

    de_DMC = calculate_delta_e_double(mutation_pair1, mutation_pair2, seq, J_dict, min_position, max_position)
    de_double_sum = 1
    for aa1 in ['A', 'B', 'C', 'D']:
        for aa2 in ['A', 'B', 'C', 'D']:
            if aa1 == mutation_pair1[0] and aa2 == mutation_pair2[0]:
                continue
            de_DMC_alt = calculate_delta_e_double(mutation_pair1[:-1] + aa1, mutation_pair2[:-1] + aa2, seq, J_dict, min_position, max_position)
            de_DMC_alt_e = math.exp(de_DMC_alt) if de_DMC_alt is not None else 0
            de_double_sum += de_DMC_alt_e
    return (math.exp(de_DMC))/de_double_sum