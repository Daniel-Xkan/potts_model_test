import re
import numpy as np

def import_test():
    return 'import succeed'

# Read the reduced dictionary from a file and apply an offset to positions
def get_redu_dict(redux_file, offset):
    redux_dict = {}
    with open(redux_file, 'r') as file:
        for line in file:
            parts = line.strip().split()
            position = int(parts[0]) + offset
            groups = parts[1:]
            for i, group in enumerate(groups):
                key = (position, chr(65 + i))  # 'A', 'B', 'C', 'D', etc.
                redux_dict[key] = list(group.replace('-', ''))  # Remove '-' from the group
    return redux_dict
#get the reduced pair in format like 'C140D' translate back to unreduced with rules that 1. look up in dictionary ('C', 140)get the last letter in list, and samely, look up in dictionary ('D', 140), get the last letter in list, and combine the new translated unreduced two letters put them back to new pairs, like 'G140S'

#usage: dict, C140D, [unreduced seq]
def reduced_to_unreduced(redux_dict, reduced_pair,in_unre_seq):
    # Parse the reduced pair
    wildtype, pos, mutate = reduced_pair[0], int(reduced_pair[1:-1]), reduced_pair[-1]
    # Look up the unreduced wildtype and mutate in the dictionary
    unreduced_wildtypes = [aa for aa in redux_dict.get((pos, wildtype), ['-'])]
    unreduced_mutates = [aa for aa in redux_dict.get((pos, mutate), ['-'])]
    # print( 'unreduced_wildtypes',unreduced_wildtypes)
    # print( 'unreduced_mutates ', unreduced_mutates)
    
    # If the most frequent amino acid is needed, calculate it
    most_freq_wildtype = max(unreduced_wildtypes, key=lambda aa: calculate_freq(aa, pos, in_unre_seq))
    most_freq_mutate = sorted(set(unreduced_mutates), key=lambda aa: calculate_freq(aa, pos, in_unre_seq), reverse=True)[0]

    
    # Combine the new translated unreduced letters with the position
    unreduced_pair = f"{most_freq_wildtype}{pos}{most_freq_mutate}"
    return unreduced_pair
#usage: dict, G140S
def unreduced_to_reduced(redux_dict, unreduced_pair):
    # Parse the unreduced pair
    wildtype, pos, mutate = unreduced_pair[0], int(unreduced_pair[1:-1]), unreduced_pair[-1]
    # Look up the reduced wildtype and mutate in the dictionary
    reduced_wildtype = next((key[1] for key, values in redux_dict.items() if key[0] == pos and wildtype in values), '-')
    reduced_mutate = next((key[1] for key, values in redux_dict.items() if key[0] == pos and mutate in values), '-')
    # Combine the reduced letters with the position
    reduced_pair = f"{reduced_wildtype}{pos}{reduced_mutate}"
    return reduced_pair

# Calculate the frequency of an amino acid at a given position in a list of sequences
def calculate_freq(aa, pos, seq_list):
    count = sum(1 for seq in seq_list if seq[pos - 1] == aa)
    frequency = count / len(seq_list)
    return frequency

# Split a string of two pairs separated by a hyphen into a tuple
def split_pairs(pair_string):
    return tuple(pair.strip() for pair in pair_string.split('-') if pair.strip())

# Join two pairs into a single string with a hyphen
def join_pairs(pair1, pair2):
    return f"{pair1}-{pair2}"

# Split a single pair like 'C140D' into its components
def split_pair(pair):
    wildtype = pair[0]
    position = int(pair[1:-1])
    mutate = pair[-1]
    return wildtype, position, mutate

# usage: kn.her2.all, C140D, True
#        kn.her2.all, G140S, False
def DE_dict(kn_file):
    de_dict = {}
    with open(kn_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            columns = re.split(r'[,\t\s]+', line)
            pairs = columns[0]
            pair1 = split_pairs(pairs)[0]  # first pair
            pair2 = split_pairs(pairs)[1]  # second pair

            ddE = float(columns[1])
            delta_e1 = float(columns[3])
            delta_e2 = float(columns[4])
            de_dict[pair1] = delta_e1
            de_dict[pair2] = delta_e2
    return de_dict

# get de from kn.her.all dictionary
def get_DE(de_dict, pair, reduced=True, redux_dict=None):
    if reduced == False:
        if redux_dict is None:
            raise ValueError("redux_dict must be provided when reduced is False")
        pair_reduced = unreduced_to_reduced(redux_dict, pair)
    else:
        pair_reduced = pair
    return de_dict.get(pair_reduced)

def calculate_delta_e(mutation_pair, seq, J_dict):
    old_amino_acid, position, new_amino_acid = split_pair(mutation_pair)
    max_position = len(seq)
    
    # E(old_amino_acid)
    energy_old = 0
    for other_pos in range(1, max_position + 1):
        if other_pos == position:
            continue
        other_aa = seq[other_pos - 1]  # Access the amino acid at other_pos
        energy_old += J_dict.get((position, other_pos, old_amino_acid, other_aa), 0)

    # E(new_amino_acid)
    energy_new = 0
    for other_pos in range(1, max_position + 1):
        if other_pos == position:
            continue
        other_aa = seq[other_pos - 1]  # Access the amino acid at other_pos
        energy_new += J_dict.get((position, other_pos, new_amino_acid, other_aa), 0)

    delta_e = energy_old - energy_new
    return delta_e

def load_J_dict(j_file, max_position):
    J_dict = {}
    J = np.load(j_file)

    row = 0

    for pos1 in range(1, max_position + 1):
        for pos2 in range(pos1 + 1, max_position + 1):
            for i, aa1 in enumerate(['A', 'B', 'C', 'D']):
                for j, aa2 in enumerate(['A', 'B', 'C', 'D']):
                    col = i * 4 + j
                    J_dict[(pos1, pos2, aa1, aa2)] = J[row, col]
                    J_dict[(pos2, pos1, aa2, aa1)] = J[row, col]
            row += 1

    return J_dict

def read_seq(seq_file):
    seq_list = []
    with open(seq_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                seq_list.append(line)
    return seq_list