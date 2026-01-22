import re
import numpy as np
import pandas as pd

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
def reduced_to_unreduced(redux_dict, reduced_pair,in_unre_seq, start_index=1):
    # Parse the reduced pair
    wildtype, pos, mutate = reduced_pair[0], int(reduced_pair[1:-1]), reduced_pair[-1]
    # Look up the unreduced wildtype and mutate in the dictionary
    unreduced_wildtypes = [aa for aa in redux_dict.get((pos, wildtype), ['-'])]
    unreduced_mutates = [aa for aa in redux_dict.get((pos, mutate), ['-'])]
    # print( 'unreduced_wildtypes',unreduced_wildtypes)
    # print( 'unreduced_mutates ', unreduced_mutates)
    
    # If the most frequent amino acid is needed, calculate it
    most_freq_wildtype = max(unreduced_wildtypes, key=lambda aa: calculate_freq(aa, pos, in_unre_seq,start_index))
    most_freq_mutate = sorted(set(unreduced_mutates), key=lambda aa: calculate_freq(aa, pos, in_unre_seq,start_index), reverse=True)[0]

    
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
def calculate_freq(aa, pos, seq_list,start_index=1):
    count = sum(1 for seq in seq_list if seq[pos - start_index] == aa)
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
    position = int(pair[1:-1])  # Ensure position is an integer
    mutate = pair[-1]
    return wildtype, position, mutate

def split_pair_multiple_mt(pair):
    m = re.match(r'^([A-Za-z])(\d+)(.+)$', pair)
    if not m:
        raise ValueError(f"Invalid mutation format: {pair}")
    wildtype, position, rest = m.group(1), int(m.group(2)), m.group(3)
    # split mutations by '/' if present, otherwise single mutation
    mutations = [seg for seg in rest.split('/') if seg]
    # return a list of (wildtype, position, mutation) tuples
    return (position, wildtype,mutations)
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

def DDE_dict(kn_file):
    dde_dict = {}
    with open(kn_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            columns = re.split(r'[,\t\s]+', line)
            pairs = columns[0]
            dde = float(columns[1])
            dde_dict[pairs] = dde
    return dde_dict
# get de from kn.her.all dictionary

#input: kn.her2.all
#output dde sorted list of tuples (pos1, pos2, ddE)
def kn_file(file_path):
    delta_delta_e_per_pair = []
    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            columns = re.split(r'[,\t\s]+', line)
            mutations = columns[0]
            ddE = float(columns[1])
            positions = [int(pos) for pos in re.findall(r'\d+', mutations)]
            if len(positions) != 2:
                continue

            pos1, pos2 = positions
            delta_delta_e_per_pair.append((pos1, pos2, ddE))

    # Sort by ddE values
    delta_delta_e_per_pair_sorted = sorted(delta_delta_e_per_pair, key=lambda x: x[-1], reverse=True)
    return delta_delta_e_per_pair_sorted

def get_pos1_pos2(two_pairs):
    pair1 = split_pairs(two_pairs)[0]
    pair1_pos = get_pos(pair1)
    pair2 = split_pairs(two_pairs)[1]
    pair2_pos = get_pos(pair2)
    return int(pair1_pos),int(pair2_pos)

def get_DE(de_dict, pair, reduced=True, redux_dict=None):
    if reduced == False:
        if redux_dict is None:
            raise ValueError("redux_dict must be provided when reduced is False")
        pair_reduced = unreduced_to_reduced(redux_dict, pair)
    else:
        pair_reduced = pair
    return de_dict.get(pair_reduced)

def get_seq_freq(mut_list,seq_list,consensus_seq):
    seq_count = 0
    ref_seq = generate_sequence_from_mutations(consensus_seq, mut_list)
    for seq in seq_list:
        if seq == ref_seq:
            seq_count += 1
    frequency = seq_count / len(seq_list)
    return frequency

def calculate_e(seq, J_dict, min_position, max_position):
    energy = 0
    for pos1 in range(min_position, max_position + 1):
        aa1 = seq[pos1 - min_position]  # Access the amino acid at pos1
        for pos2 in range(pos1 + 1, max_position + 1):
            aa2 = seq[pos2 - min_position]  # Access the amino acid at pos2
            energy += J_dict.get((pos1, pos2, aa1, aa2), 0)
            # print(f"pos1: {pos1}, aa1: {aa1}, pos2: {pos2}, aa2: {aa2}, interaction: {J_dict.get((pos1, pos2, aa1, aa2), 0)}")
    return energy

def calculate_delta_e(mutation_pair, seq, J_dict, min_position, max_position):
    old_amino_acid, position, new_amino_acid = split_pair(mutation_pair)
    # E(old_amino_acid)
    energy_old = 0
    for other_pos in range(min_position, max_position + 1):
        if other_pos == position:
            continue
        other_aa = seq[other_pos - min_position]  # Access the amino acid at other_pos
        energy_old += J_dict.get((position, other_pos, old_amino_acid, other_aa), 0)

    # E(new_amino_acid)
    energy_new = 0
    for other_pos in range(min_position, max_position + 1):
        if other_pos == position:
            continue
        other_aa = seq[other_pos - min_position]  # Access the amino acid at other_pos
        energy_new += J_dict.get((position, other_pos, new_amino_acid, other_aa), 0)

    delta_e = energy_old - energy_new
    return delta_e

def calculate_delta_e_double(mutation_pair1, mutation_pair2, seq, J_dict,min_position,max_position):
    old_amino_acid1, pos1, new_amino_acid1 = split_pair(mutation_pair1)
    old_amino_acid2, pos2 , new_amino_acid2 = split_pair(mutation_pair2)

    # Check if positions are already mutated
    current_aa1 = seq[pos1 - min_position]
    current_aa2 = seq[pos2 - min_position]

    if current_aa1 != old_amino_acid1 or current_aa2 != old_amino_acid2:
        return None

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

def calculate_delta_delta_e(mutation_pair1, mutation_pair2, seq, J_dict,min_position,max_position):
    de1 = calculate_delta_e(mutation_pair1, seq, J_dict, min_position, max_position)
    de2 = calculate_delta_e(mutation_pair2, seq, J_dict, min_position, max_position)
    de12 = calculate_delta_e_double(mutation_pair1, mutation_pair2, seq, J_dict, min_position, max_position)

    if de1 is None or de2 is None or de12 is None:
        return None

    dde =  de12- de1 - de2
    return de1,de2,de12,dde

def calculate_with_double_mutation_dde(mutation_pair1, mutation_pair2, seq, J_dict,min_position,max_position):
    wt1,pos1,mut1 = split_pair(mutation_pair1)
    wt2,pos2,mut2 = split_pair(mutation_pair2)

    if seq[pos1 - min_position] != mut1 or seq[pos2 - min_position] != mut2:
        return None
    seq_list = list(seq)
    seq_list[pos1 - min_position] = wt1
    seq_list[pos2 - min_position] = wt2
    seq_wt = ''.join(seq_list)
    
    de1 = calculate_delta_e(mutation_pair1, seq_wt, J_dict, min_position, max_position)
    de2 = calculate_delta_e(mutation_pair2, seq_wt, J_dict, min_position, max_position)
    de12 = calculate_delta_e_double(mutation_pair1, mutation_pair2, seq_wt, J_dict, min_position, max_position)


    dde =  de12- de1 - de2
    return de1,de2,de12,dde

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

def read_seq(seq_file):
    seq_list = []
    with open(seq_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                seq_list.append(line)
    return seq_list

def read_reduced_MSA(msa_file):
    msa_list = []
    with open(msa_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                msa_list.append(line)
    return msa_list

def calculate_average_DE1_DE2_DDE(mut1, mut2, MSA_seq,J_dict,min_position,max_position):
    wt1, position1, mt1 = split_pair(mut1)
    wt2, position2, mt2 = split_pair(mut2)
    valid_count = 0
    DE1_total = 0
    DE2_total = 0
    DE12_total = 0
    for seq in MSA_seq:
        if seq[position1-min_position] != wt1 or seq[position2-min_position] != wt2:
            # print(f"Expected: wt1={wt1}, wt2={wt2}, Found: seq[position1]={seq[position1-1]}, seq[position2]={seq[position2-1]}")
            continue
        valid_count += 1
        DE1_total += calculate_delta_e(mut1, seq, J_dict, min_position, max_position)
        DE2_total += calculate_delta_e(mut2, seq, J_dict, min_position, max_position)
        DE12_total += calculate_delta_e_double (mut1, mut2, seq, J_dict,min_position, max_position)
    print(f"Valid sequences count: {valid_count}")
    if valid_count == 0:
        return None, None, None
    return DE1_total/valid_count, DE2_total/valid_count, DE12_total/valid_count


def pairs_to_reduced(redux_dict, pairs):
    """
    Convert mutation pairs from unreduced to reduced format.
    
    Args:
        redux_dict: Dictionary mapping (position, reduced_letter) to list of amino acids
        pairs: List of pairs in format 'X###Y-A###B'
    
    Returns:
        List of reduced format pairs
    """
    reduced_pairs = []
    for pair in pairs:
        # Split the pair into two mutations
        mut1, mut2 = pair.split('-')
        # Convert each mutation to reduced format
        reduced_mut1 = unreduced_to_reduced(redux_dict, mut1)
        reduced_mut2 = unreduced_to_reduced(redux_dict, mut2)
        # Combine back into pair format
        reduced_pair = f"{reduced_mut1}-{reduced_mut2}"
        reduced_pairs.append(reduced_pair)
    return reduced_pairs

def analyze_sequences_mutations(consensus_file, sequences_file):
    """
    Read consensus and experimental sequences, then create a DataFrame with mutations.
    
    Args:
        consensus_file: Path to consensus sequence file
        sequences_file: Path to experimental sequences file
    
    Returns:
        pandas.DataFrame with columns: Sequence, Mutations, Mutations_count
    """
    
    # Read the consensus sequence
    with open(consensus_file, 'r') as f:
        consensus_sequence = f.read().strip()
    
    # Read the experimental sequences
    sequences = []
    with open(sequences_file, 'r') as f:
        for line in f:
            sequences.append(line.strip())
    
    # Create a DataFrame to store sequences and their mutations
    sequence_df = pd.DataFrame({'Sequence': sequences})
    
    # Function to identify mutations compared to the consensus sequence
    def find_mutations(sequence, consensus):
        mutations = []
        for i, (seq_residue, cons_residue) in enumerate(zip(sequence, consensus), start=1):
            if seq_residue != cons_residue:
                mutations.append(f"{cons_residue}{i}{seq_residue}")
        return mutations
    
    # Add columns for mutations and mutation count
    sequence_df['Mutations'] = sequence_df['Sequence'].apply(lambda seq: find_mutations(seq, consensus_sequence))
    sequence_df['Mutations_count'] = sequence_df['Mutations'].apply(len)
    
    return sequence_df

def reduced_to_unreduced_list(redux_dict, reduced_pairs, in_unre_seq, start_index=1):
    unreduced_pairs = []
    for reduced_pair in reduced_pairs:
        unreduced_pair = reduced_to_unreduced(redux_dict, reduced_pair, in_unre_seq, start_index)
        unreduced_pairs.append(unreduced_pair)
    return unreduced_pairs

def unreduced_to_reduced_pair(redux_dict,unreduced_pair):
    wildtype1, pos1, mutate1 = split_pair(unreduced_pair.split('-')[0])
    wildtype2, pos2, mutate2 = split_pair(unreduced_pair.split('-')[1])

    reduced_wildtype1 = next((key[1] for key, values in redux_dict.items() if key[0] == pos1 and wildtype1 in values), '-')
    reduced_mutate1 = next((key[1] for key, values in redux_dict.items() if key[0] == pos1 and mutate1 in values), '-')

    reduced_wildtype2 = next((key[1] for key, values in redux_dict.items() if key[0] == pos2 and wildtype2 in values), '-')
    reduced_mutate2 = next((key[1] for key, values in redux_dict.items() if key[0] == pos2 and mutate2 in values), '-')

    reduced_pair = f"{reduced_wildtype1}{pos1}{reduced_mutate1}-{reduced_wildtype2}{pos2}{reduced_mutate2}"
    return reduced_pair

def get_kn_file_dict(kn_file):
    de_dict = {}
    with open(kn_file, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) < 5:
                continue
            pair1, pair2 = split_pairs(parts[0])
            de1 = float(parts[3])
            de2 = float(parts[4])
            de_dict[pair1] = de1
            de_dict[pair2] = de2
    return de_dict

def get_wt(single_mutation):
    wt = single_mutation[:1]
    return wt

def get_mt(single_mutation):
    mt = single_mutation[-1]
    return mt

def get_pos(single_mutation):
    pos = int(single_mutation[1:-1])
    return pos

def read_list_from_file(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file if line.strip()]

##### Get mutations from consensus and experimental sequences and sequence generation######
def get_mutations_from_sequence(consensus_seq, experimental_seq):
    mutations = []
    for i, (cons_residue, exp_residue) in enumerate(zip(consensus_seq, experimental_seq), start=1):
        if cons_residue != exp_residue:
            mutations.append(f"{cons_residue}{i}{exp_residue}")
    return mutations

def generate_sequence_from_mutations(consensus_seq, mutations):
    seq_list = list(consensus_seq)
    for mutation in mutations:
        wt = mutation[0]
        pos = int(mutation[1:-1])
        mt = mutation[-1]
        if seq_list[pos - 1] != wt:
            raise ValueError(f"Mismatch at position {pos}: expected {wt}, found {seq_list[pos - 1]}")
        seq_list[pos - 1] = mt
    return ''.join(seq_list)



################################################################################################