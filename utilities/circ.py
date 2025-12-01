from pycirclize import Circos
from utilities import functions
from tqdm import tqdm
from matplotlib import cm
from matplotlib import colors
import numpy as np

global_min_dde = -3.683882
# global_min_dde = -8.50879

global_max_dde = 8.50879
global_max_DDE_sum = 33.17925
#outer rim with ticks every 10 amino acids and labels for every 10th tick
def initialize_IN_circos_w_ticks(seq_length,track):
    """
    Initialize Circos with a single sector for nl43 ΔΔE and configure the main arc and ticks.
    track: tuple (96, 98)
    """
    
    # Initialize Circos with a single sector for nl43 ΔΔE

    sectors = {"IN": seq_length}
    circos = Circos(sectors, space=0, start=0, end=340)

    # Get the sector
    sector = circos.get_sector("IN")

    # Add outer track for the main arc
    start_z, end_z = track
    mid_z = (start_z+end_z)/2
    track1 = sector.add_track(track)
    # track1.rect(0, seq_length, fc="white", ec="black", lw=1.0)
    # Draw domain boundaries with connecting lines and middle regions as "linker"
    track1.rect(0, 46, fc="lightgrey", ec="black", lw=1.0)  # NTD
    track1.rect(46, 58, fc="white", ec="black", lw=1.0)  # Linker between NTD and CCD
    track1.rect(58, 201, fc="lightgrey", ec="black", lw=1.0)  # CCD
    track1.rect(201, 222, fc="white", ec="black", lw=1.0)  # Linker between CCD and CTD
    track1.rect(222, 263, fc="lightgrey", ec="black", lw=1.0)  # CTD

    # Add ticks and labels
    major_ticks_pos = list(range(0, seq_length, 10))
    major_ticks_labels = [str(i+1) if (i+1) % 10 == 1 else "" for i in major_ticks_pos]
    all_ticks_pos = list(range(0, seq_length, 1))
    all_ticks_labels = None

    track1.xticks(major_ticks_pos, major_ticks_labels, label_size=10, tick_length=2, outer=True)
    track1.xticks(all_ticks_pos, all_ticks_labels, label_size=8, tick_length=1, outer=True)

    Segmentation_track = sector.add_track((mid_z-3, mid_z-3))
    # Define domain and linker positions and labels
    domain_positions = [23, 52, 129.5, 211.5, 242.5]  # Adjusted middle positions for each domain and linker
    domain_labels = ["NTD", "Linker", "CCD", "Linker", "CTD"]

    # Separate labels for domains and linkers
    domain_label_positions = [23, 129.5, 242.5]  # Adjusted positions for NTD, CCD, CTD
    domain_label_texts = ["NTD", "CCD", "CTD"]

    linker_label_positions = [52, 211.5]  # Adjusted positions for linkers
    # linker_label_texts = ["Linker", "Linker"]
    linker_label_texts = ["", ""]

    # Add domain labels with larger font size
    Segmentation_track.xticks(domain_label_positions, domain_label_texts, label_size=12, tick_length=0)

    # Add linker labels with smaller font size
    Segmentation_track.xticks(linker_label_positions, linker_label_texts, label_size=10, tick_length=0)
    # Return the initialized Circos and sector
    return circos, sector

def initialize_PR_circos_w_ticks(seq_length,track):
    """
    Initialize Circos with a single sector for nl43 ΔΔE and configure the main arc and ticks.
    track: tuple (96, 98)
    """
    
    # Initialize Circos with a single sector for nl43 ΔΔE

    sectors = {"IN": seq_length}
    circos = Circos(sectors, space=0, start=0, end=340)

    # Get the sector
    sector = circos.get_sector("IN")

    # Add outer track for the main arc
    start_z, end_z = track
    mid_z = (start_z+end_z)/2
    track1 = sector.add_track(track)
    # track1.rect(0, seq_length, fc="white", ec="black", lw=1.0)
    # Draw structural regions with connecting lines and middle regions as "linker"
    track1.rect(0, 4, fc="lightgrey", ec="black", lw=1.0)  # Dimer interface (residues 1–4)
    track1.rect(4, 9, fc="white", ec="black", lw=1.0)  # Linker between Dimer interface and Fulcrum
    track1.rect(9, 23, fc="lightgrey", ec="black", lw=1.0)  # Fulcrum (residues 10–23)
    track1.rect(23, 34, fc="white", ec="black", lw=1.0)  # Linker between Fulcrum and Hinge part 1
    track1.rect(34, 42, fc="lightgrey", ec="black", lw=1.0)  # Hinge part 1 (residues 35–42)
    track1.rect(42, 45, fc="white", ec="black", lw=1.0)  # Linker between Hinge part 1 and Flap
    track1.rect(45, 54, fc="lightgrey", ec="black", lw=1.0)  # Flap (residues 46–54)
    track1.rect(54, 56, fc="white", ec="black", lw=1.0)  # Linker between Flap and Hinge part 2
    track1.rect(56, 61, fc="lightgrey", ec="black", lw=1.0)  # Hinge part 2 (residues 57–61)
    track1.rect(61, 78, fc="lightgrey", ec="black", lw=1.0)  # Cantilever (residues 62–78)
    track1.rect(78, 84, fc="lightgrey", ec="black", lw=1.0)  # 80s loop (residues 79–84, boxed in black)
    track1.rect(84, 85, fc="white", ec="black", lw=1.0)  # Linker between 80s loop and α-helix
    track1.rect(85, 93, fc="lightgrey", ec="black", lw=1.0)  # α-helix (residues 86–93)
    track1.rect(93, 99, fc="lightgrey", ec="black", lw=1.0)  # Dimer interface (residues 94–99)

    # Add ticks and labels
    major_ticks_pos = list(range(0, seq_length, 5))
    major_ticks_labels = [str(i+1) if (i+1) % 5 == 1 else "" for i in major_ticks_pos]
    all_ticks_pos = list(range(0, seq_length, 1))
    all_ticks_labels = None

    track1.xticks(major_ticks_pos, major_ticks_labels, label_size=10, tick_length=2, outer=True)
    track1.xticks(all_ticks_pos, all_ticks_labels, label_size=8, tick_length=1, outer=True)

    Segmentation_track = sector.add_track((mid_z-3, mid_z-3))
    # Define structural region positions and labels
    region_positions = [2, 16.5, 50, 38.5, 59, 70, 81.5, 89.5, 96.5]  # Adjusted middle positions for each region
    region_labels = [
        "DI", "Fulcrum", "Flap", "Hinge", "Hinge", 
        "Cantilever", "80s loop", "α-helix", "DI"
    ]

    # Add region labels with appropriate font size
    Segmentation_track.xticks(region_positions, region_labels, label_size=10, tick_length=0)

    # # Add domain labels with larger font size
    # Segmentation_track.xticks(domain_label_positions, domain_label_texts, label_size=12, tick_length=0)

    # # Add linker labels with smaller font size
    # Segmentation_track.xticks(linker_label_positions, linker_label_texts, label_size=10, tick_length=0)
    # # Return the initialized Circos and sector
    return circos, sector

def initialize_RT_circos_w_ticks(seq_length,track):
    """
    Initialize Circos with a single sector for nl43 ΔΔE and configure the main arc and ticks.
    track: tuple (96, 98)
    """
    start_aa = 39
    end_aa = 226
    aa_length = end_aa - start_aa + 1
    # Initialize Circos with a single sector for nl43 ΔΔE

    sectors = {"IN": aa_length}
    circos = Circos(sectors, space=0, start=0, end=340)

    # Get the sector
    sector = circos.get_sector("IN")

    # Add outer track for the main arc
    start_z, end_z = track
    mid_z = (start_z + end_z) / 2
    track1 = sector.add_track(track)
    # Draw structural regions with connecting lines and middle regions as "linker"
    track1.rect(39 - start_aa, 86 - start_aa, fc="lightgrey", ec="black", lw=1.0)  # Fingers subdomain (residues 1–85)
    # track1.rect(85 - start_aa, 86 - start_aa, fc="white", ec="black", lw=1.0)  # Linker between Fingers and Palm
    track1.rect(86 - start_aa, 118 - start_aa, fc="lightgrey", ec="black", lw=1.0)  # Palm subdomain (residues 86–117)
    # track1.rect(117 - start_aa, 118 - start_aa, fc="white", ec="black", lw=1.0)  # Linker between Palm and Fingers
    track1.rect(118 - start_aa, 156 - start_aa, fc="lightgrey", ec="black", lw=1.0)  # Fingers subdomain (residues 118–155)
    # track1.rect(155 - start_aa, 156 - start_aa, fc="white", ec="black", lw=1.0)  # Linker between Fingers and Palm
    track1.rect(156 - start_aa, 226 - start_aa, fc="lightgrey", ec="black", lw=1.0)  # Palm subdomain (residues 156–226)

    # Add ticks and labels
    major_ticks_pos = list(range(0, seq_length, 10))
    major_ticks_labels = [str(i) for i in range(start_aa, end_aa + 1, 10)]
    all_ticks_pos = list(range(0, seq_length, 1))
    all_ticks_labels = None

    track1.xticks(major_ticks_pos, major_ticks_labels, label_size=10, tick_length=2, outer=True)
    track1.xticks(all_ticks_pos, all_ticks_labels, label_size=8, tick_length=1, outer=True)

    Segmentation_track = sector.add_track((mid_z - 3, mid_z - 3))
    # Define structural region positions and labels
    region_positions = [23, 62.5, 97.5, 157]  # Adjusted middle positions for each region
    region_labels = [
        "Fingers ", "Palm", "Fingers", "Palm"
    ]

    # Add region labels with appropriate font size
    Segmentation_track.xticks(region_positions, region_labels, label_size=10, tick_length=0)

    # Return the initialized Circos and sector
    return circos, sector


#top_mutation_pairs = ['G140S-Q148H', ....]
def add_DRM_annotation(circos, sector, track, top_mutation_pairs, start_aa = 1):

    annotation_track = sector.add_track(track)
    DRM_set = set()
    for pair in top_mutation_pairs:
        pair1 = functions.split_pairs(pair)[0]
        pair2 = functions.split_pairs(pair)[1]
        # print(pair1,pair2)

        consensus1 = functions.get_wt(pair1)
        consensus2 = functions.get_wt(pair2)

        pos1 = functions.get_pos(pair1)
        pos2 = functions.get_pos(pair2)
                                 
        DRM_set.add((consensus1,pos1))
        DRM_set.add((consensus2,pos2))
    # print(DRM_set)
    annotated_positions = set()
    for consensus, pos in DRM_set:
        circos_pos = pos - start_aa  # Adjust position to 0-based indexing for Circos
        label = f"{consensus}{pos}"
        if (consensus, pos) not in annotated_positions:
            annotation_track.annotate(circos_pos, label, label_size=13, shorten=100)
            annotated_positions.add((consensus, pos))
    
    
# def connect_background_syn_cords(kn_file, )# Example usage:
# circos, sector = initialize_IN_circos_w_ticks(seq_length, track)
# add_DRM_annotation(circos, sector, track, top_mutation_pairs)

def connect_background_syn_cords(circos,sector,kn_file, top_n,start_aa =1, chord_thickness = 0.5, color_scheme = 'IN',chord_factor = 'color',chord_comparison ='local',weight_power=2,a = 2,b = 0.1,x_min = 1.8230609, x_max = 33.17926, weight = 1):
    chord_thickness_double = chord_thickness*2
    kn_dict = functions.DDE_dict(kn_file)
    # Sort the dictionary by the highest top_n values
    sorted_DDE = sorted(kn_dict.items(), key=lambda x: max(x[1]) if isinstance(x[1], (list, tuple)) else x[1], reverse=True)[:top_n]
    max_dde = max(
        max(values) if isinstance(values, (list, tuple)) else values
        for _, values in sorted_DDE
    )
    if chord_comparison == 'global':
        max_dde = global_max_dde
    # Iterate through the top_n sorted items
    chords = []
    for key, value in sorted_DDE:
        # print(key, value)
        # Process each key-value pair as needed
        chords.extend([(functions.get_pos1_pos2(key), v) for v in value] if isinstance(value, (list, tuple)) else [(functions.get_pos1_pos2(key), value)])
        # print(f"Key: {key}, Values: {value}")
        # print(functions.get_pos1_pos2(key))

        # chord_thickness = 0.5  # Define the thickness for the background chords
        actual_length = sector.size  # Get the actual length of the sector
    

    pair_count = {}
    for (pos1, pos2), _ in chords:
        if pos1 > pos2:
            pos1, pos2 = pos2, pos1
        pair_count[(pos1, pos2)] = pair_count.get((pos1, pos2), 0) + 1

    pair_dde_sum = {}
    for (pos1, pos2), dde in chords:
        if pos1 > pos2:
            pos1, pos2 = pos2, pos1
        if (pos1, pos2) not in pair_dde_sum:
            pair_dde_sum[(pos1, pos2)] = 0
        pair_dde_sum[(pos1, pos2)] += abs(dde)

    max_pair_dde_sum = max(pair_dde_sum.values(), default=0)
    min_pair_dde_sum = min(pair_dde_sum.values(), default=0)
    print(f"Max pair DDE sum: {max_pair_dde_sum}, Key: {max(pair_dde_sum, key=pair_dde_sum.get, default=None)}")
    print(f"Min pair DDE sum: {min_pair_dde_sum}, Key: {min(pair_dde_sum, key=pair_dde_sum.get, default=None)}")

    for (pos1, pos2), dde in tqdm(chords, total=len(chords), desc="Processing background chords"):
        if pos1 > pos2:
            pos1, pos2 = pos2, pos1
        # Define regions for the chord
        if chord_factor == 'width':
            weight = 0.1 + 0.9 * (dde / max_dde)
            weight = weight**weight_power
            chord_thickness = weight * chord_thickness_double
        if chord_factor == 'cw':
            weight = 0.1 + (0.9 * (pair_count[(pos1, pos2)] / 6))**weight_power
            # weight = weight**weight_power
            chord_thickness = weight * chord_thickness_double
        
        if chord_factor == 'cw2':
            weight = 0.1 + 0.9 * (np.log(pair_dde_sum[(pos1, pos2)] + 1)/np.log(global_max_DDE_sum+1))
            chord_thickness = weight * chord_thickness_double

        if chord_factor == 'cw3':
            weight =0.2+(0.8*pair_dde_sum[(pos1, pos2)] / global_max_DDE_sum)**weight_power
            chord_thickness = weight * chord_thickness_double
        
        if chord_factor == 'cw4':
            # weight =pair_dde_sum[(pos1, pos2)]
            weight = np.log10(pair_dde_sum[(pos1, pos2)] + 1)
            chord_thickness = weight * chord_thickness_double

        if chord_factor == 'cw5':
            abs_sum_dde= abs(pair_dde_sum[(pos1, pos2)])
            chord_thickness = weight*(b+(a-b)*(abs_sum_dde-x_min)/(x_max-x_min))/2
        
        if chord_factor == 'cw6':
            abs_sum_dde= abs(pair_dde_sum[(pos1, pos2)])
            chord_thickness = weight*(b+(a-b)*(abs_sum_dde-x_min)/(x_max-x_min))/2
            
        region1 = ("IN", pos1 - chord_thickness - start_aa, pos1 + chord_thickness - start_aa)
        region2 = ("IN", pos2 - chord_thickness - start_aa, pos2 + chord_thickness - start_aa)

        # Check if the regions are within valid bounds
        if region1[1] < 0 or region1[2] > actual_length or region2[1] < 0 or region2[2] > actual_length:
            print(f"Chord out of range: region1={region1}, region2={region2}")
            continue
    
        # Use light gray color for background
        if color_scheme == 'IN':
            color = "#ffe9cf"
        elif color_scheme == 'PR':
            color = '#e6e6fa'
        elif color_scheme == 'RT':
            color = "#e2ecd4"
        else:
            color = "#d3d3d3"  # Default color if no scheme matches
        color = "#d3d3d3"  # Light grey color
        circos.link(region1, region2, color=color, alpha =1.0)
    # max_pair_count = max(pair_count.values(), default=0)
    # max_pair_key = max(pair_count, key=pair_count.get, default=None)
    # print(f"Max pair count: {max_pair_count}, Key: {max_pair_key}")

def connect_background_ant_cords(circos,sector,kn_file, bottom_n,start_aa =1, chord_thickness = 0.5, color_scheme='IN',chord_factor = 'color',chord_comparison ='local',weight_power=2,a = 2,b = 0.1,x_min = 1.8230609, x_max = 33.17926, weight = 1):
    chord_thickness_double = chord_thickness*2
    kn_dict = functions.DDE_dict(kn_file)
    # Sort the dictionary by the highest top_n values
    sorted_DDE = sorted(kn_dict.items(), key=lambda x: min(x[1]) if isinstance(x[1], (list, tuple)) else x[1])[:bottom_n]
    min_dde = min(
        min(values) if isinstance(values, (list, tuple)) else values
        for _, values in sorted_DDE
        )
    # Iterate through the top_n sorted items
    if chord_comparison == 'global':
        min_dde = global_min_dde

    chords = []
    for key, value in sorted_DDE:
        # print(key, value)
        # Process each key-value pair as needed
        chords.extend([(functions.get_pos1_pos2(key), v) for v in value] if isinstance(value, (list, tuple)) else [(functions.get_pos1_pos2(key), value)])
        # print(f"Key: {key}, Values: {value}")
        # print(functions.get_pos1_pos2(key))

        # chord_thickness = 0.5  # Define the thickness for the background chords
        actual_length = sector.size  # Get the actual length of the sector

    pair_count = {}
    for (pos1, pos2), _ in chords:
        if pos1 > pos2:
            pos1, pos2 = pos2, pos1
        pair_count[(pos1, pos2)] = pair_count.get((pos1, pos2), 0) + 1
    pair_dde_sum = {}
    for (pos1, pos2), dde in chords:
        if pos1 > pos2:
            pos1, pos2 = pos2, pos1
        if (pos1, pos2) not in pair_dde_sum:
            pair_dde_sum[(pos1, pos2)] = 0
        pair_dde_sum[(pos1, pos2)] += abs(dde)
    for (pos1, pos2), dde in tqdm(chords, total=len(chords), desc="Processing background chords"):
        # Define regions for the chord
        if pos1 > pos2:
            pos1, pos2 = pos2, pos1
        if chord_factor == 'width':
            weight = 0.1 + 0.9 * (abs(dde) / abs(min_dde))
            weight = weight**weight_power
            chord_thickness = weight * chord_thickness_double
        if chord_factor == 'cw':
            weight = 0.1 + (0.9 * (pair_count[(pos1, pos2)] / 6))**weight_power
            # weight = weight**weight_power
            chord_thickness = weight * chord_thickness_double

        if chord_factor == 'cw2':
            weight = 0.1 + 0.9 * (np.log(pair_dde_sum[(pos1, pos2)] + 1)/np.log(global_max_DDE_sum+1))
            chord_thickness = weight * chord_thickness_double

        if chord_factor == 'cw3':
            weight =0.2+(0.8*pair_dde_sum[(pos1, pos2)] / global_max_DDE_sum)**weight_power
            chord_thickness = weight * chord_thickness_double

        if chord_factor == 'cw4':
            # weight =pair_dde_sum[(pos1, pos2)]
            weight = np.log10(pair_dde_sum[(pos1, pos2)] + 1)
            chord_thickness = weight * chord_thickness_double
        if chord_factor == 'cw5':
            abs_sum_dde= abs(pair_dde_sum[(pos1, pos2)])
            chord_thickness = weight*(b+(a-b)*(abs_sum_dde-x_min)/(x_max-x_min))/2

        if chord_factor == 'cw6':
            abs_sum_dde= abs(pair_dde_sum[(pos1, pos2)])
            chord_thickness = weight*(b+(a-b)*(abs_sum_dde-x_min)/(x_max-x_min))/2
            
        region1 = ("IN", pos1 - chord_thickness - start_aa, pos1 + chord_thickness - start_aa)
        region2 = ("IN", pos2 - chord_thickness - start_aa, pos2 + chord_thickness - start_aa)

        # Check if the regions are within valid bounds
        if region1[1] < 0 or region1[2] > actual_length or region2[1] < 0 or region2[2] > actual_length:
            print(f"Chord out of range: region1={region1}, region2={region2}")
            continue

        # Use light gray color for background
        if color_scheme == 'IN':
            color = "#ffe9cf"
        elif color_scheme == 'PR':
            color = '#e6e6fa'
        elif color_scheme == 'RT':
            color = "#e2ecd4"
        else:
            color = "#d3d3d3"  # Default color if no scheme matches
        color = "#d3d3d3"  # Light grey color
        circos.link(region1, region2, color=color, alpha =1.0)

def connect_syn_cords(circos,sector,kn_file, pairs, redux_dict, start_aa = 1, chord_thickness = 0.5, chord_factor = 'color',chord_comparison ='local',weight_power=2,a = 2,b = 0.1,x_min = 1.8230609, x_max = 33.17926, weight = 1):
    chord_thickness_double = chord_thickness*2
    kn_dict = functions.DDE_dict(kn_file)
    # Get the maximum value from kn_dict
    max_dde = max(
        max(values) if isinstance(values, (list, tuple)) else values
        for values in kn_dict.values()
    )
    if chord_comparison == 'global':
        max_dde = global_max_dde
    chord_cmap = cm.get_cmap("Blues")
    # max_dde = 8.50879
    chord_norm = colors.Normalize(vmin=0, vmax=max_dde)
    print('max_dde:', max_dde)

    chords = []
    for pair in pairs: 
        pair1 = functions.split_pairs(pair)[0]
        pair2 = functions.split_pairs(pair)[1]  
        pair1_reduced = functions.unreduced_to_reduced(redux_dict,pair1)
        pair2_reduced = functions.unreduced_to_reduced(redux_dict,pair2)

        key = f"{pair1_reduced}-{pair2_reduced}" if f"{pair1_reduced}-{pair2_reduced}" in kn_dict else f"{pair2_reduced}-{pair1_reduced}"
        if key in kn_dict:
            value = kn_dict[key]
            chords.extend([(functions.get_pos1_pos2(key), v) for v in value] if isinstance(value, (list, tuple)) else [(functions.get_pos1_pos2(key), value)])
            chords = sorted(chords, key=lambda x: x[1])
    # chord_thickness = 0.5  # Define the thickness for the background chords

    actual_length = sector.size  # Get the actual length of the sector

    pair_count = {}
    for (pos1, pos2), _ in chords:
        if pos1 > pos2:
            pos1, pos2 = pos2, pos1
        pair_count[(pos1, pos2)] = pair_count.get((pos1, pos2), 0) + 1
    pair_dde_sum = {}
    for (pos1, pos2), dde in chords:
        if pos1 > pos2:
            pos1, pos2 = pos2, pos1
        if (pos1, pos2) not in pair_dde_sum:
            pair_dde_sum[(pos1, pos2)] = 0
        pair_dde_sum[(pos1, pos2)] += abs(dde)
    for (pos1, pos2), dde in tqdm(chords, total=len(chords), desc="Processing background chords"):

    # Define regions for the chord
        #CHORD THICKNESS
        if pos1 > pos2:
            pos1, pos2 = pos2, pos1
        if chord_factor == 'width':
            weight = 0.1 + 0.9 * (dde / max_dde)
            weight = weight**weight_power
            chord_thickness = weight * chord_thickness_double
            # chord_thickness = weight
        #CHORD THICKNESS
        if chord_factor == 'cw':
            weight = 0.1 + (0.9 * (pair_count[(pos1, pos2)] / 6))**weight_power

            chord_thickness = weight * chord_thickness_double

        if chord_factor == 'cw2':
            weight = 0.1 + 0.9 * (np.log(pair_dde_sum[(pos1, pos2)] + 1)/np.log(global_max_DDE_sum+1))
            chord_thickness = weight * chord_thickness_double
        
        if chord_factor == 'cw3':
            weight =0.2+(0.8*pair_dde_sum[(pos1, pos2)] / global_max_DDE_sum)**weight_power
            chord_thickness = weight * chord_thickness_double
            # if (pos1, pos2) == (140,148):
            #     print(chord_thickness)
        if chord_factor == 'cw4':
            # weight =pair_dde_sum[(pos1, pos2)]
            weight = np.log10(pair_dde_sum[(pos1, pos2)] + 1)
            chord_thickness = weight * chord_thickness_double
        
        if chord_factor == 'cw5':
            abs_sum_dde= abs(pair_dde_sum[(pos1, pos2)])
            chord_thickness = weight*(b+(a-b)*(abs_sum_dde-x_min)/(x_max-x_min))/2

        if chord_factor == 'cw6':
            abs_sum_dde= abs(pair_dde_sum[(pos1, pos2)])
            chord_thickness = weight*(b+(a-b)*(abs_sum_dde-x_min)/(x_max-x_min))/2
            
        region1 = ("IN", pos1 - chord_thickness - start_aa, pos1 + chord_thickness - start_aa)
        region2 = ("IN", pos2 - chord_thickness - start_aa, pos2 + chord_thickness - start_aa)
            
        # Check if the regions are within valid bounds
        if region1[1] < 0 or region1[2] > actual_length or region2[1] < 0 or region2[2] > actual_length:
            print(f"Chord out of range: region1={region1}, region2={region2}")
            continue

        # Use light gray color for background
        # color = chord_cmap(chord_norm(dde))
        #CHORD THICKNESS
        if chord_factor == 'width' or chord_factor == 'cw6':
            color = "#08519C"  # Simple blue color
        # elif chord_factor == 'color':
        #     color = chord_cmap(chord_norm(dde))
        elif chord_factor != 'width' and chord_factor != 'cw6':
            color = chord_cmap(chord_norm(dde))
        #CHORD THICKNESS
        
        circos.link(region1, region2, color=color)

def connect_ant_cords(circos,sector,kn_file, pairs, redux_dict, start_aa = 1, chord_thickness = 0.5, chord_factor = 'color',chord_comparison ='local',weight_power=2,a = 2,b = 0.1,x_min = 1.8230609, x_max = 33.17926, weight = 1):
    chord_thickness_double = chord_thickness*2
    kn_dict = functions.DDE_dict(kn_file)
    # Get the maximum value from kn_dict
    min_dde = min(
        max(values) if isinstance(values, (list, tuple)) else values
        for values in kn_dict.values()
    )
    print('min_dde:', min_dde)
    if chord_comparison == 'global':
        min_dde = global_min_dde
        # min_dde = -8.50879
    chord_cmap = cm.get_cmap("Reds_r")
    chord_norm = colors.Normalize(vmin=min_dde, vmax=0)

    chords = []
    for pair in pairs: 
        pair1 = functions.split_pairs(pair)[0]
        pair2 = functions.split_pairs(pair)[1]  
        pair1_reduced = functions.unreduced_to_reduced(redux_dict,pair1)
        pair2_reduced = functions.unreduced_to_reduced(redux_dict,pair2)

        key = f"{pair1_reduced}-{pair2_reduced}" if f"{pair1_reduced}-{pair2_reduced}" in kn_dict else f"{pair2_reduced}-{pair1_reduced}"
        if key in kn_dict:
            value = kn_dict[key]
            chords.extend([(functions.get_pos1_pos2(key), v) for v in value] if isinstance(value, (list, tuple)) else [(functions.get_pos1_pos2(key), value)])
    chords = sorted(chords, key=lambda x: x[1], reverse=True)
    # chord_thickness = 0.5  # Define the thickness for the background chords
    actual_length = sector.size  # Get the actual length of the sector

    pair_count = {}
    for (pos1, pos2), _ in chords:
        if pos1 > pos2:
            pos1, pos2 = pos2, pos1
        pair_count[(pos1, pos2)] = pair_count.get((pos1, pos2), 0) + 1
    pair_dde_sum = {}
    for (pos1, pos2), dde in chords:
        if pos1 > pos2:
            pos1, pos2 = pos2, pos1
        if (pos1, pos2) not in pair_dde_sum:
            pair_dde_sum[(pos1, pos2)] = 0
        pair_dde_sum[(pos1, pos2)] += abs(dde)
    for (pos1, pos2), dde in tqdm(chords, total=len(chords), desc="Processing background chords"):
    # Define regions for the chord
        if pos1 > pos2:
                pos1, pos2 = pos2, pos1
        if chord_factor == 'width':
            weight = 0.1 + 0.9 * (abs(dde) / abs(min_dde))
            weight = weight**weight_power
            chord_thickness = weight * chord_thickness_double
            # chord_thickness = weight
        if chord_factor == 'cw':
            weight = 0.1 + (0.9 * (pair_count[(pos1, pos2)] / 6))**weight_power
            # weight = weight**weight_power
            chord_thickness = weight * chord_thickness_double

        if chord_factor == 'cw2':
            weight = 0.1 + 0.9 * (np.log(pair_dde_sum[(pos1, pos2)] + 1)/np.log(global_max_DDE_sum+1))
            chord_thickness = weight * chord_thickness_double

        if chord_factor == 'cw3':
            weight =0.2+(0.8*pair_dde_sum[(pos1, pos2)] / global_max_DDE_sum)**weight_power
            chord_thickness = weight * chord_thickness_double
        if chord_factor == 'cw4':
            # weight =pair_dde_sum[(pos1, pos2)]
            weight = np.log10(pair_dde_sum[(pos1, pos2)] + 1)
            chord_thickness = weight * chord_thickness_double
            
        if chord_factor == 'cw5':
            abs_sum_dde= abs(pair_dde_sum[(pos1, pos2)])
            chord_thickness = weight*(b+(a-b)*(abs_sum_dde-x_min)/(x_max-x_min))/2

        if chord_factor == 'cw6':
            abs_sum_dde= abs(pair_dde_sum[(pos1, pos2)])
            chord_thickness = weight*(b+(a-b)*(abs_sum_dde-x_min)/(x_max-x_min))/2

        region1 = ("IN", pos1 - chord_thickness - start_aa, pos1 + chord_thickness - start_aa)
        region2 = ("IN", pos2 - chord_thickness - start_aa, pos2 + chord_thickness - start_aa)

        # Check if the regions are within valid bounds
        if region1[1] < 0 or region1[2] > actual_length or region2[1] < 0 or region2[2] > actual_length:
            print(f"Chord out of range: region1={region1}, region2={region2}")
            continue

        if chord_factor == 'width' or chord_factor == 'cw6':
            color = "#A50F15"  # Use dark red color
        # elif chord_factor == 'color':
        #     color = chord_cmap(chord_norm(dde))
        elif chord_factor != 'width' and chord_factor != 'cw6': #color and width
            color = chord_cmap(chord_norm(dde))
        #CHORD THICKNESS
        circos.link(region1, region2, color=color)

def plot_colorbar(min_v, max_v, syn=True):
    if syn == True:
        cmap = cm.get_cmap("Blues")
    else:
        cmap = cm.get_cmap("Reds_r")
    norm = colors.Normalize(vmin=min_v, vmax=max_v)
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    return sm