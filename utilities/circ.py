from pycirclize import Circos
from utilities import functions
from tqdm import tqdm
from matplotlib import cm
from matplotlib import colors

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
    linker_label_texts = ["Linker", "Linker"]

    # Add domain labels with larger font size
    Segmentation_track.xticks(domain_label_positions, domain_label_texts, label_size=12, tick_length=0)

    # Add linker labels with smaller font size
    Segmentation_track.xticks(linker_label_positions, linker_label_texts, label_size=10, tick_length=0)
    # Return the initialized Circos and sector
    return circos, sector

#top_mutation_pairs = ['G140S-Q148H', ....]
def add_DRM_annotation(circos, sector, track, top_mutation_pairs):

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
        circos_pos = pos - 1  # Adjust position to 0-based indexing for Circos
        label = f"{consensus}{pos}"
        if (consensus, pos) not in annotated_positions:
            annotation_track.annotate(circos_pos, label, label_size=13, shorten=100)
            annotated_positions.add((consensus, pos))
    
    
# def connect_background_syn_cords(kn_file, )# Example usage:
# circos, sector = initialize_IN_circos_w_ticks(seq_length, track)
# add_DRM_annotation(circos, sector, track, top_mutation_pairs)

def connect_background_syn_cords(circos,sector,kn_file, top_n,start_aa =1, chord_thickness = 0.5):
    kn_dict = functions.DDE_dict(kn_file)
    # Sort the dictionary by the highest top_n values
    sorted_DDE = sorted(kn_dict.items(), key=lambda x: max(x[1]) if isinstance(x[1], (list, tuple)) else x[1], reverse=True)[:top_n]
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
        start_aa = 1  # Define the starting amino acid position (adjust if needed)

    for (pos1, pos2), dde in tqdm(chords, total=len(chords), desc="Processing background chords"):
        # Define regions for the chord
        region1 = ("IN", pos1 - chord_thickness - start_aa, pos1 + chord_thickness - start_aa)
        region2 = ("IN", pos2 - chord_thickness - start_aa, pos2 + chord_thickness - start_aa)

        # Check if the regions are within valid bounds
        if region1[1] < 0 or region1[2] > actual_length or region2[1] < 0 or region2[2] > actual_length:
            continue

        # Use light gray color for background
        light_gray_color = "#d3d3d3"
        circos.link(region1, region2, color=light_gray_color)

def connect_background_ant_cords(circos,sector,kn_file, bottom_n,start_aa =1, chord_thickness = 0.5):
    kn_dict = functions.DDE_dict(kn_file)
    # Sort the dictionary by the highest top_n values
    sorted_DDE = sorted(kn_dict.items(), key=lambda x: min(x[1]) if isinstance(x[1], (list, tuple)) else x[1])[:bottom_n]
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
        start_aa = 1  # Define the starting amino acid position (adjust if needed)

    for (pos1, pos2), dde in tqdm(chords, total=len(chords), desc="Processing background chords"):
        # Define regions for the chord
        region1 = ("IN", pos1 - chord_thickness - start_aa, pos1 + chord_thickness - start_aa)
        region2 = ("IN", pos2 - chord_thickness - start_aa, pos2 + chord_thickness - start_aa)

        # Check if the regions are within valid bounds
        if region1[1] < 0 or region1[2] > actual_length or region2[1] < 0 or region2[2] > actual_length:
            continue

        # Use light gray color for background
        light_gray_color = "#d3d3d3"
        circos.link(region1, region2, color=light_gray_color)

def connect_syn_cords(circos,sector,kn_file, pairs, redux_dict, start_aa = 1, chord_thickness = 0.5):
    kn_dict = functions.DDE_dict(kn_file)
    # Get the maximum value from kn_dict
    max_dde = max(
        max(values) if isinstance(values, (list, tuple)) else values
        for values in kn_dict.values()
    )
    chord_cmap = cm.get_cmap("Blues")
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

    # chord_thickness = 0.5  # Define the thickness for the background chords
    actual_length = sector.size  # Get the actual length of the sector

    for (pos1, pos2), dde in tqdm(chords, total=len(chords), desc="Processing background chords"):
    # Define regions for the chord
        region1 = ("IN", pos1 - chord_thickness - start_aa, pos1 + chord_thickness - start_aa)
        region2 = ("IN", pos2 - chord_thickness - start_aa, pos2 + chord_thickness - start_aa)

        # Check if the regions are within valid bounds
        if region1[1] < 0 or region1[2] > actual_length or region2[1] < 0 or region2[2] > actual_length:
            continue

        # Use light gray color for background
        color = chord_cmap(chord_norm(dde))
        circos.link(region1, region2, color=color)

def connect_ant_cords(circos,sector,kn_file, pairs, redux_dict, start_aa = 1, chord_thickness = 0.5):
    kn_dict = functions.DDE_dict(kn_file)
    # Get the maximum value from kn_dict
    min_dde = min(
        max(values) if isinstance(values, (list, tuple)) else values
        for values in kn_dict.values()
    )
    print('max_dde:', min_dde)
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

    # chord_thickness = 0.5  # Define the thickness for the background chords
    actual_length = sector.size  # Get the actual length of the sector

    for (pos1, pos2), dde in tqdm(chords, total=len(chords), desc="Processing background chords"):
    # Define regions for the chord
        region1 = ("IN", pos1 - chord_thickness - start_aa, pos1 + chord_thickness - start_aa)
        region2 = ("IN", pos2 - chord_thickness - start_aa, pos2 + chord_thickness - start_aa)

        # Check if the regions are within valid bounds
        if region1[1] < 0 or region1[2] > actual_length or region2[1] < 0 or region2[2] > actual_length:
            continue

        # Use light gray color for background
        color = chord_cmap(chord_norm(dde))
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