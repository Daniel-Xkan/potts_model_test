from pycirclize import Circos
#outer rim with ticks every 10 amino acids and labels for every 10th tick
def initialize_circos_w_ticks(seq_length,track):
    """
    Initialize Circos with a single sector for nl43 ΔΔE and configure the main arc and ticks.
    track: tuple (96, 98)
    """
    
    # Initialize Circos with a single sector for nl43 ΔΔE

    sectors = {"nl43 ΔΔE": seq_length}
    circos = Circos(sectors, space=0, start=0, end=340)

    # Get the sector
    sector = circos.get_sector("nl43 ΔΔE")

    # Add outer track for the main arc
    track1 = sector.add_track(track)
    track1.rect(0, seq_length, fc="grey", ec="grey", lw=1.0)

    # Add ticks and labels
    major_ticks_pos = list(range(0, seq_length, 10))
    major_ticks_labels = [str(i+1) if (i+1) % 10 == 1 else "" for i in major_ticks_pos]
    all_ticks_pos = list(range(0, seq_length, 1))
    all_ticks_labels = None

    track1.xticks(major_ticks_pos, major_ticks_labels, label_size=10, tick_length=2, outer=True)
    track1.xticks(all_ticks_pos, all_ticks_labels, label_size=8, tick_length=1, outer=True)

    # Return the initialized Circos and sector
    return circos, sector


