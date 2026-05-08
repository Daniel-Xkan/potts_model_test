IN_primary_mutations = [
    # Alone reduces susceptibility meaningfully (typically ≥5-fold to at least one INSTI)
    "H51Y",        # Alone reduces EVG susceptibility 2-3-fold (borderline, but documented solo effect)
    "T66A", "T66I", "T66K",   # T66A/I reduce EVG 5-10-fold alone; T66K reduces EVG ~40-fold, RAL ~10-fold
    "E92Q",        # Alone reduces RAL ~5-fold and EVG ~30-fold
    "E92G", "E92V",# Alone reduces EVG >=10-fold
    "T97A",        # Alone reduces EVG ~3-fold (minor, but documented solo effect vs EVG)
    "G118R",       # Alone ~9-fold reduced DTG/CAB, ~3-fold BIC
    "F121Y",       # Alone >10-fold reduced RAL
    "F121C",       # High-level resistance to RAL and EVG alone
    "G140S", "G140A", "G140C",  # Alone ~5-fold reduced EVG
    "Y143C", "Y143R", "Y143H", "Y143K", "Y143S", "Y143G", "Y143A",  # Alone reduce RAL 3-20-fold
    "P145S",       # High-level EVG resistance alone
    "Q146I",       # High-level EVG resistance alone
    "Q146P", "Q146L",  # Alone 5-15-fold reduced EVG
    "S147G",       # Alone ~5-fold reduced EVG
    "Q148H", "Q148K", "Q148R",  # Alone reduce RAL/EVG 5-100-fold; Q148R/K reduce CAB ~5-fold
    "V151L",       # Alone reduces RAL/EVG 15-20-fold, CAB/DTG ~3-fold
    "V151A",       # Alone reduces RAL/EVG ~5-fold
    "S153Y", "S153F",  # Alone reduce EVG ~5-fold, DTG/CAB/BIC ~2-fold
    "N155H",       # Alone reduces RAL ~10-fold, EVG ~30-fold
    "N155S", "N155T",  # Alone reduce RAL ~5-fold, EVG >30-fold
    "R263K",       # Alone reduces DTG/BIC/CAB ~2-fold
]

IN_accessory_mutations = [
    # Alone have little/no effect; require combinations with other DRMs for impact
    "A49G",        # Few data; always seen in combination with R263K
    "M50I",        # Minimal if any solo effect
    "L74M",        # Alone does not reduce susceptibility
    "L74I",        # Not shown to reduce susceptibility alone or in combination
    "L74F",        # Only reduces susceptibility in combination (L74F + V75I + N155H or G140S/Q148H)
    "V75I",        # Rarely selected; effect alone not established
    "V75A",        # Rarely selected; solo effect unknown
    "Q95K",        # Alone has little effect
    "E138K", "E138A", "E138T",  # Alone do not reduce susceptibility
    "G140R",       # Requires specific mutational background; alone may not reduce susceptibility
    "T122N",       # Limited data; contributes only in combination
    "Q146R",       # Selected in CAB recipients but ≤2-fold solo effect on DTG/CAB/BIC
    "Q148N",       # Rare, only ~3-fold EVG reduction; likely revertant
    "G149A",       # No effect alone; only contributes in combination
    "V151I",       # Alone little or no effect
    "S153A",       # Not selected by INSTIs; no susceptibility reduction
    "E157Q",       # Minimal if any solo effect
    "G163R", "G163K",  # Accessory; phenotypic effects not well characterized
    "S230R",       # Alone minimal effect; occurs in combinations
    "S230N",       # Polymorphism; not associated with reduced susceptibility
    "D232N",       # Alone little effect
]