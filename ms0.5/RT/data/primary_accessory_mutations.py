# ============================================================
# NRTI MUTATIONS
# ============================================================

NRTI_primary_mutations = [
    # Alone causes meaningful reduction in susceptibility to ≥1 NRTI

    "M184V", "M184I",   # >200-fold reduced 3TC/FTC; ~3-fold reduced ABC
    "K65R",             # ~2-fold TFV/ABC alone; 5-10-fold 3TC/FTC
    "K65N",             # ~2-fold TFV/ABC; ~5-fold 3TC
    "T215Y",            # ~7-fold AZT alone; 1.5-fold ABC/TDF
    "T215F",            # Lesser but documented solo effect on AZT/ABC/TFV
    "K70R",             # ~5-fold reduced AZT alone
    "M41L",             # Usually with T215Y but documented combination effect (12-fold AZT)
    "L210W",            # Usually with M41L/T215Y; part of documented high-level combination
    "D67N",             # Reduces AZT susceptibility alone
    "L74V",             # With M184V reduces ABC ~5-fold; documented solo ABC effect
    "L74I",             # With M184V reduces ABC ~3-fold
    "Y115F",            # Alone reduces ABC ~3-fold
    "Q151M",            # High-level AZT/ABC resistance; intermediate 3TC/FTC/TFV
    "A62V",             # Corrects K65R replication deficit; part of Q151M complex
    "V75I",             # Reduces susceptibility to d4T/ddI; possibly ABC
    "N348I",            # Alone ~2-fold AZT reduction; also affects NVP/EFV
]

NRTI_accessory_mutations = [
    # Alone have little/no effect; contribute primarily in combination with other NRTI DRMs

    # TAM accessory mutations
    "E40F",             # Contributes only with multiple TAMs
    "E44D", "E44A",     # Contributes only with multiple TAMs
    "V118I",            # Polymorphic; contributes only with multiple TAMs
    "K43Q", "K43N",     # Poorly characterized; occurs with multiple TAMs
    "E203K",            # Poorly characterized; occurs with multiple TAMs
    "H208Y",            # Poorly characterized; occurs with multiple TAMs
    "D218E",            # Poorly characterized; occurs with multiple TAMs
    "K223Q", "K223E",   # Poorly characterized; occurs with multiple TAMs
    "L228H", "L228R",   # Poorly characterized; occurs with multiple TAMs

    # TAM position variants
    "D67G", "D67E",     # Contribute to reduced susceptibility only in combination
    "K219Q", "K219E",   # Reduce AZT susceptibility only with other TAMs
    "K219N", "K219R",   # Contribute to reduced susceptibility only in combination

    # K70 non-TAM variants
    "K70E", "K70G", "K70Q", "K70T", "K70N", "K70S",  # Minimal detectable susceptibility reduction alone

    # T215 revertants (not resistance mutations per se)
    "T215D", "T215C", "T215E", "T215I", "T215S", "T215N", "T215A", "T215V",

    # Miscellaneous accessory
    "S68G", "S68N",     # Effects on TFV susceptibility uncertain
    "T69D",             # Reduces ddI/possibly d4T; not AZT/ABC/TFV alone
    "T69S", "T69N",     # Occur with TAMs; no solo effect on AZT/ABC/TFV
    "T69G",             # Extremely rare; occurs with TAMs and position 67 deletions
    "V75T",             # Reduces d4T/ddI susceptibility; possibly ABC
    "V75M", "V75S", "V75A",  # Rare; uncertain phenotypic significance
    "M184L", "M184T",   # Rare; reduced replication fitness; not well characterized as resistance
    "K65E",             # Extremely rare; too unfit to study phenotypically
    "Q151L",            # Extremely rare; transition to Q151M; does not reduce NRTI susceptibility
]


# ============================================================
# NNRTI MUTATIONS
# ============================================================

NNRTI_primary_mutations = [
    # Alone causes meaningful reduction in susceptibility to ≥1 NNRTI

    "K101E",        # 3-10-fold NVP; ~2-fold EFV/ETR/RPV/DOR alone
    "K101P",        # >20-fold NVP/EFV/RPV alone; ~5-fold ETR
    "K103N",        # ~50-fold NVP; ~20-fold EFV alone
    "K103S",        # ~5-fold EFV; >20-fold NVP alone
    "K103H",        # ~20-fold NVP/EFV alone
    "K103T",        # >20-fold NVP alone
    "V106A",        # ~50-fold NVP; ~5-fold EFV; ~10-fold DOR alone
    "V106M",        # >30-fold NVP/EFV alone; ~3-fold DOR
    "Y181C",        # >50-fold NVP; ~5-fold ETR; ~3-fold RPV alone
    "Y181I", "Y181V",  # >50-fold NVP; 0-15-fold ETR/RPV alone
    "Y181S",        # ~30-fold NVP alone
    "Y188L",        # >50-fold NVP/EFV/DOR; ~5-fold RPV alone
    "Y188C",        # >30-fold NVP alone
    "Y188H",        # ~5-fold NVP/EFV alone
    "G190A",        # >50-fold NVP; 5-10-fold EFV alone
    "G190S",        # >50-fold NVP/EFV alone; variable DOR
    "G190E",        # >100-fold NVP/EFV; >10-fold ETR/RPV/DOR alone
    "G190Q",        # >100-fold EFV/NVP alone
    "G190C", "G190T", "G190V",  # High-level NVP/EFV resistance alone
    "M230L",        # Up to 5-fold ETR/RPV; >10-fold NVP/EFV/DOR alone
    "M230I",        # ~2-3-fold ETR/RPV; ~5-fold EFV; ~10-fold NVP alone
    "K238T",        # ~5-fold NVP alone
    "Y318F",        # ~11-fold DOR alone (median in 3 clinical isolates)
    "E138K",        # ~2-fold RPV/ETR alone; ~3-fold in combo with M184I/K101E
    "E138A",        # ~2-fold ETR/RPV alone
    "E138Q", "E138G",  # ~3-fold NVP; ~2-fold ETR/RPV alone
    "E138R",        # Slightly more than E138Q/G
    "V179D",        # ~2-fold NVP/EFV/ETR/RPV alone
    "A98G",         # ~2-fold NVP/EFV/RPV/DOR alone
    "V108I",        # ~2-fold NVP/EFV alone
    "P236L",        # ~4-fold NVP alone; high-level delavirdine resistance
    "N348I",        # ~2-fold NVP/EFV alone (also listed in NRTI primary)
]

NNRTI_accessory_mutations = [
    # Alone have little/no effect; contribute primarily in combination with other NNRTI DRMs

    "V90I",         # Little if any solo effect; contributes in ETR genotypic score
    "K101H",        # Contributes NVP/EFV reduction only in combination
    "K101Q",        # Minimal if any solo effect
    "K101N", "K101A", "K101T",  # Not well studied
    "K103R",        # No solo effect; synergizes with V179D for ~15-fold NVP/EFV
    "K103E", "K103Q",  # Not selected by NNRTIs; no reduced susceptibility
    "V106I",        # Little if any solo effect; contributes in combination
    "V106L",        # Rare; not well characterized
    "V179F",        # No solo effect; with Y181C confers >10-fold ETR/RPV
    "V179E",        # Weakly selected; effect not well defined
    "V179T",        # No solo reduced susceptibility; contributes in ETR score
    "V179L",        # Minimal solo effect; listed in RPV package insert
    "V179I",        # Highly polymorphic; no clear solo resistance effect
    "H221Y",        # Minimal solo effects; usually occurs with Y181C
    "P225H",        # Usually with K103N; alone not characterized as primary
    "F227L",        # Little solo effect; high-level resistance only with V106A
    "F227C",        # Usually selected in combination; solo effect not documented alone
    "F227I", "F227V",  # Extremely rare; not studied phenotypically
    "F227Y",        # Usually with other DRMs; few data on solo impact
    "Y232H",        # Nearly always in combination; contributes NVP/EFV reduction
    "L234I",        # Solo effect on NNRTIs not known; high-level DOR only with V106A
    "K238N",        # Minimal solo effect on NNRTI susceptibility
    "L100I",        # Rarely occurs in isolation; solo effect not well documented alone
    "L100V",        # Extremely rare; few isolates studied
]