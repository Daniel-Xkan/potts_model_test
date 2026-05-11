# ============================================================
# NRTI MUTATIONS - sorted by position
# ============================================================

NRTI_primary_mutations = [
    # Alone causes meaningful reduction in susceptibility to ≥1 NRTI

    "M41L",             # With T215Y ~12-fold AZT; primary TAM pathway
    "K65R",             # ~2-fold TFV/ABC; 5-10-fold 3TC/FTC alone
    "K65N",             # ~2-fold TFV/ABC; ~5-fold 3TC alone
    "K70R",             # ~5-fold AZT alone
    "L74V",             # With M184V ~5-fold ABC; ABC-selected
    "L74I",             # With M184V ~3-fold ABC
    "Y115F",            # ~3-fold ABC alone
    "F116Y",            # Q151M complex; high-level AZT/ABC in complex
    "Q151M",            # High-level AZT/ABC; intermediate 3TC/FTC/TFV
    "M184V", "M184I",   # >200-fold 3TC/FTC; ~3-fold ABC alone
    "L210W",            # M41L/L210W/T215Y >100-fold AZT; ~3-fold ABC/TFV
    "T215Y",            # ~7-fold AZT; 1.5-fold ABC/TDF alone
    "T215F",            # Lesser but documented solo AZT/ABC/TFV effect
    "K219Q", "K219E",   # Reduce AZT with other TAMs
    "N348I",            # ~2-fold AZT alone; also reduces NVP/EFV
]

NRTI_accessory_mutations = [
    # Alone have little/no effect; require combination with other NRTI DRMs

    # Early RT positions
    "S68G", "S68N",         # Effects on TFV susceptibility uncertain
    "T69D",                 # Reduces ddI/possibly d4T; not AZT/ABC/TFV alone
    "T69S", "T69N",         # Occur with TAMs; no solo effect on AZT/ABC/TFV
    "T69G",                 # Extremely rare; occurs with TAMs and pos.67 deletions

    # TAM accessory and position variants
    "E40F",                 # Contributes only with multiple TAMs
    "E44D", "E44A",         # Contributes only with multiple TAMs
    "A62V",                 # Corrects replication deficit of K65R/Q151M; not direct resistance
    "D67N",                 # Reduces AZT primarily but contributes more in combinations
    "D67G", "D67E",         # Contribute to reduced susceptibility only in combination
    "K70E", "K70G",         # Minimal detectable reduction alone
    "K70Q", "K70T",         # Minimal detectable reduction alone
    "K70N", "K70S",         # Minimal detectable reduction alone

    # Q151M complex accessory mutations
    "V75I",                 # Q151M complex accessory; solo effect not clearly documented
    "V75T",                 # Reduces d4T/ddI; possibly ABC
    "V75M", "V75S", "V75A", # Rare; uncertain phenotypic significance
    "F77L",                 # Q151M complex accessory mutation; solo effect not documented

    # Higher position accessory TAMs
    "V118I",                # Polymorphic; contributes only with multiple TAMs
    "K43Q", "K43N",         # Poorly characterized; occurs with multiple TAMs
    "E203K",                # Poorly characterized; occurs with multiple TAMs
    "H208Y",                # Poorly characterized; occurs with multiple TAMs
    "D218E",                # Poorly characterized; occurs with multiple TAMs
    "K219N", "K219R",       # Contribute to reduced susceptibility only in combination
    "K223Q", "K223E",       # Poorly characterized; occurs with multiple TAMs
    "L228H", "L228R",       # Poorly characterized; occurs with multiple TAMs

    # T215 revertants
    "T215D", "T215C", "T215E", "T215I",
    "T215S", "T215N", "T215A", "T215V",

    # Rare/unfit/transitional mutations
    "K65E",                 # Extremely rare; too unfit to study phenotypically
    "M184L", "M184T",       # Rare; reduced fitness; not well characterized as resistance
    "Q151L",                # Transition to Q151M; does not reduce NRTI susceptibility
]


# ============================================================
# NNRTI MUTATIONS - sorted by position
# ============================================================

NNRTI_primary_mutations = [
    # Alone causes meaningful reduction in susceptibility to ≥1 NNRTI

    "A98G",                 # ~2-fold NVP/EFV/RPV/DOR alone
    "K101E",                # 3-10-fold NVP; ~2-fold EFV/ETR/RPV/DOR alone
    "K101P",                # >20-fold NVP/EFV/RPV alone; ~5-fold ETR alone
    "K103N",                # ~50-fold NVP; ~20-fold EFV alone
    "K103S",                # ~5-fold EFV; >20-fold NVP alone
    "K103H",                # ~20-fold NVP/EFV alone
    "K103T",                # >20-fold NVP alone
    "V106A",                # ~50-fold NVP; ~5-fold EFV; ~10-fold DOR alone
    "V106M",                # >30-fold NVP/EFV alone; ~3-fold DOR
    "V108I",                # ~2-fold NVP/EFV alone
    "E138K",                # ~2-fold RPV/ETR alone
    "E138A",                # ~2-fold ETR/RPV alone
    "E138Q", "E138G",       # ~3-fold NVP; ~2-fold ETR/RPV alone
    "E138R",                # Slightly more than E138Q/G alone
    "V179D",                # ~2-fold NVP/EFV/ETR/RPV alone
    "Y181C",                # >50-fold NVP; ~5-fold ETR; ~3-fold RPV alone
    "Y181I", "Y181V",       # >50-fold NVP; 0-15-fold ETR/RPV alone
    "Y181S",                # ~30-fold NVP alone
    "Y188L",                # >50-fold NVP/EFV/DOR; ~5-fold RPV alone
    "Y188C",                # >30-fold NVP alone
    "Y188H",                # ~5-fold NVP/EFV alone
    "G190A",                # >50-fold NVP; 5-10-fold EFV alone
    "G190S",                # >50-fold NVP/EFV alone; variable DOR
    "G190E",                # >100-fold NVP/EFV; >10-fold ETR/RPV/DOR alone
    "G190Q",                # >100-fold EFV/NVP alone
    "G190C", "G190T", "G190V",  # High-level NVP/EFV resistance alone
    "M230L",                # >10-fold NVP/EFV/DOR; up to 5-fold ETR/RPV alone
    "M230I",                # ~10-fold NVP; ~5-fold EFV; ~2-3-fold ETR/RPV alone
    "K238T",                # ~5-fold NVP alone
    "Y318F",                # ~11-fold DOR alone
    "N348I",                # ~2-fold NVP/EFV alone (also in NRTI primary)
]

NNRTI_accessory_mutations = [
    # Alone have little/no effect; require combination with other NNRTI DRMs

    "V90I",                 # Little if any solo reduction; contributes in ETR score
    "L100I",                # Rarely occurs in isolation; major effects always in combination
    "L100V",                # Extremely rare; assessed in very few isolates
    "K101H",                # Contributes NVP/EFV reduction only in combination
    "K101Q",                # Minimal if any detectable solo effect
    "K101N", "K101A", "K101T",  # Not well studied
    "K101R",                # Polymorphic; not selected by NNRTIs; no effect
    "K103R",                # No solo effect; synergizes with V179D ~15-fold NVP/EFV
    "K103E", "K103Q",       # Not selected by NNRTIs; no reduced susceptibility
    "V106I",                # Little if any solo effect; contributes in combination
    "V106L",                # Rare; not well characterized
    "V179E",                # Weakly selected; effect not well defined
    "V179F",                # No solo effect; >10-fold ETR/RPV only with Y181C
    "V179I",                # Highly polymorphic; no clear solo resistance effect
    "V179L",                # Minimal solo effect; listed in RPV package insert
    "V179T",                # No solo reduced susceptibility; contributes in ETR score
    "Y181F", "Y181G",       # Rare transitional/partial revertant mutations
    "Y188F",                # Rare transitional/partial revertant mutation
    "H221Y",                # Minimal solo effects; usually occurs with Y181C
    "P225H",                # Usually with K103N; K103N+P225H >50-fold NVP/EFV
    "F227L",                # Little solo effect; high-level resistance only with V106A
    "F227C",                # Usually in combination; solo effect not documented alone
    "F227I", "F227V",       # Extremely rare; not studied phenotypically
    "F227Y",                # Usually with other DRMs; few solo data
    "Y232H",                # Nearly always in combination; contributes NVP/EFV in combination
    "L234I",                # Solo NNRTI effect unknown; >100-fold DOR only with V106A
    "P236L",                # ~4-fold NVP alone; but primarily known for delavirdine
    "K238N",                # Minimal solo effect
]