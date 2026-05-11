PI_primary_mutations = [
    # Alone or with minimal context reduces susceptibility to ≥1 PI meaningfully

    "L24I",      # Reduces LPV and ATV susceptibility
    "D30N",      # Reduces NFV susceptibility (sole PI, but clear standalone effect)
    "V32I",      # Reduces LPV, ATV, DRV susceptibility; in combo with I47V/A causes intermediate-high resistance
    "I47V",      # Reduces LPV and DRV susceptibility
    "I47A",      # High-level LPV resistance (nearly always with V32I, but the pair is a clear primary driver)
    "G48V",      # Intermediate ATV resistance alone
    "G48M",      # Similar resistance profile to G48V
    "I50V",      # Reduces LPV and DRV susceptibility
    "I50L",      # High-level ATV resistance alone; increases susceptibility to other PIs
    "I54V",      # Reduces susceptibility to each PI except DRV
    "I54M",      # Reduces LPV, ATV, DRV susceptibility
    "I54L",      # Reduces LPV, ATV, DRV susceptibility
    "I54A", "I54T", "I54S",  # Reduces susceptibility to each PI except DRV
    "L76V",      # Reduces LPV and DRV susceptibility (increases ATV susceptibility)
    "V82A",      # Reduces LPV and ATV susceptibility
    "V82T", "V82S",  # Reduces LPV and ATV susceptibility
    "V82F",      # Reduces LPV and DRV susceptibility
    "I84V",      # Reduces LPV, ATV, DRV susceptibility
    "I84A",      # Markedly reduces susceptibility to each PI
    "N88S",      # ~10-fold reduced ATV susceptibility alone
    "N88D",      # Low-level cross-resistance to ATV
    "L90M",      # Reduces ATV and to a lesser extent LPV susceptibility
    "L23I",      # Reduces NFV susceptibility (substrate-cleft mutation with clear solo effect)
]

PI_accessory_mutations = [
    # Alone have little/no effect; contribute primarily in combination with other PI DRMs

    "L10F",      # Accessory; reduced susceptibility mainly in combination
    "L10I", "L10V",  # Polymorphic accessory; increase replication of resistant viruses
    "L10R", "L10Y",  # Rare; effects not well studied
    "V11I",      # Accessory; effect only in combination (DRV genotypic score)
    "V11L",      # Weakly associated with reduced DRV susceptibility
    "K20R",      # Polymorphic accessory; increases fitness of PI-resistant viruses
    "K20I", "K20M", "K20T", "K20V",  # Accessory PI-selected; not well studied individually
    "L33F",      # Only associated with reduced susceptibility in combination
    "L33I",      # Does not appear to reduce PI susceptibility
    "L33V",      # Polymorphic; not PI-selected, does not reduce susceptibility
    "M36I",      # Increases replication fitness of PI-resistant viruses; not a direct resistance mutation
    "K43T",      # Accessory; effect obscured by co-occurring mutations
    "M46I", "M46L",  # Reduce ATV and LPV susceptibility but primarily as accessory/context-dependent
    "M46V",      # Uncommon; not well studied
    "F53L",      # Accessory; reduced ATV susceptibility only in combination
    "F53Y",      # Accessory; not well studied
    "Q58E",      # Accessory; low-level ATV resistance only in combination
    "A71V", "A71T",  # Polymorphic accessory; increase replication capacity in combination
    "A71I", "A71L",  # Accessory; occur in multi-PI-resistant viruses
    "G73S", "G73T", "G73C", "G73A",  # Accessory; minimally reduced susceptibility in combination
    "G73D", "G73V",  # Rare; not well studied
    "T74P",      # Accessory; minimally reduced ATV/DRV susceptibility in combination
    "T74S",      # Accessory polymorphic in non-B subtypes
    "V82L",      # Rare; effect on other PIs not well characterized
    "V82M",      # Reduces IDV and possibly LPV; effect not well studied
    "V82C",      # Occurs in multi-PI-resistant viruses; not well studied
    "V82I",      # Consensus in subtype G; does not reduce PI susceptibility
    "N83D",      # Possibly contributes to reduced ATV susceptibility in combination only
    "N83S",      # Extremely rare; not well studied
    "I84C",      # Less-marked effect than I84A; context-dependent
    "I85V",      # Minimal if any effect on PI susceptibility
    "N88T", "N88G",  # Extremely rare; minimal effect on PI susceptibility
    "L89V",      # Minimally associated with reduced susceptibility; effect only in combination
    "L89T",      # Effect on PI susceptibility not well studied
]