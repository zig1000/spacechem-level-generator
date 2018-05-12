#!/usr/bin/env python
# -*- coding: utf-8 -*-

import periodictable

# Default elements by bond count: He, H, O, B, C, N, S, Cl, Os
# We only really need other elements if we're making a fusion level
default_bond_count_elements = { 0:'He', 1:'H', 2:'O', 3:'B', 4:'C', 5:'N', 6:'S', 7:'Cl', 8:'Os' }

basic_elements = [default_bond_count_elements[c] for c in range(1, 9)] # Skip noble gases

max_bonds = {
    1: 1,
    2: 0,
    3: 1,
    4: 2,
    5: 3,
    6: 4,
    7: 5,
    8: 2,
    9: 1,
    10: 0,
    11: 1,
    12: 2,
    13: 4,
    14: 4,
    15: 5,
    16: 6,
    17: 7,
    18: 0,
    19: 1,
    20: 2,
    21: 3,
    22: 4,
    23: 5,
    24: 6,
    25: 7,
    26: 6,
    27: 5,
    28: 4,
    29: 4,
    30: 2,
    31: 3,
    32: 4,
    33: 5,
    34: 6,
    35: 7,
    36: 0,
    37: 1,
    38: 2,
    39: 3,
    40: 4,
    41: 5,
    42: 6,
    43: 7,
    44: 8,
    45: 6,
    46: 4,
    47: 3,
    48: 2,
    49: 3,
    50: 4,
    51: 5,
    52: 6,
    53: 7,
    54: 0,
    55: 1,
    56: 2,
    57: 3,
    58: 4,
    59: 4,
    60: 3,
    61: 3,
    62: 3,
    63: 3,
    64: 3,
    65: 4,
    66: 3,
    67: 3,
    68: 3,
    69: 3,
    70: 3,
    71: 3,
    72: 4,
    73: 5,
    74: 6,
    75: 7,
    76: 8,
    77: 6,
    78: 6,
    79: 5,
    80: 4,
    81: 3,
    82: 4,
    83: 5,
    84: 6,
    85: 7,
    86: 0,
    87: 1,
    88: 2,
    89: 3,
    90: 4,
    91: 5,
    92: 6,
    93: 7,
    94: 7,
    95: 6,
    96: 4,
    97: 4,
    98: 4,
    99: 3,
    100: 3,
    101: 3,
    102: 3,
    103: 3,
    104: 4,
    105: 5,
    106: 6,
    107: 7,
    108: 8,
    109: 6,

    # Greek elements
    200: 12,
    201: 12,
    202: 12,
    203: 12

    # Australium
    # Not included for now because it follows same fission rules as gold which will be problematic
    #204: 5
}

# Symbol -> atomic # lookups
atomic_numbers = { periodictable.elements[i].symbol: i for i in range(1, 110) }
atomic_numbers.update({ 'Θ':200, 'Ω':201, 'Σ':202, 'Δ':203 })
symbols = { atomic_numbers[s]: s for s in atomic_numbers.keys() }

# Commonly used list:
non_noble_symbols = [s for n, s in symbols.items() if max_bonds[n] != 0]

# Unused for now but could be useful
elements_with_bond_count = { k: [e for e in max_bonds.keys()
                                 if max_bonds[e] == k]
                             for k in range(9) + [12] }
