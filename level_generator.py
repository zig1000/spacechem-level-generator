#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import copy
import fractions
import random

import elements_data
from helpers import *

class RandomizationSettings(argparse.Namespace):
    '''Container class for args affecting the level randomization.'''
    pass

class ElementPicker:
    '''Class to assist with randomly picking elements from a given selection, with a chance to
    reuse a previously picked element which is independent of the number of available elements.
    Args:
        elements: The collection of elements to choose from. Default all non-noble elements.
        prob_new_element: The probability of picking a previously unpicked element. Default 0.5
        prob_decay_factor: Each time a new element is picked other than the first, prob_new_element
                           is multiplied by this factor. Default 1 (no decay).
    '''
    def __init__(self, elements=None, prob_new_element=0.5, prob_decay_factor=1):
        if elements is None:
            elements = elements_data.non_noble_elements
        self.used_elements = []
        # A shallow copy is sufficient since we never mutate Elements
        self.unused_elements = list(elements)
        self.prob_new_element = prob_new_element
        self.prob_decay_factor = prob_decay_factor

    def pick(self, min_bond_count=0):
        '''Pick an element, optionally with a minimum specified bond count.'''
        if min_bond_count > 0:
            valid_used_elements = [e for e in self.used_elements if e.max_bonds >= min_bond_count]
            valid_unused_elements = [e for e in self.unused_elements if e.max_bonds >= min_bond_count]
        else:
            valid_used_elements = self.used_elements
            valid_unused_elements = self.unused_elements

        if (valid_unused_elements and random.random() < self.prob_new_element) \
           or not valid_used_elements:
            # If this wasn't the first element picked, decay the new element probability
            if self.used_elements:
                self.prob_new_element *= self.prob_decay_factor

            element = random.choice(valid_unused_elements)
            self.used_elements.append(element)
            self.unused_elements.remove(element)
        else:
            element = random.choice(valid_used_elements)
        return element


def min_split_formula(formula, max_formulas=None):
    '''Given a formula, split it into an integer linear combination of valid (4x4 zone) formulas,
    such that the formulas have a minimum total size in atoms. Output a list of non-empty formulas
    of this allocation.
    If multiple such allocations are possible, choose one of them randomly.

    Max 1 formula: Divide the formula by its GCD.

    Max 2 formulas: Divide the formula by its GCD, then test each value of n from 1 to the highest
    count of any element. The first output formula will consist of the dividends of each element
    count divided by n, and the second formula will consist of the remainders (and divided by the
    GCD of these remainders).
    ~ O(kN) where k is the # of unique elements and N is the LCM of the highest count element.

    Max F formulas: Recurse on the remainder formula.

    Runs in ~ O(k(N^(f-1))) where:
    - k is the # of unique elements
    - N is the LCM of the highest count element
    - f is the maximum # of formulas allowed

    The pure algorithm isn't perfect, e.g. C3H4O2 should give 2*H2O + 3*C, but instead gives
    2*CH2O + 1*C. A partial fix is implemented; in cases where the second formula contains a lone
    element, that element is removed from the first formula when possible (fixing the example but
    not more complex cases). However this is probably an NP-hard ILP anyway so this approximation is
    considered 'good enough'.
    '''
    # TODO: This currently breaks in several cases of weird bond counts, e.g. if given two different
    #       noble gases it tries to put them all in the first input zone and fails.
    if not formula:
        raise Exception("min_split_formula called with empty formula")

    if max_formulas is None:
        max_formulas = 2

    # Divide the formula by its GCD
    formula = copy.copy(formula) # Don't mutate the input
    formula.divide_by_gcd()

    if max_formulas == 1:
        if not formula.is_valid():
            raise Exception('Cannot fit formula {0} into one zone'.format(formula))
        return [formula]

    best_splits = []
    min_size = 49
    for a in range(1, max(formula.values()) + 1):
        this_split = []
        # Calculate f = a*f1 + b*f2 (+ c*f3 ...) while checking that all sub-formulas are valid
        formula_1 = Formula({element: count / a for element, count in formula.items()
                             if count / a != 0})
        if not formula_1.is_valid():
            continue

        formula_2 = Formula({element: count % a for element, count in formula.items()
                             if count % a != 0})

        # Attempt to partially fix the greediness of the first formula
        if len(formula_2.elements()) == 1:
            e = formula_2.elements()[0]
            if e in formula_1:
                del formula_1[e]
        this_split.append(formula_1)

        if formula_2:
            if max_formulas == 2:
                if not formula_2.is_valid():
                    continue
                this_split.append(formula_2)
            else:
                # Recurse - this gets very costly very fast, but we only need to go up to 3
                this_split.extend(min_split_formula(formula_2, max_formulas=max_formulas-1))

        # Calculate the size of this split
        size = 0
        for f in this_split:
            size += f.num_atoms()

        if size <= min_size:
            if size < min_size:
                best_splits = []
                min_size = size
            best_splits.append(this_split)

    if not best_splits:
        raise Exception("Could not find valid split of formula {0}".format(formula))

    return random.choice(best_splits)

def randFormula(num_atoms=None,
                elements=None,
                element_diversity_factor=0.6):
    '''Generate a random molecule formula. Doesn't make nice molecules so currently unused.'''
    if num_atoms is None: # If size isn't specified keep adding atoms with 75% probability
        num_atoms = 1
        while num_atoms <= 16 and random.random() < 0.75:
            num_atoms += 1

    # ElementPicker defaults to non-nobles, but we'll allow them if num_atoms is 1
    if elements is None and num_atoms == 1:
        elements = elements_data.elements

    formula = Formula()
    tries = 0
    while not formula.is_valid():
        if tries == 500:
            raise ValueError(("Could not generate valid formula with {0} atoms from {1} in 500 "
                              + "tries.").format(num_atoms, elements))

        formula = Formula()
        element_picker = ElementPicker(elements,
                                       prob_new_element=element_diversity_factor,
                                       prob_decay_factor=element_diversity_factor)
        for _ in range(num_atoms):
            formula[element_picker.pick()] += 1
        tries += 1

    return formula

def randMoleculeFromFormula(formula):
    '''Generate a Molecule from a Formula. Raise an exception if the given formula cannot fit within
    a 4x4 zone (e.g. C3H8 is impossible to construct within a 4x4 grid).
    '''
    if not formula.is_valid():
        raise Exception("Cannot construct a molecule with formula {0}".format(formula))

    attempts = 0
    while True:
        if attempts >= 1000:
            raise Exception("No valid molecule found after 1000 retries for formula: {0}"
                            .format(formula))
            # TODO: Since we know the formula is valid, instead warn the user and then resort to the
            #       biased but guaranteed to work construction method (build middle-out using
            #       highest bond count atoms first).

        atom_elements = list(formula.atoms())
        random.shuffle(atom_elements)

        molecule = Molecule()
        while atom_elements:
            if molecule.open_bonds == 0 and len(molecule) > 0:
                # If we've accidentally blocked ourselves off, retry on a new molecule
                break

            # If it's the last atom or placing this atom can't block placement of the next atom,
            # proceed freely. It's possible to have only one open position but lots of open bonds,
            # or only one open bond but lots of open positions, so check carefully
            open_positions = molecule.open_positions()
            if not open_positions:
                raise Exception('Error in call to randMoleculeFromFormula with formula: '
                                + str(formula) + '\nConstructed partial molecule with open bonds'
                                + ' but no open positions: ' + str(molecule))
            elif (len(open_positions) > 1 and molecule.open_bonds != 1) \
                 or len(atom_elements) == 1:
                i = 0
            else:
                # Given that we only have one available bond or position to add the new atom (and it
                # isn't the last atom), avoid adding an atom that only accepts one bond, unless it's
                # the final atom to be added.
                # Walk through the list until we find a 2+ bondable atom
                # TODO: inefficient af
                retry = False # For nested break
                for j, e in enumerate(atom_elements):
                    if e.max_bonds >= 2:
                        i = j
                        break
                    elif j == len(atom_elements) - 1:
                        # If we have nothing but bond-1 atoms left, we have to retry
                        retry = True
                        break
                if retry:
                    break

            # O[1] list pop via swapping to end
            atom_elements[i], atom_elements[-1] = atom_elements[-1], atom_elements[i]
            atom = Atom(element=atom_elements.pop(), pos=random.choice(open_positions))

            # Add a single bond to the new atom, on a side facing a neighbor
            if len(molecule) > 0:
                neighbor_posns = atom.pos.neighbors()
                neighbors = [molecule[neighbor_posn] for neighbor_posn in neighbor_posns \
                             if molecule[neighbor_posn] is not None
                             and molecule[neighbor_posn].remaining_bonds() > 0]
                neighbor = random.choice(neighbors)
                if neighbor.col < atom.col:
                    atom.left_bonds = 1
                elif neighbor.col > atom.col:
                    atom.right_bonds = 1
                elif neighbor.row < atom.row:
                    atom.top_bonds = 1
                else:
                    atom.bottom_bonds = 1

            molecule.add_atom(atom)
        attempts += 1

        # If any atoms failed to be added, retry
        if atom_elements:
            continue

        # Otherwise, now that we've successfully added all the atoms, randomly add more bonds
        # For now, just walk through every atom and give a 50% chance to increase each of its bonds
        # (so max +2 bonds between any two atoms, and no need to worry about the 3-bond limit)
        added_atoms = molecule.atoms
        random.shuffle(added_atoms)
        for atom in added_atoms:
            if atom.remaining_bonds() == 0:
                continue
            neighbor_posns = atom.pos.neighbors()
            random.shuffle(neighbor_posns)
            # Identify all the neighbors we can bond to
            neighbors = [molecule[neighbor] for neighbor in neighbor_posns \
                         if molecule[neighbor] is not None
                         and molecule[neighbor].remaining_bonds() > 0]
            for i, neighbor in enumerate(neighbors):
                if random.random() < 0.5:
                    if neighbor.col < atom.col:
                        atom.left_bonds += 1
                        neighbor.right_bonds += 1
                    elif neighbor.col > atom.col:
                        atom.right_bonds += 1
                        neighbor.left_bonds += 1
                    elif neighbor.row < atom.row:
                        atom.top_bonds += 1
                        neighbor.bottom_bonds += 1
                    else:
                        atom.bottom_bonds += 1
                        neighbor.top_bonds += 1
                if atom.remaining_bonds() == 0:
                    break

        return molecule

def randMolecule(num_atoms=None, elements=None):
    '''Generate a random molecule. This doesn't make nice molecules so it is currently unused.'''
    formula = randFormula(num_atoms=num_atoms, elements=elements)
    return randMoleculeFromFormula(formula)

def randFullyBondedMolecule():
    pass # TODO

def randSymmetricalMolecule(settings):
    '''Generate a random symmetric molecule.
    It can either be rotationally symmetrical around a point (90 or 180 degree rotational symmetry),
    or else it can be reflexively symmetrical across an axis.
    '''
    elements = settings.elements
    if elements is None:
        elements = elements_data.elements
    elements = [e for e in elements if e.max_bonds >= 4] # TODO: FIX ME

    element_picker = ElementPicker(elements=elements, prob_new_element=0.5)

    # Calculate a minimum size for our molecule based on difficulty
    target_num_atoms = 2*(settings.difficulty + 1) + 4*settings.large_output

    # Randomly decide on a symmetry to determine the 'twins' of any given grid position.
    choice = random.randint(1, 6) # not 0-indexed, sue me
    if choice == 1:
        # rotational symmetry around the center of the grid
        if random.random() < 0.5:
            symmetry_type = 'rotate 180'
        else:
            symmetry_type = 'rotate 90'
        # Invalid GridPos, but convenient
        _pivot = GridPos(1.5 + 2*settings.large_output, 1.5,
                         large_output=settings.large_output)
    elif choice == 2:
        # rotational symmetry (180-degree) around the center of one of the grid's cell edges
        symmetry_type = 'rotate 180'
        # We'll randomly pick one of the four edges near the middle of the zone
        _pivot = random.choice([GridPos(1 + 2*settings.large_output, 1.5,
                                        large_output=settings.large_output),
                                GridPos(1.5 + 2*settings.large_output, 2,
                                        large_output=settings.large_output),
                                GridPos(2 + 2*settings.large_output, 1.5,
                                        large_output=settings.large_output),
                                GridPos(1.5 + 2*settings.large_output, 1,
                                        large_output=settings.large_output)])
    elif choice == 3:
        # rotational symmetry around the center of one of the grid's middle cells (180 or 90)
        _pivot = GridPos(random.randint(1, 2) + 2*settings.large_output, random.randint(1, 2),
                         large_output=settings.large_output)
        if random.random() < 0.5:
            symmetry_type = 'rotate 180'
        else:
            symmetry_type = 'rotate 90'
    elif choice == 4:
        # Horizontal reflexive symmetry (through a vertical axis)
        symmetry_type = 'horizontal'
        if random.random() < 0.5:
            # Vertical center axis of the zone
            axis = 1.5
        else:
            # Vertical center axis of either of the two middle columns
            axis = random.randint(1, 2)

        def twins_raw(pos):
            '''Given a position return its symmetrical mappings.'''
            return [GridPos(pos.row, int(2*axis - pos.col),
                            large_output=settings.large_output)]
    elif choice == 5:
        # Vertical reflexive symmetry (through a horizontal axis)
        symmetry_type = 'vertical'
        if random.random() < 0.5:
            # Horizontal center axis of the zone
            axis = 1.5
        else:
            # Horizontal center axis of either of the two middle rows
            axis = random.randint(1,2)

        def twins_raw(pos):
            '''Given a position return its symmetrical mappings.'''
            return [GridPos(int(2*axis - pos.row), pos.col,
                            large_output=settings.large_output)]
    elif choice == 6:
        # Diagonal reflexive symmetry across either diagonal of the 4x4 grid
        symmetry_type = 'diagonal'

        if random.random() < 0.5:
            # top-left to bottom-right
            def twins_raw(pos):
                '''Given a position return its symmetrical mappings.'''
                return [GridPos(pos.col, pos.row,
                                large_output=settings.large_output)]
        else:
            # bottom-left to top-right
            def twins_raw(pos):
                '''Given a position return its symmetrical mappings.'''
                return [GridPos(3 - pos.col, 3 - pos.row,
                                large_output=settings.large_output)]

    if symmetry_type == 'rotate 180':
        def twins_raw(pos):
            '''Given a position return its raw symmetrical mappings.'''
            return [GridPos(int(2*_pivot.row - pos.row),
                            int(2*_pivot.col - pos.col),
                            large_output=settings.large_output)]
    elif symmetry_type == 'rotate 90':
        # 90-degree rotational symmetry
        def twins_raw(pos):
            return [GridPos(int(_pivot.row + pos.col - _pivot.col),
                            int(_pivot.col + _pivot.row - pos.row),
                            large_output=settings.large_output),
                    GridPos(int(2*_pivot.row - pos.row),
                            int(2*_pivot.col - pos.col),
                            large_output=settings.large_output),
                    GridPos(int(_pivot.row + _pivot.col - pos.col),
                            int(_pivot.col + pos.row - _pivot.row),
                            large_output=settings.large_output)]

    def twins(pos):
        '''Return a posn's twins, cleaned of duplicates and None if any twin is out-of-bounds
        '''
        _twins = twins_raw(pos)
        for twin in _twins:
            if not twin.is_valid():
               return None
            # If pos is the pivot, it is its own twin
            if twin == pos:
                return []
        return _twins

    molecule = Molecule(large_output=settings.large_output)
    # TODO: Need proper methods for making symmetrical molecules of target sizes
    # For now just terminate as soon as the molecule has no unconnected atoms
    while not molecule.is_connected() or len(molecule) < target_num_atoms:
        available_positions = molecule.open_positions()
        # Remove any positions that have no symmetrical twins inside the 4x4 grid
        for pos in copy.copy(available_positions):
            if twins(pos) is None:
                available_positions.remove(pos)

        # If can't symmetrically add any more atoms but we didn't reach our target size, settle
        # TODO: rig symmetry choice, e.g. 90-degree rotational symmetry is limited to 9 atoms in
        #       a large output which isn't enough for difficulty 3 as is.
        if not available_positions:
            break

        # Grab a random position
        posn = random.choice(available_positions)
        posns = [posn] + twins(posn)

        # Pick a new element to add
        element = element_picker.pick()

        # Add an atom of that element to each of the symmetrical positions
        for p in posns:
            atom = Atom(element=element, pos=p)

            # TODO: For now, just going to add single bonds to all neighbors with available bonds,
            #       and use elements with 4+ max bonds
            # To be proper, I need to write a bunch of shit to make the bonds symmetrical too, as
            # well as a way to ensure I don't break symmetry or connectedness due to max bond
            # counts
            neighbors = [molecule[p] for p in atom.pos.neighbors() \
                         if molecule[p] is not None
                         and molecule[p].remaining_bonds() > 0]
            # Note that the twins might be neighbors of each other
            for neighbor in neighbors:
                if neighbor.col < atom.col:
                    atom.left_bonds = 1
                elif neighbor.col > atom.col:
                    atom.right_bonds = 1
                elif neighbor.row < atom.row:
                    atom.top_bonds = 1
                else:
                    atom.bottom_bonds = 1
            molecule.add_atom(atom)

    return molecule

def randPolymer(settings):
    '''Generate a random polymer Molecule. That is, a molecule constructed mostly from repeating
    segments. Not necessarily symmetric.

    At least one of the allowed elements must have a max bond count of 2 or more.
    '''
    if settings.elements is None:
        elements = elements_data.non_noble_elements
    else:
        # Sanitize input
        elements = [e for e in settings.elements if e.max_bonds != 0]

    element_picker = ElementPicker(elements=elements, prob_new_element=0.5)

    molecule = Molecule(large_output=settings.large_output)

    # Randomize the size and quantity of segments (at least 2 segments)
    if difficulty != 0 and random.random() < 0.25:
        segment_height = 2
    else:
        segment_height = 1
    max_height = 4 + 4*settings.large_output # Python, baby
    num_segments = random.randint(2, max_height/segment_height)

    # Determine which elements available to us allow an acceptable # of bonds to connect the
    # segments of the polymer
    if difficulty == 0:
        min_core_bonds = 2
    elif 1 <= difficulty <= 2:
        min_core_bonds = 3
    else:
        min_core_bonds = 4
    if not any(e.max_bonds >= min_core_bonds for e in elements):
        raise Exception("Cannot construct polymer of difficulty {0} from elements {1}".format(difficulty, elements))

    # Construct the top segment, with the top row randomly based on #/size of segments
    top_row = random.randint(0, max_height - num_segments*segment_height)
    # Determine the column of the 'core' that connects the segments, biased towards the center.
    core_col = random.choice([0, 1, 1, 1, 2, 2, 2, 3])

    for i in range(segment_height):
        core_atom = Atom(element_picker.pick(min_bond_count=min_core_bonds),
                         GridPos(top_row + i, core_col, large_output=settings.large_output))

        # TODO: Second core atom of the segment shouldn't have to be directly
        #       bonded to the first core atom.
        #       E.g. C-O  is a valid 2-tall polymer segment with C and B as the chainable core.
        #              |
        #            B-O
        # But for now, we'll just single-bond them
        if i > 0:
            core_atom.top_bonds = 1

        molecule.add_atom(core_atom)

    # Randomly attach atoms to the core, with the expected # of atoms added based on difficulty
    if settings.difficulty == 0:
        num_atoms = random.randint(0, 1)
    elif settings.difficulty == 1:
        num_atoms = random.randint(2, 3)
    elif settings.difficulty == 2:
        num_atoms = random.randint(4, 5)
    else:
        num_atoms = 6 # Max for 2 segments

    while num_atoms > 0:
        available_positions = [p for p in molecule.open_positions()
                               if p.row == top_row \
                                  or (segment_height == 2 and p.row == top_row + 1)]
        if not available_positions:
            break

        # TODO: Rig min bonds of element to ensure num_atoms is met
        atom = Atom(element_picker.pick(), random.choice(available_positions))

        # Randomly add bonds to neighbors, to a minimum of 1 bond
        neighbors = [p for p in atom.pos.neighbors() \
                     if molecule[p] is not None
                     and molecule[p].remaining_bonds() > 0]
        for i, neighbor in enumerate(neighbors):
            if atom.remaining_bonds() > 0 \
               and (i == 0 or random.random() < 0.5):
                # TODO: store bonds/report nighbors better so I don't have to do this rigamarole
                if neighbor.col < atom.col:
                    atom.left_bonds = 1
                elif neighbor.col > atom.col:
                    atom.right_bonds = 1
                elif neighbor.row < atom.row:
                    atom.top_bonds = 1
                else:
                    atom.bottom_bonds = 1

        molecule.add_atom(atom)
        num_atoms -= 1

    # Now that we have the first segment, copy it the needed # of times
    segment = copy.deepcopy(molecule)
    # Add a bond to connect to the previous segment
    # TODO: Allow for more than a single core bond between segments
    segment[GridPos(top_row, core_col, large_output=settings.large_output)].top_bonds = 1
    for i in range(num_segments - 1):
        segment.shift(rows=segment_height)
        molecule.add_molecule(copy.deepcopy(segment))

    # Include a chance to randomly add 1-2 atoms to any open ends of the chain
    # We'll connect them only to the spine of the polymer/each other to be kind
    bottom_row = top_row + num_segments*segment_height - 1
    open_rows = [r for r in range(max_height) if r < top_row or r > bottom_row]
    num_extra_atoms = random.randint(0, 2) # To extend past 2 the code below needs changing
    for i in range(num_extra_atoms):
        # We'll do a bit of cheaty hard-coding to make sure we only connect to the spine
        if i == 0:
            available_positions = [p for p in molecule.open_positions()
                                   if p.col == core_col]
        elif i == 1:
            # Borrow the atom from the first loop (highly unkosher) to make sure we only attach
            # to the prior extra atom or the other end of the spine
            available_positions = [p for p in molecule.open_positions()
                                   if p.col == core_col \
                                      or (atom.remaining_bonds() > 0 and p in atom.pos.neighbors())]
        if not available_positions:
            break

        atom = Atom(element_picker.pick(), random.choice(available_positions))

        neighbors = [p for p in atom.pos.neighbors() \
                     if p.col == core_col # If num_extra_atoms is allowed to be 3 this is a problem
                        and molecule[p] is not None
                        and molecule[p].remaining_bonds() > 0] # Redundant safety check
        for neighbor in neighbors:
            if neighbor.col < atom.col:
                atom.left_bonds = 1
            elif neighbor.col > atom.col:
                atom.right_bonds = 1
            elif neighbor.row < atom.row:
                atom.top_bonds = 1
            else:
                atom.bottom_bonds = 1

        molecule.add_atom(atom)
        num_extra_atoms -= 1

    return molecule

def randNiceMolecule(settings):
    '''Generate a random 'nice' molecule - currently, either symmetrical or a Polymer.'''
    # TODO: Optionally specify a number of atoms and restrict symmetry/polymers accordingly
    # TODO: When asked for large output, ensure the molecule is of a decent size
    # TODO: when asked for a large output and only 1 input zone, bias the generation algorithm
    #       accordingly
    if random.random() < 0.5:
        return randSymmetricalMolecule(settings)
    else:
        return randPolymer(settings)

def randNiceMoleculeFromFormula(formula, settings):
    '''Given a formula, try to build a 'nice' molecule.'''
    pass # TODO

def randMolecules(num_molecules, settings):
    '''Generate a specified # of random molecules. Optionally specify difficulty, available
    elements, and whether to use a large output zone.
    '''
    return [randNiceMolecule(settings) for _ in range(num_molecules)]

def randInputsFromOutputs(molecules, settings):
    '''Given a list of molecules, generate a list of minimally-sized molecules that could have
    generated them. We'll refer to the given molecules as outputs and the generated molecules as
    inputs, though how they are actually used is irrelevant.
    '''
    output_molecules = []
    balanced_reaction_formula = Formula()

    for molecule in molecules:
        # For the purposes of calculating a balanced reaction, we can effectively reduce a
        # molecule to its LCM (lowest common molecule :P)
        # e.g. C2H4 is effectively equivalent to CH2 in terms of creating a wasteless reaction
        lcm_formula = copy.copy(molecule.formula)
        lcm_formula.divide_by_gcd()
        balanced_reaction_formula += lcm_formula

    # Determine a minimally-sized set of input formulas which satisfy this reaction formula
    input_formulas = min_split_formula(balanced_reaction_formula,
                                       max_formulas=settings.num_inputs)

    # If we didn't reach the target # of molecules (if specified), we need to spread the minimal
    # formulas across the required inputs
    if settings.num_inputs is not None and len(input_formulas) != settings.num_inputs:
        # Identify the largest input formula
        largest_formula = max(input_formulas, key=lambda f: f.num_atoms())

        # If there aren't enough atoms in the min split to cover all inputs, duplicate whichever
        # contains the element with the highest count in the balanced formula
        if largest_formula.num_atoms() < settings.num_inputs:
            bottleneck_element = max(balanced_reaction_formula.elements(),
                                     key=lambda e: balanced_reaction_formula[e])
            for formula in input_formulas:
                if bottleneck_element in formula:
                    bottleneck_formula = formula
                    break
            for _ in range(settings.num_inputs - len(input_formulas)):
                input_formulas.append(copy.copy(bottleneck_formula))
        else:
            # If only N (N < max_formulas) formulas were used, there were at most N unique element
            # counts in the balanced formula. Randomly split up the largest formula to fill our quota.
            # It's possible for splitting an output up to invalidate its molecule so allow re-attempts.
            split_attempts = 0
            new_input_formulas = []
            atoms = largest_formula.atoms()
            while not (new_input_formulas
                       and all(formula.is_valid() for formula in new_input_formulas)):
                if split_attempts >= 1000:
                    raise Exception(('Could not find valid split of atoms {0} across {1} molecules'
                                     + 'in 1000 attempts').format(atoms, ))

                random.shuffle(atoms)
                # Remove chunks from the start of the atom list iteratively
                new_input_formulas = []
                remaining_atoms = copy.copy(atoms)
                for remaining_formulas in range(settings.num_inputs - len(input_formulas), 0, -1):
                    # Leave at least 1 atom per remaining formula
                    # Pivot doesn't include itself in the draw so start from idx at least 1.
                    pivot = random.randint(1, len(remaining_atoms) - remaining_formulas)
                    drawn_atoms, remaining_atoms = remaining_atoms[:pivot], remaining_atoms[pivot:]
                    new_input_formulas.append(Formula(drawn_atoms))
                new_input_formulas.append(Formula(remaining_atoms))

                split_attempts += 1
            # Delete the original formula and add the newly split formulas
            input_formulas.pop(input_formulas.index(largest_formula))
            input_formulas.extend(new_input_formulas)

    return [randMoleculeFromFormula(formula) for formula in input_formulas]

def randResearchLevel(settings, verbose=False):
    '''Generate a random Research level and return the mission code.
    If specified, also pretty-print the created molecules.'''
    level = ResearchLevel()

    # Make modifications to the settings as needed
    if settings.difficulty is None:
        settings.difficulty = random.randint(0, 2)
    level['difficulty'] = settings.difficulty

    # Determine how many output zones we're using
    level['has-large-output'] = settings.large_output
    if settings.large_output:
        settings.num_outputs = 1
        # As a quick hack, set diff of large output levels from 1-4 instead of 0-3.
        settings.difficulty += 1
    elif settings.num_outputs is None:
        settings.num_outputs = random.randint(1, 2)

    # Generate output molecules first
    output_molecules = randMolecules(num_molecules=settings.num_outputs,
                                     settings=settings)

    # Generate input molecules based on the output molecules
    input_molecules = randInputsFromOutputs(molecules=output_molecules,
                                            settings=settings)

    # Randomly swap the inputs with the outputs to reduce bias from the generating algorithm,
    # if this doesn't violate our given settings
    if not settings.large_output \
       and ((settings.num_inputs is None and settings.num_outputs is None)
            or len(input_molecules) == len(output_molecules)) \
       and random.random() < 0.5:
        input_molecules, output_molecules = output_molecules, input_molecules

    # Add input zones to the level
    for zone_idx, molecule in enumerate(input_molecules):
        if len(input_molecules) == 1:
            zone_idx = random.randint(0, 1)

        input_zone_json = {'inputs':[]}
        num_molecules_per_zone = 1 # TODO: random inputs
        for _ in range(num_molecules_per_zone):
            input_json = {}
            input_json['molecule'] = molecule.get_json_str()
            input_json['count'] = 12
            input_zone_json['inputs'].append(input_json)
        level['input-zones'][str(zone_idx)] = input_zone_json

    # Add output zones to the level
    for zone_idx, molecule in enumerate(output_molecules):
        if len(output_molecules) == 1 and not settings.large_output:
            zone_idx = random.randint(0, 1)

        output_zone_json = {}
        output_zone_json['molecule'] = molecule.get_json_str()
        output_zone_json['count'] = 10

        level['output-zones'][str(zone_idx)] = output_zone_json

    # Add 'features' to the level (bonders, fusers, etc.)
    if difficulty == 0:
        level['bonder-count'] = 4
    else:
        level['bonder-count'] = random.choice([2, 4, 4, 8]) # Bias toward 4 bonders
    level['has-sensor'] = random.random() < 0.5
    level['has-fuser'] = settings.fusion
    level['has-splitter'] = settings.fission
    level['has-teleporter'] = settings.difficulty >= 2 and random.random() < 0.5

    if verbose:
        print 'Inputs:'
        for molecule in input_molecules:
            print molecule
        print 'Outputs:'
        for molecule in output_molecules:
            print molecule

    # Return the mission code
    return level.get_code()

def randProductionLevel(settings, verbose=False):
    '''Generate a random Production level.'''
    level = ProductionLevel()

    if difficulty is None:
        difficulty = random.randint(0, 2)
    level['difficulty'] = difficulty

    # Determine how many output zones we're using
    if num_outputs is None:
        num_outputs = random.choice([1, 2, 2, 3]) # bias toward 2 outputs

    # Generate output molecules first
    output_molecules = randMolecules(num_molecules=settings.num_outputs,
                                     settings=settings)

    # Generate input molecules (up to 3 since this is a production level)
    input_molecules = randInputsFromOutputs(output_molecules,
                                            settings=settings)

    # Randomly swap the inputs with the outputs to reduce bias from the generating algorithm
    if ((settings.num_inputs is None and settings.num_outputs is None) \
        or len(input_molecules) == len(output_molecules)) \
       and random.random() < 0.5:
        input_molecules, output_molecules = output_molecules, input_molecules

    # Add the random input zone to the level
    input_zone_json = {'inputs':[]}
    input_json = {}
    input_json['molecule'] = input_molecules[0].get_json_str()
    input_json['count'] = 12
    input_zone_json['inputs'].append(input_json)
    level['random-input-zones']['0'] = input_zone_json

    # Add the fixed input zones to the level
    for zone_idx, molecule in enumerate(input_molecules[1:]):
        level['fixed-input-zones'][str(zone_idx)] = molecule.get_json_str()

    # Add the output zones
    for zone_idx, molecule in enumerate(output_molecules):
        output_zone_json = {}
        output_zone_json['molecule'] = molecule.get_json_str()
        output_zone_json['count'] = 40

        level['output-zones'][str(zone_idx)] = output_zone_json

    # Add features to level
    # I'm leaving 'max-reactors' at its default 6 as set in ProductionLevel's initialization
    level['terrain'] = random.randint(0, 4)
    # By convention, if the inputs contain greek elements the terrain is set to Flidais
    for molecule in input_molecules:
        for element in molecule.formula.elements():
            if element.symbol in ('Θ', 'Ω', 'Σ', 'Δ'):
                level['terrain'] = 5
                break

    # No point having assembly/disassembly if any other reactors are available, and there's a max of
    # 4 types of reactors available anyway, probably for this reason
    if random.random() < 0.02:
        level['has-assembly'] = True
        level['has-disassembly'] = True
    else:
        level['has-starter'] = True # Always available, in case the player wants to comment their
                                    # code by using a strictly worse reactor
        level['has-advanced'] = random.random() < 0.5 # Sensor reactor
        level['has-nuclear'] = settings.fusion or settings.fission
        level['has-superbonder'] = random.random() < 0.25

    level['has-recycler'] = random.random() < 0.5

    if verbose:
        print 'Generated research level.'
        print 'Inputs:'
        for molecule in input_molecules:
            print molecule
        print 'Outputs:'
        for molecule in output_molecules:
            print molecule

    return level.get_code()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--view', action='store_true',
                        help="Open Cearn's mission viewer in the browser" \
                             + "(http://coranac.com/spacechem/mission-viewer?mode=editor)")
    parser.add_argument('--print', dest='verbose', action='store_true', default=False,
                        help="Pretty-print the generated level.")

    level_type = parser.add_mutually_exclusive_group()
    level_type.add_argument('-r', '--research', action='store_true', default=True,
                            help="Dummy flag; Research levels are generated by default anyway.")
    level_type.add_argument('-p', '--production', action='store_true', default=False,
                            help="Generate a Production level instead of a Research level.")

    parser.add_argument('-d', '--difficulty', type=int, action="store",
                        choices=[0, 1, 2, 3], default=random.randint(0, 2),
                        help="Set the difficulty of the level.")
    parser.add_argument('--lanky', action="store_true",
                        help="Set difficulty to max")

    element_args = parser.add_mutually_exclusive_group()
    element_args.add_argument('-e', '--elements', action='store', nargs='*', type=str,
            help="Select only from the specified elements in creating the level." \
                 + " May be specified by atomic # or symbol."
                 + " Not currently guaranteed to use all of the available elements.")
    element_args.add_argument('--basic', action='store_true',
            help="Use only 'basic' elements, one for each of the non-0 bond counts:" \
                 + "[H, O, B, C, N, S, Cl, Os]")

    parser.add_argument('-i', '--inputs', action='store', type=int, default=None,
                        choices=[1, 2, 3],
                        help="The # of inputs to use.")
    parser.add_argument('-o', '--outputs', action='store', type=int, default=None,
                        choices=[1, 2, 3],
                        help="The # of outputs to use.")
    parser.add_argument('--large_output', action="store_true",
                        help="Set a research level to have a large output zone")

    parser.add_argument('--fusion', action="store_true", default=False,
                        help="Allow fusion in the level. Default none.")
    parser.add_argument('--fission', action="store_true", default=False,
                        help="Allow fission in the level. Default none.")
    parser.add_argument('--nuclear', action="store_true",
                        help="Allow both fission and fusion in the level. Default none.")
    args = parser.parse_args()

    difficulty = args.difficulty
    if args.lanky:
        difficulty = 3

    elements = None
    if args.elements:
        elements = [elements_data.elements_dict[int(e)] if e.isdigit()
                    else elements_data.elements_dict[e]
                    for e in args.elements]
    elif args.basic:
        elements = elements_data.basic_elements

    if args.nuclear:
        args.fusion = True
        args.fission = True

    settings = RandomizationSettings(difficulty=difficulty,
                                     elements=elements,
                                     num_inputs=args.inputs,
                                     num_outputs=args.outputs,
                                     large_output=args.large_output,
                                     fusion=args.fusion,
                                     fission=args.fission)

    if args.production:
        if args.large_output:
            raise argparse.ArgumentTypeError('Cannot use a large output in a production level')

        code = randProductionLevel(settings=args.settings,
                                   verbose=args.verbose)
    else:
        # Check for conflicts in the given options
        if args.inputs == 3:
            raise argparse.ArgumentTypeError('Too many inputs; did you mean to add --production ?')
        elif args.outputs == 3:
            raise argparse.ArgumentTypeError('Too many outputs; did you mean to add --production ?')
        elif args.outputs == 2 and args.large_output:
            raise argparse.ArgumentTypeError('Cannot have 2 outputs in a large output level')

        code = randResearchLevel(settings=settings,
                                 verbose=args.verbose)

    if args.view:
        import webbrowser
        webbrowser.open('http://coranac.com/spacechem/mission-viewer?mode=editor&code=' + code)

    if args.verbose:
        print 'Level code:\n'
    print code
