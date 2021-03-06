#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import copy
import random

import elements_data
from helpers import *

class RandomizationSettings(argparse.Namespace):
    '''Container class for args affecting the level randomization.'''
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # If no nuclear options were specified, add nuclear with 1/3 chance:
        # 1/6 of both + 1/12 of fusion only + 1/12 fission only.
        if self.fusion is None and self.fission is None:
            choice = random.random()
            if choice < 1/6:
                self.fusion = self.fission = True
            elif choice < 1/4:
                self.fusion = True
                self.fission = False
            elif choice < 1/3:
                self.fusion = False
                self.fission = True

        # If large output option wasn't specified, add it with 10% chance
        if self.large_output is None:
            self.large_output = random.random() < 0.1

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
        if not elements:
            raise ValueError("ElementPicker cannot be initialized with empty elements iterable")
        self.used_elements = []
        # Use list() to both copy and standardize the given iterable of elements
        self.unused_elements = list(elements)
        self.prob_new_element = prob_new_element
        self.prob_decay_factor = prob_decay_factor

    def pick(self, min_bond_count=0):
        '''Pick an element, optionally with a minimum specified bond count.'''
        if min_bond_count > 0:
            valid_used_elements = [e for e in self.used_elements if e.max_bonds >= min_bond_count]
            valid_unused_elements = [e for e in self.unused_elements if e.max_bonds >= min_bond_count]
            # Make sure at least one element with the specified bond count exists
            if not valid_unused_elements and not valid_used_elements:
                raise ValueError(f"No elements provided with minimum bond count {min_bond_count}")
        else:
            valid_used_elements = self.used_elements
            valid_unused_elements = self.unused_elements

        if (valid_unused_elements and random.random() < self.prob_new_element) \
                or not valid_used_elements:
            element = random.choice(valid_unused_elements)
            self.used_elements.append(element)
            self.unused_elements.remove(element)

            # If this wasn't the first element picked, decay the new element probability
            if self.used_elements:
                self.prob_new_element *= self.prob_decay_factor
        else:
            element = random.choice(valid_used_elements)
        return element

def mutate_molecule(molecule, settings):
    '''Given a molecule, perform an atomic mutation on it, i.e. do one of:
    * Mutate a single atom's element and update its bonds
    * Remove a single atom (without disconnecting the molecule)
    TBA:
    * Add an atom
    * Modify a bond (if possible)
    '''
    r = random.random()
    if r < 0.5 or len(molecule) == 1:
        # With 50% probability or if there are no spare atoms to remove, mutate an element
        element_picker = ElementPicker(elements=settings.elements)
        # If the molecule is a lone atom, we can only mutate its element
        atom = random.choice(list(molecule))
        atom.set_element(element_picker.pick(min_bond_count=sum(atom.bonds)))
        molecule.update_formula()
        molecule.update_open_bonds()
    else:
        # Attempt to remove random atoms until we find one that can
        # be removed without disconnecting the molecule.
        # Given that the molecule contains at least 2 atoms, we know such an atom exists
        atoms = list(molecule)
        random.shuffle(atoms)

        for atom in atoms:
            new_molecule = copy.deepcopy(molecule)
            new_molecule.remove_atom(atom)

            if new_molecule.is_connected():
                molecule.remove_atom(atom)
                return

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
    Noble gases are also pre-filtered into their own input zones.
    '''
    if max_formulas is None:
        max_formulas = 2

    if not formula:
        raise FormulaValidityError("Formula must not be empty")

    # Get the formula divided by its own GCD (but keep track of the original for logging purposes)
    given_formula = formula
    formula = formula.least_common_formula() # Also ensures we aren't mutating the given formula

    if max_formulas == 1 and formula.is_valid():
        return [formula]
    elif max_formulas <= 1:
        raise FormulaValidityError(f"Cannot fit formula {given_formula} into {max_formulas} zone")

    # If the formula contains any noble gases, filter them into their own input zones then recurse
    # on the remainder
    noble_gases = [e for e in formula.elements() if e.max_bonds == 0]
    if noble_gases:
        best_split = []
        clean_formula = copy.copy(formula)
        for e in noble_gases:
            best_split.append(Formula({e: 1}))
            del clean_formula[e]
        if clean_formula:
            try:
                best_split.extend(min_split_formula(clean_formula,
                                                    max_formulas=max_formulas - len(noble_gases)))
            except FormulaValidityError:
                raise FormulaValidityError(f"Cannot fit formula {given_formula} into {max_formulas} zones")

        random.shuffle(best_split)
        return best_split

    # Loop over possible greedy allocations to the first formula
    best_splits = []
    min_size = 49 # Initialize to larger than 3 full input zones
    for a in range(1, max(formula.values()) + 1):
        this_split = []
        # Calculate f = a*f1 + b*f2 (+ c*f3 ...) while checking that all sub-formulas are valid
        formula_1 = Formula({element: count // a for element, count in formula.items()
                             if count // a != 0})
        if not formula_1.is_valid():
            continue

        formula_2 = Formula({element: count % a for element, count in formula.items()
                             if count % a != 0})

        # Attempt to partially fix the greediness of the first formula
        if len(formula_2.elements()) == 1:
            # If the second formula is a lone element remove that element from the first formula
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
                # Recurse - this gets very costly very fast, but we only call it with up to 3
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
        raise FormulaValidityError(f"Could not find valid split of formula {formula}")

    return random.choice(best_splits)

def randFormula(num_atoms=None,
                elements=None,
                large_output=False,
                element_diversity_factor=0.6):
    '''Generate a random molecule formula. Doesn't make nice molecules so currently unused.
    Keeping it around in case randNiceMoleculeFromFormula is ever implemented and works well.
    '''
    if num_atoms is None: # If size isn't specified keep adding atoms with 75% probability
        num_atoms = 1
        while num_atoms <= 16 and random.random() < 0.75:
            num_atoms += 1

    # ElementPicker defaults to non-nobles, but we'll allow them if num_atoms is 1
    if elements is None and num_atoms == 1:
        elements = elements_data.elements

    formula = Formula()
    tries = 0
    while not formula.is_valid(large_output=large_output):
        if tries == 500:
            raise ValueError("Could not generate valid formula with"
                             + f" {num_atoms} atoms from {elements} in 500 tries.")

        formula = Formula()
        element_picker = ElementPicker(elements,
                                       prob_new_element=element_diversity_factor,
                                       prob_decay_factor=element_diversity_factor)
        for _ in range(num_atoms):
            formula[element_picker.pick()] += 1
        tries += 1

    return formula

def randMoleculeFromFormula(formula, large_output=False):
    '''Generate a Molecule from a Formula. Raise an exception if the given formula cannot fit within
    a 4x4 zone (e.g. C3H8 is impossible to construct within a 4x4 grid).
    '''
    if not formula.is_valid(large_output=large_output):
        raise FormulaValidityError(f"Cannot construct a molecule with formula {formula}")

    attempts = 0
    while True:
        if attempts >= 1000:
            raise FormulaValidityError("No valid molecule found after 1000 retries for formula: "
                                       + str(formula))
            # TODO: Since we know the formula is valid, instead warn the user and then resort to the
            #       biased but guaranteed to work construction method (build middle-out using
            #       highest bond count atoms first).

        # Do a sub-loop attempting to build a single molecule
        # First randomize the order elements are selected to be added
        atom_elements = formula.elements_collection()
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
                raise Exception("Error in call to randMoleculeFromFormula with formula: "
                                + f"{formula}\nConstructed partial molecule with open bonds"
                                + f" but no open positions: {molecule}")
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
                neighbor_dirs = [d for d, p in atom.pos.dirs_and_neighbors()
                                 if molecule[p] is not None
                                    and molecule[p].remaining_bonds() > 0]
                dir = random.choice(neighbor_dirs)
                atom.bonds[dir] = 1

            molecule.add_atom(atom)
        attempts += 1

        # If any atoms failed to be added, retry
        if not atom_elements:
            break

    # Now that we've successfully created a connected molecule, randomly add more bonds.
    # For now, just walk through every atom and give a 50% chance to increase each of its bonds
    # (so max +2 bonds between any two atoms, and no need to worry about the 3-bond limit)
    atoms = list(molecule)
    random.shuffle(atoms)
    for atom in atoms:
        if atom.remaining_bonds() == 0:
            continue
        # Identify all the neighbors we can bond to
        dirs_and_adj_atoms = [(d, molecule[pos]) for d, pos in atom.pos.dirs_and_neighbors()
                              if molecule[pos] is not None
                                 and molecule[pos].remaining_bonds() > 0]
        random.shuffle(dirs_and_adj_atoms)
        for dir, adj_atom in dirs_and_adj_atoms:
            if random.random() < 0.5:
                atom.bonds[dir] += 1
                adj_atom.bonds[opposite_dir(dir)] += 1
            if atom.remaining_bonds() == 0:
                break

    return molecule

def randSymmetricMolecule(settings):
    '''Generate a random symmetric molecule.
    It can either be rotationally symmetric around a point (90 or 180 degree rotational symmetry),
    or else it can be reflexively symmetric across an axis.
    '''
    element_picker = ElementPicker(elements=settings.elements, prob_new_element=0.5)

    # Calculate a min and max size for our molecule based on difficulty
    min_num_atoms = 2 + 3*settings.difficulty + 4*settings.large_output
    max_num_atoms = min(4*(settings.difficulty + 1 + 1*settings.large_output),
                        16 + 16*settings.large_output)

    # Randomly decide on a symmetry to determine the 'twins' of any given grid position.
    # TODO: rig symmetry choice, e.g. 90-degree rotational symmetry is limited to 9 atoms in
    #       a large output which isn't enough for difficulty 3 as is.
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

        def twins(pos):
            '''Given a position return a list of it and its twins.'''
            return [pos, GridPos(pos.row, int(2*axis - pos.col),
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

        def twins(pos):
            '''Given a position return a list of it and its twins.'''
            return [pos, GridPos(int(2*axis - pos.row), pos.col,
                                 large_output=settings.large_output)]
    elif choice == 6:
        # Diagonal reflexive symmetry across either diagonal of the 4x4 grid
        if random.random() < 0.5:
            # top-left to bottom-right axis of symmetry
            symmetry_type = 'diagonal \\'
            def twins(pos):
                '''Given a position return a list of it and its twins.'''
                return [pos, GridPos(pos.col, pos.row,
                                     large_output=settings.large_output)]
        else:
            # bottom-left to top-right
            symmetry_type = 'diagonal /'
            def twins(pos):
                '''Given a position return return a list of it and its twins.'''
                return [pos, GridPos(3 - pos.col, 3 - pos.row,
                                     large_output=settings.large_output)]

    if symmetry_type == 'rotate 180':
        def twins(pos):
            '''Given a position return a list of it and its twins.'''
            return [pos, GridPos(int(2*_pivot.row - pos.row),
                                 int(2*_pivot.col - pos.col),
                                 large_output=settings.large_output)]
    elif symmetry_type == 'rotate 90':
        # 90-degree rotational symmetry
        def twins(pos):
            '''Given a position return a list of it and its twins in clockwise order.'''
            return [pos, GridPos(int(_pivot.row + pos.col - _pivot.col), # 90
                                 int(_pivot.col + _pivot.row - pos.row),
                                 large_output=settings.large_output),
                         GridPos(int(2*_pivot.row - pos.row), # 180
                                 int(2*_pivot.col - pos.col),
                                 large_output=settings.large_output),
                         GridPos(int(_pivot.row + _pivot.col - pos.col), # 270
                                 int(_pivot.col + pos.row - _pivot.row),
                                 large_output=settings.large_output)]

    tries, max_tries = 0, 20
    molecule = Molecule(large_output=settings.large_output)
    while not (min_num_atoms <= len(molecule) <= max_num_atoms):
        tries += 1
        if tries > max_tries:
            raise MoleculeValidityError(f"Could not create Molecule with {symmetry_type} symmetry"
                                        + f" between {min_num_atoms} and {max_num_atoms} atoms size"
                                        + f" within {max_tries} tries.")
        molecule = Molecule(large_output=settings.large_output)

        # Loop until we've got a connected molecule of minimum size, then retry if it got too large
        while len(molecule) < min_num_atoms or not molecule.is_connected():
            # Get all available positions in the molecule, ignoring those which have invalid
            # symmetric twins
            available_positions = [p for p in molecule.open_positions()
                                   if all(twin.is_valid() for twin in twins(p))]
            if not available_positions:
                break

            # For now we just need any element that allows 4 bonds
            element = element_picker.pick(min_bond_count=4)

            # Grab a random position and add atoms of the selected element to it and each of its
            # unique twins
            for pos in set(twins(random.choice(available_positions))):
                atom = Atom(element=element, pos=pos)

                # Add fully-connected bonds to all neighbors
                neighbor_dirs = [d for d, p in atom.pos.dirs_and_neighbors()
                                 if molecule[p] is not None
                                    and molecule[p].remaining_bonds() > 0]
                for dir in neighbor_dirs:
                    atom.bonds[dir] = 1

                molecule.add_atom(atom)

    # Once we have a fully-connected molecule of acceptable size, symmetrically remove bonds without
    # leaving any parts of the molecule unconnected, then update to new elements based on the new
    # bond counts
    element_picker = ElementPicker(elements=settings.elements, prob_new_element=0.5)

    def update_symmetric_bonds(twin_posns, bond_dir, bond_count):
        '''Modify bonds on molecule symmetrically.
        Args:
            twins: List of symmetric twin positions in clockwise order. Note that a position should
                   be listed twice if it is its own twin, or 4 times, etc (makes the code cleaner).
            bond_dir: The direction of the bond to update on the first twin in the list
            bond_count: the new bond count to be set
        '''
        for twin_idx, pos in enumerate(twin_posns):
            if twin_idx == 0:
                new_bond_dir = bond_dir
            elif symmetry_type == 'rotate 90':
                new_bond_dir = (bond_dir + twin_idx) % 4
            elif symmetry_type == 'rotate 180':
                new_bond_dir = (bond_dir + 2*twin_idx) % 4
            elif symmetry_type == 'vertical' and bond_dir in (UP, DOWN):
                new_bond_dir = (bond_dir + 2*twin_idx) % 4
            elif symmetry_type == 'horizontal' and bond_dir in (RIGHT, LEFT):
                new_bond_dir = (bond_dir + 2*twin_idx) % 4
            elif symmetry_type == 'diagonal /':
                if bond_dir in (UP, RIGHT):
                    new_bond_dir = UP + RIGHT - bond_dir # UP <-> RIGHT
                else:
                    new_bond_dir = DOWN + LEFT - bond_dir # DOWN <-> LEFT
            elif symmetry_type == 'diagonal \\':
                new_bond_dir = 3 - bond_dir # UP <-> LEFT and DOWN <-> RIGHT
            else:
                new_bond_dir = bond_dir

            molecule[pos].bonds[new_bond_dir] = bond_count

            # Also update the bond on the neighbor
            neighbor_atom = molecule[pos.neighbor(new_bond_dir)]
            if neighbor_atom is not None:
                neighbor_atom.bonds[opposite_dir(new_bond_dir)] = bond_count


    checked_posns = set()
    for atom in molecule:
        if atom.pos in checked_posns:
            continue

        # Identify the directions this atom has bonds in. If it already has only one bond, we can
        # skip it since removing its bond will disconnect it from the molecule
        atom_bond_dirs = [bond_dir for bond_dir, bond_count in enumerate(atom.bonds)
                          if bond_count != 0]
        if len(atom_bond_dirs) <= 1:
            continue

        # Take it in turn removing a bond from this atom and the corresponding bond on each
        # of its twins (note that when it is its own twin(s) this will neatly take care of symmetry
        # across its own bonds)
        twin_posns = twins(atom.pos)

        # Randomly shuffle the order bonds are attempted to be removed
        random.shuffle(atom_bond_dirs)

        # Don't remove the last bond regardless of the success of the other bonds
        for original_bond_dir in atom_bond_dirs[:-1]:
            update_symmetric_bonds(twin_posns, original_bond_dir, 0)

            # Undo this bond removal if we disconnected the molecule
            if not molecule.is_connected():
                update_symmetric_bonds(twin_posns, original_bond_dir, 1)

        # Once we've done all the bond removals allowed on these twins, update their element to
        # more appropriately match their new bond count
        new_element = element_picker.pick(min_bond_count=sum(atom.bonds))

        for twin_posn in set(twin_posns): # Remove duplicate positions
            molecule[twin_posn].set_element(new_element)
            # Also mark each twin position as visited
            checked_posns.add(twin_posn)

    # Fix the molecule's formula and open bond count after our changes
    molecule.update_formula()
    molecule.update_open_bonds()

    assert molecule.is_connected(), f"Created molecule is disconnected:\n{molecule}"

    return molecule

def randPolymer(settings):
    '''Generate a random polymer Molecule. That is, a molecule constructed mostly from repeating
    segments. Not necessarily symmetric.

    At least one of the allowed elements must have a max bond count of 2 or more.
    '''
    # Sanitize the input a little
    element_picker = ElementPicker(elements=settings.elements, prob_new_element=0.5)

    molecule = Molecule(large_output=settings.large_output)

    # Randomize the size and quantity of segments (at least 2 segments)
    if settings.difficulty != 0 and random.random() < 0.25:
        segment_height = 2
    else:
        segment_height = 1
    max_height = 4 + 4*settings.large_output
    num_segments = random.randint(2, max_height // segment_height)

    # Determine which elements available to us allow an acceptable # of bonds to connect the
    # segments of the polymer
    if settings.difficulty == 0:
        min_core_bonds = 2
    elif 1 <= settings.difficulty <= 2:
        min_core_bonds = 3
    else:
        min_core_bonds = 4
    if not any(e.max_bonds >= min_core_bonds for e in elements):
        raise Exception(f"Cannot construct polymer of difficulty {settings.difficulty}"
                        + f" from elements {elements}")

    # Construct the top segment, with the top row randomly based on #/size of segments
    top_row = random.randint(0, max_height - num_segments*segment_height)
    # Determine the column of the 'core' that connects the segments, biased towards the center.
    core_col = random.choice([0, 1, 1, 1, 2, 2, 2, 3])

    for i in range(segment_height):
        core_atom = Atom(element=element_picker.pick(min_bond_count=min_core_bonds),
                         pos=GridPos(top_row + i, core_col, large_output=settings.large_output))

        # TODO: Second core atom of the segment shouldn't have to be directly
        #       bonded to the first core atom.
        #       E.g. C-O  is a valid 2-tall polymer segment with C and B as the chainable core.
        #              |
        #            B-O
        # But for now, we'll just single-bond them
        if i > 0:
            core_atom.bonds[UP] = 1

        molecule.add_atom(core_atom)

    # Randomly attach atoms to the core, with the expected # of atoms added based on difficulty
    if settings.difficulty == 0:
        num_atoms = random.randint(0, 1)
    else:
        num_atoms = random.randint(2*settings.difficulty - 1, 2*settings.difficulty)

    while num_atoms > 0:
        available_positions = [p for p in molecule.open_positions()
                               if p.row == top_row \
                                  or (segment_height == 2 and p.row == top_row + 1)]
        if not available_positions:
            break

        # TODO: Rig min bonds of element to ensure num_atoms is met
        atom = Atom(element=element_picker.pick(min_bond_count=1),
                    pos=random.choice(available_positions))

        # Randomly add bonds to the atom for its neighbors, to a minimum of 1 bond
        neighbor_dirs = [d for d, p in atom.pos.dirs_and_neighbors()
                         if molecule[p] is not None
                            and molecule[p].remaining_bonds() > 0]
        random.shuffle(neighbor_dirs)
        for i, dir in enumerate(neighbor_dirs):
            if (atom.remaining_bonds() > 0
                    and (i == 0 or random.random() < 0.5)):
                atom.bonds[dir] = 1

        molecule.add_atom(atom)
        num_atoms -= 1

    # Now that we have the first segment, copy it the needed # of times
    segment = copy.deepcopy(molecule)
    # Add a bond to connect to the previous segment
    # TODO: Allow for more than a single core bond between segments
    segment[GridPos(top_row, core_col, large_output=settings.large_output)].bonds[UP] = 1
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
            # Borrow the atom added in the first loop (highly unkosher) to make sure we only attach
            # to the prior extra atom or the other end of the spine
            available_positions = [p for p in molecule.open_positions()
                                   if p.col == core_col \
                                      or (atom.remaining_bonds() > 0 and p in atom.pos.neighbors())]
        if not available_positions:
            break

        atom = Atom(element=element_picker.pick(min_bond_count=1),
                    pos=random.choice(available_positions))

        # Identify neighbouring atoms and add bonds to them
        neighbor_dirs = [d for d, p in atom.pos.dirs_and_neighbors() \
                         if p.col == core_col # If num_extra_atoms is allowed to be 3 this is a problem
                            and molecule[p] is not None
                            and molecule[p].remaining_bonds() > 0] # Redundant safety check
        for dir in neighbor_dirs:
            atom.bonds[dir] = 1

        molecule.add_atom(atom)
        num_extra_atoms -= 1

    return molecule

def randMolecule(settings):
    '''Generate a random 'nice' molecule - currently, either symmetric or a polymer.
    To spice things up, add a chance to slightly mutate an otherwise 'nice' molecule.
    '''
    if settings.symmetric and not settings.polymer:
        molecule = randSymmetricMolecule(settings)
    elif settings.polymer and not settings.symmetric:
        molecule = randPolymer(settings)
    elif random.random() < 0.5:
        molecule = randSymmetricMolecule(settings)
    else:
        molecule = randPolymer(settings)

    # Add a chance to slightly mutate a 'too perfect' molecule
    if random.random() < 0.5:
        mutate_molecule(molecule, settings=settings)

    return molecule

def randNiceMoleculeFromFormula(formula, settings):
    '''Given a formula, try to build a 'nice' molecule.'''
    pass # TODO

def randMolecules(num_molecules, settings):
    '''Generate a specified # of random molecules with the given randomization settings.
    '''
    return [randMolecule(settings) for _ in range(num_molecules)]


def rand_inverse_fission(formula):
    '''Mutate the given formula using a random series of valid inverse fission operations which add
    up to a single element. Return the element which was added (may not be a new element).
    '''
    # Walk through the formula and build up the list of elements that fission could
    # successfully be used on to produce some of the outputs (without waste).
    valid_fissile_inputs_dict = formula.fission_sources() # <- bulk of the work hidden here
    if not valid_fissile_inputs_dict:
        raise FormulaValidityError(f"Given formula unobtainable via fission:\n{formula}")

    new_element_mass, new_element_count = random.choice(list(valid_fissile_inputs_dict.items()))
    new_element = elements_data.elements_dict[new_element_mass]

    # Calculate and remove the element's dependencies, then add the new element
    formula.remove_fissile_element(new_element, new_element_count)
    formula[new_element] = new_element_count

    return new_element

def rand_inverse_fusion(formula, target_element=None):
    '''Mutate the given formula using any series of valid inverse fusion operations on a single
    element. The target element can optionally be specified.
    '''
    # Pick an element(s) and split it up randomly.
    # Bias toward splitting into elements already in the formula, and toward perfect divisors
    # TODO: For now will split all copies of an element that were given, but could change that
    #       such that some of the atoms to be fused up to also occur in the input.
    if target_element is None:
        fusable_elements = [e for e in formula.elements()
                            if 2 <= e.atomic_num <= 109]
        if not fusable_elements:
            raise FormulaValidityError(f"Given formula unobtainable via fusion:\n{formula}")
        target_element = random.choice(fusable_elements)

    # Randomly split the element up, biasing toward a split that produces an element already
    # in the formula, or one that divides evenly into the output element and isn't noble
    # (fuck Helium).
    # It should be impossible for nice_split_elements to end up empty since we always have
    # Hydrogen as a fallback divisor
    nice_split_elements = []
    ugly_split_elements = []
    for i in range(1, target_element.atomic_num):
        element = elements_data.elements_dict[i]
        if (element in formula or (target_element.atomic_num % i == 0
                                   and element.max_bonds != 0)):
            nice_split_elements.append(element)
        else:
            ugly_split_elements.append(element)
    # Go for clean splits 80% of the time if possible
    if not ugly_split_elements or random.random() < 0.8:
        divisor_element = random.choice(nice_split_elements)
    else:
        divisor_element = random.choice(ugly_split_elements)
    dividend = target_element.atomic_num // divisor_element.atomic_num
    formula[divisor_element] += dividend * formula[target_element]

    # Add leftover element if we didn't pick a perfectly splitting element
    remainder = target_element.atomic_num % divisor_element.atomic_num
    if remainder != 0:
        remainder_element = elements_data.elements_dict[remainder]
        formula[remainder_element] += formula[target_element]

    # Remove the transformed element from the input side
    del formula[target_element]

    return target_element


def randInputsFromOutputs(molecules, settings):
    '''Given a list of molecules, generate a list of minimally-sized molecules that could have
    generated them based on the available level features. We'll refer to the given molecules as
    outputs and the generated molecules as inputs, though how they are actually used is irrelevant.
    '''
    if not molecules:
        raise ValueError("Given molecules must not be empty.")

    outputs_formula = Formula()

    # Sum up the given output molecules' formulas
    for molecule in molecules:
        if not molecule:
            raise ValueError("Received empty output molecule.")
        # For the purposes of calculating a balanced reaction, we can effectively divide a
        # molecule by its GCD
        # e.g. C2H4 is effectively equivalent to CH2 in terms of creating a wasteless reaction
        outputs_formula += molecule.formula.least_common_formula()

    inputs_formula = copy.deepcopy(outputs_formula)

    try:
        if settings.fission and settings.fusion:
            tries, max_tries = 0, 10
            while set(inputs_formula.elements()) == set(outputs_formula.elements()):
                if tries >= max_tries:
                    raise MoleculeValidityError("Could not find non-trivial nuclear operations with given molecule(s).")
                # Fission then fusion
                # TODO: More complex order of nuclear operations.
                #       Fusion then fission is more involved as it can easily end up not requiring
                #       fission anyway.
                nuclear_element = rand_inverse_fission(inputs_formula)
                if nuclear_element.atomic_num <= 109:
                    _ = rand_inverse_fusion(inputs_formula, target_element=nuclear_element)
                tries += 1
        elif settings.fission:
            _ = rand_inverse_fission(inputs_formula)
        elif settings.fusion:
            _ = rand_inverse_fusion(inputs_formula)
    except FormulaValidityError as e:
        raise MoleculeValidityError("Could not complete nuclear operations with given molecule(s).") from e

    try:
        # Determine a minimally-sized set of input formulas which satisfy this reaction formula
        input_formulas = min_split_formula(inputs_formula,
                                           max_formulas=settings.num_inputs)
    except FormulaValidityError as e:
        raise MoleculeValidityError("Failed to split formula across inputs.") from e

    # If we didn't reach the target # of molecules (if specified), we need to spread the minimal
    # formulas across the required inputs
    if settings.num_inputs is not None and len(input_formulas) != settings.num_inputs:
        # Identify the largest input formula
        largest_formula = max(input_formulas, key=lambda f: f.num_atoms())

        # If there aren't enough atoms in the min split to cover all inputs, duplicate whichever
        # contains the element with the highest count in the balanced formula
        if largest_formula.num_atoms() < settings.num_inputs:
            bottleneck_element = inputs_formula.most_common(1)[0][0]
            for formula in input_formulas:
                if bottleneck_element in formula:
                    bottleneck_formula = formula
                    break
            for _ in range(settings.num_inputs - len(input_formulas)):
                input_formulas.append(copy.deepcopy(bottleneck_formula))
        else:
            # If only N (N < max_formulas) formulas were used, there were at most N unique element
            # counts in the balanced formula. Randomly split up the largest formula to fill our quota.
            # It's possible for splitting an output up to invalidate its molecule so allow re-attempts.
            split_attempts = 0
            new_input_formulas = []
            elements_collection = largest_formula.elements_collection()
            while not (new_input_formulas
                       and all(formula.is_valid() for formula in new_input_formulas)):
                if split_attempts >= 1000:
                    raise Exception(f"Could not find valid split of atoms {elements_collection}"
                                    + f" across {settings.num_inputs} molecules in 1000 attempts")

                remaining_elements = copy.deepcopy(elements_collection)
                random.shuffle(remaining_elements)
                # Remove chunks from the start of the atom list iteratively
                new_input_formulas = []
                for remaining_formulas in range(settings.num_inputs - len(input_formulas), 0, -1):
                    # Leave at least 1 atom per remaining formula
                    # Pivot doesn't include itself in the draw so start from idx at least 1.
                    pivot = random.randint(1, len(remaining_elements) - remaining_formulas)
                    drawn_elements, remaining_elements = remaining_elements[:pivot], remaining_elements[pivot:]
                    new_input_formulas.append(Formula(drawn_elements))
                new_input_formulas.append(Formula(remaining_elements))

                split_attempts += 1
            # Delete the original formula and add the newly split formulas
            input_formulas.pop(input_formulas.index(largest_formula))
            input_formulas.extend(new_input_formulas)

    return [randMoleculeFromFormula(formula) for formula in input_formulas]

def randResearchLevel(settings, verbose=False):
    '''Generate a random Research level and return the mission code.
    If specified, also pretty-print the created molecules.'''
    level = ResearchLevel()
    level['difficulty'] = settings.difficulty

    # Determine how many output zones we're using
    level['has-large-output'] = settings.large_output
    if settings.large_output:
        settings.num_outputs = 1
        # As a quick hack, internally set difficulty of large output levels from 1-4 instead of 0-3.
        settings.difficulty += 1
    elif settings.num_outputs is None:
        settings.num_outputs = random.randint(1, 2)

    # Generate molecules. Our molecule generation isn't sophisticated enough to avoid occasionally
    # making molecules that would violate the given settings (e.g. generating output molecules with
    # only greek elements in a nuclear level), so allow some loops until success
    tries, max_tries = 0, 50
    while True:
        tries += 1
        try:
            # Generate output molecules first
            output_molecules = randMolecules(num_molecules=1, settings=settings)
            # If we ended up generating 2 output zones, we'll decrease the difficulty of the second
            # by 1, to make it slightly more reasonable
            if settings.num_outputs == 2:
                lower_settings = copy.copy(settings)
                lower_settings.difficulty = max(settings.difficulty - 1, 0)
                output_molecules.extend(randMolecules(num_molecules=1, settings=lower_settings))

            # Generate input molecules based on the output molecules
            input_molecules = randInputsFromOutputs(molecules=output_molecules,
                                                    settings=settings)
            break
        except MoleculeValidityError as e:
            if tries >= max_tries:
                raise Exception("Could not generate valid molecules with the given settings within"
                                + f" {max_tries} attempts") from e

    # Randomly swap the inputs with the outputs to reduce bias from the generating algorithm,
    # if this doesn't violate our given settings
    if (not settings.large_output
            and not settings.fusion
            and not settings.fission
            and ((settings.num_inputs is None and settings.num_outputs is None)
                 or len(input_molecules) == len(output_molecules))
            and random.random() < 0.5):
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
        level['input-zones'][zone_idx] = input_zone_json

    # Add output zones to the level
    for zone_idx, molecule in enumerate(output_molecules):
        if len(output_molecules) == 1 and not settings.large_output:
            zone_idx = random.randint(0, 1)

        output_zone_json = {}
        output_zone_json['molecule'] = molecule.get_json_str()
        output_zone_json['count'] = 10

        level['output-zones'][zone_idx] = output_zone_json

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
        print('Inputs:')
        for molecule in input_molecules:
            print(molecule)
        print('Outputs:')
        for molecule in output_molecules:
            print(molecule)

    # Return the mission code
    return level.get_code()

def randProductionLevel(settings, verbose=False):
    '''Generate a random Production level.'''
    level = ProductionLevel()
    level['difficulty'] = settings.difficulty

    # Determine how many output zones we're using
    if settings.num_outputs is None:
        if settings.difficulty == 0:
            settings.num_outputs = random.randint(1, 2)
        else:
            settings.num_outputs = random.choice([1, 2, 2, 3]) # bias toward 2 outputs

    # Generate molecules. Our molecule generation isn't sophisticated enough to avoid occasionally
    # making molecules that would violate the given settings (e.g. generating output molecules with
    # only greek elements in a nuclear level), so allow some loops until success
    tries, max_tries = 0, 50
    while True:
        tries += 1
        try:
            # Generate output molecules first
            output_molecules = randMolecules(num_molecules=settings.num_outputs,
                                             settings=settings)
            # Generate input molecules based on the output molecules
            input_molecules = randInputsFromOutputs(molecules=output_molecules,
                                                    settings=settings)
            break
        except MoleculeValidityError as e:
            if tries >= max_tries:
                raise Exception("Could not generate valid molecules with the given settings within"
                                + f" {max_tries} attempts") from e

    # Randomly swap the inputs with the outputs to reduce bias from the generating algorithm
    # Must still follow # input / output constraints and fission / fusion results
    if (not settings.fusion
            and not settings.fission
            and ((settings.num_inputs is None and settings.num_outputs is None)
                 or len(input_molecules) == len(output_molecules))
            and random.random() < 0.5):
        input_molecules, output_molecules = output_molecules, input_molecules

    # Add the random input zone to the level
    input_zone_json = {'inputs':[]}
    input_json = {}
    input_json['molecule'] = input_molecules[0].get_json_str()
    input_json['count'] = 12
    input_zone_json['inputs'].append(input_json)
    level['random-input-zones'][0] = input_zone_json

    # Add the fixed input zones to the level
    for zone_idx, molecule in enumerate(input_molecules[1:]):
        level['fixed-input-zones'][zone_idx] = molecule.get_json_str()

    # Add the output zones
    for zone_idx, molecule in enumerate(output_molecules):
        output_zone_json = {}
        output_zone_json['molecule'] = molecule.get_json_str()
        output_zone_json['count'] = 40

        level['output-zones'][zone_idx] = output_zone_json

    # Add features to level
    # I'm leaving 'max-reactors' at its default 6 as set in ProductionLevel's initialization
    level['terrain'] = random.randint(0, 4)
    # By convention, if the inputs contain greek elements the terrain is set to Flidais
    for molecule in input_molecules:
        if any(e.symbol in ('Θ', 'Ω', 'Σ', 'Δ') for e in molecule.formula.elements()):
            level['terrain'] = 5
            break

    # No point having assembly/disassembly if any other reactors are available, and there's a max of
    # 4 types of reactors available anyway, probably for this reason
    if (not (settings.fusion or settings.fission)
            and random.random() < 0.02):
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
        print('Generated research level.')
        print('Inputs:')
        for molecule in input_molecules:
            print(molecule)
        print('Outputs:')
        for molecule in output_molecules:
            print(molecule)

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
    parser.add_argument('--large_output', action="store_true", default=None,
                        help="Set a research level to have a large output zone." \
                             + " Default has a small chance.")

    parser.add_argument('--fusion', action="store_true", default=None,
                        help="Use fusion in the level. Default has a small chance.")
    parser.add_argument('--fission', action="store_true", default=None,
                        help="Allow fission in the level. Default has a small chance.")
    parser.add_argument('--nuclear', action="store_true", default=None,
                        help="Allow both fission/fusion in the level. Default has a small chance.")

    molecule_topology_args = parser.add_mutually_exclusive_group()
    molecule_topology_args.add_argument('--polymer', action='store_true',
            help="Include a polymer molecule.")
    molecule_topology_args.add_argument('--symmetric', '--symmetrical', action='store_true',
            help="Include a symmetric molecule.")
    args = parser.parse_args()

    difficulty = args.difficulty
    if args.lanky:
        difficulty = 3

    if args.elements:
        elements = [elements_data.elements_dict[int(e)] if e.isdigit()
                    else elements_data.elements_dict[e]
                    for e in args.elements]
    elif args.basic:
        elements = elements_data.basic_elements
    else:
        elements = elements_data.elements

    if args.nuclear:
        args.fusion = True
        args.fission = True

    settings = RandomizationSettings(difficulty=difficulty,
                                     elements=elements,
                                     num_inputs=args.inputs,
                                     num_outputs=args.outputs,
                                     large_output=args.large_output,
                                     fusion=args.fusion,
                                     fission=args.fission,
                                     polymer=args.polymer,
                                     symmetric=args.symmetric)

    if args.production:
        if args.large_output:
            raise argparse.ArgumentTypeError('Cannot use a large output in a production level')

        code = randProductionLevel(settings=settings,
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
        print('Level code:\n')
    print(code)
