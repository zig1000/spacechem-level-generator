#!/usr/bin/env python

'''Maybe allow as args which tools the reactor contains (fusion laser, fission
laser, etc) and have the generated molecules be based on those.

TODO:
    - add fusion/fission
    - add support for production levels
    - add support for 'difficulty' (0 to 3)
    - do some more thinking on what makes for nice bond randomization rules.
      How homogenous should the bond levels be? How frequent should triple bonds be? etc.
    - add support for generating fully-bonded molecules
    - add support for waste-generating levels (need to calculate waste and limit it to ~40 at most)
    - add some sense of 'symmetry' to molecule generation preferences
    - Support greek/australium
'''

import argparse
import base64
import copy
import fractions
import gzip
import json
import random
import StringIO

import elements_data
from helpers import Formula, GridPos, Atom, Molecule

def randAtom(bond_count=None):
    '''Get a random atom. Optionally specify the desired bond count.
    '''
    if bond_count is None:
        return Atom(random.choice(elements_data.all_elements))
    else:
        return Atom(random.choice(elements_data.elements_with_bond_count[bond_count]))

def randBasicAtom():
    return Atom(random.choice(elements_data.basic_elements))

def randMolecule(num_atoms=None,
                 elements=elements_data.all_elements,
                 difficulty=random.randint(0, 3)):
    '''Generate a random molecule by randomy generating and adding one atom at a time. Bias the
    elements new atoms are drawn from to increase the chances of drawing elements already present
    in the molecule.
    Args:
        num_atoms: The target number of atoms, 1 to 16. Note that the current algorithm can fall short of the
                   target # of atoms if it runs out of available bonds.
        elements: The pool of elements that atoms may be drawn from (with replacement)
        difficulty: From 0 to 3 with 3 being the hardest, as displayed in ResearchNet. If num_atoms
                    is not specified, the difficulty affects the randomly selected # of atoms.
    '''
    molecule = Molecule()

    # Generate size exponentially, with the probability dropping sharply past ~6
    if num_atoms is None:
        num_atoms = 0
        if difficulty == 0:
            num_atoms = random.randint(1, 3)
        elif difficulty == 1:
            for i in range(16):
                num_atoms += 1
                if num_atoms < 5 and random.randint(1, 6) == 1 \
                   or num_atoms >= 5 and random.randint(1, 3) == 1:
                    break
        elif difficulty == 2:
            for i in range(16):
                num_atoms += 1
                if num_atoms < 6 and random.randint(1, 8) == 1 \
                   or num_atoms >= 6 and random.randint(1, 3) == 1:
                    break
        else:
            num_atoms = random.randint(8, 16)

    # Track which elements have already appeared in the molecule - adding new
    # unique elements within a single atom should become progressively less
    # likely
    used_elements = []
    unused_elements = elements

    for i in range(num_atoms):
        # Figure out where we're inserting the new atom
        open_positions = molecule.open_positions()
        if not open_positions:
            return molecule

        # Determine the element of the new atom
        prob_new_element = 0.5**len(used_elements)

        # Only the first element is allowed to be noble
        if i == 1:
            used_elements = [e for e in used_elements if elements_data.max_bonds[e] != 0]
            unused_elements = [e for e in unused_elements if elements_data.max_bonds[e] != 0]

        if random.uniform(0,1) < prob_new_element and len(unused_elements) != 0:
            # Use a new element -
            # Some finagling to make this an O(1) random pop
            i = random.randint(0, len(unused_elements) - 1)
            unused_elements[i], unused_elements[-1] = unused_elements[-1], unused_elements[i]
            elem = unused_elements.pop(-1)
            atom = Atom(elem)
            used_elements.append(elem)
        else:
            # Used element
            atom = Atom(random.choice(used_elements))

        atom.set_pos(random.choice(open_positions))

        # Assign bonds to the neighbors who still have available bonds, to a
        # minimum of 1 bond
        total_bonds = 0
        neighbor_posns = atom.pos.neighbors()
        random.shuffle(neighbor_posns)
        # Identify all the neighbors we can bond to
        neighbors = [molecule[neighbor] for neighbor in neighbor_posns \
                     if molecule[neighbor] is not None
                     and molecule[neighbor].remaining_bonds() > 0]
        for i, neighbor in enumerate(neighbors):
            tmp_max_bonds = min(neighbor.remaining_bonds(),
                                atom.remaining_bonds(),
                                3)
            # If we have yet to add bonds and this is the last bondable neighbor,
            # set the min bond count to 1. Else, 0.
            if total_bonds == 0 and i == len(neighbors) - 1:
                bond_count = random.randint(1, tmp_max_bonds)
            else:
                bond_count = random.randint(0, tmp_max_bonds)
            total_bonds += bond_count

            if neighbor.col < atom.col:
                atom.left_bonds = bond_count
            elif neighbor.col > atom.col:
                atom.right_bonds = bond_count
            elif neighbor.row < atom.row:
                atom.top_bonds = bond_count
            else:
                atom.bottom_bonds = bond_count

        molecule.add_atom(atom)

        # Use an exponential distribution to decide when to stop adding atoms
        # TODO: Rig this to specifically reduce single-atom rate
        if random.randint(1, 4) == 1:
            break

    return molecule

def randMoleculeFromFormula(formula):
    '''Generate a molecule containing exactly the atoms indicated (by atomic numbers).
    Here I call 'atoms' what is actually just a list of atomic symbols, whoops.
    It must be possible to construct such a molecule, e.g. [1, 1, 1] is an invalid atom list.
    C3H8 is similarly impossible to construct without hitting the reactor wall, though this is less
    obvious.
    '''
    if not formula.isValid():
        raise Exception("Cannot construct a molecule with formula {0}".format(formula))

    attempts = 0
    while True:
        if attempts >= 1000:
            raise Exception("No valid molecule found after 1000 retries for formula: {0}"
                            .format(list(formula.elements())))

        atoms = list(elements_data.atomic_numbers[s] for s in formula.elements())
        random.shuffle(atoms)

        molecule = Molecule()
        while atoms:
            if molecule.open_bonds == 0 and len(molecule) > 0:
                # If we've accidentally blocked ourselves off, retry on a new molecule
                break

            # If it's the last atom or placing this atom can't block placement of the next atom,
            # proceed freely. It's possible to have only one open position but lots of open bonds, or
            # only one open bond but lots of open positions.
            open_positions = molecule.open_positions()
            if (len(open_positions) > 1 and molecule.open_bonds != 1) \
               or len(atoms) == 1:
                i = 0
            else:
                # Given that we only have one available bond or position to add the new atom (and it
                # isn't the last atom), we have to make
                # sure we don't add an atom that only accepts one bond, unless it's the final atom to be added
                # Walk through the list until we find a 2+ bondable atom
                # TODO: inefficient af
                retry = False # For nested break
                for j, e in enumerate(atoms):
                    if elements_data.max_bonds[e] >= 2:
                        i = j
                        break
                    elif j == len(atoms) - 1:
                        # If we have nothing but bond-1 atoms left, we have to retry
                        retry = True
                        break
                if retry:
                    break

            atom = Atom(atoms[i])
            atoms[i], atoms[-1] = atoms[-1], atoms[i] # O[1] list pop via swapping to end
            atom = Atom(atoms.pop(-1))
            atom.set_pos(random.choice(open_positions))

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

        # If we failed to add any atoms, restore the original list and retry
        if atoms:
            atoms = list(elements_data.atomic_numbers[s] for s in formula.elements())
            continue

        # Otherwise, now that we've successfully added all the molecules, randomly add more bonds
        # For now, just walk through every atom and give a 50% chance to increase each of its bonds
        # (so max +2 bonds between any two atoms, and no need to worry about the 3-bond limit)
        added_atoms = molecule.atoms
        if added_atoms: # Handle the case when we were given no atoms
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
                if random.randint(0, 1):
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


def randBasicMolecule(num_atoms=None, difficulty=None):
    '''Generate a random non-noble molecule using only the 'basic' elements:
    [O, B, C, N, S, Cl, Os]
    '''
    return randMolecule(num_atoms=num_atoms,
                        elements=elements_data.basic_elements,
                        difficulty=difficulty)

def randFullyBondedMolecule():
    pass # TODO

def randFullyBondedBasicMolecule():
    pass # TODO

def randOutputMolecules(input_molecules,
                        difficulty,
                        bonders=True,
                        fusion=True,
                        fission=True):
    '''Given a list of input molecules, generate a valid list of output molecules.
    '''
    if type(input_molecules) == type(Molecule):
        input_molecules = [input_molecules]

    # Randomize the ratio of inputs in the balanced equation
    balanced_inputs_formula = Formula()
    for input_molecule in input_molecules:
        # For the purposes of calculating possible outputs, we can effectively reduce an input
        # molecule to its LCM (lowest common molecule :P)
        # e.g. C2H4 is effectively equivalent to CH2 in terms of creating a wasteless output
        lcm_formula = copy.copy(input_molecule.formula)
        lcm_formula.divide_by_gcd()

        # TODO: This is a kludgy workaround for now to make sure that we don't
        # have to be too smart about not exceeding the output grid size
        max_input_use_factor = 16 / lcm_formula.total_atoms()
        # Exponentially bias the input group created so its not usually too large
        input_use_factor = 1 # Leave it here for difficulty 0
        if difficulty == 1:
            for i in range(max_input_use_factor - 1):
                if random.randint(1, 2) == 1:
                    break
                input_use_factor += 1
        elif difficulty == 2:
            for i in range(max_input_use_factor - 1):
                if random.randint(1, 4) == 1:
                    break
                input_use_factor += 1
        else:
            input_use_factor = max_input_use_factor

        balanced_inputs_formula += input_use_factor * lcm_formula

    # Now we need to determine what ratio of outputs to match to this input ratio
    balanced_outputs_formula = balanced_inputs_formula

    # Determine which output zones we're using
    # First, we need to check if it's impossible to fit the balanced reaction formula into a single
    # input zone (keeping in mind that we can reduce the formula to its LCM to help)
    balanced_outputs_lcm_formula = copy.copy(balanced_outputs_formula)
    balanced_outputs_lcm_formula.divide_by_gcd()
    if not balanced_outputs_lcm_formula.isValid():
        output_zones = [0, 1]
    else:
        # 33% chance of double output
        output_zones = random.choice([[0], [1], [0, 1]])

    # Now attempt to allocate the outputs among the output zones and reduce each to a valid LCM
    output_lcm_formulas = []
    if len(output_zones) == 1:
        # In this case we've already checked that the full formula's LCM fits in one input zone
        output_lcm_formulas.append(balanced_outputs_lcm_formula)
    if len(output_zones) == 2:
        split_attempts = 0
        while not (len(output_lcm_formulas) == 2
                   and output_lcm_formulas[0].isValid()
                   and output_lcm_formulas[1].isValid()):
            if split_attempts >= 1000:
                raise Exception("Could not find valid split of output atoms {0} in 1000 attempts"
                                .format(balanced_outputs_formula))

            # Randomly split the formula between the two output zones
            # TODO: use numpy.random.choice(p=[1,2,3]) for performance
            num_atoms = balanced_outputs_formula.total_atoms()
            atoms = list(balanced_outputs_formula.elements())

            # If we only had one atom and want two output zones, just duplicate it
            # we also know in this case that the LCM is valid
            if num_atoms == 1:
                output_lcm_formulas = [Formula(atoms), Formula(atoms)]
                break
            else: # Randomly split the formula
                random.shuffle(atoms)
                pivot = random.randint(max(1, num_atoms - 16),
                                       min(16, num_atoms - 1))
                output_lcm_formulas = [Formula(atoms[:pivot]),
                                       Formula(atoms[pivot:])]
                # Take the LCM of the split lists
                output_lcm_formulas[0].divide_by_gcd()
                output_lcm_formulas[1].divide_by_gcd()
            split_attempts += 1

    # Now take the LCM formulas that we just generated for the output zones and multiply them by
    # some reasonable value, while keeping them valid
    result = []
    for output_lcm_formula in output_lcm_formulas:
        # Figure out the most we can multiply the formula by while still being constructable
        # Inefficient, but this isn't the part of the code that could run 1000 times
        max_output_use_factor = 16 / output_lcm_formula.total_atoms()
        for i in range(1, max_output_use_factor+1):
            if not (i*output_lcm_formula).isValid():
                max_output_use_factor = i-1
                break

        # Exponentially weight the multiplier toward lower values
        # Scale difficulty the same way as the input molecules, I guess
        output_use_factor = 1 # Leave this alone for difficulty 0
        if difficulty == 1:
            for i in range(max_output_use_factor - 1):
                if random.randint(1, 2) == 1:
                    break
                output_use_factor += 1
        elif difficulty == 2:
            for i in range(max_output_use_factor - 1):
                if random.randint(1, 4) == 1:
                    break
                output_use_factor += 1
        else:
            output_use_factor = max_output_use_factor

        output_formula = output_use_factor * output_lcm_formula
        result.append(randMoleculeFromFormula(formula=output_formula))
    return result


def randResearchLevel(difficulty,
                      fusion=True,
                      fission=True,
                      balanced=True,
                      rand_inputs=False,
                      large_output=False):
    '''Return a string of the mission code.
    '''
    level_json = {} # Still haven't figured out what hashing function is used
                    # for mission codes
    # 33% chance for each of alpha input, beta input, or both
    input_zones = random.choice([[0], [1], [0, 1]])
    input_molecules = []

    # Generate Input Zones
    input_zones_json = {}
    for zone_idx in input_zones:
        input_zone_json = {}
        num_molecules = 1 # TODO: random inputs
        molecules_list = []
        for _ in range(num_molecules):
            molecule_json = {}
            molecule = randBasicMolecule(difficulty=difficulty)
            molecule_json['molecule'] = molecule.get_json_str()
            molecule_json['count'] = 12
            molecules_list.append(molecule_json)
            input_molecules.append(molecule)
        input_zone_json['inputs'] = molecules_list
        input_zones_json[str(zone_idx)] = input_zone_json
    level_json['input-zones'] = input_zones_json

    # Generate Output Zones
    output_zones_json = {}
    output_molecules = randOutputMolecules(input_molecules, difficulty=difficulty)
    for i, molecule in enumerate(output_molecules):
        if len(output_molecules) == 1:
            zone_idx = random.randint(0, 1)
        else:
            zone_idx = i

        molecule_json = {}
        molecule_json['molecule'] = molecule.get_json_str()
        molecule_json['count'] = 10

        output_zones_json[str(zone_idx)] = molecule_json

    level_json['output-zones'] = output_zones_json

    # Set features of the level (bonders, fusers, etc.)
    level_json['has-large-output'] = False
    if difficulty == 0:
        level_json['bonder-count'] = 4
    else:
        level_json['bonder-count'] = 2*random.randint(1, 4)
    level_json['has-sensor'] = False
    level_json['has-fuser'] = False
    level_json['has-splitter'] = False
    level_json['has-teleporter'] = difficulty >= 2 and random.randint(0,1) == 1

    # Level meta-data
    level_json['type'] = 'research'
    level_json['name'] = 'RandomlyGenerated'
    level_json['author'] = 'Zig'
    level_json['difficulty'] = difficulty

    raw_json = json.dumps(level_json)

    # Convert to mission code - gzip then b64 the result
    out = StringIO.StringIO()
    with gzip.GzipFile(fileobj=out, mode="w") as f:
        f.write(raw_json)
    mission_code = base64.b64encode(out.getvalue())

    return mission_code

def randProductionLevel():
    pass # TODO

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--lanky', action="store_true",
                          help="Set difficulty to max")
    parser.add_argument('--difficulty', type=int, action="store", default=2,
                          help="Set difficulty to max")
    args = parser.parse_args()

    difficulty = args.difficulty
    if args.lanky:
        difficult = 3

    print randResearchLevel(difficulty=difficulty)
