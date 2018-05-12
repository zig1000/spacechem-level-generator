#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
import copy
import fractions
import random

import elements_data
from helpers import *

def randTotalInputsSize(difficulty, min_size=1):
    '''Generate a combined size for a level's inputs, in # of atoms.
    '''
    if difficulty == 0:
        return random.randint(min_size, 3)
    elif difficulty == 1:
        return random.randint(max(4, min_size), 6)
    elif difficulty == 2:
        return random.randint(max(7, min_size), 10)
    else:
        return random.randint(max(11, min_size), 16)


def randFormula(num_atoms=None,
                elements=None,
                element_reuse_factor=0.6,
                nobles=False):
    '''Generate a random molecule formula, optionally specifying a # of atoms, and drawn with replacement
    from the specified list of elements (default all elements).
    '''
    if num_atoms is None: # If size isn't specified keep adding atoms with 75% probability
        num_atoms = 1
        while num_atoms <= 16 and random.randint(1, 4) != 1:
            num_atoms += 1

    # Restrict noble gases
    if elements is None:
        if num_atoms == 1 and nobles:
            elements = elements_data.symbols.values()
        else:
            elements = elements_data.non_noble_symbols

    formula = Formula()
    tries = 0
    while not formula.isValid():
        if tries == 500:
            raise ValueError(("Could not generate valid formula with {0} atoms from {1} in 500 "
                             + "tries.").format(num_atoms, elements))

        formula = Formula()
        # Bias the algorithm to increase the chance of reusing elements within a molecule
        used_elements = []
        unused_elements = copy.copy(elements)
        prob_new_element = 1
        for _ in range(num_atoms):
            if random.uniform(0, 1) < prob_new_element and len(unused_elements) != 0:
                # Use a new element and move it to the used list
                i = random.randint(0, len(unused_elements) - 1)
                # Some finagling to make this an O(1) random pop
                unused_elements[i], unused_elements[-1] = unused_elements[-1], unused_elements[i]
                element = unused_elements.pop(-1)
                used_elements.append(element)
                # Exponentially decrease the probability of getting a new element
                prob_new_element *= element_reuse_factor
            else:
                # Used element
                element = random.choice(used_elements)
            formula[element] += 1
        tries += 1

    return formula

def randMolecule(num_atoms=None,
                 elements=None):
    '''Generate a random molecule.
    Args:
        num_atoms: The target number of atoms, 1 to 16. Note that the current algorithm can fall short of the
                   target # of atoms if it runs out of available bonds.
        elements: The pool of elements that atoms may be drawn from (with replacement).
    '''
    formula = randFormula(num_atoms=num_atoms, elements=elements)
    return randMoleculeFromFormula(formula)

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
                            .format(list(formula.atoms())))

        atoms = list(elements_data.atomic_numbers[s] for s in formula.atoms())
        random.shuffle(atoms)

        molecule = Molecule()
        while atoms:
            if molecule.open_bonds == 0 and len(molecule) > 0:
                # If we've accidentally blocked ourselves off, retry on a new molecule
                break

            # If it's the last atom or placing this atom can't block placement of the next atom,
            # proceed freely. It's possible to have only one open position but lots of open bonds,
            # or only one open bond but lots of open positions, so check carefully
            open_positions = molecule.open_positions()
            if not open_positions:
                # In theory this should never occur
                print 'Error in call to randMoleculeFromFormula with formula: {0}'.format(formula) \
                      + '\nConstructed partial molecule with open bonds but no open positions: ' \
                      + str(molecule)
            elif (len(open_positions) > 1 and molecule.open_bonds != 1) \
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

        # If any atoms failed to be added, restore the original list and retry
        if atoms:
            atoms = list(elements_data.atomic_numbers[s] for s in formula.atoms())
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

def randBasicMolecule(num_atoms=None):
    '''Generate a random non-noble molecule using only the 'basic' elements:
    [O, B, C, N, S, Cl, Os]
    '''
    return randMolecule(num_atoms=num_atoms,
                        elements=elements_data.basic_elements)

def randFullyBondedMolecule():
    pass # TODO

def randFullyBondedBasicMolecule():
    pass # TODO

def randOutputMolecules(input_molecules,
                        difficulty,
                        max_outputs=2, # Set to 3 for a production level
                        bonders=True,
                        fusion=True,
                        fission=True):
    '''Given a list of input molecules, generate a valid list of output molecules.
    '''
    if type(input_molecules) == type(Molecule): # Allow passing in a lone molecule
        input_molecules = [input_molecules]

    # Randomize the relative ratio of inputs in the balanced equation
    balanced_inputs_formula = Formula()
    for input_molecule in input_molecules:
        # For the purposes of calculating possible outputs, we can effectively reduce an input
        # molecule to its LCM (lowest common molecule :P)
        # e.g. C2H4 is effectively equivalent to CH2 in terms of creating a wasteless output
        lcm_formula = copy.copy(input_molecule.formula)
        lcm_formula.divide_by_gcd()
        # TODO: I was multiplying the GCD molecule by some factor before, but the outputs were
        #       getting way too complex so for now I'll just pass the gcd over, and let the
        #       randomness of the input molecules' GCDs decide their relative ratios
        balanced_inputs_formula += lcm_formula

        # TODO: This is a kludgy workaround for now to make sure that we don't
        # have to be too smart about not exceeding the output grid size
        #max_input_use_factor = 16 / lcm_formula.num_atoms()
        # Exponentially bias the input group created so its not usually too large
        #input_use_factor = 1 # Leave it here for difficulty 0
        #if difficulty == 1:
        #    for i in range(max_input_use_factor - 1):
        #        if random.randint(1, 2) == 1:
        #            break
        #        input_use_factor += 1
        #elif difficulty == 2:
        #    for i in range(max_input_use_factor - 1):
        #        if random.randint(1, 4) == 1:
        #            break
        #        input_use_factor += 1
        #else:
        #    input_use_factor = max_input_use_factor

        #balanced_inputs_formula += input_use_factor * lcm_formula

    # Now we need to determine what ratio of outputs to match to this input ratio
    balanced_outputs_formula = balanced_inputs_formula

    # Determine which output zones we're using
    # First, we need to check if it's impossible to fit the balanced reaction formula into a single
    # input zone (keeping in mind that we can reduce the formula to its LCM to help)
    balanced_outputs_lcm_formula = copy.copy(balanced_outputs_formula)
    balanced_outputs_lcm_formula.divide_by_gcd()
    if not balanced_outputs_lcm_formula.isValid():
        if balanced_outputs_lcm_formula.num_atoms() <= max_outputs:
            # In the rare case where we have noble gas(es) among 2 or 3 inputs, just restrict the
            # # of outputs to match them.
            # TODO: This biases cases where 3 atoms fit in 2 outputs (but not 1)
            num_outputs = balanced_outputs_lcm_formula.num_atoms()
        else:
            # We'll assume the balanced inputs fit within 2 output zones
            num_outputs = random.randint(2, max_outputs)
    else:
        # Equal chance for any # of outputs (research = 1 or 2, production = 1, 2, or 3)
        # TODO: 50% chance of 2 outputs in productions?
        num_outputs = random.randint(1, max_outputs)

    # Now attempt to allocate the outputs among the output zones and reduce each to a valid LCM
    # In case we have less atoms to allocate than output zones, just duplicate
    # the balanced inputs across all outputs.
    # TODO: fine for 1:2 or 1:3, but heavy-handed for 2:3
    if num_outputs == 1 or balanced_outputs_formula.num_atoms() < num_outputs:
        output_lcm_formulas = [copy.copy(balanced_outputs_formula)
                               for i in range(num_outputs)]
    else:
        # Otherwise, randomly split the balanced formula between output zones
        split_attempts = 0
        output_lcm_formulas = []
        atoms = list(balanced_outputs_formula.atoms())
        num_atoms = len(atoms)
        while not (output_lcm_formulas
                   and all(formula.isValid() for formula in output_lcm_formulas)):
            if split_attempts >= 1000:
                raise Exception("Could not find valid split of output atoms {0} in 1000 attempts"
                                .format(balanced_outputs_formula))

            # Remove chunks from the start of the atom list iteratively
            # TODO: Current implementation is biased toward initial inputs in the way it splits
            output_lcm_formulas = []
            random.shuffle(atoms)
            remaining_atoms = copy.copy(atoms)
            num_remaining_atoms = num_atoms
            for i in range(num_outputs - 1):
                remaining_outputs = num_outputs - i - 1
                # Leave at least 1 and at most 16 atoms per remaining output zone
                # Pivot doesn't include itself in the draw so start from idx at least 1.
                pivot = random.randint(max(1, num_remaining_atoms - 16*remaining_outputs),
                                       min(16, num_remaining_atoms - remaining_outputs))
                drawn_atoms, remaining_atoms = remaining_atoms[:pivot], remaining_atoms[pivot:]
                output_lcm_formulas.append(Formula(drawn_atoms))
                num_remaining_atoms -= pivot
            output_lcm_formulas.append(Formula(remaining_atoms))
            # Take the LCM of the split lists before checking validity
            for formula in output_lcm_formulas:
                formula.divide_by_gcd()
            split_attempts += 1

    # Now take the LCM formulas that we just generated for the output zones and multiply them by
    # some reasonable value, while keeping them valid
    result = []
    for output_lcm_formula in output_lcm_formulas:
        # Multiply the LCM molecule by some amount, with the total molecule size capped by the
        # difficulty (unless that cap has already been exceeded somehow, in which case do nothing).
        if difficulty == 0:
            max_output_size = 3
        elif difficulty == 1:
            max_output_size = 6
        elif difficulty == 2:
            max_output_size = 10
        else:
            max_output_size = 16
        max_output_use_factor = max(1, max_output_size / output_lcm_formula.num_atoms())

        # Figure out the most we can multiply the formula by while still being constructable
        # Inefficient, but this isn't the part of the code that could run 1000 times
        for i in range(1, max_output_use_factor):
            if (i*output_lcm_formula).isValid():
                max_output_use_factor = i
            else:
                break

        output_formula = random.randint(1, max_output_use_factor) * output_lcm_formula
        result.append(randMoleculeFromFormula(formula=output_formula))
    return result

def randMolecules(num_molecules, total_size=None, elements=None):
    '''Generate a specified # of random molecules. Optionally specify the elements to use and the
    total size in # of atoms of all molecules combined.
    '''
    molecules = []
    for i in range(num_molecules):
        if total_size is not None:
            remaining_molecules = num_molecules - i - 1
            if i < num_molecules - 1:
                # Leave at least 1 and at most 16 atoms per remaining input
                size = random.randint(max(1, total_size - 16*remaining_molecules),
                                      min(16, total_size - remaining_molecules))
            else:
                size = total_size
            molecules.append(randMolecule(num_atoms=size, elements=elements))
            total_size -= size
        else:
            molecules.append(randMolecule(elements=elements))
    return molecules

def randResearchLevel(difficulty=random.randint(0, 2),
                      elements=None,
                      num_inputs=None,
                      fusion=True,
                      fission=True,
                      balanced=True,
                      rand_inputs=False,
                      large_output=False):
    '''Generate a random Research level and return the mission code.
    '''
    level = ResearchLevel()
    level['difficulty'] = difficulty

    # Determine how many input zones we're using
    if num_inputs is None:
        num_inputs = random.randint(1, 2)

    # Generate input molecules
    total_inputs_size = randTotalInputsSize(difficulty=difficulty, min_size=num_inputs)
    input_molecules = randMolecules(num_inputs, total_size=total_inputs_size, elements=elements)

    # Generate output molecules
    output_molecules = randOutputMolecules(input_molecules, difficulty=difficulty)

    # Randomly swap the inputs with the outputs to reduce bias from the inputs -> outputs algorithm
    # TODO: This will need fixing when --outputs gets added
    if ((num_inputs is None or len(input_molecules) == len(output_molecules))
        and random.randint(0, 1) == 1):
        input_molecules, output_molecules = output_molecules, input_molecules

    # Add input zones to the level
    for zone_idx, molecule in enumerate(input_molecules):
        if len(input_molecules) == 1:
            zone_idx = random.randint(0, 1)

        input_zone_json = {'inputs':[]}
        num_inputs = 1 # TODO: random inputs
        for _ in range(num_inputs):
            input_json = {}
            input_json['molecule'] = molecule.get_json_str()
            input_json['count'] = 12
            input_zone_json['inputs'].append(input_json)
        level['input-zones'][str(zone_idx)] = input_zone_json

    # Add output zones to the level
    for zone_idx, molecule in enumerate(output_molecules):
        if len(output_molecules) == 1:
            zone_idx = random.randint(0, 1)

        output_zone_json = {}
        output_zone_json['molecule'] = molecule.get_json_str()
        output_zone_json['count'] = 10

        level['output-zones'][str(zone_idx)] = output_zone_json

    # Add 'features' to the level (bonders, fusers, etc.)
    level['has-large-output'] = False
    if difficulty == 0:
        level['bonder-count'] = 4
    else:
        level['bonder-count'] = random.choice([2, 4, 4, 8]) # Bias toward 4 bonders
    level['has-sensor'] = random.randint(0, 1) == 1 # Sure, why not
    level['has-fuser'] = False
    level['has-splitter'] = False
    level['has-teleporter'] = difficulty >= 2 and random.randint(0, 1) == 1

    # Return the mission code
    return level.getCode()

def randProductionLevel(difficulty=random.randint(0, 2),
                        elements=None,
                        num_inputs=None):
    '''Create a random Production level.
    Args:
        difficulty: Max 3, defaults to 0-2. Controls total size of molecules.
        elements: The elements to select from for the level.
        num_inputs: The number of input zones to use (1-3).
    '''
    level = ProductionLevel()

    level['difficulty'] = difficulty

    # Determine how many input zones we're using
    if num_inputs is None:
        num_inputs = random.choice([1, 2, 2, 3]) # bias toward 2 inputs

    # Generate input molecules
    total_inputs_size = randTotalInputsSize(difficulty=difficulty, min_size=num_inputs)
    input_molecules = randMolecules(num_inputs, total_size=total_inputs_size, elements=elements)

    # Generate output molecules (up to 3 since this is a production level)
    output_molecules = randOutputMolecules(input_molecules, difficulty=difficulty, max_outputs=3)

    # Randomly swap the inputs with the outputs to reduce bias from the inputs -> outputs algorithm
    # TODO: This will need fixing when --outputs gets added
    if ((num_inputs is None or len(input_molecules) == len(output_molecules))
        and random.randint(0, 1) == 1):
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
        for element in molecule.formula.keys():
            if element in ('Θ', 'Ω', 'Σ', 'Δ'):
                level['terrain'] = 5
                break

    # No point having assembly/disassembly if any other reactors are available, and there's a max of
    # 4 types of reactors available anyway, probably for this reason
    if random.randint(1, 50) == 1:
        level['has-assembly'] = True
        level['has-disassembly'] = True
    else:
        level['has-starter'] = True # In case anyone wants to comment their code by using a
                                        # strictly worse reactor
        level['has-advanced'] = random.randint(0, 1) == 1 # Sensor reactor; 50%
        level['has-nuclear'] = False # TODO
        level['has-superbonder'] = random.randint(1, 4) == 1 # 25%

    level['has-recycler'] = random.randint(0, 1) == 1

    return level.getCode()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--view', action='store_true',
                        help="Open Cearn's mission viewer in the browser")
    parser.add_argument('--research', action='store_true',
                        help="Dummy flag; Research levels are generated by default anyway.")
    parser.add_argument('--production', action='store_true',
                        help="Generate a Production level instead of a Research level.")
    parser.add_argument('--difficulty', '-d', type=int, action="store",
                        choices=[0, 1, 2, 3], default=random.randint(0, 2),
                        help="Set the difficulty of the level.")
    parser.add_argument('--lanky', action="store_true",
                        help="Set difficulty to max")
    parser.add_argument('--basic', action='store_true',
                        help="Use only 'basic' elements, one for each of the non-0 bond counts:" \
                             + "[H, O, B, C, N, S, Cl, Os]")
    parser.add_argument('--elements', action='store', nargs='*', type=str,
                        help="Select only from the specified elements in creating the level." \
                             + "Not currently guaranteed to use all elements.")
    parser.add_argument('--inputs', action='store', type=int,
                        choices=[1, 2, 3],
                        help="The # of inputs to use.")
    args = parser.parse_args()

    difficulty = args.difficulty
    if args.lanky:
        difficulty = 3
    elements = None
    if args.elements:
        # Allow specifying elements by atomic #, because why not
        elements = [elements_data.symbols[int(e)] if e.isdigit() else e for e in args.elements]
    elif args.basic:
        elements = elements_data.basic_elements

    if args.inputs == 3 and not args.production:
        raise argparse.ArgumentTypeError('Too many inputs; did you mean to add --production ?')

    # Generate the level
    if args.production:
        code = randProductionLevel(difficulty=difficulty,
                                   elements=elements,
                                   num_inputs=args.inputs)
    else:
        code = randResearchLevel(difficulty=difficulty,
                                 elements=elements,
                                 num_inputs=args.inputs)

    if args.view:
        import webbrowser
        webbrowser.open('http://coranac.com/spacechem/mission-viewer?mode=editor&code=' + code)
    print code
