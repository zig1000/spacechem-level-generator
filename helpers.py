#!/usr/bin/env python
import collections
import fractions

import elements_data

# Confusingly, Counter's list-ifier is already called 'elements()'
class Formula(collections.Counter):
    # Have to override Counter's add method or else adding two Formula's will make a Counter
    def __add__(self, other):
        result = Formula()
        for k in self.keys():
            result[k] += self[k]
        for k in other.keys():
            result[k] += other[k]
        return result

    def __mul__(self, other):
        return Formula({i: other * self[i] for i in self.keys()})
    __rmul__ = __mul__

    def total_atoms(self):
        return sum(self.values())

    def divide_by_gcd(self):
        '''Divide this formula's numeric values by their GCD.
        '''
        gcd = reduce(fractions.gcd, self.values())
        for i in self.keys():
            self[i] = self[i] / gcd

    def get_json_str(self):
        # Use Hill System for chemical formulas; C, then H, then alphabetized
        result = ''
        elements = ['C', 'H'] + sorted(k for k in self.keys() if k not in ['C', 'H'])
        for element in elements:
            if self[element] != 0:
                result += element
                if self[element] != 1:
                    result += '~' + str(self[element]).rjust(2, '0')
        return result

    def isValid(self):
        '''Check if it's possible to form a molecule with this formula within a 4x4 grid.
        '''
        # Verify size constraints for a basic input/output
        if not (1 <= self.total_atoms() <= 16):
            return False

        # We'll calculate validity based on whether there are enough bonds in the atom list to form a
        # minimally-connected molecule (no cycles in the graph)
        # In order not to be fooled by extraneous max bond counts that can't be used to connect the
        # graph, consider a minimally-connected 16-atom graph, constructed by adding middle atoms first
        # and always using an element with exactly the bond count needed for its new connections.
        #  H   H   H   H
        #  |   |   |   |
        #  O - C - C - O
        #      |   |
        #  O - B   B - O
        #  |   |   |   |
        #  H   H   H   H
        # It is left as an exercise to the reader to prove that constructing a molecule in this way,
        # with higher bond count atoms placed first and preferring the middle, will fail iff it is
        # impossible to construct any molecule from the given atoms.
        # Thus, we can count the 'useful' bonds added by each successive atom, and if it fails to
        # equal or exceed the minimum # of bonds in a fully-connected graph of the required size,
        # then a molecule cannot be constructed.
        # It is also required that each atom be able to form at least 1 bond if the molecule is not
        # size 1.
        total_atoms = self.total_atoms()
        total_connections = 0
        max_connections_by_idx = [4, 4, 3, 3, 2, 2, 2, 2] + [1]*8
        i = 0
        for element in sorted(self.keys(), key=lambda x: elements_data.max_bonds[elements_data.atomic_numbers[x]], reverse=True):
            element_max_bonds = elements_data.max_bonds[elements_data.atomic_numbers[element]]
            # If there's a noble gas and more than 1 atom, it's invalid
            if total_atoms > 1 and element_max_bonds == 0:
                return False
            # Count 'useful' connections
            for i in range(self[element]):
                total_connections += min(element_max_bonds, max_connections_by_idx[i])
                i += 1
        # For N atoms there are min N-1 connections, each contributing min 1 bond to each of two atoms
        return total_connections >= 2*(total_atoms - 1)


class GridPos:
    def __init__(self, row, col):
        self.row = row
        self.col = col

    def __str__(self):
        return '({0}, {1})'.format(self.row, self.col)

    def __repr__(self):
        return str(self)

    def neighbors(self):
        result = []
        for row, col in [(self.row, self.col - 1), # left
                         (self.row, self.col + 1), # right
                         (self.row - 1, self.col), # top
                         (self.row + 1, self.col)]: # bottom
            if 0 <= row <= 3 and 0 <= col <= 3:
                result.append(GridPos(row, col))
        return result


class Atom:
    '''Represent an atom and its neighbouring bonds.
    '''
    def __init__(self, id):
        if type(id) == int:
            self.atomic_num = id
            self.symbol = elements_data.symbols[id]
        else:
            self.symbol = id
            self.atomic_num = elements_data.atomic_numbers[id]
        self.max_bonds = elements_data.max_bonds[self.atomic_num]

        self.row = None
        self.col = None
        self.pos = None

        self.left_bonds = 0
        self.right_bonds = 0
        self.top_bonds = 0
        self.bottom_bonds = 0

    def __str__(self):
        return self.symbol.rjust(2) # Pad element symbol to two chars

    def __repr__(self):
        return 'Atom({0}, {1}, {2}, {3}, {4}, {5})'.format(self.symbol, self.pos,
                                                                self.left_bonds, self.top_bonds,
                                                                self.right_bonds, self.bottom_bonds)

    def get_json_str(self):
        return '{0}{1}{2}{3}{4}'.format(self.col, self.row,
                                        str(self.atomic_num).rjust(3, '0'),
                                        self.right_bonds,
                                        self.bottom_bonds)

    def set_pos(self, pos):
        self.row = pos.row
        self.col = pos.col
        self.pos = pos

    def remaining_bonds(self):
        return self.max_bonds - self.left_bonds \
                              - self.right_bonds \
                              - self.top_bonds \
                              - self.bottom_bonds


class Molecule:
    '''For convenience a molecule will be postioned in a grid rather than
    abstracting the grid away.
    '''
    def __init__(self):
        self.atoms = []
        self.grid = [[None, None, None, None] for i in range(4)]
        self.formula = Formula() # a Counter object for the constituent elements
        self.open_bonds = 0 # A rough meausure of the # of open bonds available
                            # Only needs to be correct insofar as 0 vs 1 vs many

    def __str__(self):
        result = ''
        for r in range(4):
            for c in range(4):
                atom = self.grid[r][c]
                if atom is None:
                    result += 2*' '
                    if c < 3:
                        result += 3*' '
                else:
                    result += str(atom)
                    if c < 3:
                        if atom.right_bonds > 0:
                            result += ' ' + str(atom.right_bonds) + ' '
                        else:
                            result += 3*' '
            result += '\n'
            if r < 3:
                for c in range(4):
                    atom = self.grid[r][c]
                    if atom is None or atom.bottom_bonds == 0:
                        result += 5*' '
                    else:
                        result += ' ' + str(atom.bottom_bonds) + 3*' '
                result += '\n'
        return result

    def __getitem__(self, pos):
        return self.grid[pos.row][pos.col]

    def __len__(self):
        return len(self.atoms)

    def get_json_str(self):
        result = "Randite;"
        result += self.formula.get_json_str() + ';'
        for atom in self.atoms:
            result += atom.get_json_str() + ';'
        return result

    def open_positions(self):
        # For an empty molecule, all positions are open
        if len(self) == 0:
            return [GridPos(r, c) for r in range(4) for c in range(4)]
        # If there are no remaining bonds, we can skip the overhead of walking through the atoms
        elif self.open_bonds == 0:
            return []

        result_dict = {} # For O(1) checks on whether a position has already been added
        for atom in self.atoms:
            for pos in atom.pos.neighbors():
                if atom.remaining_bonds() > 0  \
                   and self[pos] is None  \
                   and pos not in result_dict:
                    result_dict[pos] = None # None takes up less memory
        return result_dict.keys() # We need a list back though

    def add_atom(self, atom):
        self.grid[atom.row][atom.col] = atom

        # Quick helper to check if an atom within this molecule's grid has at least 1 open side
        def has_open_side(atom):
            for neighbor_posn in atom.pos.neighbors():
                if self[neighbor_posn] is None:
                    return True
            return False

        # Add any open bonds donated by the new atom to the molecule
        if has_open_side(atom):
            self.open_bonds += atom.remaining_bonds()

        for pos in atom.pos.neighbors():
            adj_atom = self[pos]

            # If this neighbor is an atom, update its bonds to match this atom's
            if adj_atom is not None:
                if adj_atom.col < atom.col:
                    adj_atom.right_bonds = added_bonds = atom.left_bonds
                elif adj_atom.col > atom.col:
                    adj_atom.left_bonds = added_bonds = atom.right_bonds
                elif adj_atom.row < atom.row:
                    adj_atom.bottom_bonds = added_bonds = atom.top_bonds
                else:
                    adj_atom.top_bonds = added_bonds = atom.bottom_bonds

                # Assuming correct input, we've at least removed added_bonds worth of bonds from
                # the molecule since we know any neighbor used to have at least one open side
                self.open_bonds -= added_bonds
                # If we closed off the neighbor's last open face, we've additionally removed
                # however many bonds it now has left from the molecule's 'open' bonds
                if not has_open_side(adj_atom):
                    self.open_bonds -= adj_atom.remaining_bonds()

        # Update this molecule's formula
        self.formula[atom.symbol] += 1

        self.atoms.append(atom)
