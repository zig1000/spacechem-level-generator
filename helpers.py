#!/usr/bin/env python
# -*- coding: utf-8 -*-

import base64
import collections
import fractions
import gzip
import json
import StringIO

class Element:
    '''Class representing a SpaceChemical element.'''
    def __init__(self, atomic_num, symbol, max_bonds):
        self.atomic_num = atomic_num
        self.symbol = symbol
        self.max_bonds = max_bonds

    def __str__(self):
        return self.symbol
    __repr__ = __str__

    # These are necessary so we can use Elements as dict keys
    def __eq__(self, other):
        return type(self) == type(other) and self.atomic_num == other.atomic_num

    def __hash__(self):
        return hash(self.atomic_num)

class Formula(collections.Counter):
    '''Represent a chemical formula, as a Counter of elements.'''
    def __init__(self, *args, **kwargs):
        '''For checking validity, need to know the size of the zone'''
        if 'large_output' in kwargs:
            self.large_output = kwargs['large_output']
            del kwargs['large_output']
        else:
            self.large_output = False
        super(Formula, self).__init__(*args, **kwargs)

    # Redefine Counter-s built-in 'elements()' method to return the list of unique ELements in the
    # formula, and move its original functionality to 'atoms()'.
    def elements(self):
        '''Return a list of unique elements in this formula.'''
        # Make sure not to include any 0-counts that Counter leaves hanging around
        return [e for e in self.keys() if self[e] != 0]

    def atoms(self):
        '''Return a list containing each element as many times as its count.'''
        return list(super(Formula, self).elements())

    # Have to override Counter's add method or else adding two Formulas will make a Counter
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

    def num_atoms(self):
        return sum(self.values())

    def divide_by_gcd(self):
        '''Divide this formula's numeric values by their GCD. Return the GCD or 0 if this formula
        is empty.'''
        if self:
            gcd = reduce(fractions.gcd, self.values())
            for e in self.elements():
                self[e] = self[e] / gcd
            return gcd
        return 0

    def get_json_str(self):
        '''Return a string representing this formula using the Hill System (C, then H, then
        alphabetized), in the game's accepted format. E.g. Glycine would be 'C~02H~05NO~02'.
        '''
        result = ''
        # Sort Carbon and Hydrogen to the front and alphabetize the rest
        elements = sorted(self.elements(),
                          key=lambda e: 0 if e.symbol == 'C'
                                        else 1 if e.symbol == 'H'
                                        else e.symbol)
        for element in elements:
            result += element.symbol
            if self[element] != 1:
                result += '~' + str(self[element]).rjust(2, '0')
        return result
    __str__ = get_json_str # For debugging convenience

    def is_valid(self):
        '''Check if it's possible to form a molecule with this formula within a 4x4 grid.
        Empty formulas are considered invalid.
        '''
        # Verify size constraints
        if not (1 <= self.num_atoms() <= 16 + 16*self.large_output):
            return False

        # We'll calculate validity based on whether there are enough bonds in the atom list to form
        # a minimally-connected molecule. To check this, we start with a simple linearly connected
        # molecule with 2 endpoints. An endpoint element has max bonds 1, while elements with 3 or
        # 4+ max bonds each allow for 1 or 2 additional endpoints in the molecule, respectively.
        # Elements of max bonds 2 do not affect the count and thus do not affect the validity of
        # the formula, apart from the 16 or 32-atom limit of the zone.
        #
        # Though not formally proven, it appears that the maximum number of endpoints is in each
        # case a constant equal to half the number of cells in the zone. By trial and error with the
        # 4x4 zone it appears that this maximum can be reached with any composition of 3 vs 4 bond
        # atoms, and a simplifying assumption is made that this holds for the 8x4 case.
        # The cases for which an incorrect return value by this method could cause an exception
        # would in any case be prohibitively rare and of little concern, and well worth the tradeoff
        # for O(k) runtime (where k is the # of unique elements in the formula).
        num_atoms = self.num_atoms()
        usable_connections = 0
        if self.large_output:
            extra_endpoints_dict = {3:14, 4:6}
        else:
            extra_endpoints_dict = {3:6, 4:2}

        allowed_endpoint_count = 2 # H-H base case
        for element in self.elements():
            if element.max_bonds == 0:
                # A noble gas is only valid as the only atom in a molecule
                return self.num_atoms() == 1
            elif element.max_bonds == 1:
                allowed_endpoint_count -= self[element]
            # Count one extra endpoint per atom with bond count 3 or more (up to limit)
            elif element.max_bonds >= 3:
                # As long as our formula has no negatives this should be fine
                extra_endpoints = min(self[element],
                                      extra_endpoints_dict[3])
                allowed_endpoint_count += extra_endpoints
                extra_endpoints_dict[3] -= extra_endpoints

                # Count an additional extra endpoint per atom with bond count 4 or more (up to limit)
                if element.max_bonds >= 4:
                    extra_endpoints = min(self[element],
                                          extra_endpoints_dict[3],
                                          extra_endpoints_dict[4])
                    allowed_endpoint_count += extra_endpoints
                    extra_endpoints_dict[3] -= extra_endpoints
                    extra_endpoints_dict[4] -= extra_endpoints

        return allowed_endpoint_count >= 0


class GridPos:
    '''Represent a 0-indexed (row, col) position within a 4x4 input/output zone.'''
    def __init__(self, row, col, large_output=False):
        self.row = row
        self.col = col
        self.large_output = large_output
        self.num_rows = 4 + 4*large_output
        self.num_cols = 4

    def __str__(self):
        return '({0}, {1})'.format(self.row, self.col)
    __repr__ = __str__

    # __eq__ and __hash__ so we can use GridPos as dictionary keys
    def __eq__(self, other):
        return type(self) == type(other) and (self.row, self.col) == (other.row, other.col)

    def __hash__(self):
        return hash((self.row, self.col))

    def is_valid(self):
        '''Check that this position consists of integer positions within the zone's grid.
        '''
        return isinstance(self.row, int) and isinstance(self.col, int) \
               and (0 <= self.row < self.num_rows) and (0 <= self.col < self.num_cols)

    def neighbors(self):
        '''Return all orthogonally adjacent positions within the zone's grid.
        '''
        result = []
        for row, col in [(self.row, self.col - 1), # left
                         (self.row, self.col + 1), # right
                         (self.row - 1, self.col), # top
                         (self.row + 1, self.col)]: # bottom
            if 0 <= row < self.num_rows and 0 <= col < self.num_cols:
                result.append(GridPos(row, col, large_output=self.large_output))
        return result


class Atom:
    '''Represent an Atom, including its element, grid position, and attached bonds.
    '''
    def __init__(self, element, pos):
        self.element = element
        self.pos = pos
        self.left_bonds = 0
        self.right_bonds = 0
        self.top_bonds = 0
        self.bottom_bonds = 0

        # Exposing some sub-attributes for convenience
        self.atomic_num = element.atomic_num
        self.symbol = element.symbol
        self.max_bonds = element.max_bonds
        self.row = pos.row
        self.col = pos.col

    def __str__(self):
        return self.symbol.rjust(2) # Pad element symbol to two chars

    def __repr__(self):
        return 'Atom({0}, {1}, {2}, {3}, {4}, {5})'.format(self.symbol, self.pos,
                                                           self.left_bonds, self.top_bonds,
                                                           self.right_bonds, self.bottom_bonds)

    def get_json_str(self):
        '''Return a string representing this atom in the level json's format.'''
        return '{0}{1}{2}{3}{4}'.format(self.col, self.row,
                                        str(self.atomic_num).rjust(3, '0'),
                                        self.right_bonds,
                                        self.bottom_bonds)

    def remaining_bonds(self):
        '''Return the # of remaining bonds this atom is allowed.'''
        return self.max_bonds - self.left_bonds \
                              - self.right_bonds \
                              - self.top_bonds \
                              - self.bottom_bonds

    def set_pos(self, pos):
        '''Change this atom's position in the grid.'''
        self.pos = pos
        self.row = self.pos.row
        self.col = self.pos.col

    def set_element(self, element):
        '''Mutate this atom to a different element. That element must not have a lower max_bonds
        than the current number of bonds on this atom
        '''
        if element.max_bonds < (self.left_bonds + self.right_bonds
                                + self.top_bonds + self.bottom_bonds):
            raise Exception("Cannot set atom {0} to element {1} with insufficient bond count".format(self, element))
        self.element = element
        self.atomic_num = element.atomic_num
        self.symbol = element.symbol
        self.max_bonds = element.max_bonds


class Molecule:
    '''Represents an input/output zone and the molecule constructed therein.
    '''
    def __init__(self, large_output=False):
        self.name = 'Randite'
        self.atoms = [] # Order of atoms in this list is not guaranteed
        self.large_output = large_output
        self.num_rows = 4 + 4*large_output
        self.num_cols = 4
        self.grid = [[None, None, None, None] for r in range(self.num_rows)]
        self.formula = Formula(large_output=large_output)
        self.open_bonds = 0 # A meausure of the # of open bonds available
                            # An atom with no open adjacencies contributes 0 to this count.

    def __getitem__(self, pos):
        '''Return the atom at the specified grid position or None.'''
        return self.grid[pos.row][pos.col]

    def __setitem__(self, pos, item):
        '''Set the specified grid position (item should be None or an Atom).'''
        self.grid[pos.row][pos.col] = item

    def __len__(self):
        '''Return the # of atoms in this molecule.'''
        return len(self.atoms)

    def __str__(self):
        '''Pretty-print this molecule.'''
        result = ' _________________ \n' # Border of the input/output zone
        for r in range(self.num_rows):
            result += '|'
            for c in range(self.num_cols):
                atom = self.grid[r][c]
                # Represent any atoms here
                if atom is None:
                    result += 2*' '
                else:
                    result += str(atom)
                # Represent any bonds to the right of the atom
                if c < self.num_cols - 1:
                    if atom is not None and atom.right_bonds > 0:
                        result += ' ' + str(atom.right_bonds) + ' '
                    else:
                        result += 3*' '
            result += '|\n'
            # Add a row of vertical bonds
            if r < self.num_rows - 1:
                result += '|'
                for c in range(self.num_cols):
                    top_atom = self.grid[r][c]
                    if r + 1 < self.num_rows:
                        bottom_atom = self.grid[r + 1][c]
                    else:
                        bottom_atom = None
                    bond_str = '  '
                    if top_atom is not None and bottom_atom is not None \
                       and top_atom.bottom_bonds != bottom_atom.top_bonds:
                        bond_str = '??'
                    elif top_atom is not None and top_atom.bottom_bonds != 0:
                        bond_str = ' ' + str(top_atom.bottom_bonds)
                    elif bottom_atom is not None and bottom_atom.top_bonds != 0:
                        bond_str = ' ' + str(bottom_atom.top_bonds)
                    result += bond_str
                    if c < self.num_cols - 1:
                        result += 3*' '
                result += '|\n'
        result += '|_________________|\n'
        return result

    __repr__ = __str__

    def get_json_str(self):
        '''Return a string representing this molecule in the level json's format.'''
        result = self.name + ';' + self.formula.get_json_str() + ';'
        for atom in self.atoms:
            result += atom.get_json_str() + ';'
        return result

    def open_positions(self):
        '''Return a list of valid grid positions where an atom could be added to this molecule.'''
        # For an empty molecule, all positions are open
        if len(self) == 0:
            return [GridPos(r, c, large_output=self.large_output)
                    for r in range(self.num_rows) for c in range(self.num_cols)]
        # If there are no remaining bonds, we can skip the overhead of walking through the atoms
        elif self.open_bonds == 0:
            return []

        result_dict = {} # For O(1) checks on whether a position has already been added
        for atom in self.atoms:
            if atom.remaining_bonds() > 0:
                for pos in atom.pos.neighbors():
                   if self[pos] is None and pos not in result_dict:
                      result_dict[pos] = True
        return result_dict.keys()

    def add_atom(self, new_atom):
        if self[new_atom.pos] is not None:
            raise Exception("Conflict with existing atom; cannot add {0} to \n{1}".format(repr(atom), self))
        self[new_atom.pos] = new_atom

        # Quick helper to check if an atom within this molecule's grid has at least 1 open side
        def has_open_side(atom):
            for neighbor_posn in atom.pos.neighbors():
                if self[neighbor_posn] is None:
                    return True
            return False

        # Partial update of the number of open bonds this molecule has
        if has_open_side(new_atom):
            self.open_bonds += new_atom.remaining_bonds()

        for pos in new_atom.pos.neighbors():
            adj_atom = self[pos]

            # Add bonds to any neighbours of the added atom
            if adj_atom is not None:
                if adj_atom.col < new_atom.col:
                    adj_atom.right_bonds = added_bonds = new_atom.left_bonds
                elif adj_atom.col > new_atom.col:
                    adj_atom.left_bonds = added_bonds = new_atom.right_bonds
                elif adj_atom.row < new_atom.row:
                    adj_atom.bottom_bonds = added_bonds = new_atom.top_bonds
                else:
                    adj_atom.top_bonds = added_bonds = new_atom.bottom_bonds

                # Subtract the bond we just added from the molecule's 'open bonds'
                self.open_bonds -= added_bonds
                # If we closed off the neighbor's last open face, we've additionally removed
                # however many bonds it now has left from the molecule's 'open' bonds
                if not has_open_side(adj_atom):
                    self.open_bonds -= adj_atom.remaining_bonds()

        # Update this molecule's formula
        self.formula[new_atom.element] += 1

        self.atoms.append(new_atom)

    def is_connected(self):
        '''For more advanced construction algorithms we're going to allow adding atoms in
        unconnected cells. This checks if the molecule is currently 'connected' and thus valid.
        We'll count empty molecules as unconnected.
        '''
        if len(self) == 0:
            return False

        # Do a DFS starting from one atom and following the bonds of the molecule. If we don't
        # find every atom, it's not connected
        stack = [self.atoms[0]]
        visited = {} # Track the grid positions of visited atoms
        while True:
            if len(visited) == len(self):
                return True
            elif not stack:
                return False

            atom = stack.pop()
            if atom.pos not in visited:
                visited[atom.pos] = True

            # Check for connected neighbors
            for pos in [p for p in atom.pos.neighbors() if p not in visited]:
                adj_atom = self[pos]
                if adj_atom is not None:
                    if (adj_atom.col < atom.col and atom.left_bonds != 0) \
                         or (adj_atom.col > atom.col and atom.right_bonds != 0) \
                         or (adj_atom.row < atom.row and atom.top_bonds != 0) \
                         or (adj_atom.row > atom.row and atom.bottom_bonds != 0):
                       stack.append(adj_atom)

    def shift(self, rows=0, cols=0):
        '''Shift the current contents of this molecule downward/rightward by the specified number
        of rows/columns. Negative numbers shift upward/leftward.

        Raise an exception if this would place atoms out-of-bounds.
        '''
        # Make sure this is a legal shift
        for atom in self.atoms:
            if (atom.row + rows < 0 or atom.row + rows > self.num_rows) \
               or (atom.col + cols < 0 or atom.col + cols > self.num_cols):
               raise Exception('Cannot shift molecule {0} by {1} rows and {2} cols'.format(self, rows, cols))

        # Wipe the grid clean and re-add the atoms in the new positions
        self.grid = [[None, None, None, None] for r in range(self.num_rows)]
        for atom in self.atoms:
            atom.set_pos(GridPos(atom.row + rows, atom.col + cols, large_output=self.large_output))
            self[atom.pos] = atom

        # Recount open bonds once we're done since some atoms may no longer have open sides
        self.open_bonds = 0
        for atom in self.atoms:
            if [True for p in atom.pos.neighbors() if self[p] is None]:
                self.open_bonds += atom.remaining_bonds()

    def add_molecule(self, other):
        '''Add the specified molecule to this molecule. Must not have any atoms in conflicting
        positions.
        '''
        # Check for conflicts
        if any([self[atom.pos] is not None for atom in other.atoms]):
            raise Exception('Cannot add molecule \n{0} to molecule \n{1}; conflicting atoms'.format(other, self))

        # Add the new atoms
        for atom in other.atoms:
            self.add_atom(atom)

class Level:
    '''Parent class for research and production levels.'''
    def __init__(self):
        self.dict = {}
        self.dict['name'] = 'RandomlyGenerated'
        self.dict['author'] = "Zig's Random"
        self.dict['difficulty'] = 0

    def __getitem__(self, item):
        return self.dict[item]

    def __setitem__(self, item, val):
        self.dict[item] = val

    def __str__(self):
        return json.dumps(self.dict)

    def get_code(self):
        '''Get the mission code - gzip then b64 the level json.'''
        out = StringIO.StringIO()
        with gzip.GzipFile(fileobj=out, mode="w") as f:
            f.write(json.dumps(self.dict))
        return base64.b64encode(out.getvalue())

class ResearchLevel(Level):
    def __init__(self):
        Level.__init__(self)
        self.dict['type'] = 'research'
        self.dict['input-zones'] = {}
        self.dict['output-zones'] = {}

        # Features of the level
        self.dict['bonder-count'] = 0
        self.dict['has-sensor'] = False
        self.dict['has-fuser'] = False
        self.dict['has-splitter'] = False
        self.dict['has-teleporter'] = False

        self.dict['has-large-output'] = False

class ProductionLevel(Level):
    def __init__(self):
        Level.__init__(self)
        self.dict['type'] = 'production'
        self.dict['terrain'] = 0
        self.dict['random-input-zones'] = {} # Max 1 random
        self.dict['fixed-input-zones'] = {} # Max 2 fixed
        self.dict['output-zones'] = {} # Max 3 outputs

        self.dict['max-reactors'] = 6 # Default maximum allowed

        self.dict['has-starter'] = False
        self.dict['has-assembly'] = False
        self.dict['has-disassembly'] = False
        self.dict['has-advanced'] = False # Sensor reactor
        self.dict['has-nuclear'] = False
        self.dict['has-superbonder'] = False

        self.dict['has-recycler'] = False
