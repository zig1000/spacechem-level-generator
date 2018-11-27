#!/usr/bin/env python
# -*- coding: utf-8 -*-

import base64
import collections
import copy
import elements_data
import fractions
import gzip
import json
import StringIO

# Element class is defined in elements_data.py to avoid circular dependencies

'''Classes For more apropros errors'''
class FormulaValidityError(ValueError):
    pass
class MoleculeValidityError(ValueError):
    pass

class Formula(collections.Counter):
    '''Represent a chemical formula, as a Counter of elements.'''
    # Redefine Counter's built-in elements() method to return the list of unique ELements in the
    # formula, and move its original functionality to 'atoms()'.
    def elements(self):
        '''Return a list of unique elements in this formula.'''
        # Make sure not to include any 0-counts that Counter leaves hanging around
        return [e for e in self.keys() if self[e] != 0]

    def elements_collection(self):
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

    def least_common_formula(self):
        '''Return a new formula which is this formula divided by the GCD of its element counts'''
        new_formula = copy.copy(self)
        if new_formula:
            gcd = reduce(fractions.gcd, new_formula.values())
            for e in new_formula.elements():
                new_formula[e] = new_formula[e] / gcd
        return new_formula

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

    def is_valid(self, large_output=False):
        '''Check if it's possible to form a molecule with this formula within an input/output zone.
        Empty formulas are considered invalid. Default 4x4 zone, optionally large (8x4) output zone.
        '''
        # Verify size constraints
        if not (1 <= self.num_atoms() <= 16 + 16*large_output):
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

        # Due to the limited sizes of the zones, each max bound count element can only contribute
        # extra connections so many times before it reaches its limit. For example, in a 4x4 zone,
        # only two atoms with max bond count 4 will provide two extra endpoints each. Subsequent
        # max bond count 4 atoms will only allow for 1 extra endpoint each due to space
        # constraints. E.g. C3H8 is impossible to construct in a 4x4 zone.
        #
        if large_output:
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

    def fission_sources(self):
        '''Return a dict of atomic masses and their counts for all elements that could have fission
        performed on them to obtain part of this formula.
        The heavy lifting of this method is tucked away at the bottom of this file in
        splittable_sources since it's a monster of a function.
        '''
        output_masses = collections.Counter({e.atomic_num: count for e, count in self.items()
                                             if count > 0})
        return splittable_sources(output_masses)

    def remove_fissile_element(self, element, count):
        '''Used while converting an output formula to an input formula via inverse fission.
        Given a fissile element and count, remove the target count of the element from
        this formula, drilling down into its 'fission tree' as needed.
        Raise FormulaValidityError if the element / its fission tree doesn't add up to the count.
        '''
        # Remove as much as we can of this element without fission
        direct_removals = min(count, self[element])
        count -= direct_removals
        self[element] -= direct_removals
        # Keep this object clean
        if self[element] == 0:
            del self[element]

        if count != 0:
            # If we hit the bottom of the fission tree and aren't done, raise an exception
            if element.atomic_num == 1:
                raise FormulaValidityError("Couldn't remove {0} of {1} from formula".format(count, element))

            try:
                if element.atomic_num % 2 == 0:
                    child = elements_data.elements_dict[element.atomic_num / 2]
                    self.remove_fissile_element(child, 2*count)
                else:
                    child_A, child_B = (elements_data.elements_dict[element.atomic_num / 2 + 1],
                                        elements_data.elements_dict[element.atomic_num / 2])
                    self.remove_fissile_element(child_A, count)
                    self.remove_fissile_element(child_B, count)
            except FormulaValidityError:
                raise FormulaValidityError("Couldn't remove {0} of {1} from formula".format(count, element))

# Enum-esque directional vars for convenience
DIRECTIONS = UP, RIGHT, DOWN, LEFT = (0, 1, 2, 3) # Python, baby

# Doing a proper Direction(Int) subclass requires too much boilerplate for my taste
def opposite_dir(dir):
    '''Given an Int representing a direction return its opposite direction.'''
    if dir % 2 == 0:
        return 2 - dir # UP <-> DOWN
    else:
        return 4 - dir # LEFT <-> RIGHT

class GridPos:
    '''Represent a 0-indexed (row, col) position within an input/output zone.
    Indices increase from left to right and top to bottom.
    '''
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
        '''Check that this position consists of integer positions within the zone's grid.'''
        return isinstance(self.row, int) and isinstance(self.col, int) \
               and (0 <= self.row < self.num_rows) and (0 <= self.col < self.num_cols)

    def neighbor(self, dir):
        '''Return the neighbor GridPos in the indicated direction, or None if out-of-bounds.'''
        if dir == UP:
            r, c = self.row - 1, self.col
        elif dir == RIGHT:
            r, c = self.row, self.col + 1
        elif dir == DOWN:
            r, c = self.row + 1, self.col
        elif dir == LEFT:
            r, c = self.row, self.col - 1
        else:
            raise ValueError("Invalid direction: " + str(dir))

        if 0 <= r < self.num_rows and 0 <= c < self.num_cols:
            return GridPos(r, c, self.large_output)
        return None

    def neighbors(self):
        '''Return all orthogonally adjacent positions within the zone's grid.'''
        return [p for p in (self.neighbor(dir) for dir in DIRECTIONS) if p is not None]

    def dirs_and_neighbors(self):
        '''Return a list of (dir, pos) pairs for each neighboring position within the grid.'''
        return [(d, p) for d, p in ((_d, self.neighbor(_d)) for _d in DIRECTIONS) if p is not None]


class Atom:
    '''Represent an Atom, including its element, grid position, and attached bonds.
    '''
    def __init__(self, element, pos):
        self.bonds = [0, 0, 0, 0]
        self.set_element(element)
        self.set_pos(pos)

    def __str__(self):
        return self.symbol.rjust(2) # Pad element symbol to two chars

    def __repr__(self):
        return 'Atom({0}, {1}, {2})'.format(self.symbol, self.pos, self.bonds)

    def get_json_str(self):
        '''Return a string representing this atom in the level json's format.'''
        return '{0}{1}{2}{3}{4}'.format(self.col, self.row,
                                        str(self.atomic_num).rjust(3, '0'),
                                        self.bonds[RIGHT],
                                        self.bonds[DOWN])

    def remaining_bonds(self):
        '''Return the # of remaining bonds this atom is allowed.'''
        return self.max_bonds - sum(self.bonds)

    def set_pos(self, pos):
        '''Change this atom's position in the grid.'''
        self.pos = pos
        self.row = self.pos.row
        self.col = self.pos.col

    def set_element(self, element):
        if sum(self.bonds) > element.max_bonds:
            raise ValueError("Too many bonds to change atom {} to element {}".format(self, element))

        self.element = element

        # Exposing some sub-attributes for convenience
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
        self.formula = Formula()
        # To optimize the performance of available_positions(), we'll roughly track the # of open
        # bonds available on this molecule.
        # An atom with no open adjacencies in the grid contributes 0 to this count.
        self.open_bonds = 0

    def __getitem__(self, pos):
        '''Return the atom at the specified grid position or None.'''
        return self.grid[pos.row][pos.col]

    def __setitem__(self, pos, item):
        '''Set the specified grid position (item should be None or an Atom).'''
        self.grid[pos.row][pos.col] = item

    def __iter__(self):
        for atom in self.atoms:
            yield atom

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
                left_atom = atom
                right_atom = self.grid[r][c + 1] if c + 1 < self.num_cols else None

                bond_str = ' '
                if left_atom is not None and right_atom is not None \
                   and left_atom.bonds[RIGHT] != right_atom.bonds[LEFT]:
                    bond_str = '?'
                elif left_atom is not None and left_atom.bonds[RIGHT] != 0:
                    bond_str = str(left_atom.bonds[RIGHT])
                elif right_atom is not None and right_atom.bonds[LEFT] != 0:
                    bond_str = str(right_atom.bonds[LEFT])
                if c < self.num_cols - 1:
                    result += ' ' + bond_str + ' '
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
                       and top_atom.bonds[DOWN] != bottom_atom.bonds[UP]:
                        bond_str = '??'
                    elif top_atom is not None and top_atom.bonds[DOWN] != 0:
                        bond_str = ' ' + str(top_atom.bonds[DOWN])
                    elif bottom_atom is not None and bottom_atom.bonds[UP] != 0:
                        bond_str = ' ' + str(bottom_atom.bonds[UP])
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

    def update_formula(self):
        '''To be called after mutating any atom elements. Update the formula of this molecule.'''
        self.formula = Formula()
        for atom in self.atoms:
            self.formula[atom.element] += 1

    def update_open_bonds(self):
        '''Update the count of open bonds. Since we only care about updating it well
        enough to know when it's 0, we'll ignore the triple bond limit, and count any open side of
        an atom as adding the remainder of its max bond count to the open bonds.
        '''
        self.open_bonds = 0
        for atom in self.atoms:
            if any(self[pos] is None for pos in atom.pos.neighbors()):
                self.open_bonds += atom.remaining_bonds() # Not exact but we don't need it to be

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
        '''Adds the given Atom to this molecule. The Atom's position must be open in this molecule.
        Also adds any bonds specified by the incoming atom to its neighboring atoms.
        For convenience of more complex operations, it is allowable to add an atom with unfulfilled
        bonds or which is not connected to the rest of the molecule.
        '''
        if self[new_atom.pos] is not None:
            raise Exception("Conflict with existing atom; cannot add {0} to \n{1}".format(repr(new_atom), self))
        self[new_atom.pos] = new_atom

        # Quick helper to check if an atom within this molecule's grid has at least 1 open side
        def has_open_side(atom):
            return any(self[pos] is None for pos in atom.pos.neighbors())

        # Partial update of the number of open bonds this molecule has
        if has_open_side(new_atom):
            self.open_bonds += new_atom.remaining_bonds()

        # Add bonds to neighbours matching the bonds indicated on this atom
        for dir, pos in new_atom.pos.dirs_and_neighbors():
            adj_atom = self[pos]
            if adj_atom is not None:
                adj_atom.bonds[opposite_dir(dir)] = new_atom.bonds[dir]
                # Subtract the bond we just added from the molecule's 'open bonds'
                self.open_bonds -= new_atom.bonds[dir]

                # If we closed off the neighbor's last open face, we've additionally removed
                # however many bonds it now has left from the molecule's 'open' bonds
                if not has_open_side(adj_atom):
                    self.open_bonds -= adj_atom.remaining_bonds()

        # Update this molecule's formula
        self.formula[new_atom.element] += 1

        self.atoms.append(new_atom)

    def is_connected(self):
        '''For the purposes of more advanced construction algorithms we allow adding atoms in
        unconnected cells. This checks if the molecule is currently 'connected' and thus valid.
        We'll count empty molecules as unconnected.
        '''
        if len(self) == 0:
            return False

        # Do a DFS starting from one atom and following the bonds of the molecule. If we don't
        # find every atom, it's not connected
        stack = [self.atoms[0]]
        # We don't have to actually 'visit' every atom, seeing them as neighbors is sufficient
        seen = {self.atoms[0].pos: True} # Track the grid positions of seen connected atoms
        while stack:
            if len(seen) == len(self):
                return True

            atom = stack.pop()
            # Check for connected neighbors. When we see an unseen connected atom, add it to the
            # stack
            for dir, adj_pos in atom.pos.dirs_and_neighbors():
                if atom.bonds[dir] != 0 and adj_pos not in seen:
                    seen[adj_pos] = True
                    adj_atom = self[adj_pos]
                    stack.append(adj_atom)
        return False

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
        self.update_open_bonds()

    def add_molecule(self, other):
        '''Add the specified molecule to this molecule. Must not have any atoms in conflicting
        positions.
        '''
        # Check for conflicts
        if any(self[atom.pos] is not None for atom in other.atoms):
            raise Exception('Cannot add molecule \n{0} to molecule \n{1}; conflicting atoms'.format(other, self))

        # Add the new atoms
        for atom in other.atoms:
            self.add_atom(atom)


class Level:
    '''Parent class for Research and Production levels.'''
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


def splittable_sources(given):
    '''Given a Counter of ints, return a dict of ints and their total counts that can be
    non-trivially constructed from the given integers, using only series' of addition of integers
    within 1 of each other, and where no given int is used more times than its count in the Counter.
    In other words, we're using the inverse fission operation to calculate viable input elements
    that could have been split any (non-0) # of times to create part of the given output.
    '''
    # NOTE: We ask for a Counter as input because unlike dicts, they implicitly return 0 if
    #       asked for a key they don't contain, which simplifies our code

    # Tally tracking ints we were given/discover, the max of each we can create at once, and
    # some additional helper values - as needed we'll create dicts that track how many N-1 ints
    # we can create at the same time as any possible count of N's. We'll also track the
    # 'most balanced' such allocation of N vs N-1, which will allow us to properly get counts for
    # odd ints, as well as assisting in the creation of higher dicts.
    tally = {}
    # Dict of ints that were constructable (not just given), and their max counts
    constructed = {}

    # Keep a running sum of what we were given so we don't waste time on clearly impossible sums
    givens_sum = 0

    # To avoid the overhead of a priority queue, use one queue for the given ints,
    # and one queue for ints we obtained by addition (we'll add them in numerical order).
    # Loop on whichever's next value is lower
    given_queue = collections.deque(sorted(given.keys()))
    added_queue = collections.deque()
    while given_queue or added_queue:
        # Pop the element we're iterating on - pop from both queues at once if they match
        if (not added_queue
            or (given_queue and given_queue[0] < added_queue[0])):
            n = given_queue.popleft()
        else:
            if given_queue and given_queue[0] == added_queue[0]:
                given_queue.popleft()
            n = added_queue.popleft()

        # Calculate how many copies of n we can obtain at once
        if n % 2 == 0:
            # If n is even, we only need to update its count, based on
            # the count of n / 2 and however much of n we were given to start
            component_count = tally[n / 2]['count'] if n / 2 in tally else 0
            this_count = component_count / 2 + given[n]
        else:
            # To count odd n, we must first make a dict that pairs the max # of n / 2 that can
            # be created for any given count of n / 2 + 1. We can do this recursively off
            # a previous such dict. When creating this dict we will also store the count
            # of n / 2 + 1 for which there can simultaneously be created a most closely balanced
            # count of n / 2. This can be used directly to count n and also to speed up
            # the recursive creation of dicts.
            # However we can skip this if either of n / 2 or n / 2 + 1 are unavailable.
            # Note that even if both are available they may not be addable so our count could
            # still come out to 0
            upper_child, lower_child = n / 2 + 1, n / 2
            if upper_child in tally and lower_child in tally:
                # In this case, calculate and store upper_child's neighbour dict
                tally[upper_child]['neighbour_counts'] = {}
                if upper_child == 2: # Calc 2->1 dict
                    balanced_upper_count, min_count_in_best_balance = -1, -1
                    for upper_count in range(1, tally[upper_child]['count'] + 1):
                        lower_count = (tally[lower_child]['count']
                                       - 2*max(upper_count - given[upper_child], 0))
                        tally[upper_child]['neighbour_counts'][upper_count] = lower_count

                        # Check how balanced this allocation is
                        worst_count = min(upper_count, lower_count)
                        if worst_count >= min_count_in_best_balance:
                            balanced_upper_count = upper_count
                            min_count_in_best_balance = worst_count
                    tally[upper_child]['balanced_count'] = balanced_upper_count
                elif upper_child == 3: # Calc 3->2 dict
                    balanced_upper_count, min_count_in_best_balance = -1, -1
                    for upper_count in range(1, tally[upper_child]['count'] + 1):
                        nongiven_upper_count = max(upper_count - given[upper_child], 0)
                        # 2s count = (2s constructable given 1s used in 3s) - (2s used in 3s)
                        lower_count = (given[2] + (given[1] - nongiven_upper_count) / 2
                                       - nongiven_upper_count)
                        tally[upper_child]['neighbour_counts'][upper_count] = lower_count

                        # Check how balanced this allocation is
                        worst_count = min(upper_count, lower_count)
                        if worst_count >= min_count_in_best_balance:
                            balanced_upper_count = upper_count
                            min_count_in_best_balance = worst_count
                    # Store the most balanced count for upper_child
                    tally[upper_child]['balanced_count'] = balanced_upper_count
                # If either the upper child or the lower child had no compound components,
                # the upper_child's neighbour_counts dict is just the max count of lower_child,
                # regardless of the count of upper_child
                elif (tally[upper_child]['count'] == given[upper_child]
                      or tally[lower_child]['count'] == given[lower_child]):
                    tally[upper_child]['neighbour_counts'] = {
                            i: tally[lower_child]['count']
                            for i in range(1, tally[upper_child]['count'] + 1) }
                    # Since the lower_child gets the same count no matter what, just maximize
                    # upper_child's count for the 'balanced' allocation
                    tally[upper_child]['balanced_count'] = tally[upper_child]['count']
                # Otherwise, based on our recursion principle, the upper child's upper
                # child must already have its neighbour_counts dict set. Use that to calculate
                # the upper child's neighbour_counts. The algorithm for this depends on which of
                # upper_child/lower_child is even.
                # We also have a couple of base cases to handle when building the neighbour dict
                # dict for 3->2 and 2-> 1, since in those cases lower_child is also a component
                # of upper_child
                elif upper_child % 2 == 0:
                    # If the upper child is even, calculate how much of lower_child's components
                    # are used up by any valid count of upper_child, and thus the max
                    # lower_child count for that count of upper_child.
                    # Call A upper_child / 2 and B the other (lower) component of lower_child
                    A = upper_child / 2
                    balanced_upper_count, min_count_in_best_balance = -1, -1
                    for upper_count in range(1, tally[upper_child]['count'] + 1):
                        if upper_count <= given[upper_child]:
                            lower_count = tally[lower_child]['count']
                        else:
                            A_used_in_upper = 2*(upper_count - given[upper_child])
                            if A_used_in_upper == tally[A]['count']:
                                lower_count = given[lower_child]
                            else:
                                # Search to the right of the original balance point and/or our
                                # new A limit, to find a balance given the unusable As:
                                start_idx = max(tally[A]['balanced_count'], A_used_in_upper + 1)
                                built_lower_count = 0
                                for used_A in range(start_idx, tally[A]['count'] + 1):
                                    worst_count = min(used_A - A_used_in_upper,
                                                      tally[A]['neighbour_counts'][used_A])
                                    built_lower_count = max(built_lower_count, worst_count)
                                lower_count = built_lower_count + given[lower_child]
                        tally[upper_child]['neighbour_counts'][upper_count] = lower_count

                        # Check how balanced this allocation is
                        worst_count = min(upper_count, lower_count)
                        if worst_count >= min_count_in_best_balance:
                            balanced_upper_count = upper_count
                            min_count_in_best_balance = worst_count
                    # Store the most balanced count for upper child
                    tally[upper_child]['balanced_count'] = balanced_upper_count
                else:
                    # If the upper child is odd, call its subchildren A and B, where B =
                    # lower_child / 2. Using A's neighbour_counts dict, calculate how much B and
                    # from that how much lower_child we can make for any valid count of
                    # upper_child
                    A = upper_child / 2 + 1
                    # For each possible count of upper_child, calculate how many copies of
                    # lower_child can be simultaneously constructed from the leftovers
                    # Also track the 'most balanced' count we can assign to upper_child
                    balanced_upper_count, min_count_in_best_balance = -1, -1
                    for upper_count in range(1, tally[upper_child]['count'] + 1):
                        if upper_count <= given[upper_child]:
                            lower_count = tally[lower_child]['count']
                        else:
                            used_A = used_B = upper_count - given[upper_child]
                            available_B = (tally[A]['neighbour_counts'][used_A] - used_A)
                            lower_count = available_B / 2 + given[lower_child]
                        tally[upper_child]['neighbour_counts'][upper_count] = lower_count

                        # Check how balanced this allocation is
                        worst_count = min(upper_count, lower_count)
                        if worst_count >= min_count_in_best_balance:
                            balanced_upper_count = upper_count
                            min_count_in_best_balance = worst_count
                    # Store the most balanced count for upper child
                    tally[upper_child]['balanced_count'] = balanced_upper_count

                # Calculate the count of n based on upper_child's most balanced count
                balanced_upper_count = tally[upper_child]['balanced_count']
                this_count = (min(balanced_upper_count,
                                  tally[upper_child]['neighbour_counts'][balanced_upper_count])
                              + given[n])
            else:
                # If n only occurred as an input and not a compound, set it to the given count
                # The n = 1 case is handled here since 1 can never be compound
                # We don't need to calculate its neighbour dict in this case.
                this_count = given[n]

        # If the count came out to 0, ignore this int
        if this_count == 0:
            continue
        # Update the tally with the discovered count
        tally[n] = {'count': this_count}
        #  Add this int to the output dict if it was possible to construct
        # (not just obtained from the givens)
        if this_count != given[n]:
            constructed[n] = this_count

        # Add any viable sums (restricted to valid atomic masses) obtained from n to the queue
        # As a mini-optimization, we won't add odd numbers to the queue that exceed the sum of
        # the givens up to n
        givens_sum += n*given[n]
        # If n - 1 is in the tally, add 2n - 1 to the queue
        if (n - 1 in tally
            and (2*n - 1 <= 109 or 2*n - 1 in (201, 203))
            and 2*n - 1 <= givens_sum):
            added_queue.append(2*n - 1)
        # If the count for n was at least 2, add 2n to the queue
        if tally[n]['count'] >= 2 and (2*n <= 109 or 2*n in (200, 202)):
            added_queue.append(2*n)

    # Once we've looped over all possible sums, return a dict of the relevant ints and their counts
    return constructed
