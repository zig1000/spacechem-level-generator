# Spacechem Random Level Generator

# Usage

`python level_generator.py`

Prints the mission code for a research level with randomly assigned input zone and output zones,
with molecules sizes probabilistically based on the difficulty.

To quickly view the level, try Cearn's handy mission viewer:
http://www.coranac.com/spacechem/mission-viewer?mode=editor

*Note*: there is a small chance for the tool to fail upon attempting to construct too complex an
        output molecule (e.g. a 16-atom molecule, half of which are Hydrogens), but the probability
        of this is pretty statistically insignificant even at high difficulties.

Developed with python 2.7 (.10)

## Args

--difficulty [int] : Int from 0 to 3 approximately controlling the level complexity
                     (basically just molecule sizes, at the moment). Default 2.

--lanky : Set the difficulty to the maximum value

TBA:
 --production
 --fusion
 --fission
 --inputs ('alpha', 'beta', 'both')
 --outputs ('phi', 'omega', 'both')
 --elements (list of elements to exclusively use, e.g. ['C','H'])
 --num_atoms
 --large-output (make the output zone large and scale up the output molecule's size)

# Behavior

The process for generating a level is as follows:

* Randomly draw # of input zones - currently 33% chance of both input zones being used.

* Randomly generate input molecules with sizes based on the difficulty. Currently just chooses from
  what I consider the 'basic' elements of non-0 max-bond counts - H, O, B, C, N, S, Cl, Os.
  Added atoms are biased toward reusing elements that have already appeared in the molecule.

* Randomly generate balanced ratios for the reaction. To achieve this, first divide each input
  molecule's formula by its GCD (e.g. C2H4 / 2 -> CH2 or O3 / 3 -> O) to get that input's atom
  ratio, then multiply each of these simplified molecules by a random number (higher for higher
  difficulty, not exceeding 16 atoms total per input), representing how many are needed for a
  wasteless reaction.
  Adding those molecules up gives a chemical formula for the atoms to be used by the outputs.

* Randomly select what output zones to use - again, 33% chance of both.

* Randomly allocate the input atoms between the output zones (currently uses a uniform split, i.e.
  no attempt to equalize how many atoms are assigned to each zone), making sure not to give an
  un-constructable set of atoms. Then once again divide each output by its GCD and multiply each by
  a random number, weighted higher for higher difficulties.
  This creates the list of atoms to be used by the output zone's molecule.

* Construct each output molecule by drawing without replacement from the atoms assigned to it (retry
  the construction until all have been placed successfully).

* Randomly add level features (currently just bonders and teleporter)

* Print mission code
