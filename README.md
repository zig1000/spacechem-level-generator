# Spacechem Random Level Generator

## Installs

Built in python 2.7 (.10)

### periodictable

A python package used to get the atomic symbols of all elements

 `pip install periodictable`

## Usage

`python level_generator.py`

Default creates a Research level. See args for options.

*Note*: there is a small chance for the tool to fail upon attempting to construct too complex an
        output molecule (e.g. a 16-atom molecule, half of which are Hydrogens), but the probability
        of this is pretty statistically insignificant even at high difficulties.

## Args

--view : Auto-open Cearn's mission viewing website with the generated code

--production : Create a production level instead of a research level

--difficulty [int] : Int from 0 to 3 approximately controlling the level complexity
                     (basically just molecule sizes, at the moment).
                     Default random from 0 to 2

--lanky : Set the difficulty to the maximum value

--inputs : # of input zones to use; 1-2 for research, 1-3 for production.

--elements : Set of elements to exclusively select atoms from. Can be specified by atomic symbol or
             number - Greek elements can only be specified by number (200-203).

--basic : Select atoms exclusively from the set of 'basic' bond count atoms
          equivalent to --elements H O B C N S Cl Os

e.g. `python level_generator.py --view -d 1 --inputs 1 --elements C O H`

TBA:

 --outputs

 --fusion

 --fission

 --sorting

 --large-output (make the output zone large and scale up the output molecule's size)


## Behaviour

The process for generating a level is as follows:

* Randomly generate # of input zones

* Randomly generate a total # of atoms across all input zones based on difficulty:
  0: 1-3, 1: 4-6, 2: 7-10, 3: 11-16

* Generate input molecules of specified sizes (and from specified elements if given).
  These will be our final input molecules.

* Divide each input molecule by its 'GCD', and combine all the resulting atoms. This sum creates the
  total formula of the inputs in the 'balanced' reaction equation.
  I'm no longer multiplying these reduced values, in order to make it easier to prevent output
  molecules from becoming bloated - so the relative ratios of the two inputs are only randomized
  based on which get GCD reductions if any.

* Randomly choose a # of output zones, ensuring the given atoms will fit within that many zones.

* Randomly allocate the input atoms between the output zones (currently uses a uniform split, i.e.
  no attempt to equalize how many atoms are assigned to each zone), making sure not to give an
  un-constructable set of atoms.

* Divide the atoms allocated to each output zone by their GCD, then multiply them by a uniform
  random # from 1 to either 3, 6, 10, or 16, depending on the 0-3 difficulty level.

* Construct each output molecule from these atoms (randomly drawn without replacement until a
  successful construction is found).

* To reduce bias from the input -> output process, randomly swap inputs with outputs with 50%
  probability.

* Randomly add level features - anything except fusion, fission, or nuclear reactors.
  Production levels have a 1/50 chance to use assemblers/disassemblers instead of regular reactors.

* Convert to and display the mission code (and open coranac viewer if specified).
