# Spacechem Random Level Generator

## Installs

Built in python 2.7 (.10)

### periodictable

A python package used to get the atomic symbols of all elements

 `pip install periodictable`

## Usage

`python level_generator.py --imlazy`

Prints the mission code for a research level with randomly assigned input zone and output zones,
with molecules sizes probabilistically based on the difficulty.

*Note*: there is a small chance for the tool to fail upon attempting to construct too complex an
        output molecule (e.g. a 16-atom molecule, half of which are Hydrogens), but the probability
        of this is pretty statistically insignificant even at high difficulties.

## Args

--imlazy : Auto-open Cearn's website with the generated code

--difficulty [int] : Int from 0 to 3 approximately controlling the level complexity
                     (basically just molecule sizes, at the moment). Default random from 0 to 3.

--lanky : Set the difficulty to the maximum value

--inputs : # of inputs to use; 1 or 2

--elements : Set of elements to exclusively select atoms from

--basic : Select atoms exclusively from the set of 'basic' bond count atoms
          equivalent to --elements H O B C N S Cl Os

e.g. `python level_generator.py --imlazy -d 1 --inputs 1 --elements C O H`

TBA:

 --outputs

 --production

 --fusion

 --fission

 --sorting

 --large-output (make the output zone large and scale up the output molecule's size)


## Behavior

The process for generating a level is as follows:

* Randomly generate a total # of atoms between both input zones based on difficulty:
  0: 1-3, 1: 4-6, 2: 7-10, 3: 11-16

* Randomly generate # of input zones if not specified - 50% chance of two input zones.
  Allocate the total # of atoms randomly between the zones.

* Generate input molecules of specified sizes (and from specified elements if given).
  These will be our final molecules.

* Divide each molecule by its 'GCD', and combine the resulting input atoms. This sum creates the
  total formula of the inputs in the 'balanced' reaction equation.
  I'm no longer multiplying these reduced values in order to make it easier to prevent output
  molecules from becoming bloated - so the relative ratios of the two inputs are only randomized
  based on which get GCD reductions if any.

* Randomly choose a # of outputs (force 2 outputs if the total set of input atoms was too large and
  couldn't be simplified). Barring the forced double-output, this is another 50% chance for 2 outputs.

* Randomly allocate the input atoms between the output zones (currently uses a uniform split, i.e.
  no attempt to equalize how many atoms are assigned to each zone), making sure not to give an
  un-constructable set of atoms.

* Divide the atoms allocated to each output zone by their GCD, then multiply them by a uniform
  random # from 1 to either 3, 6, 10, or 16, depending on the difficulty level.

* Construct each output molecule from these atoms (randomly drawn without replacement until a
  successful construction is found).

* To reduce bias from the input -> output process, randomly swap inputs with outputs at 50% chance


* Randomly add level features (currently just bonders, sensor, and teleporter)

* Convert to mission code
