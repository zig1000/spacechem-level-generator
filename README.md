# Spacechem Random Level Generator

Built in python 2.7 (.10)

## Usage

`python level_generator.py`

Default creates a Research level, printing the importable level code to cmd line. See [Args](#args) for options.

## Args

--view : Auto-open Cearn's mission-viewing website with the generated code
         (http://coranac.com/spacechem/mission-viewer?mode=editor).

--print : Pretty-print the level's molecules in addition to the level code.

--production : Create a production level instead of a research level.

--difficulty [int] : Int from 0 to 3 approximately controlling the level complexity
                     (basically just molecule sizes, at the moment).
                     Default random from 0 to 2.

--lanky : Set the difficulty to the maximum value.

--inputs : # of input zones to use; 1-2 for research, 1-3 for production.

--outputs : # of output zones to use; 1-2 for research, 1-3 for production.

--elements : Set of elements to exclusively select atoms from. Can be specified by atomic symbol or
             number - Greek elements can only be specified by number (200-203).

--basic : Select atoms exclusively from the set of 'basic' bond count atoms;
          equivalent to `--elements H O B C N S Cl Os`.

--large_output: Make a large output (8x4) research level.

--fusion: Add a fusion laser to the level and adjust random generation accordingly.

--fission: Add a fission laser to the level and adjust random generation accordingly.

--nuclear: Shortcut for combining `--fission` and `--fusion`.

--symmetric: Use symmetric molecules

--polymer: Use polymer molecules

e.g. `python level_generator.py --view -d 0 --inputs 1 --elements C O H --fission`

TBA:

 --random (random inputs)

 --sorting

 --waste (imbalanced reaction that requires storing waste atoms in the reactor)

## Behaviour

The process for generating a level is roughly as follows:

* Randomly generate level features (e.g. fuser) and # of output zones (if unspecified in args)

* Randomly generate output molecules which are either symmetrical (rotationally, reflexively, etc)
  or polymers (mostly composed of repeating segments).

* Divide each output molecule by its 'GCD', and combine all the resulting atoms. This sum creates
  the total formula of the outputs in the 'balanced' reaction equation.

* Randomly perform inverse nuclear transformations on the balanced output formula

* Calculate a set of input formulas of minimum total size, such that they span the balanced reaction
  formula. Ensure these formulas are valid based on element bond counts and input zone borders

* Construct each input molecule from the generated formulas (currently ugly construction algorithm)

* To reduce bias from the input -> output process, randomly swap inputs with outputs with 50%
  probability, if input/output specifications allow this.

* Randomly add level features (bonders, recycler, etc).
  Production levels have a 1/50 chance to use assemblers/disassemblers instead of regular reactors.

* Display the mission code (and prettyprint molecules or open cearn's mission viewer website, as specified).
