# Chemical Pathways analysis Package

## Introduction

This is a Chemical Pathway Analysis package for Python.
Based on the algorithm described in Lehmann, 2004.

### Basic install and running example

At the moment, the package is not implement in a way to be installed via ```pip```.
Nonetheless, the user can run the model via command line taking as example the test_cpa.py file.

In the directory where the github repository is installed, run

``` bash
python3 -m test_cpa
```

This should run the model and provide you with results files.

Some information will be displayed in your terminal and some files will be created:

>simple_output.txt

This is the list of relevant pathways in a text format

>chronicles.txt

If the option ```chronicle_writing=True``` in ```run_cp()``` then ```chronicles.txt``` will be created.
The purpose of this file is to write down each steps of the model with much more details.

### test_cpa.py use

The file ```test_cpa.py``` has multiple example with input files already provided.

* Typical O3 example from Lehmann 2004
* Typical O3 example from Lehmann 2004 with P4 deleted
* Example from Chempath module with corrected reaction
* Example from Uranus 1D model results

## Status of the development

This is still an early stage. It's not production ready.

Here's a short list of functionnality:

* [x] Selection of min rate threshold
* [x] Selection of min lifetime of species (discard long-lived species as branching points)
* [x] Subpathways analysis implementation
* [x] Checking the conservation of atoms in the User's chemical reaction system
* [x] Basic visual outputs
* [ ] Security checks of the correct distribution of each reaction rate through pathways/sub-pathways
* [x] Selection of species of interest
* [ ] Packaging (install via ```pip``` or ```conda```)
* [x] Tools for visualization: Pie charts
* [ ] Rate evaluation for steady-state (info)
* [ ] Tools for visualization: Links, Step Progression
* [ ] Check the conservation of atoms in the User Inputs files
* [ ] Flagging the unused reactions
* [ ] Flagging the long-lived species according to the timestep

### Notes

There is still some work to do and the documentation is thin for the moment, hang in there!
