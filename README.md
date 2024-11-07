# Chemical Pathways analysis Package

## Introduction

This is a Chemical Pathway Analysis package for Python.
Reproduction of the method described in Lehmann, 2004.

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

If the option ```chronicle_writing=True``` in ```run_cp()``` then _chronicles.txt_ will be created. The purpose of this file is to write down each steps of the model with much more details.

## Status of the development

This is still an early stage. It's not production ready.

Here's a short list of functionnality:

- [x] Selection of min rate threshold
- [x] Selection of min lifetime of species (discard long-lived species as branching points)
- [x] Subpathways analysis implementation
- [x] Basic visual outputs
- [ ] Security checks of the correct distribution of each reaction rate through pathways/sub-pathways
- [x] Selection of species of interest (discarding automatically species that have a longer lifetime)
- [ ] Packaging (install via ```pip``` or ```conda```)
- [ ] Tools for visualization (Pie charts, links, etc)

### Notes

There is still some work to do and the documentation is thin for the moment, hang in there!
