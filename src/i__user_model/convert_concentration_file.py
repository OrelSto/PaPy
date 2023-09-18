"""
Module Name
convert_reaction_system_file

A brief description of what this module does or provides.

Detailed Description:
- Provide additional details about the module's functionality.
- Mention any key classes, functions, or variables defined in the module.
- Explain the module's role in the larger project or system.

Usage:
- Describe how to import this module.
- Provide examples of how to use the module's features or classes.
- Mention any common use cases.

Dependencies:
- List any external libraries or modules that this module depends on.
- Include version requirements if necessary.

Author:
- Your name or the name of the module's author.

License:
- Specify the module's licensing information if applicable.

Note:
- Include any important notes, warnings, or considerations.
- This module is part of the sub-package i__user_model from CPA

"""

import json

def format_line(reaction_equation:str):
    """format_line _summary_

    _extended_summary_

    Parameters
    ----------
    reaction_equation : str
        _description_

    Returns
    -------
    _type_
        _description_
    """
    # Split the reaction equation into compound and concentration
    compound, concentration = reaction_equation.split(' = ')

    # Create a JSON structure
    chemical_species_data = {
        "name": compound,
        "concentration": float(concentration),
    }

    return chemical_species_data

def convert_concentration_file(filename:str):
    """convert_concentration_file _summary_

    _extended_summary_

    Parameters
    ----------
    filename : str
        _description_
    """
    # Empty list of all reactions
    chemical_system = []
    
    # Read the content of the input text file line by line
    with open(filename, 'r') as file:
        while line := file.readline():
            reaction_equation = line.rstrip()

            chemical_system.append(format_line(reaction_equation))

    # Write the JSON data to an output file
    with open('chemical_species.json', 'w') as output_file:
        json.dump(chemical_system, output_file, indent=2)

    print("Conversion of",filename,"to JSON chemical species format complete.")