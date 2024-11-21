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

from ..p__data_management import global_var
from ..o__cpap_output import output_tools as o_tools

def format_line(reaction_equation:str,reaction_system:list,timestep:float):
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

    # adding all the prod/destruction rates to evaluate the change in concentration of the compound
    D_conc = 0.0
    d_compound = 0.0
    prod = 0.0
    destr = 0.0
    for item in reaction_system:
        for r in item["results"]:
            if compound == r["compound"]:
                if r["stoichiometry"] > 0:
                    prod += r["stoichiometry"] * item["rate"]
                else:
                    destr += -r["stoichiometry"] * item["rate"]

                D_conc += r["stoichiometry"] * item["rate"] * timestep
                d_compound += r["stoichiometry"] * item["rate"]

    # Check for low values
    if (D_conc < 1.0e-16) and (D_conc > -1.0e-16):
        D_conc = 0.0
    if D_conc == 0.0:
        d_compound = 0.0
    
    # Create a JSON structure
    chemical_species_data = {
        "name": compound,
        "concentration": float(concentration),
        "production rate":{
            "initial":prod,
            "active pathways":0.0,
            "deleted pathways":0.0,
            },
        "destruction rate":{
            "initial":destr,
            "active pathways":0.0,
            "deleted pathways":0.0,
            },
        "lifetime":0.0,
        "Delta concentration":D_conc,
        "delta":d_compound,
        "used_as_BP":False
    }

    return chemical_species_data

def convert_concentration_file(filename:str,timestep:float):
    """convert_concentration_file _summary_

    _extended_summary_

    Parameters
    ----------
    filename : str
        _description_
    """
    # Opening JSON file
    reactions = open('chemical_reaction_system.json')

    # returns JSON object as a dictionary
    reaction_system = json.load(reactions)

    # closing the file
    reactions.close()

    # Empty list of all reactions
    chemical_species = []
    
    if global_var.chronicle_writing:
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('Starting the convertion to JSON of the user chemical species concentration file')
        o_tools.write_line_chronicle('Adding the following species')

    # Read the content of the input text file line by line
    with open(filename, 'r') as file:
        while line := file.readline():
            reaction_equation = line.rstrip()
            chemical_species.append(format_line(reaction_equation,reaction_system,timestep))
            if global_var.chronicle_writing:
                o_tools.write_line_chronicle(chemical_species[-1]["name"])

    # Write the JSON data to an output file
    with open('chemical_species.json', 'w') as output_file:
        json.dump(chemical_species, output_file, indent=2)

    if global_var.chronicle_writing:
        o_tools.write_line_chronicle('Conversion of '+filename+' to JSON  chemical species format complete.')
        o_tools.write_line_chronicle('Saved as chemical_species.json')
    
    return chemical_species
