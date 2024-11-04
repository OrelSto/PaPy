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
import re

from ..p__data_management import data_tools as d_tools
from ..p__data_management import global_var
from ..o__cpap_output import output_tools as o_tools

# Function to extract the compound and stoichiometry from a term
def extract_compound_and_stoichiometry(term:str):
    """extract_compound_and_stoichiometry _summary_

    _extended_summary_

    Parameters
    ----------
    term : str
        _description_

    Returns
    -------
    _type_
        _description_
    """

    # Define a regular expression pattern to match stoichiometry numbers
    # stoichiometry_pattern = re.compile(r'(\d*)\s*(\w+)')
    stoichiometry_pattern = re.compile(r'(\d*\s)(\w+)')

    match = stoichiometry_pattern.match(term)
    split = stoichiometry_pattern.split(term)
    if match:
        stoichiometry = float(match.group(1)) if match.group(1) else 1
        compound = match.group(2)
        return {"compound": compound, "stoichiometry": stoichiometry}
    else:
        if split:
            stoichiometry = 1
            compound = split[0]
            return {"compound": compound, "stoichiometry": stoichiometry}
        else:
            exit()


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
    
    # list of string that can be present in a chemical reaction but does not represent an actual chemical species of interest, or with stoichiometry = 0, or hv
    exluded_str =["M","hv","e","e-"]

    # Split the reaction equation into reactants and products
    reactants, products = reaction_equation.split(' => ')
    products, rate = products.split(' R ',1)

    # Split reactants and products into individual compounds
    reactants = reactants.split(' + ')
    products = products.split(' + ')

    # Create a list to store reactant, product and result data
    reactant_data = []
    product_data = []
    result_data = []

    # Process and merge reactant terms
    for term in reactants:
        compound_data = extract_compound_and_stoichiometry(term.strip())
        if compound_data:
            # Check if the compound already exists in reactant_data
            existing_compound = next((c for c in reactant_data if c["compound"] == compound_data["compound"]), False)
            if existing_compound:
                # If the compound exists, add its stoichiometry
                existing_compound["stoichiometry"] += compound_data["stoichiometry"]
            else:
                # If the compound does not exist, add it to reactant_data
                reactant_data.append(compound_data)
    
    # Process and merge products terms
    for term in products:
        compound_data = extract_compound_and_stoichiometry(term.strip())
        if compound_data:
            # Check if the compound already exists in reactant_data
            existing_compound = next((c for c in product_data if c["compound"] == compound_data["compound"]), False)
            if existing_compound:
                # If the compound exists, add its stoichiometry
                existing_compound["stoichiometry"] += compound_data["stoichiometry"]
            else:
                # If the compound does not exist, add it to product_data
                product_data.append(compound_data)
    
    # Process the results
    # First do the reactants
    for reactant in reactant_data:
        result = {}
        result["compound"] = reactant["compound"]
        # Check if the compound already exists in product_data
        exist_in_products = next((c for c in product_data if c["compound"] == reactant["compound"]), False)
        if exist_in_products:
            # If the compound exists, add its stoichiometry
            result["stoichiometry"] = -reactant["stoichiometry"] + exist_in_products["stoichiometry"]
        else:
            # If the compound does not exist, just the reactant stoichiometry
            result["stoichiometry"] = -reactant["stoichiometry"]
        # Appending to the results if this is an actual chemical species
        if not (result["compound"] in exluded_str):
            result_data.append(result)

    # Second do the products
    for product in product_data:
        # Check if the compound already exists in product_data
        exist_in_reactants = next((c for c in reactant_data if c["compound"] == product["compound"]), False)
        if exist_in_reactants:
            # If the compound exists,do nothing
            pass
        else:
            # If it does not exist, appending to the results
            if not (product["compound"] in exluded_str):
                result_data.append(product)

    # Create a JSON structure
    reaction_data = {
        "reactants": reactant_data,
        "products": product_data,
        "initial rate": float(rate),
        "rate": float(rate),
        "deleted rate":0.0,
        "results": result_data,
        "is_pseudo":False
    }

    return reaction_data

def convert_chemical_reaction_file(filename:str):
    """convert_chemical_reaction_file _summary_

    _extended_summary_

    Parameters
    ----------
    filename : str
        _description_
    """

    if global_var.chronicle_writing:
        o_tools.write_line_chronicle('Starting the convertion to JSON of the user chemical reaction system file')
        o_tools.write_line_chronicle('Adding the following reactions')

    # Empty list of all reactions
    chemical_system = []
    
    # Read the content of the input text file line by line
    with open(filename, 'r') as file:
        while line := file.readline():
            reaction_equation = line.rstrip()
            chemical_system.append(format_line(reaction_equation))
            if global_var.chronicle_writing:
                o_tools.write_line_chronicle(o_tools.reaction_to_str(
                    reaction={"index":chemical_system.index(chemical_system[-1]),
                               "multiplicity":1},
                    chem_system_data=chemical_system
                ))

    # Write the JSON data to an output file
    with open('chemical_reaction_system.json', 'w') as output_file:
        json.dump(chemical_system, output_file, indent=2)

    print("Conversion of",filename,"to JSON chemical system format complete.")
    if global_var.chronicle_writing:
        o_tools.write_line_chronicle('Conversion of '+filename+' to JSON chemical system format complete.')
        o_tools.write_line_chronicle('Saved as chemical_reaction_system.json')


def adding_pseudo_reactions():
    if global_var.chronicle_writing:
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('Addition of pseudo-reactions')
    # This routine adds 2 pseudo-reactions (prod and destruction) for each species in the system
    # Opening JSON file
    cs = open('chemical_species.json')
    crs = open('chemical_reaction_system.json')
    # returns JSON object as a dictionary
    chemical_species = json.load(cs)
    chemical_system = json.load(crs)

    chemical_pseudo_reaction_system = []
    for item in chemical_species:
        chemical_pseudo_reaction_system.append(d_tools.format_pseudo_reaction(species=item["name"],flag='prod'))
        chemical_pseudo_reaction_system.append(d_tools.format_pseudo_reaction(species=item["name"],flag='destroy'))
    
    chemical_system = chemical_system + chemical_pseudo_reaction_system
    # Write the JSON data to an output file
    with open('chemical_reaction_system.json', 'w') as output_file:
        json.dump(chemical_system, output_file, indent=2)
    
    if global_var.chronicle_writing:
        o_tools.write_line_chronicle('Addition complete')
