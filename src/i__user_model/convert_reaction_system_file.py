import json
import re

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
- This module is part of the sub-package i__user_model from CPAP

"""

# Function to extract the compound and stoichiometry from a term
def extract_compound_and_stoichiometry(term):
    
    # Define a regular expression pattern to match stoichiometry numbers
    stoichiometry_pattern = re.compile(r'(\d*)\s*(\w+)')

    match = stoichiometry_pattern.match(term)
    if match:
        stoichiometry = int(match.group(1)) if match.group(1) else 1
        compound = match.group(2)
        return {"compound": compound, "stoichiometry": stoichiometry}
    else:
        return None

def format_line(reaction_equation):
    # Split the reaction equation into reactants and products
    reactants, products = reaction_equation.split(' => ')
    products, rate = products.split(' R ',1)

    # Split reactants and products into individual compounds
    reactants = reactants.split(' + ')
    products = products.split(' + ')

    # Create a list to store reactant, product data
    reactant_data = []
    product_data = []

    # Process and merge reactant terms
    for term in reactants:
        compound_data = extract_compound_and_stoichiometry(term.strip())
        if compound_data:
            # Check if the compound already exists in reactant_data
            existing_compound = next((c for c in reactant_data if c["compound"] == compound_data["compound"]), None)
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
            existing_compound = next((c for c in product_data if c["compound"] == compound_data["compound"]), None)
            if existing_compound:
                # If the compound exists, add its stoichiometry
                existing_compound["stoichiometry"] += compound_data["stoichiometry"]
            else:
                # If the compound does not exist, add it to product_data
                product_data.append(compound_data)
    
    # Create a JSON structure
    reaction_data = {
        "reactants": reactant_data,
        "products": product_data,
        "rate": float(rate)
    }
    return reaction_data

def convert_chemical_reaction_file():

    # # Read the content of the input text file
    # with open('user_model_example.txt', 'r') as file:
    #     reaction_equation = file.read()

    # Empty list of all reactions
    chemical_system = []
    
    # Read the content of the input text file line by line
    with open('user_model_example.txt', 'r') as file:
        while line := file.readline():
            reaction_equation = line.rstrip()

            chemical_system.append(format_line(reaction_equation))

    # Write the JSON data to an output file
    with open('chemical_reaction_system.json', 'w') as output_file:
        json.dump(chemical_system, output_file, indent=2)

    print("Conversion to JSON format complete.")

