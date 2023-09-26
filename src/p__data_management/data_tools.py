import json

def find_compound_in_merged_list(listing:list,compound:str):
    # Initialize a variable to store the index
    index_found = []

    # Iterate through the list
    for index, item in enumerate(listing):
        if "compound" in item and item["compound"] == compound:
            index_found.append(index)

    if index_found is []:
        print("No dictionary with 'compound' equals ",compound," was found in the list.")
        exit()

    return index_found

def save_pathways_to_JSON(pathways:list,filename:str):
    # need a doc here?
    # Write the JSON data to an output file
    with open(filename, 'w') as output_file:
        json.dump(pathways, output_file, indent=2)
    print("Pathways saved as",filename)

def D_compound(compound:dict):
    # It's the contribution of prod or consum for the pathways of a specific compound
    p = compound["production rate"]["active pathways"] + compound["production rate"]["deleted pathways"]
    d = compound["destruction rate"]["active pathways"] + compound["destruction rate"]["deleted pathways"]
    return max(p,d)

def get_compound_dict(compound:str):
    # return the dict of a specified compound
    # Opening JSON file
    cs = open('chemical_species.json')

    # returns JSON object as a dictionary
    chemical_species = json.load(cs)

    # iterate and return the correct dict
    for s in chemical_species:
        if s["name"] == compound:
            return s

def format_pseudo_reaction():
    # this is used in sub_pathway analysis.
    # when a pathway is not a steady-state for any branching point species
    # we need to add some pseudo-reaction to enforce it!

    return