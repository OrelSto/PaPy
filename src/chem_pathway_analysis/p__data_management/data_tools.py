import json
from itertools import compress

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


def format_pseudo_reaction(species:str,flag:str):
    # this is used in sub_pathway analysis.
    # when a pathway is not a steady-state for any branching point species
    # we need to add some pseudo-reaction to enforce it!
    if flag == 'destroy':
        return {
        "reactants": [{"compound":species,"stoichiometry":1}],
        "products": [{"compound":"...","stoichiometry":1}],
        "rate": 0.0,
        "deleted rate":0.0,
        "results": [{"compound":species,"stoichiometry":-1},
                    {"compound":"...","stoichiometry":1}],
        "is_pseudo":True
        }
    if flag == 'prod':
        return {
        "reactants": [{"compound":"...","stoichiometry":1}],
        "products": [{"compound":species,"stoichiometry":1}],
        "rate": 0.0,
        "deleted rate":0.0,
        "results": [{"compound":"...","stoichiometry":-1},
                    {"compound":species,"stoichiometry":1}],
        "is_pseudo":True
        }


def format_pseudo_pathway(species:str,multiplicity:int,flag:str,chemical_system:list):
    # we format a pathway containing the pseudo-reaction (prod or destr) of the species provided
    # First, we have to find the reaction
    if flag == 'prod':
        reaction = {
        "reactants": [
        {
            "compound": "...",
            "stoichiometry": 1
        }
        ],
        "products": [
        {
            "compound": species,
            "stoichiometry": 1
        }
        ],
        "rate": 0.0,
        "deleted rate": 0.0,
        "results": [
        {
            "compound": "...",
            "stoichiometry": -1
        },
        {
            "compound": species,
            "stoichiometry": 1
        }
        ],
        "is_pseudo":True
        }
    elif flag == 'destroy':
        reaction = {
        "reactants": [
        {
            "compound": species,
            "stoichiometry": 1
        }
        ],
        "products": [
        {
            "compound": "...",
            "stoichiometry": 1
        }
        ],
        "rate": 0.0,
        "deleted rate": 0.0,
        "results": [
        {
            "compound": species,
            "stoichiometry": -1
        },
        {
            "compound": "...",
            "stoichiometry": 1
        }
        ],
        "is_pseudo":True
        }
    else:
        print('You need to set flag as prod or destroy in data_tools.py')
        exit()
    
    ind = chemical_system.index(reaction)

    return {
        "reactions": [
            {
            "index": ind,
            "multiplicity": multiplicity
            }
        ],
        "branching points": [
            {
            "compound": species,
            "stoichiometry": multiplicity
            }
        ],
        "list branching points used": [],
        "rate":0.0
        }


def list_connecting_pathways(set_of_pathways:list,species:str):
    # we connect pathways that are producing species to pathways that are consuming species
    # list of pathways that produce species
    list_pathways_prod = []
    list_pathways_destroy = []
    mask = []
    for item in set_of_pathways:
        no_compound = True
        for bp in item["branching points"]:
            # checking if species is in the pathway
            if bp["compound"] == species:
                if bp["stoichiometry"] > 0:
                    # adding the index number of the pathway to the list
                    list_pathways_prod.append(set_of_pathways.index(item))
                    no_compound = False
                elif bp["stoichiometry"] < 0:
                    # adding the index number of the pathway to the list
                    list_pathways_destroy.append(set_of_pathways.index(item))
                    no_compound = False
        if no_compound:
            mask.append(True)
        else:
            mask.append(False)

    # Now that we have our lists, we can update set_of_pathways
    # Basically we have the untouched pathways set_of_pathways - (lists)
    # And a new list of merged pathways
    pathways_non_affected = list(compress(set_of_pathways, mask))

    return list_pathways_prod,list_pathways_destroy,pathways_non_affected


def is_there_a_pseudo_reaction(pathway:list,chemical_system_data:list):
    for r in pathway:
        if chemical_system_data[r["index"]]["is_pseudo"]:
            return True
    
    return False


def find_index_pseudo_reaction(species:str,flag:str):
    # We want to find the index of a pseudo_reaction and return it
    # first we open the JSON chemical system file
    cs = open('chemical_reaction_system.json')
    # returns JSON object as a dictionary
    chemical_system_data = json.load(cs)

    # loop over it
    # Iterate through the list
    for index, r in enumerate(chemical_system_data):
        # if it is a pseudo
        if r["is_pseudo"]:
            if flag == "prod":
                if (r["results"][1]["compound"]==species) and (r["results"][1]["stoichiometry"]>0):
                    return index
            elif flag == "destr":
                if (r["results"][0]["compound"]==species) and (r["results"][0]["stoichiometry"]<0):
                    return index
            else :
                print('flag should be prod or destr in find_index_pseudo_reaction')
                exit()
    
    # if we exit the loop, smth is wrong
    print('no pseudo reaction of ',species,' found?')
    print('error in find_index_pseudo_reaction')
    exit()
