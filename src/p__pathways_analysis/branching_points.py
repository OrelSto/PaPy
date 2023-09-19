import json
from itertools import compress
from p__data_management import data_update as data

def list_next_branching_points(t_min:float):
    # here we determine the list of species by their lifetime
    # loading the chemical species
    # Opening JSON file
    cs = open('chemical_species.json')

    # returns JSON object as a dictionary
    chemical_species = json.load(cs)

    # empty dict
    tmp_dict = {}

    # updating the rates
    for item in chemical_species:
        tmp_dict[item["name"]]=item["lifetime"]
    
    # Filter the keys based on values lower than 2.5
    filtered_keys = [key for key, value in tmp_dict.items() if value < t_min]

    # Sort the filtered keys by their corresponding values
    sorted_keys = sorted(filtered_keys, key=lambda key: tmp_dict[key])

    # return the list of specie with the lowest lifetime
    # return list(sorted(tmp_dict, key=tmp_dict.__getitem__))
    return sorted_keys

def connecting_pathways(active_pathways:list,species:str):
    # we connect pathways that are producing species to pathways that are consuming species
    # list of pathways that produce species
    list_pathways_prod = []
    list_pathways_destroy = []
    mask = []
    for item in active_pathways:
        no_compound = True
        for bp in item["branching points"]:
            # checking if species is in the pathway
            if bp["compound"] == species:
                if bp["stoichiometry"] > 0:
                    # adding the index number of the pathway to the list
                    list_pathways_prod.append(active_pathways.index(item))
                    no_compound = False
                elif bp["stoichiometry"] < 0:
                    # adding the index number of the pathway to the list
                    list_pathways_destroy.append(active_pathways.index(item))
                    no_compound = False
        if no_compound:
            mask.append(True)
        else:
            mask.append(False)

    # Now that we have our lists, we can update active_pathways
    # Basically we have the untouched pathways active_pathways - (lists)
    # And a new list of merged pathways
    pathways_non_affected = list(compress(active_pathways, mask))

    print(pathways_non_affected)

    # New list of new pathways
    new_pathways = []
    # We connect each prod to each destroy
    for p_from in list_pathways_prod:
        for p_to in list_pathways_destroy:
            print('connecting', p_from, 'to', p_to)
            new_pathways.append(data.connect_two_pathway(active_pathways[p_from],active_pathways[p_to]))

    return

