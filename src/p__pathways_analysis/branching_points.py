import json
from itertools import compress
from p__data_management import data_update as data
from p__data_management import data_tools as d_tools

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
    
    # Filter the keys based on values lower than t_min
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

    print()
    print('Here are the list of pathways for', species)
    print('Productive pathways:',list_pathways_prod)
    print('destructive pathways:',list_pathways_destroy)
    print('unaffected pathways:',pathways_non_affected)
    print()

    # New list of new pathways
    new_pathways = []
    new_pathways_Dbp = []

    # we checked if one of the list is not empty
    # if so, the species, even if "short-lived" according to t_min is flagged
    # This prevent the case where some species are stil long-term but the user setup a very large t_min that do not math the destruction rate of species.
    # In that case, if there is no available pathways to produce ... nothing to do
    # And the main_loop can go forever
    flagged_species = False
    if list_pathways_prod:
        cond_prod = True
    else:
        cond_prod = False
    if list_pathways_destroy:
        cond_destroy = True
    else:
        cond_destroy = False
    
    print('Condition Prod/Destroy:',cond_prod,cond_destroy)

    if (cond_prod and cond_destroy):
        # We connect each prod pathway to each destroy pathway
        for p_from in list_pathways_prod:
            for p_to in list_pathways_destroy:
                n_from = [n["index"] for n in active_pathways[p_from]["reactions"]]
                n_to = [n["index"] for n in active_pathways[p_to]["reactions"]]
                print('connecting', str(n_from), 'to', str(n_to))
                new_pathways.append(data.connect_two_pathway(active_pathways[p_from],active_pathways[p_to],species))
                print('with rate of:',new_pathways[-1]["rate"])

        # Now we have new_pathways. We need the same pathways that end up explaining the contribution of Delta concentration of branching poing species
        species_dict = d_tools.get_compound_dict(species)
        if species_dict["Delta concentration"] > 0:
            # we have the productive case
            for p_prod in list_pathways_prod:
                n_from = [n["index"] for n in active_pathways[p_prod]["reactions"]]
                print('connecting prod',n_from,' to D[',species,']')
                new_pathways_Dbp.append(data.connect_pathway_to_Dbp(active_pathways[p_prod],species,flag_update='production'))
                print('with rate of:',new_pathways_Dbp[-1]["rate"])
        else:
            # we have the destructive case
            for p_destruct in list_pathways_destroy:
                n_to = [n["index"] for n in active_pathways[p_destruct]["reactions"]]
                print('connecting destr',n_to,' to D[',species,']')
                new_pathways_Dbp.append(data.connect_pathway_to_Dbp(active_pathways[p_destruct],species,flag_update='destruction'))
                print('with rate of:',new_pathways_Dbp[-1]["rate"])

        # We return the unaffected + new pathways . It means that the old productive and destructive pathways of species are deleted from active_p
        return flagged_species, pathways_non_affected + new_pathways + new_pathways_Dbp
    else:
        flagged_species = True

        return flagged_species, active_pathways


def cleaning_slow_pathways(active_pathways:list,deleted_pathways:list,f_min:float):
    # we iterate through active pathways to check if rate < f_min
    list_to_remove = []
    for item in active_pathways:
        # print('item rate:',item["rate"],'f_min:',f_min)
        if item["rate"] < f_min:
            deleted_pathways.append(item)
            # print("item deleted",item)
            list_to_remove.append(item)
    
    for item in list_to_remove:
        active_pathways.remove(item)
    
    # Now we have to update the rates.
    
    return active_pathways,deleted_pathways

