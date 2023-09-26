import json
import numpy as np
import copy

from p__data_management import data_tools as d_tools

def evaluate_production_rate_active_pathways(species:str):
    # Evaluates the production rate of a given species according to the active pathways
    # Opening JSON file
    ap = open('active_pathways.json')

    # returns JSON object as a dictionary
    active_p = json.load(ap)

    # setting destruction rate to 0
    production_rate = 0.0

    # Check the JSON dict for the presence of species
    for item in active_p:
        for bp in item["branching points"]:
            # checking of right specie and actual destruction
            if species == bp["compound"] and bp["stoichiometry"] > 0:
                # production_rate is formatted positive
                production_rate += bp["stoichiometry"] * item["rate"]

    return production_rate

def evaluate_production_rate_deleted_pathways(species:str):
    # Evaluates the production rate of a given species according to the deleted pathways
    # Opening JSON file
    dp = open('deleted_pathways.json')

    # returns JSON object as a dictionary
    deleted_p = json.load(dp)

    # setting destruction rate to 0
    production_rate = 0.0

    # Check the JSON dict for the presence of species
    for item in deleted_p:
        for bp in item["branching points"]:
            # checking of right specie and actual destruction
            if species == bp["compound"] and bp["stoichiometry"] > 0:
                # production_rate is formatted positive
                production_rate += bp["stoichiometry"] * item["rate"]
    
    return production_rate


def evaluate_destruction_rate_active_pathways(species:str):
    # Evaluates the destruction rate of a given species according to the active pathways
    # Opening JSON file
    ap = open('active_pathways.json')
    
    # returns JSON object asa dictionary
    active_p = json.load(ap)

    # setting destruction rate to 0
    destruction_rate = 0.0

    # Check the JSON dict for the presence of species
    for item in active_p:
        for bp in item["branching points"]:
            # checking of right specie and actual destruction
            if species == bp["compound"] and bp["stoichiometry"] < 0:
                # destruction_rate is formatted positive
                destruction_rate += -bp["stoichiometry"] * item["rate"]

    return destruction_rate

def evaluate_destruction_rate_deleted_pathways(species:str):
    # Evaluates the destruction rate of a given species according to deleted pathways
    # Opening JSON file
    dp = open('deleted_pathways.json')
    
    # returns JSON object asa dictionary
    deleted_p = json.load(dp)

    # setting destruction rate to 0
    destruction_rate = 0.0

    # Check the JSON dict for the presence of species
    for item in deleted_p:
        for bp in item["branching points"]:
            # checking of right specie and actual destruction
            if species == bp["compound"] and bp["stoichiometry"] < 0:
                # destruction_rate is formatted positive
                destruction_rate += -bp["stoichiometry"] * item["rate"]
    
    return destruction_rate

def update_rates_chemical_species():
    # Opening JSON file
    cs = open('chemical_species.json')

    # returns JSON object as a dictionary
    chemical_species = json.load(cs)

    # updating the rates
    for item in chemical_species:
        print("We are updating ",item["name"])
        item["production rate"]["active pathways"] = evaluate_production_rate_active_pathways(species=item["name"])
        item["production rate"]["deleted pathways"] = evaluate_production_rate_deleted_pathways(species=item["name"])
        item["destruction rate"]["active pathways"] = evaluate_destruction_rate_active_pathways(species=item["name"])
        item["destruction rate"]["deleted pathways"] = evaluate_destruction_rate_deleted_pathways(species=item["name"])
        
        # updating its lifetime
        item["lifetime"] = item["concentration"] / (item["destruction rate"]["active pathways"] + item["destruction rate"]["deleted pathways"])

    # Now we save it
    # Write the JSON data to an output file
    with open('chemical_species.json', 'w') as output_file:
        json.dump(chemical_species, output_file, indent=2)

    print("Updating rates in chemical_species.json complete.")


def connect_two_pathway(pathway_prod:dict,pathway_destr:dict,species:str):
    # self explanatory ... connecting two pathways
    # I means we merge two pathways that are connected via a species that the one produce and the latter destroys

    # 1. Updating the pathways to a new multiplicity if necessary
    # making deepcopies because we don't want to modify the pathways
    # we want ones with modified multiplicity to construct a new one
    pathway_prod_tmp = copy.deepcopy(pathway_prod)
    pathway_destr_tmp = copy.deepcopy(pathway_destr)
    list_bp_used = pathway_prod_tmp["list branching points used"] + pathway_destr_tmp["list branching points used"]
    list_bp_used.append(species)
    list_bp_used = list(set(list_bp_used))

    bp_stoichiometry,pathway_prod_tmp,pathway_destr_tmp = update_pathway_multiplicity(pathway_prod_tmp,pathway_destr_tmp,species)
    # print('pathway_prod_tmp',pathway_prod_tmp)
    # print('pathway_prod',pathway_prod)

    # 2. We check the deleted pathways effects on branching point species
    pathway_prod_tmp = update_pathway_rate_from_deleted_p(pathway_prod_tmp,species,case_flag=1)
    pathway_destr_tmp = update_pathway_rate_from_deleted_p(pathway_destr_tmp,species,case_flag=2)

    # 2.1 Now that we have the according rates from the deleted pathways
    # We have to update the reaction system file
    # DO WE ??? the actual rates of the chemical reaction are uased as a conservation property
    # DO IT AT THE END OF THE WHILE IN MAIN LOOP
    # OR NOT because we'll loose the f_del_prod,f_del_destr
    # update_chemical_system(f_prod=f_del_prod,p_prod=pathway_prod,f_destr=f_del_destr,p_destr=pathway_destr)

    # 3. Merging the pathways into a new one
    new_pathway = {"reactions":pathway_prod_tmp["reactions"] + pathway_destr_tmp["reactions"],
        "branching points":clean_branching_points(pathway_prod_tmp["branching points"],pathway_destr_tmp["branching points"]),
        "list branching points used":list_bp_used,
        "rate":bp_stoichiometry * (pathway_prod_tmp["rate"] * pathway_destr_tmp["rate"]) / d_tools.D_compound(d_tools.get_compound_dict(species))
    }
    
    # We can clean for redundancy in the reactions
    new_pathway = clean_reactions(new_pathway)
    new_pathway = clean_multiplicity(new_pathway)
    # You might think that cleaning the redundancy of reactions is bad practice since we will look at
    # sub-pathways later on. BUT IT'S NOT !!!
    # Sub-pathways are constructed starting with each individual reactions!
    # it does not start with the actual suit of event of the pathway like:
    # R1 => R3 => R0
    # So, a suit of event like R0 => R3 => R1 => R3 => R7
    # might be deconstructed in R0 => R3 => R7 and R3 <=> R1
    # it does not matters if in the pathway reactions are listed as R0, R1, R3, R7.
    # Because it is the connection between Branching points and Reactions that is important here
    # So, we clean here, better sooner than later where it can be messy ^^

    return new_pathway

def clean_multiplicity(pathway:dict):
    # Now we're looking for a common divisor
    list_multiplicity = []
    for r in pathway["reactions"]:
        list_multiplicity.append(r["multiplicity"])
    common_divisor = np.gcd.reduce(list_multiplicity)

    # Updating the reactions if necessary
    if common_divisor > 1:
        print('There is a common divisor')
        print('previous rate:',pathway["rate"])
        pathway["rate"] = pathway["rate"] * float(common_divisor)
        for r in pathway["reactions"]:
            r["multiplicity"] = r["multiplicity"] / common_divisor
        for bp in pathway["branching points"]:
            bp["stoichiometry"] = bp["stoichiometry"] / common_divisor

    return pathway

def clean_reactions(pathway:dict):
    # we add the same reactions appearing and check for common divisor
    reactions = pathway["reactions"]
    new_reactions = []
    list_ind_r = [r["index"] for r in reactions]
    print('we clean the list of reactions', list_ind_r)

    # Create a dictionary to store indexes of values
    index_dict = {}

    # Iterate through the list and record the indexes of each value
    for idx, value in enumerate(list_ind_r):
        if value in index_dict:
            index_dict[value].append(idx)
        else:
            index_dict[value] = [idx]
    
    # Find values with more than one index (duplicates)
    duplicates = {value: indexes for value, indexes in index_dict.items() if len(indexes) > 1}
    # Find values not duplicates
    not_duplicates = {value: indexes for value, indexes in index_dict.items() if len(indexes) == 1}

    # duplicate reactions
    total_multiplicity = 0
    print('duplicates:',duplicates)
    # Print the indexes of the duplicates
    for value, indexes in duplicates.items():
        print(f"Value {value} appears at indexes: {indexes}")
        for i in indexes:
            total_multiplicity += reactions[i]["multiplicity"]
        new_reactions.append({"index":value,"multiplicity":total_multiplicity})

    # same loop but we remove the old duplicates
    for value, indexes in not_duplicates.items():
        print(f"Value {value} appears at indexes: {indexes}")
        for i in indexes:
            new_reactions.append({"index":value,"multiplicity":reactions[i]["multiplicity"]})
    
    # updating pathway
    pathway["reactions"] = new_reactions
    return pathway


def connect_pathway_to_Dbp(pathway:dict,species:str,flag_update:str):
    # self explanatory ... connecting pathway to the change in species concentration
    # I means we merge two pathways that are connected via a species that the one produce and the latter destroys
    pathway_tmp = copy.deepcopy(pathway)
    species_dict = d_tools.get_compound_dict(species)
    list_bp_used = pathway_tmp["list branching points used"]
    list_bp_used.append(species)
    list_bp_used = list(set(list_bp_used))

    match flag_update:
        case 'production':
            # we remove the destructive rate from deleted pathways because we are producing species
            pathway_tmp = update_pathway_rate_from_deleted_p(pathway_tmp,species,case_flag=1)
        case 'destruction':
            # we remove the productive rate from deleted pathways because we are destroying species
            pathway_tmp = update_pathway_rate_from_deleted_p(pathway_tmp,species,case_flag=2)

    # Updating the pathway into a new one
    new_pathway = {"reactions":pathway_tmp["reactions"] ,
        "branching points":pathway_tmp["branching points"],
        "list branching points used":list_bp_used,
        "rate":(pathway_tmp["rate"] * abs(species_dict["delta"])) / d_tools.D_compound(species_dict)
    }

    return new_pathway

def update_pathway_rate_from_deleted_p(pathway_prod:dict,species:str,case_flag:int):
    # This routine update the rate of pathway_prod according to the destruction rate of species via the deleted pathways
    rate_p_prod = pathway_prod["rate"]
    species_dict = d_tools.get_compound_dict(species)
    rate_deleted_destr_species = species_dict["destruction rate"]["deleted pathways"]
    rate_deleted_prod_species = species_dict["production rate"]["deleted pathways"]
    D_species = d_tools.D_compound(species_dict)

    # associated rate of pathway_prod from the deleted pathways destr species
    f_deleted_destr = (rate_p_prod * rate_deleted_destr_species) / D_species

    # associated rate of pathway_prod from the deleted pathways prod species
    f_deleted_prod = (rate_p_prod * rate_deleted_prod_species) / D_species

    match case_flag:
        case 1:
            # updating the rate
            pathway_prod["rate"] -= f_deleted_destr
        case 2:
            # updating the rate
            pathway_prod["rate"] -= f_deleted_prod
    
    return pathway_prod

def update_pathway_multiplicity(pathway_prod:dict,pathway_destroy:dict,species:dict):
    # do smth with the multipicty to get the lowest denominator

    # First check where we have the branching point species
    bp_tmp = pathway_prod["branching points"] + pathway_destroy["branching points"]
    ind_prod,ind_destroy=d_tools.find_compound_in_merged_list(bp_tmp,species)

    # Second set the branching point species is totaly comsumed, i.e. its stoichiometry = 0
    bp_prod_m = bp_tmp[ind_prod]["stoichiometry"]
    bp_dest_m = -bp_tmp[ind_destroy]["stoichiometry"]

    # We check if the multiplicity are different
    if bp_prod_m != bp_dest_m:
        new_multiplicty = max(bp_prod_m,bp_dest_m)
        print('new multiplicity:', new_multiplicty,'for species:',species)
        print('with:', bp_prod_m,'and :',bp_dest_m)
        
        # Third, updating the pathway prod or pathway destroy
        if bp_prod_m < bp_dest_m:
            print('updating bp_prod_m', bp_prod_m, 'to',new_multiplicty)
            pathway_prod["rate"] = pathway_prod["rate"] / new_multiplicty
            for r in pathway_prod["reactions"]:
                r["multiplicity"] = int(new_multiplicty) * r["multiplicity"]
            for bp in pathway_prod["branching points"]:
                bp["stoichiometry"] = int(new_multiplicty) * bp["stoichiometry"]
        else:
            print('updating bp_dest_m', bp_dest_m, 'to',new_multiplicty)
            pathway_destroy["rate"] = pathway_destroy["rate"] /new_multiplicty
            for r in pathway_destroy["reactions"]:
                r["multiplicity"] = int(new_multiplicty) * r["multiplicity"]
            for bp in pathway_destroy["branching points"]:
                bp["stoichiometry"] = int(new_multiplicty) * bp["stoichiometry"]
    
    # We return the two pathways updated so that merged together the species stoichiometry = 0
    # Last we check again where we have the branching point species with the updated pathways
    bp_tmp = pathway_prod["branching points"] + pathway_destroy["branching points"]
    ind_prod,ind_destroy=d_tools.find_compound_in_merged_list(bp_tmp,species)
    # We return the correct stoichiometry coeff with the pathways
    print('the stoichiometry for',species,'is:',bp_tmp[ind_prod]["stoichiometry"])
    return bp_tmp[ind_prod]["stoichiometry"],pathway_prod,pathway_destroy

def update_chemical_system(f_prod:float,p_prod:dict,f_destr:float,p_destr:dict):
    # We have to update the rate associated to deleted pathways
    return

def clean_branching_points(pathway_one_bp:list,pathway_two_bp:list):
    # init the new_branching points list
    new_branching_points = []
    # We get two list of freshly merged branching points and we clean it
    bp_tmp = pathway_one_bp + pathway_two_bp
    # list of unique chemical species in both pathways
    bp_list = list(set(item["compound"] for item in bp_tmp))
    
    # we iterate on bp_list to merge the branching points
    for species in bp_list:
        ind=d_tools.find_compound_in_merged_list(bp_tmp,compound=species)
        new_st = 0
        for i in ind:
            new_st += bp_tmp[i]["stoichiometry"]
        new_branching_points.append({
            "compound":species,
            "stoichiometry":new_st
        })
    
    return new_branching_points
