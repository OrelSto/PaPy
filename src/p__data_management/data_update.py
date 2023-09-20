import json
import numpy as np

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


def connect_two_pathway(pathway_one:dict,pathway_two:dict,species:str):
    # self explanatory ... connecting two pathways
    # I means we merge two pathways that are connected via a species that the one produce and the latter destroys

    # First Updating the pathways to a new multiplicity if necessary
    pathway_one_tmp,pathway_two_tmp = update_pathway_multiplicity(pathway_one,pathway_two,species)

    # Second Merging the pathways into a new one
    new_pathway = {"reactions":pathway_one_tmp["reactions"] + pathway_two_tmp["reactions"],
        "branching points":clean_branching_points(pathway_one_tmp["branching points"],pathway_two_tmp["branching points"]),
        "pathway multiplicity":1,
        "rate":(pathway_one_tmp["rate"] * pathway_two_tmp["rate"]) / d_tools.D_compound(d_tools.get_compound_dict(species))
    }

    return new_pathway

def connect_pathway_to_Dbp(pathway:dict,species:str):
    # self explanatory ... connecting pathway to the change in species concentration
    # I means we merge two pathways that are connected via a species that the one produce and the latter destroys

    species_dict = d_tools.get_compound_dict(species)

    # Updating the pathway into a new one
    new_pathway = {"reactions":pathway["reactions"] ,
        "branching points":pathway["branching points"],
        "pathway multiplicity":1,
        "rate":(pathway["rate"] * abs(species_dict["delta"])) / d_tools.D_compound(species_dict)
    }

    return new_pathway

def update_pathway_multiplicity(pathway_prod:dict,pathway_destroy:dict,species:dict):
    # do smth with the multipicty to get the lowest denominator

    # First check where we have the branching point species
    bp_tmp = pathway_prod["branching points"] + pathway_destroy["branching points"]
    ind_prod,ind_destroy=d_tools.find_compound_in_merged_list(bp_tmp,species)

    # Second set the new pathway multiplicity as the branching point species is totaaly comsumed, i.e. its stoichiometry = 0
    bp_prod_m = bp_tmp[ind_prod]["stoichiometry"]
    bp_dest_m = -bp_tmp[ind_destroy]["stoichiometry"]
    # We check if the multiplicity are different, if not just return
    if bp_prod_m != bp_dest_m:
        new_multiplicty = max(bp_prod_m,bp_dest_m)

        print('new multiplicity:', new_multiplicty,'for species:',species)
        print('with:', bp_prod_m,'and :',bp_dest_m)
        
        # Third, updating the pathway prod or pathway destroy
        if bp_prod_m < bp_dest_m:
            pathway_prod["pathway multiplicity"] = 1
            pathway_prod["rate"] = new_multiplicty * pathway_prod["rate"]
            for r in pathway_prod["reactions"]:
                r["multiplicity"] = int(new_multiplicty) * r["multiplicity"]
            for bp in pathway_prod["branching points"]:
                bp["stoichiometry"] = int(new_multiplicty) * bp["stoichiometry"]
        else:
            pathway_destroy["pathway multiplicity"] = 1
            pathway_destroy["rate"] = new_multiplicty * pathway_destroy["rate"]
            for r in pathway_destroy["reactions"]:
                r["multiplicity"] = int(new_multiplicty) * r["multiplicity"]
            for bp in pathway_destroy["branching points"]:
                bp["stoichiometry"] = int(new_multiplicty) * bp["stoichiometry"]
        
    # Now we're looking for a common divisor
    list_multiplicity = []
    for r in pathway_prod["reactions"]:
        list_multiplicity.append(r["multiplicity"])
    for r in pathway_destroy["reactions"]:
        list_multiplicity.append(r["multiplicity"])
    common_divisor = np.gcd.reduce(list_multiplicity)

    # Updating the reactions if necessary
    if common_divisor > 1:
        pathway_prod["rate"] = pathway_prod["rate"] / float(common_divisor)
        pathway_destroy["rate"] = pathway_destroy["rate"] / float(common_divisor)
        for r in pathway_prod["reactions"]:
            r["multiplicity"] = r["multiplicity"] / common_divisor
        for bp in pathway_prod["branching points"]:
            bp["stoichiometry"] = bp["stoichiometry"] / common_divisor
        for r in pathway_destroy["reactions"]:
            r["multiplicity"] = r["multiplicity"] / common_divisor
        for bp in pathway_destroy["branching points"]:
            bp["stoichiometry"] = bp["stoichiometry"] / common_divisor

    # We return the two pathways updated so that merged together the species stoichiometry = 0
    return pathway_prod,pathway_destroy


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