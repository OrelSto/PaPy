import json

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

def connect_two_pathway(pathway_one:dict,pathway_two:dict):
    # self explanatory
    return []