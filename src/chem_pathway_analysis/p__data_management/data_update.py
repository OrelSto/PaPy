import json
import numpy as np
import copy

from ..p__data_management import data_tools as d_tools
from ..p__data_management import global_var
from ..o__cpap_output import output_tools as o_tools

def eval_specie_prod_rate_active_p(species:str):
    # Evaluates the production rate of a given species according to the active pathways
    # Opening JSON file
    ap = open('active_pathways.json')

    # returns JSON object as a dictionary
    active_p = json.load(ap)

    # closing file
    ap.close()

    # setting production rate to 0
    production_rate = 0.0

    # Check the JSON dict for the presence of species
    for item in active_p:
        for bp in item["branching points"]:
            # checking of right specie and actual destruction
            if species == bp["compound"] and bp["stoichiometry"] > 0:
                # production_rate is formatted positive
                production_rate += bp["stoichiometry"] * item["rate"]

    return production_rate


def eval_specie_prod_rate_deleted_p(species:str):
    # Evaluates the production rate of a given species according to the deleted pathways
    # Opening JSON file
    dp = open('deleted_pathways.json')

    # returns JSON object as a dictionary
    deleted_p = json.load(dp)

    # closing file
    dp.close()

    # setting production rate to 0
    production_rate = 0.0

    # Check the JSON dict for the presence of species
    for item in deleted_p:
        for bp in item["branching points"]:
            # checking of right specie and actual destruction
            if species == bp["compound"] and bp["stoichiometry"] > 0:
                # production_rate is formatted positive
                production_rate += bp["stoichiometry"] * item["rate"]
    
    return production_rate


def eval_specie_destr_rate_active_p(species:str):
    # Evaluates the destruction rate of a given species according to the active pathways
    # Opening JSON file
    ap = open('active_pathways.json')
    
    # returns JSON object asa dictionary
    active_p = json.load(ap)

    # setting destruction rate to 0
    destruction_rate = 0.0

    # closing file
    ap.close()

    # Check the JSON dict for the presence of species
    for item in active_p:
        for bp in item["branching points"]:
            # checking of right specie and actual destruction
            if species == bp["compound"] and bp["stoichiometry"] < 0:
                # destruction_rate is formatted positive
                destruction_rate += -bp["stoichiometry"] * item["rate"]

    return destruction_rate


def eval_specie_destr_rate_deleted_p(species:str):
    # Evaluates the destruction rate of a given species according to deleted pathways
    # Opening JSON file
    dp = open('deleted_pathways.json')
    
    # returns JSON object asa dictionary
    deleted_p = json.load(dp)

    # closing file
    dp.close()

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


def eval_reaction_rate_deleted_p(reaction_index:int):
    # Evaluates the rate of a given reaction according to the deleted pathways
    # Opening JSON file
    dp = open('deleted_pathways.json')

    # returns JSON object as a dictionary
    deleted_p = json.load(dp)

    # closing file
    dp.close()

    # setting deleted rate to 0
    deleted_rate = 0.0

    # Check the JSON dict for the presence of species
    for item in deleted_p:
        for r in item["reactions"]:
            # checking of right specie and actual destruction
            if reaction_index == r["index"]:
                # deleted_rate is formatted positive
                deleted_rate += item["rate"] * r["multiplicity"]
    
    return deleted_rate


def update_reaction_rate_delete_p(reactions:list,rate_del:float):
    # Update the rate of a given reaction according to its connection of a deleted pathway and an active one
    # Opening JSON file
    crs = open('chemical_reaction_system.json')

    # returns JSON object as a dictionary
    chemical_system = json.load(crs)

    # closing file
    crs.close()

    # for each reaction r in the pathway we update the rates in the reactions file
    for r in reactions:
        reaction = chemical_system[r["index"]]
        reaction["rate"] -= r["multiplicity"]*rate_del
        reaction["deleted rate"] += r["multiplicity"]*rate_del

    # Now we save it
    # Write the JSON data to an output file
    with open('chemical_reaction_system.json', 'w') as output_file:
        json.dump(chemical_system, output_file, indent=2)


def update_rate_chemical_species_delete_p(branching_points:list,rate_del:float):
    # Update the rate of a given species according to its connection of a deleted pathway and an active one
    # We update every species present in the pathway
    for bp in branching_points:
        s = bp["compound"]
        m = bp["stoichiometry"]
        species_dict = d_tools.get_compound_dict(s)

        # we check the stoichiometry of the species result in the pathway
        # if m=0, we don't update obviously
        if m > 0:
            p_rate_dict=species_dict["production rate"]
            p_rate_dict["active pathways"] -= m * rate_del
            p_rate_dict["deleted pathways"] += m * rate_del
        elif m < 0:
            d_rate_dict=species_dict["destruction rate"]
            d_rate_dict["active pathways"] -= -m * rate_del
            d_rate_dict["deleted pathways"] += -m * rate_del


def update_rates_chemical_species(species:str):

    if global_var.chronicle_writing:
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('**************************************')
        o_tools.write_line_chronicle('Updating the rates of chemical species')
        o_tools.write_line_chronicle('\n')

    # Opening JSON file
    cs = open('chemical_species.json')

    # returns JSON object as a dictionary
    chemical_species = json.load(cs)

    # closing file
    cs.close()

    # updating the rates
    for item in chemical_species:
        print("We are updating ",item["name"])
        if global_var.chronicle_writing:
            o_tools.write_line_chronicle('Updating:'+item["name"])
            o_tools.write_line_chronicle('    prod:'+'{:0.3e}'.format(item["production rate"]["active pathways"])+' to '+'{:0.3e}'.format((eval_specie_prod_rate_active_p(species=item["name"]))))
            o_tools.write_line_chronicle('   destr:'+'{:0.3e}'.format(item["destruction rate"]["active pathways"])+' to '+'{:0.3e}'.format((eval_specie_destr_rate_active_p(species=item["name"]))))
        item["production rate"]["active pathways"] = eval_specie_prod_rate_active_p(species=item["name"])
        item["production rate"]["deleted pathways"] = eval_specie_prod_rate_deleted_p(species=item["name"])
        item["destruction rate"]["active pathways"] = eval_specie_destr_rate_active_p(species=item["name"])
        item["destruction rate"]["deleted pathways"] = eval_specie_destr_rate_deleted_p(species=item["name"])
        
        # updating its lifetime
        if global_var.chronicle_writing:
            sav_lt = item["lifetime"]
        
        destruction_rate = (item["destruction rate"]["active pathways"] + item["destruction rate"]["deleted pathways"])
        if destruction_rate == 0.0:
            item["lifetime"] = 1e99
        else:
            item["lifetime"] = item["concentration"] / (item["destruction rate"]["active pathways"] + item["destruction rate"]["deleted pathways"])
        
        if global_var.chronicle_writing:
             o_tools.write_line_chronicle('lifetime:'+'{:0.3e}'.format(sav_lt)+' to '+'{:0.3e}'.format(item["lifetime"]))

    # Now we save it
    # Write the JSON data to an output file
    with open('chemical_species.json', 'w') as output_file:
        json.dump(chemical_species, output_file, indent=2)
    with open('chemical_species_'+species+'.json', 'w') as output_file:
        json.dump(chemical_species, output_file, indent=2)

    print("Updating rates in chemical_species.json complete.")
    if global_var.chronicle_writing:
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('Updating rates in chemical_species.json complete.')
        o_tools.write_line_chronicle('**************************************')


def update_rates_reaction_system():

    if global_var.chronicle_writing:
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('**************************************')
        o_tools.write_line_chronicle('Updating the rates of reaction system')
        o_tools.write_line_chronicle('\n')

    # Opening JSON file
    crs = open('chemical_reaction_system.json')

    # returns JSON object as a dictionary
    chemical_system = json.load(crs)

    # closing file
    crs.close()

    # TODO CHANGE THE LOOP INSIDE CODE TO ADAPT IT TO THE REACTION SYSTEM
    # updating the rates
    for item in chemical_system:
        if not item["is_pseudo"]:
            print("We are updating ",o_tools.reaction_to_str({"index":chemical_system.index(item),
                               "multiplicity":1},chemical_system))
            # updating the rate and deleted rate associated with the reaction
            item["deleted rate"] = eval_reaction_rate_deleted_p(reaction_index=chemical_system.index(item))
            item["rate"] = item["initial rate"] - item["deleted rate"]
            # Writing in chronicles
            if global_var.chronicle_writing:
                o_tools.write_line_chronicle('Updating:'+o_tools.reaction_to_str({"index":chemical_system.index(item),
                               "multiplicity":1},chemical_system))
                o_tools.write_line_chronicle('    rate:'+'{:0.3e}'.format(item["rate"]))
                o_tools.write_line_chronicle('del rate:'+'{:0.3e}'.format(item["deleted rate"]))

    # Now we save it
    # Write the JSON data to an output file
    with open('chemical_reaction_system.json', 'w') as output_file:
        json.dump(chemical_system, output_file, indent=2)

    print("Updating rates in chemical_reaction_system.json complete.")
    if global_var.chronicle_writing:
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('Updating rates in chemical_reaction_system.json complete.')
        o_tools.write_line_chronicle('**************************************')


def connect_two_pathway(pathway_prod:dict,pathway_destr:dict,species:str,list_species_done:list,building_SP:bool):
    # self explanatory ... connecting two pathways
    # I means we merge two pathways that are connected via a species that the one produce and the latter destroys

    # 1. Updating the pathways to a new multiplicity if necessary
    # making deepcopies because we don't want to modify the pathways
    # we want ones with modified multiplicity to construct a new one
    pathway_prod_tmp = copy.deepcopy(pathway_prod)
    pathway_destr_tmp = copy.deepcopy(pathway_destr)
    list_bp_used = pathway_prod_tmp["list branching points used"] + pathway_destr_tmp["list branching points used"]
    list_bp_used.append(species)
    # list_bp_used = list(set(list_bp_used))
    list_bp_used = sorted(set(list_bp_used),key=list_species_done.index)

    bp_stoichiometry,pathway_prod_tmp,pathway_destr_tmp = update_pathway_multiplicity(pathway_prod_tmp,pathway_destr_tmp,species)
    # print('pathway_prod_tmp',pathway_prod_tmp)
    # print('pathway_prod',pathway_prod)

    # 2. We check the deleted pathways effects on branching point species
    pathway_prod_tmp["rate"] = update_pathway_rate_from_deleted_p(pathway_prod_tmp,species,case_flag='prod')
    pathway_destr_tmp["rate"] = update_pathway_rate_from_deleted_p(pathway_destr_tmp,species,case_flag='destr')

    # 3. Merging the pathways into a new one
    # we check 0 values for D_compound
    if d_tools.D_compound(d_tools.get_compound_dict(species)) != 0.0:
        # prod and/or destr != 0, hence rate to be distributed
        print('with BP stoechio :',bp_stoichiometry)
        print('with rate prod pw:',pathway_prod_tmp["rate"])
        print('with rate dest pw:',pathway_destr_tmp["rate"])
        print('with Db          :', d_tools.D_compound(d_tools.get_compound_dict(species)))
        if (global_var.chronicle_writing and not building_SP):
                o_tools.write_line_chronicle('with BP stoechio :'+'{:0.3e}'.format(bp_stoichiometry))
                o_tools.write_line_chronicle('with rate prod pw:'+'{:0.3e}'.format(pathway_prod_tmp["rate"]))
                o_tools.write_line_chronicle('with rate dest pw:'+'{:0.3e}'.format(pathway_destr_tmp["rate"]))
                o_tools.write_line_chronicle('with Db          :'+'{:0.3e}'.format(d_tools.D_compound(d_tools.get_compound_dict(species))))
        
        # Evaluation of the rate of the new pathway
        rate = bp_stoichiometry*(pathway_prod_tmp["rate"] * pathway_destr_tmp["rate"]) / d_tools.D_compound(d_tools.get_compound_dict(species))
        # we multipliply by bp_stoichiometry in order to conserve molecule number
        # with DSb which is for 1.Sb where the constructed pathway is for 
        # bp_stoichiometry.Sb Hence we divide DSb by bp_stoichiometry
    else:
        # it means that prod and destr = 0, hence no rate to distribute
        rate = 0.0
    
    new_pathway = {"reactions":pathway_prod_tmp["reactions"] + pathway_destr_tmp["reactions"],
        "branching points":clean_branching_points(pathway_prod_tmp["branching points"],pathway_destr_tmp["branching points"]),
        "list branching points used":list_bp_used,
        "rate":rate
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

    # 4. We check the deleted pathways if some could have connected, prod and destr
    # those deleted connected pathways must be removed from the final rate
    # of the new pathway
    # new_pathway = update_pathway_rate_from_deleted_p(new_pathway,species,case_flag='new_pathway')

    return new_pathway,pathway_prod_tmp["rate"],pathway_destr_tmp["rate"]


def clean_multiplicity(pathway:dict):
    # Now we're looking for a common divisor
    list_multiplicity = []
    for r in pathway["reactions"]:
        list_multiplicity.append(r["multiplicity"])
    common_divisor = np.gcd.reduce(list_multiplicity)

    # Updating the reactions if necessary
    if common_divisor > 1:
        print('There is a common divisor: ',common_divisor)
        # print('previous rate:',pathway["rate"])
        pathway["rate"] = pathway["rate"] * float(common_divisor)
        for r in pathway["reactions"]:
            r["multiplicity"] = int(r["multiplicity"] / common_divisor)
        for bp in pathway["branching points"]:
            bp["stoichiometry"] = int(bp["stoichiometry"] / common_divisor)
        print('returning: ',pathway)

    return pathway


def clean_reactions(pathway:dict):
    # we add the same reactions appearing
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
    # print('duplicates:',duplicates)
    # Print the indexes of the duplicates
    for value, indexes in duplicates.items():
        # init total_multiplicity for each duplicated reaction
        total_multiplicity = 0
        # print(f"Value {value} appears at indexes: {indexes}")
        for i in indexes:
            total_multiplicity += int(reactions[i]["multiplicity"])
        new_reactions.append({"index":value,"multiplicity":total_multiplicity})

    # same loop but we remove the old duplicates
    for value, indexes in not_duplicates.items():
        # print(f"Value {value} appears at indexes: {indexes}")
        for i in indexes:
            new_reactions.append({"index":value,"multiplicity":int(reactions[i]["multiplicity"])})
    
    # updating pathway
    print('updating reactions to',new_reactions)
    pathway["reactions"] = new_reactions
    return pathway


def connect_pathway_to_Dbp(pathway:dict,species:str,flag_update:str):
    # self explanatory ... connecting pathway to the change in species concentration
    # I means we merge two pathways that are connected via a species that the one produce and the latter destroys
    pathway_tmp = copy.deepcopy(pathway)
    species_dict = d_tools.get_compound_dict(species)
    list_bp_used = pathway_tmp["list branching points used"]
    # we don't add a species since this is not a "full" connection
    # Meaning that there is only a prod/destr part of the connection.
    # species cannot be list as a "branching point"
    # list_bp_used.append(species)
    # list_bp_used = list(set(list_bp_used))

    match flag_update:
        case 'production':
            # we remove the destructive rate from deleted pathways because we are producing species
            # pathway_tmp = update_pathway_rate_from_deleted_p(pathway_tmp,species,case_flag='prod')
            pass
        case 'destruction':
            # we remove the productive rate from deleted pathways because we are destroying species
            # pathway_tmp = update_pathway_rate_from_deleted_p(pathway_tmp,species,case_flag='destr')
            pass

    # Updating the pathway into a new one
    # We don't need bp_stoichiometry here because we have a ratio of dSb/DSb
    # Even if the pathway produce N.Sb, we multiply it by a ratio hence, no problem
    new_pathway = {"reactions":pathway_tmp["reactions"] ,
        "branching points":pathway_tmp["branching points"],
        "list branching points used":list_bp_used,
        "rate":(pathway_tmp["rate"] * abs(species_dict["delta"])) / d_tools.D_compound(species_dict)
    }

    return new_pathway


def update_pathway_rate_from_deleted_p(pathway:dict,species:str,case_flag:int):
    # This routine update the rate of pathway according to the rate of the deleted pathways affecting species
    pathway_tmp = copy.deepcopy(pathway)
    rate_p = pathway_tmp["rate"]
    species_dict = d_tools.get_compound_dict(species)
    rate_deleted_destr_species = species_dict["destruction rate"]["deleted pathways"]
    rate_deleted_prod_species = species_dict["production rate"]["deleted pathways"]
    D_species = d_tools.D_compound(species_dict)
    # index = d_tools.find_compound_in_pathway_BP(pathway_tmp["branching points"],species)
    # stoichiometry = pathway_tmp["branching points"][index]["stoichiometry"]

    # check for 0 value
    if D_species != 0.0 :
        # associated rate of pathway from the deleted pathways destr species
        f_deleted_destr = (rate_p * rate_deleted_destr_species) / D_species
        # f_deleted_destr = stoichiometry * (rate_p * rate_deleted_destr_species) / D_species

        # associated rate of pathway from the deleted pathways prod species
        f_deleted_prod = (rate_p * rate_deleted_prod_species) / D_species
        # f_deleted_prod = stoichiometry * (rate_p * rate_deleted_prod_species) / D_species
        
        # associated rate of pathway from the deleted pathways prod and destr species
        f_deleted_prod_destr = (rate_deleted_prod_species * rate_deleted_prod_species) / D_species
        # f_deleted_prod_destr = stoichiometry * (rate_deleted_destr_species * rate_deleted_prod_species) / D_species



        match case_flag:
            case 'prod':
                if f_deleted_destr > 0.0:
                    # We open deleted_pathways in order to add it some
                    # Opening JSON file
                    dp = open('deleted_pathways.json')

                    # returns JSON object as a list of dict
                    deleted_p = json.load(dp)

                    # closing file
                    dp.close()

                    print('Updating pathway prod rate from connecting with deleted one')
                    print('with rate deleted: ',f_deleted_destr,' from: ',rate_p)
                    # updating the rate
                    pathway_tmp["rate"] -= f_deleted_destr
                    # Now we need to update the deleted rate of reactions that destroyed species and species affected by the pathway
                    # We do that by adding the same pathway in deleted_pathways, but with the are f_deleted
                    if d_tools.is_pathway_in_list(pathway_to_be_checked=pathway_tmp,list_of_pathways=deleted_p):
                        # We have the pathway already in deleted_p
                        # we need to modify its rate then in deleted_p
                        index = d_tools.find_pathway_in_list(pathway_to_be_found=pathway_tmp,list_of_pathways=deleted_p)
                        deleted_p[index]["rate"] += f_deleted_destr
                    else:
                        # pathway is NOT in deleted_p
                        # Hence, we add it with the deleted rate associated
                        deleted_p.append(copy.deepcopy(pathway_tmp))
                        deleted_p[-1]["rate"] = f_deleted_destr
                    print('the new pathway rate is: ',pathway_tmp["rate"])

                    # update_reaction_rate_delete_p(reactions=pathway["reactions"],rate_del=f_deleted_destr)
                    # Now we update the species affected by the pathway
                    # update_rate_chemical_species_delete_p(branching_points=pathway["branching points"],rate_del=f_deleted_destr)
                    # Before we leave we save the json file
                    d_tools.save_pathways_to_JSON(pathways=deleted_p,filename='deleted_pathways.json')
                
            case 'destr':
                if f_deleted_prod > 0.0:
                    # We open deleted_pathways in order to add it some
                    # Opening JSON file
                    dp = open('deleted_pathways.json')

                    # returns JSON object as a list of dict
                    deleted_p = json.load(dp)

                    # closing file
                    dp.close()

                    print('Updating pathway destr rate from connecting with deleted one')
                    print('with rate deleted: ',f_deleted_prod,' from: ',rate_p)
                    # updating the rate
                    pathway_tmp["rate"] -= f_deleted_prod
                    print('the new pathway rate is: ',pathway_tmp["rate"])
                    ## Now we need to update the deleted rate of reactions that destroyed species and species affected by the pathway
                    # We do that by adding the same pathway in deleted_pathways, but with the are f_deleted
                    if d_tools.is_pathway_in_list(pathway_to_be_checked=pathway_tmp,list_of_pathways=deleted_p):
                        # We have the pathway already in deleted_p
                        # we need to modify its rate then in deleted_p
                        index = d_tools.find_pathway_in_list(pathway_to_be_found=pathway_tmp,list_of_pathways=deleted_p)
                        deleted_p[index]["rate"] += f_deleted_prod
                    else:
                        # pathway is NOT in deleted_p
                        # Hence, we add it with the deleted rate associated
                        deleted_p.append(copy.deepcopy(pathway))
                        deleted_p[-1]["rate"] = f_deleted_prod
                    # update_reaction_rate_delete_p(reactions=pathway["reactions"],rate_del=f_deleted_prod)
                    # Now we update the species affected by the pathway
                    # update_rate_chemical_species_delete_p(branching_points=pathway["branching points"],rate_del=f_deleted_prod)
                    # Before we leave we save the json file
                    d_tools.save_pathways_to_JSON(pathways=deleted_p,filename='deleted_pathways.json')
            case 'new_pathway':
                if f_deleted_prod_destr > 0.0:
                    # updating the rate
                    # pathway["rate"] -= f_deleted_prod_destr
                    pass
    
    return pathway_tmp["rate"]


def update_pathway_multiplicity(pathway_prod:dict,pathway_destroy:dict,species:dict):
    # do smth with the multipicty of the branching point to when you connect
    # pathway_prod to pathway_destroy, it gets to 0

    # First: check where we have the branching point species
    bp_tmp = pathway_prod["branching points"] + pathway_destroy["branching points"]
    ind_prod,ind_destroy=d_tools.find_compound_in_merged_list(bp_tmp,species)

    # Second: set the branching point species is totaly comsumed, i.e. its stoichiometry = 0
    bp_prod_m = bp_tmp[ind_prod]["stoichiometry"]
    bp_dest_m = -bp_tmp[ind_destroy]["stoichiometry"]

    # We check if the multiplicity are different
    if bp_prod_m != bp_dest_m:
        # The new multiplicity is the product of the BP stoichiometry from prod and destr pathway
        new_multiplicty = int(bp_prod_m * bp_dest_m)
        m_prod = int(bp_dest_m)
        m_dest = int(bp_prod_m)
        print('new multiplicity:', new_multiplicty,'for species:',species)
        print('with m_prod:', m_prod,'and m_dest:',m_dest)
        
        # Third: updating the pathway prod or pathway destroy
        # print('updating bp_prod_m', bp_prod_m, 'to',new_multiplicty)
        # we divide "rate" by new_multiplicty because it's a multiplicator
        # we need to conserve the molecule flux.
        pathway_prod["rate"] = pathway_prod["rate"] / m_prod
        for r in pathway_prod["reactions"]:
            r["multiplicity"] = int(m_prod * r["multiplicity"])
        for bp in pathway_prod["branching points"]:
            bp["stoichiometry"] = int(m_prod * bp["stoichiometry"])
        # print('updating bp_dest_m', bp_dest_m, 'to',new_multiplicty)
        # we divide "rate" by new_multiplicty because it's a multiplicator
        # we need to conserve the molecule flux.
        pathway_destroy["rate"] = pathway_destroy["rate"] / m_dest
        for r in pathway_destroy["reactions"]:
            r["multiplicity"] = int(m_dest * r["multiplicity"])
        for bp in pathway_destroy["branching points"]:
            bp["stoichiometry"] = int(m_dest * bp["stoichiometry"])
    else:
        new_multiplicty = 1
    
    # If we don't go inside the previous IF, we still return the correct stoechiometry
    # We return the two pathways updated so that merged together the species stoichiometry = 0
    # Last we check again where we have the branching point species with the updated pathways
    bp_tmp = pathway_prod["branching points"] + pathway_destroy["branching points"]
    ind_prod,ind_destroy=d_tools.find_compound_in_merged_list(bp_tmp,species)
    # We return the correct stoichiometry coeff with the pathways
    # print('the stoichiometry prod for',species,'is:',bp_tmp[ind_prod]["stoichiometry"])
    # print('the stoichiometry dest for',species,'is:',bp_tmp[ind_destroy]["stoichiometry"])
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


def clean_pathway_of_PR(pathway:dict,chemical_system_data:list):
    list_r_to_remove = []
    for r in pathway["reactions"]:
        reaction_data = chemical_system_data[r["index"]]
        # If it is a pseudo-reaction CLEAN IT
        if reaction_data["is_pseudo"]:
            print('We are cleaning the pseudo-reaction',r["index"])
            list_r_to_remove.append(r)
            for s in reaction_data["results"]:
                # we re looking in the pseudo-reaction in chemical_system_data
                # for the compound that is associated to it
                if s["compound"] != '...':
                    multiplicity = r["multiplicity"]
                    compound = s["compound"]
                    stoichiometry =s["stoichiometry"]
                    # updating the branching points in the pathway
                    for bp in pathway["branching points"]:
                        # finding the BP == to the compound of the pseudo reaction
                        if bp["compound"] == compound:
                            bp["stoichiometry"] += - stoichiometry * multiplicity
                            # # if the BP has a Null Net, hence it does not appear as a BP
                            # if bp["stoichiometry"] == 0:
                            #     pathway["branching points"].remove(bp)
    
    for r in list_r_to_remove:
        pathway["reactions"].remove(r)
    
    for bp in pathway["branching points"]:
        if bp["compound"] == '...':
            pathway["branching points"].remove(bp)

    return pathway


def clean_pathways_of_pseudo_reaction(set_pathways:list,chemical_system_data:list):
    # We need to check if there is a reaction that is a pseudo-reaction.
    # If so, remove it and update the stoichiometry

    # print('set_pathways',set_pathways)
    # print()
    for pathway in set_pathways:
        # print('pathway in set_pathways',pathway)
        pathway = clean_pathway_of_PR(pathway=pathway,chemical_system_data=chemical_system_data)

    return set_pathways


def add_pseudo_reaction_to_pathway_to_0NET(pathway:dict,species:str):
    # we want to add a pseudo_reaction to an existing pathway

    # we want the stoich of the species in the pathway
    i_comp = d_tools.find_compound_in_merged_list(listing=pathway["branching points"],compound=species)
    # we want the stoechiometry of the compond to define if it needs a prod or destr to balance it out.
    mult = pathway["branching points"][i_comp[0]]["stoichiometry"]

    # we define the flag
    if mult > 0:
        # looking for a destr to balance a prod
        flag = "destr"
    elif mult < 0:
        # looking for a prod to balance a destr
        flag = "prod"
    else:
        flag = 0
        # print('mult == 0 in add_pseudo_reaction_to_pathway_to_0NET')
        print('mult == 0')
        print('No addition of pseudo-reactions for ',species)
        # print('NOT POSSIBLE')
        # exit()

    if flag != 0:
        print('Adding a ',flag,' pseudo-reaction for specie ',species)
        # Then, get the index of the pseudo reaction to add
        index = d_tools.find_index_pseudo_reaction(species=species,flag=flag)

        # Now we add an entry to the reaction index
        new_react = {
            "index": index,
            "multiplicity": abs(mult)
        }
        pathway["reactions"].append(new_react)

        # Updating the species stoichiometry
        pathway["branching points"][i_comp[0]]["stoichiometry"] = 0

    return pathway


def update_active_p_rates_from_deleted_p(active_pathways:list,deleted_pathways:list):
    """update_active_p_rates_from_deleted_p Updating active_p rates from same pathways present in deleted_p

    _extended_summary_

    Parameters
    ----------
    active_pathways : list
        _description_
    deleted_pathways : list
        _description_
    """
    for pathway in active_pathways:
        if d_tools.is_pathway_in_list(pathway_to_be_checked=pathway,list_of_pathways=deleted_pathways):
            index = d_tools.find_pathway_in_list(pathway_to_be_found=pathway,list_of_pathways=deleted_pathways)
            # We have the rate that is deleted
            p_del_rate = deleted_pathways[index]["rate"]
            # Now we update the rate
            pathway["rate"] -= p_del_rate

    return active_pathways

