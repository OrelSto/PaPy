import json

from p__initialization import init_pathways as p_init
from p__data_management import data_tools as d_tools
from p__data_management import data_update as data

# This is where we analyze the pathway to determine if it is a combination of subpathways

def subpathway_analysis(pathway:dict,active_pathways:list):

    # This is the "main loop" for the subpathway analysis.
    # init of the final set of Sub Pathways as an empty list
    final_set_SP = []
    # initialization of rank

    # We run through the used branching point species and build potential subpathways
    for s in pathway["list branching points used"]:
        # we initliaze the set of sub-pathways to the individual reactions present in the pathway
        # set_SP is our base to work with the sub-pathways. We'll construct the final set of SP building combination of the initial reactions defined in set_SP
        set_SP = sub_pathway_set_init(pathway=pathway,species=s)
        
        # Now we are going to connect subpathways inside the set_SP
        connecting_subpathways(set_SP=set_SP,final_set_SP=final_set_SP,species=s)

        # Ranking the subpathways according to a few rules:
        # From best to worst:
        # 1. Present in the active_pathways
        # 1.1 If present, order with best is higher rate
        # if not in the active pathways list:
        # 2. Simplier is better! Less number of reactions is better
        # What we need to do is reorganize

        # After ranking, minimize!


    # saving
    save_subpathways_to_JSON(set_SP=final_set_SP,filename='subpathways_tmp.json')
    return

def sub_pathway_set_init(pathway:dict,species:str):
    # we open the chemical reaction file to retrieve the actual reaction
    # Opening JSON file
    cs = open('chemical_reaction_system.json')
    # returns JSON object as a dictionary
    chemical_system = json.load(cs)

    set_SP = []
    for r in pathway["reactions"]:
        reaction = chemical_system[r["index"]]
        print('reaction index',r["index"])
        # we need to record each individual reaction into a new sub-pathway set
        sub_pathway = {}
        sub_pathway = p_init.format_first_pathway(reaction=reaction,index=r["index"])

        # Now we have the sub-pathway as an individual reaction, it is the set P of Lehmann
        # # its rate
        # # its multiplicity
        # # watch out here, reaction = r and not variable reaction ! We want the multiplicity of the reaction r in pathway !!!
        # we are just setting the rate to 0
        update_subpathway(sub_pathway=sub_pathway,reaction=r,pathway=pathway)
        set_SP.append(sub_pathway)
    
    # We need to check if the pathway is not at steady-state for all former branching point species
    index = d_tools.find_compound_in_merged_list(listing=pathway["branching points"],compound=species)[0]
    stoichiometry_s = pathway["branching points"][index]["stoichiometry"]
    if stoichiometry_s < 0:
        print('adding a pseudo reaction that produce',s)
        pseudo_rection_pathway = d_tools.format_pseudo_reaction(species=species,multiplicity=-stoichiometry_s)
        # adding the pseudo reaction
        set_SP.append(pseudo_rection_pathway)
    elif stoichiometry_s > 0:
        print('adding a pseudo reaction that destroy',species)
        pseudo_rection_pathway = d_tools.format_pseudo_reaction(species=species,multiplicity=-stoichiometry_s)
        # adding the pseudo reaction
        set_SP.append(pseudo_rection_pathway)

    return set_SP

def update_subpathway(sub_pathway:dict,reaction:dict,pathway:dict):
    sub_pathway_reaction = sub_pathway["reactions"][0]
    # setting the multiplicity
    # sub_pathway_reaction["multiplicity"] = reaction["multiplicity"]
    # doing stuff with pathway rate
    sub_pathway["rate"] = 0.0

def connecting_subpathways(set_SP:list,final_set_SP:list,species:str):
    # We connect subpathways as we did for pathways in branching_points.py
    # This is the same abstract idea.
    list_pathways_prod,list_pathways_destroy,_ = d_tools.list_connecting_pathways(set_of_pathways=set_SP,species=species)
    # We connect each prod pathway to each destroy pathway
    for p_from in list_pathways_prod:
        for p_to in list_pathways_destroy:
            n_from = [n["index"] for n in set_SP[p_from]["reactions"]]
            n_to = [n["index"] for n in set_SP[p_to]["reactions"]]
            print('connecting', str(n_from), 'to', str(n_to))
            # new_SP is at stoichiometry 0, so i fulfills the zero net production of species s condition for a subpathway to be added
            new_SP = data.connect_two_pathway(set_SP[p_from], set_SP[p_to],species)
            print('with rate of:',new_SP["rate"])
            # checking if new_SP is already present in final set of SP
            adding_SP(final_set_SP=final_set_SP,pathway_to_be_checked=new_SP)


    return final_set_SP

def adding_SP(final_set_SP:list,pathway_to_be_checked:dict):
    # we check if the pathway_to_be_checked is already in final_set_SP 
    if pathway_to_be_checked in final_set_SP:
        return
    else:
        final_set_SP.append(pathway_to_be_checked)
        return

def save_subpathways_to_JSON(set_SP:list,filename:str):
    # need a doc here?
    # Write the JSON data to an output file
    with open(filename, 'w') as output_file:
        json.dump(set_SP, output_file, indent=2)
    print("Pathways saved as",filename)