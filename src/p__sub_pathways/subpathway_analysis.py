import json
import scipy.optimize as opt
import numpy as np
import copy
from p__initialization import init_pathways as p_init
from p__data_management import data_tools as d_tools
from p__data_management import data_update as data

# This is where we analyze the pathway to determine if it is a combination of subpathways

def subpathway_analysis(pathway:dict,active_pathways:list,ind:int,species_done:list):

    # This is the "main loop" for the subpathway analysis.
    # initialization of rank

    # we initliaze the set of sub-pathways to the individual reactions present in the pathway
    # set_SP is our base to work with the sub-pathways. We'll construct the final set of SP building combination of the initial reactions defined in set_SP
    final_set_SP = sub_pathway_set_init(pathway=pathway)
    print('we have final_set_SP init ',final_set_SP)

    # We run through the used branching point species and build potential subpathways
    # for s in pathway["list branching points used"]:
    for s in species_done:
        set_SP_tmp = []
        print('connecting sub-pathways to',s)

        # We check if any sub-pathway in set_SP has a zero net production of species s
        set_SP_tmp = checking_zero_net_SP(set_SP=final_set_SP,species=s)
        print('we have set_SP_tmp after ckeck net zero',set_SP_tmp)
        
        # Now we are going to connect subpathways inside the set_SP
        set_SP = connecting_subpathways(set_SP=set_SP,species=s)
        print('we have set_SP after connecting SP',set_SP)

    final_set_SP = copy.deepcopy(set_SP)
    print('we have final_set_SP after connecting SP',final_set_SP)
    # Now that we are done, we have final_set_SP that is a collection 
    # of elementary sub pathways. With that collection, we do the ranking part
    # Ranking the subpathways according to a few rules:
    # From best to worst:
    # 1. Present in the active_pathways
    # 1.1 If present, order with best is higher rate
    # if not in the active pathways list:
    # 2. Simplier is better! Less number of reactions is better
    # We have a list index_list_ranked which is the 
    index_list_ranked = ranked_list(final_set_SP=final_set_SP,active_pathways=active_pathways)
    print('We have the ranked list of index: ',index_list_ranked)

    # After ranking, minimize!
    # scipy.optimize.linprog()
    # c is the rank (so the index of index_list_ranked) attributed to the elementary pathways
    # so c has the length of number_elementary_pathways == len(final_set_SP)
    c = [(index_list_ranked.index(i)+1)**2 for i in index_list_ranked]
    # b_eq is the equality condition of our linear problem
    # b_eq is the array of multiplicities for each chemical reaction present in the pathway
    b_eq = [r["multiplicity"] for r in pathway["reactions"]]
    # b_eq_reaction_ind is the index of reaction in order to build A_eq
    b_eq_reaction_ind = [r["index"] for r in pathway["reactions"]]
    # Building A_eq matrix
    A_eq = []
    print('looping over reactions:',len(b_eq_reaction_ind))
    print('looping over Sub-Pathways:',len(index_list_ranked))
    for r in b_eq_reaction_ind:
        A_eq_row = []
        for i in index_list_ranked:
            sum_m = 0
            for r_sp in final_set_SP[i]["reactions"]:
                if r == r_sp["index"]:
                    sum_m += r_sp["multiplicity"]
                else:
                    pass
            A_eq_row.append(sum_m)
        A_eq.append(A_eq_row)

    A_eq = np.array(A_eq)
    c = np.array(c)
    b_eq = np.array(b_eq)
    print('A_eq',A_eq)
    print('len A_eq',np.shape(A_eq))
    print('c',c)
    print('len c',len(c))
    print('b_eq',b_eq)
    print('len b_eq',len(b_eq))

    result = opt.linprog(c=c,b_eq=b_eq,A_eq=A_eq)
    result = list(result.x)
    print()
    print('THIS IS THE RESULTS AFTER LINPROG')
    print(result)
    print()

    # Now that we have the result for each sub-pathway, we update the final set of SP
    # by updating the rates!
    # ind_sp is the indice for index_list_ranked where we store the actual subpathway indice
    ind_sp = 0
    for res in result:
        print('res is',res,'in result.x',result,'at index',index_list_ranked[ind_sp])
        # if we have a weight > 0.0
        if res > 0.0:
            final_set_SP[index_list_ranked[ind_sp]]["rate"] = pathway["rate"] * res
        ind_sp += 1
    
    # saving
    save_subpathways_to_JSON(set_SP=final_set_SP,filename='subpathways_tmp_'+str(ind)+'.json')

    # All the Sub-Pathways are updated. Now we can return the final_set of sub pathways by only sending the SP with rates > 0
    returned_set_SP = []
    for item in final_set_SP:
        if item["rate"] > 0:
            # we keep track of the species used as BP from the MAIN MAIN loop
            item["list branching points used"] = pathway["list branching points used"]
            returned_set_SP.append(item)
    
    if len(returned_set_SP) == 1:
        print('the sub pathways are actually the pathway itself')
        print('DO NOTHING')
        flag_continue = False
    else:
        flag_continue = True

    return returned_set_SP,flag_continue

def sub_pathway_set_init(pathway:dict):
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
    for species in pathway["list branching points used"]:
        index = d_tools.find_compound_in_merged_list(listing=pathway["branching points"],compound=species)[0]
        stoichiometry_s = pathway["branching points"][index]["stoichiometry"]
        if stoichiometry_s < 0:
            print('adding a pseudo reaction that produce',species)
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

def checking_zero_net_SP(set_SP:list,species:str):
    # checking if any pathway in set_SP has a zero net prod of species
    # If so, adding it to final_set_SP
    final_set_SP = []
    for item in set_SP:
        for bp in item["branching points"]:
            if bp["compound"] == species and bp["stoichiometry"] == 0:
                adding_SP(final_set_SP=final_set_SP,pathway_to_be_checked=item)
    return set_SP

def connecting_subpathways(set_SP:list,species:str):
    # We connect subpathways as we did for pathways in branching_points.py
    # This is the same abstract idea.
    # final_set_SP = copy.deepcopy(set_SP)
    new_pathways = []
    list_pathways_prod,list_pathways_destroy,pathways_non_affected = d_tools.list_connecting_pathways(set_of_pathways=set_SP,species=species)
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
            adding_SP(final_set_SP=set_SP,pathway_to_be_checked=new_SP)
            # if we added new_SP, we have to delete 

    return pathways_non_affected + 

def adding_SP(final_set_SP:list,pathway_to_be_checked:dict):
    # we check if the pathway_to_be_checked is already in final_set_SP 
    if pathway_to_be_checked in final_set_SP:
        return
    else:
        final_set_SP.append(pathway_to_be_checked)
        return

def ranked_list(final_set_SP:list,active_pathways:list):
    # what we do is to create a dictionnary where we put the sub-pathways present in active pathways
    # rank 1 is a dict {Position index in final set SP : Rate of the pathway in active patways}
    rank1_SP = {}
    # rank 2 is a dict {Position index in final set SP : Sum of reactions multiplicities}
    rank2_SP = {}

    # we check the reaction list
    active_pathways_reactions = [i["reactions"] for i in active_pathways]
    final_set_SP_reactions = [i["reactions"] for i in final_set_SP]

    for item in final_set_SP_reactions:
        # This is rank 1
        if item in active_pathways_reactions:
            ind_active_pathway = active_pathways_reactions.index(item)
            rank1_SP.update({final_set_SP_reactions.index(item):active_pathways[ind_active_pathway]["rate"]})
        # This is rank 2
        else:
            n = 0
            for r in item:
                n += r["multiplicity"]
            rank2_SP.update({final_set_SP_reactions.index(item):n})

    # print('ranking the SP present in active pathways: ',rank1_SP)
    # for r in rank1_SP.keys():
    #     print(final_set_SP_reactions[r])
    # print('sorted list of rank 1',sorted(rank1_SP,key=rank1_SP.get,reverse=True))

    # We sorted in a descending order of rate for rank 1
    list_of_rank1_SP = sorted(rank1_SP,key=rank1_SP.get,reverse=True)
    # We sorted in a ascending order of summed multiplicities for rank 2
    list_of_rank2_SP = sorted(rank2_SP,key=rank2_SP.get)
    # We return the merged list
    return list_of_rank1_SP + list_of_rank2_SP


def save_subpathways_to_JSON(set_SP:list,filename:str):
    # need a doc here?
    # Write the JSON data to an output file
    with open(filename, 'w') as output_file:
        json.dump(set_SP, output_file, indent=2)
    print("Pathways saved as",filename)