import json
import scipy.optimize as opt
import numpy as np
import copy
from ..p__initialization import init_pathways as p_init
from ..p__data_management import data_tools as d_tools
from ..p__data_management import data_update as data
from ..o__cpap_output import output_tools as o_tools
from ..p__data_management import global_var

# This is where we analyze the pathway to determine if it is a combination of subpathways

def subpathway_analysis(pathway:dict,active_pathways:list,ind:int,list_species_done:list):

    # print('The PATHWAY STUDIED IS:',pathway['reactions'])
    # we open the chemical reaction file to clean pseudo-reaction later
    # Opening JSON file
    cs = open('chemical_reaction_system.json')
    # returns JSON object as a dictionary
    chemical_system_data = json.load(cs)
    # closing file
    cs.close()

    # Just trying to check if we can just work with the species done by being the species listed as BP used
    species_done = pathway['list branching points used']
    # NOOOOPE we need the same ORDER of construction of BP. It is not conserved, at the moment in the connecting process of pathways
    # Basically pathway['list branching points used'] do not "retain" the timescale order of BPs

    if global_var.chronicle_writing:
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('species_done for the subpathway analysis '+' '.join(species_done))

    # we initliaze the set of sub-pathways to the individual reactions present in the pathway
    # set_SP_init is our base to work with the sub-pathways. We'll construct the final set of SP building combination of the initial reactions defined in set_SP_init
    # set_SP_init = sub_pathway_set_init(pathway=pathway,first_specie=species_done[0])
    set_SP_init,pathway = sub_pathway_set_init(pathway=pathway,list_species=species_done)
    # print('we have set_SP_init',set_SP_init)
    # print('The PATHWAY STUDIED AFTER INIT IS:',pathway['reactions'])
    if global_var.chronicle_writing:
        for sp in set_SP_init:
            o_tools.write_line_chronicle('\n')
            o_tools.write_pathway_chronicle(pathway=sp,chem_system_data=chemical_system_data)
            o_tools.write_line_chronicle('\n')

    # We can have a pathway that is a single entry R1 => Delta Sb
    # This means that we have to check for this possibility later on

    # final_set_SP is the final set of elementary sub pathways
    final_set_SP = []

    # We run through the used branching point species and build potential subpathways
    # for s in pathway["list branching points used"]:
    # init_SP is set to true because for the first pass through the loop we need set_SP_init. After, we need the updated final_set_SP.
    init_SP = True
    for s in species_done:
        if global_var.chronicle_writing:
            o_tools.write_line_chronicle('\n')
            o_tools.write_line_chronicle('Building subpathways with BP '+s)
        # set_SP_tmp is the temporary set of pathway, it is equivalent to P^ in Lehmann, 2004
        set_SP_tmp = []
        
        # We look at if we are at the first step of the loop
        if init_SP:
            # print()
            # print('We are at init of SP:')
            # print()
            # print('connecting sub-pathways to',s)
            # We check if any sub-pathway in set_SP has a zero net production of species s
            # At init, there is only singular reaction but this helps "save" reaction where Sb do not appear
            # Meaning, equivalent to 0 net prod and likely containing other Sb used as BP
            set_SP_tmp = checking_zero_net_SP_v2(set_SP=set_SP_init,set_SP_tmp=set_SP_tmp,species=s)
            # Now we are going to connect subpathways inside the set_SP
            set_SP_tmp = connecting_subpathways(set_SP=set_SP_init,set_SP_tmp=set_SP_tmp,species=s,list_species_done=list_species_done)
            # setting to False for next loop step
            # cleaning pseudo reactions
            # print('cleaning pseudo-reactions in set_SP_tmp')
            # set_SP_tmp = data.clean_pathways_of_pseudo_reaction(set_pathways=set_SP_tmp,chemical_system_data=chemical_system_data)
            init_SP = False
            if global_var.chronicle_writing:
                for sp in set_SP_tmp:
                    o_tools.write_line_chronicle('\n')
                    o_tools.write_pathway_chronicle(pathway=sp,chem_system_data=chemical_system_data)
                    o_tools.write_line_chronicle('\n')
        else:
            # print()
            # print('connecting sub-pathways to',s)
            # We check if any sub-pathway in set_SP has a zero net production of species s
            set_SP_tmp = checking_zero_net_SP_v2(set_SP=final_set_SP,set_SP_tmp=set_SP_tmp,species=s)
            # Now we are going to connect subpathways inside the set_SP
            set_SP_tmp = connecting_subpathways(set_SP=final_set_SP,set_SP_tmp=set_SP_tmp,species=s,list_species_done=list_species_done)
            # cleaning pseudo reactions
            # print('cleaning pseudo-reactions in set_SP_tmp')
            # set_SP_tmp = data.clean_pathways_of_pseudo_reaction(set_pathways=set_SP_tmp,chemical_system_data=chemical_system_data)
            if global_var.chronicle_writing:
                for sp in set_SP_tmp:
                    o_tools.write_line_chronicle('\n')
                    o_tools.write_pathway_chronicle(pathway=sp,chem_system_data=chemical_system_data)
                    o_tools.write_line_chronicle('\n')
        
        # We set up final as empty and copy set_SP_tmp
        final_set_SP = []
        for sp in set_SP_tmp:
            # print('adding SP to final SP',sp["reactions"],'for species',s)
            final_set_SP.append(sp)
        
        # This is what we did previously
        # # this loop is to get a list, and not a nested list
        # for sp in set_SP_tmp:
        #     print('adding SP to final SP',sp)
        #     final_set_SP.append(sp)
        # print('We have',len(final_set_SP),'number of SP in final_set_SP')
        # for p in final_set_SP:
        #     print(p["reactions"])

    # Now that we are done, we have final_set_SP that is a collection 
    # of elementary sub pathways. With that collection, we do the ranking part
    # Ranking the subpathways according to a few rules:
    # From best to worst:
    # 1. Present in the active_pathways
    # 1.1 If present, order with best is higher rate
    # if not in the active pathways list:
    # 2. Simplier is better! Less number of reactions is better
    # We have a list index_list_ranked which is the 

    # print('cleaning final_set_SP')
    final_set_SP = data.clean_pathways_of_pseudo_reaction(set_pathways=final_set_SP,chemical_system_data=chemical_system_data)
    # print('cleaning pathway')
    pathway = data.clean_pathway_of_PR(pathway=pathway,chemical_system_data=chemical_system_data)
    # print('pathway after cleaning',pathway)

    # if len(final_set_SP) > 1:
    if len(final_set_SP) > 0:
        # Chronicles
        if global_var.chronicle_writing:
            o_tools.write_line_chronicle('\n')
            o_tools.write_line_chronicle('The elementary subpathways is/are:')
            for elem_sp in final_set_SP:
                o_tools.write_line_chronicle(o_tools.pathway_to_str(pathway=elem_sp,chem_system_data=chemical_system_data))
                o_tools.write_line_chronicle('\n')
        
        index_list_ranked = ranked_list(final_set_SP=final_set_SP,active_pathways=active_pathways,chemical_system_data=chemical_system_data)
        # print('We have the ranked list of index: ',index_list_ranked)

        # After ranking, minimize!
        # scipy.optimize.linprog()
        # c is the rank (so the index of index_list_ranked) attributed to the elementary pathways
        # so c has the length of number_elementary_pathways == len(final_set_SP)
        # c = [(index_list_ranked.index(i)+1)**2 for i in index_list_ranked]
        c = [(index_list_ranked.index(i)+1)**2 for i in index_list_ranked]
        # b_eq is the equality condition of our linear problem
        # b_eq is the array of multiplicities for each chemical reaction present in the pathway
        b_eq = [r["multiplicity"] for r in pathway["reactions"]]
        # b_eq_reaction_ind is the index of reaction in order to build A_eq
        b_eq_reaction_ind = [r["index"] for r in pathway["reactions"]]
        # Building A_eq matrix
        A_eq = []
        # print('looping over reactions:',len(b_eq_reaction_ind))
        # print('looping over Sub-Pathways:',len(index_list_ranked))
        for r in b_eq_reaction_ind:
            A_eq_row = []
            # print('NEW ROW')
            for i in index_list_ranked:
                sum_m = 0
                for r_sp in final_set_SP[i]["reactions"]:
                    # print('comparing',r,'and',r_sp["index"])
                    if r == r_sp["index"]:
                        sum_m += r_sp["multiplicity"]
                    else:
                        # print('weird stuff in A_eq_row')
                        pass
                A_eq_row.append(sum_m)
            A_eq.append(A_eq_row)

        A_eq = np.array(A_eq)
        c = np.array(c)
        b_eq = np.array(b_eq)
        # print('A_eq',A_eq)
        # print('len A_eq',np.shape(A_eq))
        # print('c',c)
        # print('len c',len(c))
        # print('b_eq',b_eq)
        # print('len b_eq',len(b_eq))

        result = opt.linprog(c=c,b_eq=b_eq,A_eq=A_eq,method='highs-ds',options={'disp':False})
        # print()
        # print('THIS IS THE RESULTS AFTER LINPROG')
        # print(result)
        # print()
        
        # we can have moment where the result is bad
        if not result.success:
            print('WE FAILED DURING THE LINPROG')
            print('b_eq_reaction_ind:',b_eq_reaction_ind)
            print()
            print('The pathway is the succession of reactions')
            for r in pathway["reactions"]:
                print(r)
            print('The final_set_SP is')
            for p in final_set_SP:
                print(p["reactions"])
            exit()

        result = list(result.x)

        # Now that we have the result for each sub-pathway, we update the final set of SP
        # by updating the rates!
        # ind_sp is the indice for index_list_ranked where we store the actual subpathway indice
        ind_sp = 0
        save_solo_index = 0
        for res in result:
            # print('res is',res,'in result.x',result,'at index',index_list_ranked[ind_sp])
            # if we have a weight > 0.0
            if res > 0.0:
                # print('updating rate')
                save_solo_index = ind_sp
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
        # print('the sub pathways are actually the pathway itself')
        # print('DO NOTHING')
        # print('Comparing unique SP at index',index_list_ranked[save_solo_index],' with result Linear Prob. >0:', final_set_SP[index_list_ranked[save_solo_index]]["reactions"])
        # print('with the initial pathway', pathway["reactions"])
        flag_continue = False
    elif len(returned_set_SP) == 0:
        # print('the sub pathways list is empty')
        # print('DO NOTHING')
        flag_continue = False
    else:
        # print('returned_set_SP length is',len(returned_set_SP))
        flag_continue = True

    return returned_set_SP,flag_continue

# def sub_pathway_set_init(pathway:dict,first_specie:str):
def sub_pathway_set_init(pathway:dict,list_species:list) -> list:
    # we open the chemical reaction file to retrieve the actual reaction
    # Opening JSON file
    cs = open('chemical_reaction_system.json')
    # returns JSON object as a dictionary
    chemical_system = json.load(cs)
    # closing file
    cs.close()

    # We want to add to pathway, which is the initial pathway studied, the pseudo_react
    for s in list_species:
        # print('adding pseudo-reactions for species ',s,' to the pathway studied for sub-pathways')
        pathway = data.add_pseudo_reaction_to_pathway_to_0NET(pathway=pathway,species=s)

    set_SP = []
    for r in pathway["reactions"]:
        reaction = chemical_system[r["index"]]
        # print('reaction index',r["index"])
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
    

    return set_SP,pathway

def update_subpathway(sub_pathway:dict,reaction:dict,pathway:dict):
    sub_pathway_reaction = sub_pathway["reactions"][0]
    # setting the multiplicity
    # sub_pathway_reaction["multiplicity"] = reaction["multiplicity"]
    # doing stuff with pathway rate
    sub_pathway["rate"] = 0.0


def checking_zero_net_SP_v2(set_SP:list,set_SP_tmp:list,species:str):
    # checking if any pathway in set_SP has a zero net prod of species
    # If so, adding it to set_SP_tmp
    # 28/06/2024:
    # This is v2, I changed things but I do not recall what I did actually
    no_net = True
    for item in set_SP:
        # print('We re looking at SP:',item["reactions"])
        bp_present = False
        for bp in item["branching points"]:
            if bp["compound"] == species:
                bp_present = True
                # print('Specie',species,'is present')
            if bp["compound"] == species and bp["stoichiometry"] == 0:
                # print('We have a zero net prod pathway for',species)
                # print('SP:',item["reactions"])
                no_net = False
                set_SP_tmp.append(item)
            else:
                # we are in the case where species is present but no net
                # we have to do smth, add pseudo-reactions to get 0 NET?
                pass
            
        # if the species is not present in the pathway, then this is equivalent to a 0 net prod
        if not bp_present:
            # print('Specie',species,'is not present')
            # print('This is an equivalent zero net prod pathway for',species)
            # print('SP:',item["reactions"])
            set_SP_tmp.append(item)

    if no_net:
        # print('NO zero net prod for',species)
        pass

    return set_SP_tmp


def checking_zero_net_SP(set_SP:list,set_SP_tmp:list,species:str,init_SP:bool):
    # checking if any pathway in set_SP has a zero net prod of species
    # If so, adding it to set_SP_tmp
    no_net = True
    for item in set_SP:
        bp_present = False
        for bp in item["branching points"]:
            if bp["compound"] == species:
                bp_present = True
            if bp["compound"] == species and bp["stoichiometry"] == 0:
                print('We have a zero net prod pathway for',species)
                print('SP:',item["reactions"])
                no_net = False
                adding_SP(set_SP=set_SP_tmp,final_set_SP=set_SP,pathway_to_be_checked=item)
        # if the species is not present in the pathway, then this is equivalent to a 0 net prod
        if not bp_present:
            print('We have an equivalent zero net prod pathway for',species)
            if init_SP:
                print('adding without check')
                print('SP:',item["reactions"])
                set_SP_tmp.append(item)
            else:
                print('SP:',item["reactions"])
                adding_SP(set_SP=set_SP_tmp,final_set_SP=set_SP,pathway_to_be_checked=item)

    if no_net:
        print('NO net prod for',species)
    return set_SP_tmp


def connecting_subpathways(set_SP:list,set_SP_tmp:list,species:str,list_species_done:list):
    # We connect subpathways as we did for pathways in branching_points.py
    # This is the same abstract idea.
    # final_set_SP = copy.deepcopy(set_SP)
    list_pathways_prod,list_pathways_destroy,_ = d_tools.list_connecting_pathways(set_of_pathways=set_SP,species=species)
    # print()
    # print('We will connect prod',list_pathways_prod,'to',list_pathways_destroy)

    # Checking the prod/destr cond
    if list_pathways_prod:
        cond_prod = True
    else:
        cond_prod = False
    if list_pathways_destroy:
        cond_destroy = True
    else:
        cond_destroy = False
    
    # print('Condition Prod/Destroy:',cond_prod,cond_destroy)

    if (cond_prod and cond_destroy):
        # We connect each prod pathway to each destroy pathway
        for p_from in list_pathways_prod:
            for p_to in list_pathways_destroy:
                n_from = [n["index"] for n in set_SP[p_from]["reactions"]]
                n_to = [n["index"] for n in set_SP[p_to]["reactions"]]
                # print('connecting', str(n_from), 'to', str(n_to))
                # new_SP is at stoichiometry 0, so it fulfills the zero net production of species s condition for a subpathway to be added
                new_SP,_,_ = data.connect_two_pathway(set_SP[p_from], set_SP[p_to],species,list_species_done,True)
                # print('with rate of:',new_SP["rate"])
                # checking if new_SP is already present in set_SP_tmp and if it is an actual 
                # elementary pathway!
                adding_SP(set_SP=set_SP_tmp,final_set_SP=set_SP,pathway_to_be_checked=new_SP)
                # if we added new_SP, we have to delete
    # only prod of Sb
    elif (cond_prod and not cond_destroy):
        print('Only prod ',species)
        for p_from in list_pathways_prod:
            n_from = [n["index"] for n in set_SP[p_from]["reactions"]]
            # print('adding prod pseudo_reaction to', str(n_from))
            # adding the pseudo reaction
            # new_SP = data.add_pseudo_reaction_to_pathway_to_0NET(pathway=set_SP[p_from],species=species)
            new_SP =set_SP[p_from]
            # print('set_SP',set_SP)
            # print('set_SP_tmp',set_SP_tmp)
            # adding_SP(set_SP=set_SP_tmp,final_set_SP=set_SP,pathway_to_be_checked=new_SP)
            set_SP_tmp.append(new_SP)

    # only destr of Sb
    elif (not cond_prod and cond_destroy):
        print('Only destr ',species)
        for p_to in list_pathways_destroy:
            n_to = [n["index"] for n in set_SP[p_to]["reactions"]]
            # print('adding destr pseudo_reaction to', str(n_to))
            # adding the pseudo reaction
            # new_SP = data.add_pseudo_reaction_to_pathway_to_0NET(pathway=set_SP[p_to],species=species)
            new_SP = set_SP[p_to]
            # print('set_SP',set_SP)
            # print('set_SP_tmp',set_SP_tmp)
            # adding_SP(set_SP=set_SP_tmp,final_set_SP=set_SP,pathway_to_be_checked=new_SP)
            set_SP_tmp.append(new_SP)
    
    # No prod 
    # nor destr, 0 NET hence the pathways are already added
    elif (not cond_prod and not cond_destroy):
        print('0 Net ',species)
        pass

    # Are we returning the unaffected pathways ?
    # Not really since normally they are already added from the check_zero_prod before
    # And weirdly enough, the next if does not add anything to set_SP_tmp
    # if pathways_non_affected:
    #     print('We return also the unaffected pathways:')
    #     for p_unaffected in pathways_non_affected:
    #         n_unaffected = [n["index"] for n in p_unaffected["reactions"]]
    #         print(n_unaffected)
    
    return set_SP_tmp


def adding_SP(set_SP:list,final_set_SP:list,pathway_to_be_checked:dict):
    # we check if the pathway_to_be_checked is already in set_SP 
    p_2b_checked_reactions = [r["index"] for r in pathway_to_be_checked["reactions"]]
    p_2b_checked_reactions = sorted(p_2b_checked_reactions)
    r_2b_checked_against = [[i["index"] for i in r["reactions"]] for r in set_SP]
    r_2b_checked_against = sorted(r_2b_checked_against)
    final_r_2b_checked_against = []
    # print('final_set_SP',final_set_SP)
    for item in final_set_SP:
        list_item = []
        for r in item["reactions"]:
            list_item.append(r["index"])
        # list_item = sorted(list_item)
        # print('list_item',sorted(list_item))
        final_r_2b_checked_against.append(sorted(list_item))
    # final_r_2b_checked_against = sorted(final_r_2b_checked_against)

    # final_r_2b_checked_against = [[i["index"] for i in r["reactions"]] for r in final_set_SP]
    # final_r_2b_checked_against = sorted(final_r_2b_checked_against)

    # print('p_2b_checked_reactions',p_2b_checked_reactions)
    # print('r_2b_checked_against',r_2b_checked_against)
    # print('final_r_2b_checked_against',final_r_2b_checked_against)

    # if pathway_to_be_checked in set_SP:
    if (p_2b_checked_reactions in r_2b_checked_against) or (p_2b_checked_reactions in final_r_2b_checked_against):
        # print('not connected: Already present')
        return
    else:
        # pathway_to_be_checked is not already in set_SP
        # check if it is an elementary pathway of set_SP
        # if cond_elementary_pathway(set_SP=set_SP,pathway_to_be_checked=pathway_to_be_checked):
        if cond_elementary_pathway(set_SP=set_SP,pathway_to_be_checked=pathway_to_be_checked):
            # print('connected: All cond checked')
            set_SP.append(pathway_to_be_checked)
        else:
            # print('not connected: subset exist')
            pass
        return


def cond_elementary_pathway(set_SP:list,pathway_to_be_checked:dict):
    # It's here that we check if that pathway_to_be_checked is not a subset of other pathways in set_SP.
    # WARNING: in Lehmann 2004, the condition is erroneous put on set P where it should be on set P^ !! == to set_SP_tmp !!
    # sorted list of reaction numbers
    p_2b_checked_reactions = [r["index"] for r in pathway_to_be_checked["reactions"]]
    p_2b_checked_reactions = sorted(p_2b_checked_reactions,reverse=True)

    for item in set_SP:
        r_2b_checked_against = [r["index"] for r in item["reactions"]]
        r_2b_checked_against = sorted(r_2b_checked_against,reverse=True)
        is_subset = set(r_2b_checked_against).issubset(set(p_2b_checked_reactions))
        if is_subset:
            return False
    # If we never return False ... then good job it is a new elementary subpathway
    return True


def ranked_list(final_set_SP:list,active_pathways:list,chemical_system_data:list):
    # what we do is to create a dictionnary where we put the sub-pathways present in active pathways
    # rank 1 is a dict {Position index in final set SP : Rate of the pathway in active patways}
    rank1_SP = {}
    # rank 2 is a dict {Position index in final set SP : Sum of reactions multiplicities}
    rank2_SP = {}
    # rank 3 is all SP with pseudo-reaction
    rank3_SP = {}

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
            if d_tools.is_there_a_pseudo_reaction(pathway=item,chemical_system_data=chemical_system_data):
                rank3_SP.update({final_set_SP_reactions.index(item):n})
            else:
                rank2_SP.update({final_set_SP_reactions.index(item):n})

    # print('ranking the SP present in active pathways: ',rank1_SP)
    # for r in rank1_SP.keys():
    #     print(final_set_SP_reactions[r])
    # print('sorted list of rank 1',sorted(rank1_SP,key=rank1_SP.get,reverse=True))

    # We sorted in a descending order of rate for rank 1
    list_of_rank1_SP = sorted(rank1_SP,key=rank1_SP.get,reverse=True)
    # We sorted in a ascending order of summed multiplicities for rank 2
    list_of_rank2_SP = sorted(rank2_SP,key=rank2_SP.get)
    # Do we move the SP containing a pseudo_reaction to the end of the list?
    list_of_rank3_SP = sorted(rank3_SP,key=rank3_SP.get)
    # We return the merged list
    list_ranked_merged = list_of_rank1_SP + list_of_rank2_SP + list_of_rank3_SP

    return list_ranked_merged


def save_subpathways_to_JSON(set_SP:list,filename:str):
    # need a doc here?
    # Write the JSON data to an output file
    with open(filename, 'w') as output_file:
        json.dump(set_SP, output_file, indent=2)
    # print("Pathways saved as",filename)