import json
import copy
from ..p__sub_pathways import subpathway_analysis as sub
from ..p__data_management import data_tools as d_tools
from ..o__cpap_output import output_tools as o_tools
from ..p__data_management import global_var

# This routine is set to initialize the subpathway set for each active pathway defined early on
# Steps:
# 1. reading the active_pathways.json
# 2. looping through the pathways and initialze the set of sub-pathways

def main_subpathways(pathways:list,species_done:list):
    # Opening JSON file
    crs = open('chemical_reaction_system.json')
    # returns JSON object as a dictionary
    chemical_system = json.load(crs)

    active_pathways_data_tmp = copy.deepcopy(pathways)
    ind=0

    print('We are working with these active pathways')
    for item in pathways:
        print(item["reactions"])
    
    print()

    for pathway in active_pathways_data_tmp:
        # for each pathway we run the subpathway analysis if the length of the pathway is > 1
        # Obviously, a single reaction does not have subpathways!
        # And a pathway with two reactions also!
        if global_var.chronicle_writing:
            o_tools.write_line_chronicle('\n')
            o_tools.write_line_chronicle('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
            o_tools.write_line_chronicle('We are looking at the pathway')
            o_tools.write_line_chronicle(o_tools.pathway_to_str(pathway=pathway,chem_system_data=chemical_system))
        if len(pathway["reactions"]) > 2:
            print()
            print('We start the initialization for a new pathway of active_pathways')

            if global_var.chronicle_writing:
                o_tools.write_line_chronicle('The pathway has at least than 3 reactions. Looking for subpathways!')

            returned_set_SP,flag_update = sub.subpathway_analysis(pathway=pathway,active_pathways=active_pathways_data_tmp,ind=ind,species_done=species_done)
            # with that returned_set_SP, we have to clean the active pathways data.
            # First we remove the actual pathway since it is no longer needed
            # we check the list of reactions and not the entire item because "rate" key is changing !
            if flag_update:
                list_reactions_active_pathways = [item["reactions"] for item in pathways]
                print('We are deleting',pathway["reactions"],'from activ pathways')
                del pathways[list_reactions_active_pathways.index(pathway["reactions"])]
                # Then,
                # adding the rate if some sub-pathways are already in active_pathways
                # And adding the sub-pathway if it is not
                # new list_reactions_active_pathways because we deleted some entry
                list_reactions_active_pathways = [item["reactions"] for item in pathways]
                list_reactions_SP = [item["reactions"] for item in returned_set_SP]
                # Chronicles
                if global_var.chronicle_writing:
                    o_tools.write_line_chronicle('\n')
                    o_tools.write_line_chronicle('We have '+str(len(returned_set_SP))+' sub-pathways:')
                    for sub_p in returned_set_SP:
                        o_tools.write_line_chronicle('\n')
                        o_tools.write_line_chronicle(o_tools.pathway_to_str(pathway=sub_p,chem_system_data=chemical_system))
                    o_tools.write_line_chronicle('\n')
                    o_tools.write_line_chronicle('explaining the pathway:')
                    o_tools.write_line_chronicle(o_tools.pathway_to_str(pathway=pathway,chem_system_data=chemical_system))
                    o_tools.write_line_chronicle('\n')
                    o_tools.write_line_chronicle('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')

                print('and we are adding/updating:')
                for sub_p in list_reactions_SP:
                    print(sub_p)
                for sub_p in list_reactions_SP:
                    i_sp = list_reactions_SP.index(sub_p)
                    if sub_p in list_reactions_active_pathways:
                        print('upadting',sub_p,'already present')
                        i_ap = list_reactions_active_pathways.index(sub_p)
                        print('from rate',pathways[i_ap]["rate"],' to',pathways[i_ap]["rate"] +returned_set_SP[i_sp]["rate"] )
                        pathways[i_ap]["rate"] += returned_set_SP[i_sp]["rate"]
                    else:
                        print('adding',sub_p)
                        pathways.append(returned_set_SP[i_sp])

            else:# Chronicles
                if global_var.chronicle_writing:
                    o_tools.write_line_chronicle('\n')
                    o_tools.write_line_chronicle('No subpathways found')
                    o_tools.write_line_chronicle('\n')
                    o_tools.write_line_chronicle('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')

        else:
            print('One or Two reactions in pathway => NO SUBPATHWAYS')
            if global_var.chronicle_writing:
                o_tools.write_line_chronicle('The pathway has less than 3 reactions. No need to look at subpathways')
                o_tools.write_line_chronicle('\n')
                o_tools.write_line_chronicle('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
        # incrementing the indice
        ind += 1
    
    return pathways


