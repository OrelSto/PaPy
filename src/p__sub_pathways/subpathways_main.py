import json
import copy
from p__sub_pathways import subpathway_analysis as sub
from p__data_management import data_tools as d_tools

# This routine is set to initialize the subpathway set for each active pathway defined early on
# Steps:
# 1. reading the active_pathways.json
# 2. looping through the pathways and initialze the set of sub-pathways

def main_subpathways(pathways:list,species_done:list):
    active_pathways_data_tmp = copy.deepcopy(pathways)
    ind=0

    for pathway in active_pathways_data_tmp:
        # for each pathway we run the subpathway analysis
        print('We start the initialization for a new pathway of active_pathways')
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
        ind += 1
    
    return pathways


