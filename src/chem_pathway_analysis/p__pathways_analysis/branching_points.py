import json
from ..p__data_management import data_update as data
from ..p__data_management import data_tools as d_tools
from ..p__data_management import global_var
from ..o__cpap_output import output_tools as o_tools

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
        if not item["used_as_BP"]:
            tmp_dict[item["name"]]=item["lifetime"]
    
    # Filter the keys based on values lower than t_min
    filtered_keys = [key for key, value in tmp_dict.items() if value < t_min]

    # Sort the filtered keys by their corresponding values
    sorted_keys = sorted(filtered_keys, key=lambda key: tmp_dict[key])

    # return the list of specie with the lowest lifetime
    # return list(sorted(tmp_dict, key=tmp_dict.__getitem__))
    return sorted_keys


def connecting_pathways(active_pathways:list,species:str):
    # Opening JSON file
    crs = open('chemical_reaction_system.json')
    # returns JSON object as a dictionary
    chemical_system = json.load(crs)
    # we connect pathways that are producing species to pathways that are consuming species
    # list of pathways that produce species

    list_pathways_prod,list_pathways_destroy,pathways_non_affected = d_tools.list_connecting_pathways(set_of_pathways=active_pathways,species=species)

    print()
    print('Here are the list of pathways for', species)
    print('Productive pathways:',list_pathways_prod)
    print('destructive pathways:',list_pathways_destroy)
    print('unaffected pathways:',pathways_non_affected)
    print()

    # Writing some stuff in chronicles
    if global_var.chronicle_writing:
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('Here are the list of pathways for: '+species)
        o_tools.write_line_chronicle('Productive pathways : ')
        for i in list_pathways_prod:
            o_tools.write_line_chronicle(o_tools.pathway_to_str(pathway=active_pathways[i],chem_system_data=chemical_system))
            o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('Destructive pathways: ')
        for i in list_pathways_destroy:
            o_tools.write_line_chronicle(o_tools.pathway_to_str(pathway=active_pathways[i],chem_system_data=chemical_system))
            o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('Unaffected pathways : ')
        for p in pathways_non_affected:
            o_tools.write_line_chronicle(o_tools.pathway_to_str(pathway=p,chem_system_data=chemical_system))
            o_tools.write_line_chronicle('\n')

    # New list of new pathways
    new_pathways = []
    new_pathways_Dbp = []

    # In that case, if there is no available pathways to produce ... nothing to do
    # And the main_loop can go forever
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

        if global_var.chronicle_writing:
            o_tools.write_line_chronicle('\n')
            o_tools.write_line_chronicle('Starting to connect Pathways:')
        for p_from in list_pathways_prod:
            for p_to in list_pathways_destroy:
                n_from = [n["index"] for n in active_pathways[p_from]["reactions"]]
                n_to = [n["index"] for n in active_pathways[p_to]["reactions"]]
                print('connecting', str(n_from), 'to', str(n_to))
                # Chronicles
                if global_var.chronicle_writing:
                    o_tools.write_line_chronicle('\n')
                    o_tools.write_line_chronicle('Connecting with '+species+' as branching point:\n'+
                                                o_tools.pathway_to_str(pathway=active_pathways[p_from],chem_system_data=chemical_system)+
                                                '\n to \n'+
                                                o_tools.pathway_to_str(pathway=active_pathways[p_to],chem_system_data=chemical_system))
                    o_tools.write_line_chronicle('\n')
                
                # we append the mew pathway
                new_pathways.append(data.connect_two_pathway(active_pathways[p_from],active_pathways[p_to],species))
                print('with rate of:',new_pathways[-1]["rate"])
                # Chronicles
                if global_var.chronicle_writing:
                    o_tools.write_line_chronicle('leads to:\n'+
                                                o_tools.pathway_to_str(pathway=new_pathways[-1],chem_system_data=chemical_system)+
                                                '\n with rate of: '+'{:0.3e}'.format(new_pathways[-1]["rate"]))


        # Now we have new_pathways. We need the pathways (prod or destr) that end up explaining the contribution of Delta concentration of branching poing species
        species_dict = d_tools.get_compound_dict(species)
        if species_dict["Delta concentration"] != 0.0:
            if species_dict["Delta concentration"] > 0.0:
                # we have the productive case
                for p_prod in list_pathways_prod:
                    n_from = [n["index"] for n in active_pathways[p_prod]["reactions"]]
                    print('connecting prod',n_from,' to D[',species,']')
                    # Chronicles
                    if global_var.chronicle_writing:
                        o_tools.write_line_chronicle('\n')
                        o_tools.write_line_chronicle('Connecting production of '+species+' to D['+species+']:\n'+
                                                    o_tools.pathway_to_str(pathway=active_pathways[p_prod],chem_system_data=chemical_system))
                    
                    new_pathways_Dbp.append(data.connect_pathway_to_Dbp(active_pathways[p_prod],species,flag_update='production'))
                    print('with rate of:',new_pathways_Dbp[-1]["rate"])

                    # Chronicles
                    if global_var.chronicle_writing:
                        o_tools.write_line_chronicle('with rate of: '+'{:0.3e}'.format(new_pathways_Dbp[-1]["rate"]))
            else:
                # we have the destructive case
                for p_destruct in list_pathways_destroy:
                    n_to = [n["index"] for n in active_pathways[p_destruct]["reactions"]]
                    print('connecting destr',n_to,' to D[',species,']')
                    # Chronicles
                    if global_var.chronicle_writing:
                        o_tools.write_line_chronicle('\n')
                        o_tools.write_line_chronicle('Connecting destruction of '+species+' to D['+species+']:\n'+
                                                    o_tools.pathway_to_str(pathway=active_pathways[p_destruct],chem_system_data=chemical_system))
                    new_pathways_Dbp.append(data.connect_pathway_to_Dbp(active_pathways[p_destruct],species,flag_update='destruction'))
                    print('with rate of:',new_pathways_Dbp[-1]["rate"])

                    # Chronicles
                    if global_var.chronicle_writing:
                        o_tools.write_line_chronicle('with rate of: '+'{:0.3e}'.format(new_pathways_Dbp[-1]["rate"]))
        else:
            # Chronicles
            if global_var.chronicle_writing:
                o_tools.write_line_chronicle('\n')
                o_tools.write_line_chronicle('No connection: to D['+species+'] with Delta ='+'{:0.3e}'.format(species_dict["Delta concentration"]))
            
        
        # We return the unaffected + new pathways . It means that the old productive and destructive pathways of species are deleted from active_p
        active_pathways = pathways_non_affected + new_pathways + new_pathways_Dbp
        for p in active_pathways:
            print(p["reactions"])

        return active_pathways
    else:
        # No prod and destr at the same time
        return active_pathways


def cleaning_slow_pathways(active_pathways:list,deleted_pathways:list,f_min:float):
    if global_var.chronicle_writing:
        # Opening JSON file
        crs = open('chemical_reaction_system.json')
        # returns JSON object as a dictionary
        chemical_system = json.load(crs)
    # we iterate through active pathways to check if rate < f_min
    list_to_remove = []
    for item in active_pathways:
        # print('item rate:',item["rate"],'f_min:',f_min)
        if item["rate"] < f_min:
            deleted_pathways.append(item)
            # print("item deleted",item)
            list_to_remove.append(item)
            # Chronicles
            if global_var.chronicle_writing:
                o_tools.write_line_chronicle('\n')
                o_tools.write_line_chronicle('Deleting pathway (too slow):\n'+o_tools.pathway_to_str(pathway=item,chem_system_data=chemical_system))
    
    for item in list_to_remove:
        active_pathways.remove(item)
    
    # Now we have to update the rates.
    
    return active_pathways,deleted_pathways

