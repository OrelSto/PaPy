import json

from ..o__cpap_output import output_tools as o_tools
from ..p__data_management import global_var

# we make a moche output ^^'
def moche(target_specie:str) -> None:
    with open('active_pathways.json', 'r') as active_pathways_file:
        # Parse the JSON data and store it in a variable
        active_pathways_data = json.load(active_pathways_file)

    with open('deleted_pathways.json', 'r') as deleted_pathways_file:
        # Parse the JSON data and store it in a variable
        deleted_pathways_data = json.load(deleted_pathways_file)

    with open('chemical_reaction_system.json', 'r') as chem_system_file:
        # Parse the JSON data and store it in a variable
        chem_system_data = json.load(chem_system_file)

    with open('output_moche.txt', 'w') as output_moche_file:
        output_moche_file.write('**************************')
        output_moche_file.write('\n')

        rate_sum = 0.0
        for pathway in active_pathways_data:
            rate_sum += pathway["rate"]
        # we re gonna sort the pathways in order of importance
        pathway_sorted = {}
        i = 0
        for pathway in active_pathways_data:
            pathway_sorted.update({i:pathway["rate"]/rate_sum * 100})
            i += 1
        ind_pathway_sorted = sorted(pathway_sorted,key=pathway_sorted.get,reverse=True)

        rate_deleted = 0.0
        for pathway in deleted_pathways_data:
            rate_sum += pathway["rate"]
            rate_deleted += pathway["rate"]
        
        # for pathway in active_pathways_data:
        for i in ind_pathway_sorted:
            moche_writing_pathway(pathway=active_pathways_data[i],output_moche_file=output_moche_file,chem_system_data=chem_system_data,rate_sum=rate_sum)
        
        # Now the rate from deleted pathways
        output_moche_file.write(' RATE DELETED  : ' + '{:0.3e}'.format(rate_deleted))
        output_moche_file.write(' \n')
        output_moche_file.write(' RATE DELETED %: ' + '{:0.3f}'.format(rate_deleted/rate_sum * 100))
        output_moche_file.write(' \n')

        # Now we call the output function for a target specie
        if target_specie != 'None':
            moche_target_specie_output(target_specie)

def moche_writing_pathway(pathway:list,output_moche_file,chem_system_data,rate_sum:float) -> None:

    o_tools.writing_pathway(pathway=pathway,output_file=output_moche_file,chem_system_data=chem_system_data)

    output_moche_file.write(' \n')
    
    # Now the rate
    output_moche_file.write(' \n')
    # output_moche_file.write(' RATE  : ' + str(round(pathway["rate"],3)))
    output_moche_file.write(' RATE  : ' + '{:0.3e}'.format(pathway["rate"]))
    output_moche_file.write(' \n')
    output_moche_file.write(' RATE %: ' + '{:0.3f}'.format(pathway["rate"]/rate_sum * 100))
    output_moche_file.write(' \n')
    

    output_moche_file.write(' \n')
    output_moche_file.write('**************************')
    output_moche_file.write('\n')


def moche_target_specie_output(target_specie:str) -> None:
    # The idea is that the user wants a specific output for a specified chemical species target_specie
    # We ll go through every active pathways and check if target_specie is present and list them.
    # Then we express their rate in a ratio over the prod/destr rate of the target_specie
    with open('active_pathways.json', 'r') as active_pathways_file:
        # Parse the JSON data and store it in a variable
        active_pathways_data = json.load(active_pathways_file)

    with open('deleted_pathways.json', 'r') as deleted_pathways_file:
        # Parse the JSON data and store it in a variable
        deleted_pathways_data = json.load(deleted_pathways_file)

    with open('chemical_reaction_system.json', 'r') as chem_system_file:
        # Parse the JSON data and store it in a variable
        chem_system_data = json.load(chem_system_file)
    
    with open('output_moche_'+target_specie+'.txt', 'w') as output_moche_file:
        output_moche_file.write('**************************')
        output_moche_file.write('\n')

        rate_sum = 0.0
        for pathway in active_pathways_data:
            rate_sum += pathway["rate"]
        # we re gonna sort the pathways in order of importance
        pathway_sorted = {}
        i = 0
        for pathway in active_pathways_data:
            pathway_sorted.update({i:pathway["rate"]/rate_sum * 100})
            i += 1
        ind_pathway_sorted = sorted(pathway_sorted,key=pathway_sorted.get,reverse=True)

        rate_deleted = 0.0
        for pathway in deleted_pathways_data:
            rate_sum += pathway["rate"]
            rate_deleted += pathway["rate"]
        
        # for pathway in active_pathways_data:
        for i in ind_pathway_sorted:
            moche_writing_pathway(pathway=active_pathways_data[i],output_moche_file=output_moche_file,chem_system_data=chem_system_data,rate_sum=rate_sum)
        
        # Now the rate from deleted pathways
        output_moche_file.write(' RATE DELETED  : ' + '{:0.3e}'.format(rate_deleted))
        output_moche_file.write(' \n')
        output_moche_file.write(' RATE DELETED %: ' + '{:0.3f}'.format(rate_deleted/rate_sum * 100))
        output_moche_file.write(' \n')

    pass
