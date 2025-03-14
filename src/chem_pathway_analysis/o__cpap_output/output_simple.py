import json
import matplotlib.pyplot as plt

from . import output_tools as o_tools
from ..p__data_management import global_var
from ..p__data_management import data_tools as d_tools

# we make a moche output ^^'
def text_output(target_species:list) -> None:
    with open('active_pathways.json', 'r') as active_pathways_file:
        # Parse the JSON data and store it in a variable
        active_pathways_data = json.load(active_pathways_file)

    with open('deleted_pathways.json', 'r') as deleted_pathways_file:
        # Parse the JSON data and store it in a variable
        deleted_pathways_data = json.load(deleted_pathways_file)

    with open('chemical_reaction_system.json', 'r') as chem_system_file:
        # Parse the JSON data and store it in a variable
        chem_system_data = json.load(chem_system_file)

    with open('simple_output.txt', 'w') as simple_output_file:
        simple_output_file.write('**************************')
        simple_output_file.write('\n')

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
            simple_output_writing_pathway(pathway=active_pathways_data[i],simple_output_file=simple_output_file,chem_system_data=chem_system_data,rate_sum=rate_sum)
        
        # Now the rate from deleted pathways
        simple_output_file.write(' RATE DELETED  : ' + '{:0.3e}'.format(rate_deleted))
        simple_output_file.write(' \n')
        simple_output_file.write(' RATE DELETED %: ' + '{:0.3f}'.format(rate_deleted/rate_sum * 100))
        simple_output_file.write(' \n')

        # Now we call the output function for a target specie
        if target_species != ['None']:
            for s in target_species:
                print('Writing outputs for branching point species: ',s)
                target_species_output(s)
                print('\n')


def simple_output_writing_pathway(pathway:list,simple_output_file,chem_system_data,rate_sum:float,stoich_coeff=1.0) -> None:
    # stoich_coeff is used when you want an output for a target specie

    o_tools.writing_pathway(pathway=pathway,output_file=simple_output_file,chem_system_data=chem_system_data)

    simple_output_file.write(' \n')
    
    # Now the rate
    simple_output_file.write(' \n')
    # simple_output_file.write(' RATE  : ' + str(round(pathway["rate"],3)))
    simple_output_file.write(' RATE  : ' + '{:0.3e}'.format(pathway["rate"]*stoich_coeff))
    simple_output_file.write(' \n')
    simple_output_file.write(' RATE %: ' + '{:0.3f}'.format(pathway["rate"]*stoich_coeff/rate_sum * 100))
    simple_output_file.write(' \n')
    

    simple_output_file.write(' \n')
    simple_output_file.write('**************************')
    simple_output_file.write('\n')


def target_species_output(target_specie:str) -> None:
    # The idea is that the user wants a specific output for a specified chemical species target_specie
    # We ll go through every active pathways and check if target_specie is present and list them.
    # Then we express their rate in a ratio over the prod/destr rate of the target_specie
    with open('active_pathways.json', 'r') as active_pathways_file:
        # Parse the JSON data and store it in a variable
        active_pathways_data = json.load(active_pathways_file)
    # We set up the list of pathways acting on target_specie
    act_pathways_data_t_specie = []
    for pathway in active_pathways_data:
        # if (target_specie in pathway["list branching points used"]) and (d_tools.find_compound_in_merged_list(pathway["branching points"],target_specie)):
        # species does not have to be in "list branching points used"
        if (d_tools.find_compound_in_merged_list(pathway["branching points"],target_specie)):
            ind = d_tools.find_compound_in_merged_list(pathway["branching points"],target_specie)[0]
            # Then we check if this is a pathway with prod or destr of target specie
            if pathway["branching points"][ind]["stoichiometry"] != 0:
                act_pathways_data_t_specie.append(pathway)
        # for species in pathway["list branching points used"]:
        #     list_ind = d_tools.find_compound_in_merged_list(listing=pathway["branching points"],compound=species)
        # if list_ind:
        #     act_pathways_data_t_specie.append(pathway)

    with open('deleted_pathways.json', 'r') as deleted_pathways_file:
        # Parse the JSON data and store it in a variable
        deleted_pathways_data = json.load(deleted_pathways_file)
    # We set up the list of pathways acting on target_specie
    del_pathways_data_t_specie = []
    for pathway in deleted_pathways_data:
        if (d_tools.find_compound_in_merged_list(pathway["branching points"],target_specie)):
            ind = d_tools.find_compound_in_merged_list(pathway["branching points"],target_specie)[0]
            # Then we check if this is a pathway with prod or destr of target specie
            if pathway["branching points"][ind]["stoichiometry"] != 0:
                del_pathways_data_t_specie.append(pathway)

    with open('chemical_reaction_system.json', 'r') as chem_system_file:
        # Parse the JSON data and store it in a variable
        chem_system_data = json.load(chem_system_file)
    
    t_species = d_tools.get_compound_dict(compound=target_specie)
    
    with open('simple_output_'+target_specie+'.txt', 'w') as simple_output_file:

        # The rate_sum for a species is the sum of its prod and destr
        rate_sum = t_species["production rate"]["active pathways"] + t_species["production rate"]["deleted pathways"] + t_species["destruction rate"]["active pathways"] + t_species["destruction rate"]["deleted pathways"]
        rate_sum_prod = t_species["production rate"]["active pathways"] + t_species["production rate"]["deleted pathways"] 
        rate_sum_dest = t_species["destruction rate"]["active pathways"] + t_species["destruction rate"]["deleted pathways"]
        # if rate_sum == 0.0 there is no prod or destr of target_specie
        if rate_sum == 0.0:
            simple_output_file.write('No prod/destr of '+target_specie)
            simple_output_file.write('\n')
            simple_output_file.write('rate_sum = 0')
            simple_output_file.write('\n')
            simple_output_file.write('\n')
            # To avoid /0.0 error
            rate_sum = 1.0
        elif rate_sum_prod > rate_sum_dest:
            simple_output_file.write('Production of '+target_specie)
            simple_output_file.write('\n')
            simple_output_file.write('rate_sum   = '+ '{:0.3e}'.format(rate_sum))
            simple_output_file.write('\n')
            simple_output_file.write('rate_prod  = '+ '{:0.3e}'.format(rate_sum_prod))
            simple_output_file.write('\n')
            simple_output_file.write('rate_destr = '+ '{:0.3e}'.format(rate_sum_dest))
            simple_output_file.write('\n')
            simple_output_file.write('\n')
        else:
            simple_output_file.write('Destruction of '+target_specie)
            simple_output_file.write('\n')
            simple_output_file.write('rate_sum   = '+ '{:0.3e}'.format(abs(rate_sum)))
            simple_output_file.write('\n')
            simple_output_file.write('rate_prod  = '+ '{:0.3e}'.format(abs(rate_sum_prod)))
            simple_output_file.write('\n')
            simple_output_file.write('rate_destr = '+ '{:0.3e}'.format(abs(rate_sum_dest)))
            simple_output_file.write('\n')
            simple_output_file.write('\n')

        simple_output_file.write('**************************')
        simple_output_file.write('\n')

        # we are working with absolute values
        rate_sum = abs(rate_sum)

        # we re gonna sort the pathways in order of importance
        pathway_sorted = {}
        i = 0
        for pathway in act_pathways_data_t_specie:
            stoich = abs(pathway["rate"]*pathway["branching points"][d_tools.find_compound_in_merged_list(listing=pathway["branching points"],compound=target_specie)[0]]["stoichiometry"])/rate_sum * 100
            # We check that the pathway account for more that 0.0001% of the total rate of the species
            # Meaning we explain the last 0.001% with pathways
            if stoich >= 0.0001:
                print('stoich',stoich)
                pathway_sorted.update({i:stoich})
                i += 1
            # If not, then we we will not print it
            else:
                print('Not printing pathway with stoich',stoich)
                # we advance the indice also
                i += 1
        
        # Now we sort the indices
        ind_pathway_sorted = sorted(pathway_sorted,key=pathway_sorted.get,reverse=True)

        rate_deleted = 0.0
        for pathway in del_pathways_data_t_specie:
            rate_deleted += abs(pathway["rate"] *pathway["branching points"][d_tools.find_compound_in_merged_list(listing=pathway["branching points"],compound=target_specie)[0]]["stoichiometry"])
        
        # for pathway in active_pathways_data:
        for ind in ind_pathway_sorted:
            pathway = act_pathways_data_t_specie[ind]
            simple_output_writing_pathway(pathway=act_pathways_data_t_specie[ind],simple_output_file=simple_output_file,chem_system_data=chem_system_data,rate_sum=rate_sum,stoich_coeff=abs(act_pathways_data_t_specie[ind]["branching points"][d_tools.find_compound_in_merged_list(listing=pathway["branching points"],compound=target_specie)[0]]["stoichiometry"]))
        
        # Now the rate from deleted pathways
        simple_output_file.write(' RATE DELETED  : ' + '{:0.3e}'.format(rate_deleted))
        simple_output_file.write(' \n')
        simple_output_file.write(' RATE DELETED %: ' + '{:0.3f}'.format(rate_deleted/rate_sum * 100))
        simple_output_file.write(' \n')


def pie_output(target_species:str,act_P_json:str,del_P_json:str,chem_R_json:str):
    """pie_output _summary_

    _extended_summary_

    Parameters
    ----------
    target_species : str
        _description_
    """
    from matplotlib import rcParams
    rcParams['text.usetex'] = True

    for s in target_species:
        # Before everything we init the text to print the pathways
        text = ''
        # The idea is that the user wants a specific output for a specified chemical species target_specie
        # We ll go through every active pathways and check if target_specie is present and list them.
        # Then we express their rate in a ratio over the prod/destr rate of the target_specie
        with open(act_P_json, 'r') as active_pathways_file:
            # Parse the JSON data and store it in a variable
            active_pathways_data = json.load(active_pathways_file)
        # We set up the list of pathways acting on target_specie
        act_pathways_prod_data_t_specie = []
        act_pathways_dest_data_t_specie = []
        for pathway in active_pathways_data:
            # if (target_specie in pathway["list branching points used"]) and (d_tools.find_compound_in_merged_list(pathway["branching points"],target_specie)):
            # species does not have to be in "list branching points used"
            if (d_tools.find_compound_in_merged_list(pathway["branching points"],s)):
                ind = d_tools.find_compound_in_merged_list(pathway["branching points"],s)[0]
                # Then we check if this is a pathway with prod or destr of target specie
                k = pathway["branching points"][ind]["stoichiometry"]
                if k > 0:
                    act_pathways_prod_data_t_specie.append(pathway)
                elif k < 0:
                    act_pathways_dest_data_t_specie.append(pathway)

        with open(del_P_json, 'r') as deleted_pathways_file:
            # Parse the JSON data and store it in a variable
            deleted_pathways_data = json.load(deleted_pathways_file)
        # We set up the list of pathways acting on s
        del_pathways_prod_data_t_specie = []
        del_pathways_dest_data_t_specie = []
        for pathway in deleted_pathways_data:
            if (d_tools.find_compound_in_merged_list(pathway["branching points"],s)):
                ind = d_tools.find_compound_in_merged_list(pathway["branching points"],s)[0]
                # Then we check if this is a pathway with prod or destr of target specie
                k = pathway["branching points"][ind]["stoichiometry"]
                if k > 0:
                    del_pathways_prod_data_t_specie.append(pathway)
                elif k < 0:
                    del_pathways_dest_data_t_specie.append(pathway)

        with open(chem_R_json, 'r') as chem_system_file:
            # Parse the JSON data and store it in a variable
            chem_system_data = json.load(chem_system_file)
        
        t_species = d_tools.get_compound_dict(compound=s)

        # The rate_sum for a species is the sum of its prod and destr
        rate_sum = t_species["production rate"]["active pathways"] + t_species["production rate"]["deleted pathways"] + t_species["destruction rate"]["active pathways"] + t_species["destruction rate"]["deleted pathways"]
        rate_sum_prod = t_species["production rate"]["active pathways"] + t_species["production rate"]["deleted pathways"] 
        rate_sum_dest = t_species["destruction rate"]["active pathways"] + t_species["destruction rate"]["deleted pathways"]
        # if rate_sum == 0.0 there is no prod or destr of s

        # Cases for text of the title
        if rate_sum == 0.0:
            rate_case = ''
        elif rate_sum_prod > rate_sum_dest:
            rate_case = r' Production '
        else:
            rate_case = r' Destruction '
        # we are working with absolute values
        rate_sum = abs(rate_sum)

        # Init Shits for plotting
        labels = []
        sizes = []

        print('For species: '+s)

        # FOR DESTRUCTION PATHWAYS
        # we re gonna sort the pathways in order of importance
        pathway_sorted = {}
        i = 0
        for pathway in act_pathways_prod_data_t_specie:
            stoich = abs(pathway["rate"]*pathway["branching points"][d_tools.find_compound_in_merged_list(listing=pathway["branching points"],compound=s)[0]]["stoichiometry"])/rate_sum * 100
            # We check that the pathway account for more that 0.0001% of the total rate of the species
            # Meaning we explain the last 0.1% with pathways
            # FOR THE PIE CHART:
            # We only save Pathways with % above 0.1%
            if stoich >= 0.1:
                print('stoich',stoich)
                pathway_sorted.update({i:stoich})
                i += 1
            # If not, then we we will not print it
            else:
                print('Not adding pathway with stoich',stoich,' to the Pie Chart')
                # we advance the indice also
                i += 1
        
        # Now we sort the indices
        ind_pathway_sorted = sorted(pathway_sorted,key=pathway_sorted.get,reverse=True)
        print('We have the ind dict sorted:',ind_pathway_sorted)

        if ind_pathway_sorted:
            text = text + ' \n'
            text = text + r'Production pathways of '+s
            text = text + ' \n'
        # for pathway in active_pathways_data:
        for ind in ind_pathway_sorted:
            pathway = act_pathways_prod_data_t_specie[ind]
            ind_p_in_AP = d_tools.find_pathway_in_list(pathway_to_be_found=pathway,list_of_pathways=active_pathways_data)
            labels.append('P'+str(ind_p_in_AP))
            rate_P_tmp = abs(pathway["rate"]*pathway["branching points"][d_tools.find_compound_in_merged_list(listing=pathway["branching points"],compound=s)[0]]["stoichiometry"])
            sizes.append(rate_P_tmp/rate_sum * 100)
            # Adding the text of the pathway
            text = text + ' \n'
            text = text + r'\textbf{P'+str(ind_p_in_AP)+r'} rate: '+r'{:0.3e}'.format(rate_P_tmp)
            text = text + ' \n'
            text = text + o_tools.pathway_to_latex_str(pathway=pathway,chem_system_data=chem_system_data)
            text = text + ' \n'

        # we re gonna sort the pathways in order of importance
        # FOR DESTRUCTION PATHWAYS
        pathway_sorted = {}
        i = 0
        for pathway in act_pathways_dest_data_t_specie:
            stoich = abs(pathway["rate"]*pathway["branching points"][d_tools.find_compound_in_merged_list(listing=pathway["branching points"],compound=s)[0]]["stoichiometry"])/rate_sum * 100
            # We check that the pathway account for more that 0.0001% of the total rate of the species
            # Meaning we explain the last 0.1% with pathways
            # FOR THE PIE CHART:
            # We only save Pathways with % above 0.1%
            if stoich >= 0.1:
                print('stoich',stoich)
                pathway_sorted.update({i:stoich})
                i += 1
            # If not, then we we will not print it
            else:
                print('Not adding pathway with stoich',stoich,' to the Pie Chart')
                # we advance the indice also
                i += 1
        
        # Now we sort the indices
        ind_pathway_sorted = sorted(pathway_sorted,key=pathway_sorted.get,reverse=True)
        print('We have the ind dict sorted:',ind_pathway_sorted)

        if ind_pathway_sorted:
            text = text + ' \n'
            text = text + r'Destruction pathways of '+s
            text = text + ' \n'
        # for pathway in active_pathways_data:
        for ind in ind_pathway_sorted:
            pathway = act_pathways_dest_data_t_specie[ind]
            ind_p_in_AP = d_tools.find_pathway_in_list(pathway_to_be_found=pathway,list_of_pathways=active_pathways_data)
            labels.append('P'+str(ind_p_in_AP))
            rate_P_tmp = abs(pathway["rate"]*pathway["branching points"][d_tools.find_compound_in_merged_list(listing=pathway["branching points"],compound=s)[0]]["stoichiometry"])
            sizes.append(rate_P_tmp/rate_sum * 100)
            # Adding the text of the pathway
            text = text + ' \n'
            text = text + r'\textbf{P'+str(ind_p_in_AP)+r'} rate: '+r'{:0.3e}'.format(rate_P_tmp)
            text = text + ' \n'
            text = text + o_tools.pathway_to_latex_str(pathway=pathway,chem_system_data=chem_system_data)
            text = text + ' \n'

        rate_deleted = 0.0
        for pathway in del_pathways_prod_data_t_specie:
            rate_deleted += abs(pathway["rate"]*pathway["branching points"][d_tools.find_compound_in_merged_list(listing=pathway["branching points"],compound=s)[0]]["stoichiometry"])/rate_sum * 100
        if rate_deleted >= 0.1:
            labels.append('P_del prod')
            sizes.append(rate_deleted)
            # Adding the text of the pathway
            text = text + ' \n'
            text = text + r'\textbf{P_del} Production of '+s+r' rate: '+r'{:0.3e}'.format(rate_deleted)
            text = text + ' \n'
        
        rate_deleted = 0.0
        for pathway in del_pathways_dest_data_t_specie:
            rate_deleted += abs(pathway["rate"]*pathway["branching points"][d_tools.find_compound_in_merged_list(listing=pathway["branching points"],compound=s)[0]]["stoichiometry"])/rate_sum * 100
        if rate_deleted >= 0.1:
            labels.append('P_del destr')
            sizes.append(rate_deleted)
            # Adding the text of the pathway
            text = text + ' \n'
            text = text + r'\textbf{P_del} Destruction of '+s+r' rate: '+r'{:0.3e}'.format(rate_deleted)
            text = text + ' \n'
        

        fig, ax = plt.subplots(nrows=1,ncols=2)
        plt.rc('text', usetex=True)
        plt.rc('text.latex', preamble=r'\usepackage[version=4]{mhchem}')
        fig.suptitle(rate_case+r'for species: \ce{'+s+r'} with rate '+r'{:0.3e}'.format(abs(rate_sum_prod-rate_sum_dest)))
        ax[0].pie(sizes, labels=labels, autopct='%1.1f%%')

        # Hide the axes
        ax[1].axis("off")

        # Display the text in the plot
        ax[1].text(0,0,text, fontsize=10)

        # Display the plot
        plt.tight_layout()
        plt.show()
