import json
import os
import matplotlib.pyplot as plt
from matplotlib import gridspec

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


def pie_output(target_species:list,act_P_json:str,del_P_json:str,chem_R_json:str,spec_L_json:str,slow_percent:float):
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
        text_prod = ''
        text_destr = ''
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

        # P_trash collect all the small Pathways
        is_P_trash_prod = False
        P_trash_p = {'rate':0.0,'stoich':0.0}
        is_P_trash_destr = False
        P_trash_d = {'rate':0.0,'stoich':0.0}

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
        
        t_species,_ = d_tools.get_compound_dict_from_results(compound=s,SpecL=spec_L_json)

        # The rate_sum for a species is the sum of its prod and destr
        rate_sum = t_species["production rate"]["active pathways"] + t_species["production rate"]["deleted pathways"] + t_species["destruction rate"]["active pathways"] + t_species["destruction rate"]["deleted pathways"]

        # Check of the total rate
        if rate_sum < 1e-30:
            print('The species ',s,' has no prod or destr rate')
            break

        rate_sum_prod = t_species["production rate"]["active pathways"] + t_species["production rate"]["deleted pathways"] 
        rate_sum_dest = t_species["destruction rate"]["active pathways"] + t_species["destruction rate"]["deleted pathways"]
        # if rate_sum == 0.0 there is no prod or destr of s

        # Cases for text of the title
        if abs(rate_sum_prod - rate_sum_dest) <= 1e-20:
            rate_case = r' \textbf{Steady State} '
            rate_diff = 0.0
        elif rate_sum_prod > rate_sum_dest:
            rate_case = r' \textbf{Production} '
            rate_diff = abs(rate_sum_prod - rate_sum_dest)
        else:
            rate_case = r' \textbf{Destruction} '
            rate_diff = abs(rate_sum_prod - rate_sum_dest)
        # we are working with absolute values
        rate_sum = abs(rate_sum)

        # Init Shits for plotting
        labels = []
        sizes = []

        print('For species: '+s)

        # FOR PRODUCTION PATHWAYS
        # we re gonna sort the pathways in order of importance
        pathway_sorted = {}
        i = 0
        i_slow = 0
        for pathway in act_pathways_prod_data_t_specie:
            rate_P_tmp = abs(pathway["rate"]*pathway["branching points"][d_tools.find_compound_in_merged_list(listing=pathway["branching points"],compound=s)[0]]["stoichiometry"])
            stoich = rate_P_tmp/rate_sum * 100.0
            # We check that the pathway account for more that slow_percent of the total rate of the species
            # Meaning we explain the last slow_percent% with pathways
            # FOR THE PIE CHART:
            # We only save Pathways with % above slow_percent%
            if stoich >= slow_percent:
                print('stoich',stoich,'rate',rate_P_tmp)
                pathway_sorted.update({i:stoich})
            # If not, then we we will not print it 
                i += 1
            else:
                print('adding P with stoich',stoich,'rate',rate_P_tmp,' to P_trash prod')
                is_P_trash_prod = True
                P_trash_p['rate'] += rate_P_tmp
                P_trash_p['stoich'] += stoich
                # we advance the indice also
                i += 1
                i_slow += 1
        
        # Now we sort the indices
        ind_pathway_sorted = sorted(pathway_sorted,key=pathway_sorted.get,reverse=True)
        print('We have the ind dict sorted:',ind_pathway_sorted)

        if ind_pathway_sorted:
            text_prod = text_prod + ' \n'
            text_prod = text_prod + r'\textbf{Production pathways} of '+s
            text_prod = text_prod + ' \n'
        # for pathway in active_pathways_data:
        for ind in ind_pathway_sorted:
            pathway = act_pathways_prod_data_t_specie[ind]
            ind_p_in_AP = d_tools.find_pathway_in_list(pathway_to_be_found=pathway,list_of_pathways=active_pathways_data)
            labels.append('P'+str(ind_p_in_AP))
            rate_P_tmp = abs(pathway["rate"]*pathway["branching points"][d_tools.find_compound_in_merged_list(listing=pathway["branching points"],compound=s)[0]]["stoichiometry"])
            stoich=rate_P_tmp/rate_sum * 100.0
            sizes.append(stoich)
            # Adding the text_prod of the pathway
            text_prod = text_prod + ' \n'
            text_prod = text_prod + r'\textbf{P'+str(ind_p_in_AP)+r'} rate: '+r'{:0.2e}'.format(rate_P_tmp)+r' or 'r'{:0.2f}'.format(stoich)+r'\%'
            text_prod = text_prod + ' \n'
            text_prod = text_prod + o_tools.pathway_to_latex_str(pathway=pathway,chem_system_data=chem_system_data)
            text_prod = text_prod + ' \n'

        # we re gonna sort the pathways in order of importance
        # FOR DESTRUCTION PATHWAYS
        pathway_sorted = {}
        j = 0
        j_slow = 0
        for pathway in act_pathways_dest_data_t_specie:
            rate_P_tmp = abs(pathway["rate"]*pathway["branching points"][d_tools.find_compound_in_merged_list(listing=pathway["branching points"],compound=s)[0]]["stoichiometry"])
            stoich = rate_P_tmp/rate_sum * 100.0
            # We check that the pathway account for more that 0.0001% of the total rate of the species
            # Meaning we explain the last 0.1% with pathways
            # FOR THE PIE CHART:
            # We only save Pathways with % above 0.1%
            if stoich >= slow_percent:
                print('stoich',stoich,'rate',rate_P_tmp)
                pathway_sorted.update({j:stoich})
                j += 1
            # If not, then we we will not print it
            else:
                print('adding P with stoich',stoich,'rate',rate_P_tmp,' to P_trash destr')
                is_P_trash_destr = True
                P_trash_d['rate'] += rate_P_tmp
                P_trash_d['stoich'] += stoich
                # we advance the indice also
                j += 1
                j_slow += 1
        
        # Now we sort the indices
        ind_pathway_sorted = sorted(pathway_sorted,key=pathway_sorted.get,reverse=True)
        print('We have the ind dict sorted:',ind_pathway_sorted)

        if ind_pathway_sorted:
            text_destr = text_destr + ' \n'
            text_destr = text_destr + r'\textbf{Destruction pathways} of '+s
            text_destr = text_destr + ' \n'
        # for pathway in active_pathways_data:
        for ind in ind_pathway_sorted:
            pathway = act_pathways_dest_data_t_specie[ind]
            ind_p_in_AP = d_tools.find_pathway_in_list(pathway_to_be_found=pathway,list_of_pathways=active_pathways_data)
            labels.append('P'+str(ind_p_in_AP))
            rate_P_tmp = abs(pathway["rate"]*pathway["branching points"][d_tools.find_compound_in_merged_list(listing=pathway["branching points"],compound=s)[0]]["stoichiometry"])
            stoich=rate_P_tmp/rate_sum * 100.0
            sizes.append(stoich)
            # Adding the text_destr of the pathway
            text_destr = text_destr + ' \n'
            text_destr = text_destr + r'\textbf{P'+str(ind_p_in_AP)+r'} rate: '+r'{:0.2e}'.format(rate_P_tmp)+r' or 'r'{:0.2f}'.format(stoich)+r'\%'
            text_destr = text_destr + ' \n'
            text_destr = text_destr + o_tools.pathway_to_latex_str(pathway=pathway,chem_system_data=chem_system_data)
            text_destr = text_destr + ' \n'

        rate_deleted = 0.0
        for pathway in del_pathways_prod_data_t_specie:
            rate_deleted += abs(pathway["rate"]*pathway["branching points"][d_tools.find_compound_in_merged_list(listing=pathway["branching points"],compound=s)[0]]["stoichiometry"])
        stoich = rate_deleted/rate_sum * 100.0
        print('stoich deleted for destr: ',stoich,' rate ',rate_deleted,)
        if stoich >= slow_percent:
            labels.append('P_del prod')
            sizes.append(stoich)
            # Adding the text_prod of the pathway
            text_prod = text_prod + ' \n'
            text_prod = text_prod + r'\textbf{P$_{\mathbf{\mathrm{del}}}$ Production} of '+s+r' rate: '+r'{:0.2e}'.format(rate_deleted)+r' or 'r'{:0.2f}'.format(stoich)+r'\%'
            text_prod = text_prod + ' \n'
            text_prod = text_prod + r'\ce{... ->'+s+r'}'
            text_prod = text_prod + ' \n'
        elif stoich < slow_percent and rate_deleted != 0:
            print('adding P del with total stoich',stoich,'rate',rate_deleted,' to P_trash prod')
            is_P_trash_prod = True
            P_trash_p['rate'] += rate_deleted
            P_trash_p['stoich'] += stoich
            # we need to add the number of del pathways added to P_trash
        
        rate_deleted = 0.0
        for pathway in del_pathways_dest_data_t_specie:
            rate_deleted += abs(pathway["rate"]*pathway["branching points"][d_tools.find_compound_in_merged_list(listing=pathway["branching points"],compound=s)[0]]["stoichiometry"])
        stoich = rate_deleted/rate_sum * 100.0
        print('stoich deleted for prod: ',stoich,' rate ',rate_deleted,)
        if stoich >= slow_percent:
            labels.append('P_del destr')
            sizes.append(stoich)
            # Adding the text_destr of the pathway
            text_destr = text_destr + ' \n'
            text_destr = text_destr + r'\textbf{P$_{\mathbf{\mathrm{del}}}$ Destruction} of '+s+r' rate: '+r'{:0.2e}'.format(rate_deleted)+r' or 'r'{:0.2f}'.format(stoich)+r'\%'
            text_destr = text_destr + ' \n'
            text_destr = text_destr + r'\ce{'+s+r' -> ...}'
            text_destr = text_destr + ' \n'
        elif stoich < slow_percent and rate_deleted != 0:
            print('adding P del with total stoich',stoich,'rate',rate_deleted,' to P_trash destr')
            is_P_trash_destr = True
            P_trash_d['rate'] += rate_deleted
            P_trash_d['stoich'] += stoich
            # we need to add the number of del pathways added to P_trash

        # if is_P_trash and P_trash['stoich']>=slow_percent:
        if is_P_trash_destr:
            print('We have P_slow destr stoich: ',P_trash_d['stoich'])
            labels.append('P_slow d')
            sizes.append(P_trash_d['stoich'])
            # Adding the text of the pathway
            text_destr = text_destr + ' \n'
            text_destr = text_destr + r'\textbf{P$_{\mathbf{\mathrm{slow}}}$ destr} of '+s+r' rate: '+r'{:0.2e}'.format(P_trash_d['rate'])+r' or 'r'{:0.2f}'.format(P_trash_d['stoich'])+r'\%'
            text_destr = text_destr + ' \n'
        
        # if is_P_trash and P_trash['stoich']>=slow_percent:
        if is_P_trash_prod:
            print('We have P_slow destr stoich: ',P_trash_d['stoich'])
            labels.append('P_slow p')
            sizes.append(P_trash_p['stoich'])
            # Adding the text_prod of the pathway
            text_prod = text_prod + ' \n'
            text_prod = text_prod + r'\textbf{P$_{\mathbf{\mathrm{slow}}}$ prod} of '+s+r' rate: '+r'{:0.2e}'.format(P_trash_p['rate'])+r' or 'r'{:0.2f}'.format(P_trash_p['stoich'])+r'\%'
            text_prod = text_prod + ' \n'

        if is_P_trash_prod and not is_P_trash_destr:
            text_prod = text_prod + ' \n'
            text_prod = text_prod + str(i_slow+j_slow) + r' pathways in \textbf{P$_{\mathbf{\mathrm{slow}}}$}'
            text_prod = text_prod + ' \n'
        elif is_P_trash_destr and not is_P_trash_prod:
            text_destr = text_destr + ' \n'
            text_destr = text_destr + str(i_slow+j_slow) + r' pathways in \textbf{P$_{\mathbf{\mathrm{slow}}}$}'
            text_destr = text_destr + ' \n'
        elif is_P_trash_destr and is_P_trash_prod:
            text_prod = text_prod + ' \n'
            text_prod = text_prod + str(i_slow+j_slow) + r' pathways in \textbf{P$_{\mathbf{\mathrm{slow}}}$}'
            text_prod = text_prod + ' \n'

        # fig, ax = plt.subplots(nrows=1,ncols=3)
        fig = plt.figure()
        gs = gridspec.GridSpec(1, 3, width_ratios=[2, 1, 1])
        ax0 = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1])
        ax2 = plt.subplot(gs[2])
        plt.rc('text', usetex=True)
        plt.rc('text.latex', preamble=r'\usepackage[version=4]{mhchem}')
        fig.suptitle(rate_case+r'for species: \ce{'+s+r'} with rate '+r'{:0.2e}'.format(rate_diff))
        ax0.pie(sizes, labels=labels, autopct='%1.2f%%')

        # Hide the axes
        ax1.axis("off")
        ax2.axis("off")

        # Display the text prod in the plot
        ax1.text(0,0,text_prod, fontsize=10)

        # Display the text destr in the plot
        ax2.text(0,0,text_destr, fontsize=10)

        # Display the plot
        # plt.tight_layout()
        plt.show()


def table_Tex(target_species:list,unit:str,act_P_json:str,del_P_json:str,chem_R_json:str,spec_L_json:str,slow_percent:float):

    header = r"""\documentclass{article}
\usepackage{makecell} % Required for multiple lines cell
\usepackage[version=4]{mhchem}
\usepackage{siunitx}
"""
    for s in target_species:
        # We start the text
        text = ''
        text = text + header
        # Before everything we init the text to print the pathways
        text = text + r'\begin{document}'
        text = text + ' \n'
        text = text + r'\begin{table}'
        text = text + ' \n'
        text = text + r'\centering'
        text = text + ' \n'
        text = text + r'\begin{tabular}{ | c | c | c | c | }'
        text = text + ' \n'
        text = text + r' \hline'
        text = text + ' \n'
        text = text + r'Pathway & Cycle & rate \unit{'+unit+r'} & \% \\'
        text = text + ' \n'
        text = text + r' \hline'
        text = text + ' \n'

        # NOW we start to deal with the results!
        with open(act_P_json, 'r') as active_pathways_file:
            # Parse the JSON data and store it in a variable
            active_pathways_data = json.load(active_pathways_file)
        # We set up the list of pathways acting on target_specie
        act_pathways_prod_data_t_specie = []
        act_pathways_dest_data_t_specie = []

        # P_trash collect all the small Pathways
        is_P_trash = False
        P_trash = {'rate':0.0,'stoich':0.0}

        for pathway in active_pathways_data:
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
        
        t_species,_ = d_tools.get_compound_dict_from_results(compound=s,SpecL=spec_L_json)
        print(t_species)

        # The rate_sum for a species is the sum of its prod and destr
        rate_sum = t_species["production rate"]["active pathways"] + t_species["production rate"]["deleted pathways"] + t_species["destruction rate"]["active pathways"] + t_species["destruction rate"]["deleted pathways"]

        # Check of the total rate
        if rate_sum < 1e-20:
            print('The species ',s,' has no prod or destr rate')
            break
        
        # FOR PRODUCTION PATHWAYS
        # we re gonna sort the pathways in order of importance
        pathway_sorted = {}
        i = 0
        for pathway in act_pathways_prod_data_t_specie:
            rate_P_tmp = abs(pathway["rate"]*pathway["branching points"][d_tools.find_compound_in_merged_list(listing=pathway["branching points"],compound=s)[0]]["stoichiometry"])
            stoich = rate_P_tmp/rate_sum * 100.0
            # We check that the pathway account for more that 0.0001% of the total rate of the species
            # Meaning we explain the last 0.1% with pathways
            # FOR THE PIE CHART:
            # We only save Pathways with % above 0.1%
            if stoich >= slow_percent:
                print('stoich',stoich)
                pathway_sorted.update({i:stoich})
                i += 1
            # If not, then we we will not print it
            else:
                print('adding pathway with stoich',stoich,' to the P_trash')
                is_P_trash = True
                P_trash['rate'] += rate_P_tmp
                P_trash['stoich'] += stoich
                # we advance the indice also
                i += 1
        
        # Now we sort the indices
        ind_pathway_sorted = sorted(pathway_sorted,key=pathway_sorted.get,reverse=True)
        print('We have the ind dict sorted:',ind_pathway_sorted)

        for ind in ind_pathway_sorted:
            pathway = act_pathways_prod_data_t_specie[ind]
            ind_p_in_AP = d_tools.find_pathway_in_list(pathway_to_be_found=pathway,list_of_pathways=active_pathways_data)
            cell1 = r'P'+str(ind_p_in_AP)
            rate_P_tmp = abs(pathway["rate"]*pathway["branching points"][d_tools.find_compound_in_merged_list(listing=pathway["branching points"],compound=s)[0]]["stoichiometry"])
            cell3 = r'{:0.3e}'.format(rate_P_tmp)
            cell4 = r'{:.2f}'.format(rate_P_tmp/rate_sum * 100.0)
            cell2 = o_tools.pathway_to_latex_cell(pathway=pathway,chem_system_data=chem_system_data)
            text = text + cell1 + r' & ' + cell2 + r' & ' + cell3 + r' & ' + cell4 + r' \\' 
            text = text + ' \n'
            text = text + r' \hline'
            text = text + ' \n'

        # FOR DESTRUCTION PATHWAYS
        # we re gonna sort the pathways in order of importance
        pathway_sorted = {}
        i = 0
        for pathway in act_pathways_dest_data_t_specie:
            rate_P_tmp = abs(pathway["rate"]*pathway["branching points"][d_tools.find_compound_in_merged_list(listing=pathway["branching points"],compound=s)[0]]["stoichiometry"])
            stoich = rate_P_tmp/rate_sum * 100.0
            # We check that the pathway account for more that 0.0001% of the total rate of the species
            # Meaning we explain the last 0.1% with pathways
            # FOR THE PIE CHART:
            # We only save Pathways with % above 0.1%
            if stoich >= slow_percent:
                print('stoich',stoich)
                pathway_sorted.update({i:stoich})
                i += 1
            # If not, then we we will not print it
            else:
                print('adding pathway with stoich',stoich,' to the P_trash')
                is_P_trash = True
                P_trash['rate'] += rate_P_tmp
                P_trash['stoich'] += stoich
                # we advance the indice also
                i += 1
        
        # Now we sort the indices
        ind_pathway_sorted = sorted(pathway_sorted,key=pathway_sorted.get,reverse=True)
        print('We have the ind dict sorted:',ind_pathway_sorted)

        for ind in ind_pathway_sorted:
            pathway = act_pathways_dest_data_t_specie[ind]
            ind_p_in_AP = d_tools.find_pathway_in_list(pathway_to_be_found=pathway,list_of_pathways=active_pathways_data)
            cell1 = r'P'+str(ind_p_in_AP)
            rate_P_tmp = abs(pathway["rate"]*pathway["branching points"][d_tools.find_compound_in_merged_list(listing=pathway["branching points"],compound=s)[0]]["stoichiometry"])
            cell2 = o_tools.pathway_to_latex_cell(pathway=pathway,chem_system_data=chem_system_data)
            cell3 = r'{:0.3e}'.format(rate_P_tmp)
            cell4 = r'{:.2f}'.format(rate_P_tmp/rate_sum * 100.0)
            text = text + cell1 + r' & ' + cell2 + r' & ' + cell3 + r' & ' + cell4 + r' \\' 
            text = text + ' \n'
            text = text + r' \hline'
            text = text + ' \n'

        # PRODUCTION DELETED RATES
        rate_deleted = 0.0
        for pathway in del_pathways_prod_data_t_specie:
            rate_deleted += abs(pathway["rate"]*pathway["branching points"][d_tools.find_compound_in_merged_list(listing=pathway["branching points"],compound=s)[0]]["stoichiometry"])
        stoich = rate_deleted/rate_sum * 100.0
        if stoich >= slow_percent:
            cell1 = r'P$_{del}$ prod'
            cell2 = r'\ce{ ... -> '+s+r'}'
            cell3 = r'{:0.3e}'.format(rate_deleted)
            cell4 = r'{:.2f}'.format(stoich)
            text = text + cell1 + r' & ' + cell2 + r' & ' + cell3 + r' & ' + cell4 + r' \\' 
            text = text + ' \n'
            text = text + r' \hline'
            text = text + ' \n'
        elif stoich < slow_percent and rate_deleted != 0:
            # Add to P_slow
            print('adding del pathway with stoich',stoich,' to the P_trash')
            is_P_trash = True
            P_trash['rate'] += rate_deleted
            P_trash['stoich'] += stoich

        
        # DESTRUCTION DELETED RATES
        rate_deleted = 0.0
        for pathway in del_pathways_dest_data_t_specie:
            rate_deleted += abs(pathway["rate"]*pathway["branching points"][d_tools.find_compound_in_merged_list(listing=pathway["branching points"],compound=s)[0]]["stoichiometry"])
        stoich = rate_deleted/rate_sum * 100.0
        if stoich >= slow_percent:
            cell1 = r'P$_{del}$ destr'
            cell2 = r'\ce{ '+s+r' -> ...}'
            cell3 = r'{:0.3e}'.format(rate_deleted)
            cell4 = r'{:.2f}'.format(stoich)
            text = text + cell1 + r' & ' + cell2 + r' & ' + cell3 + r' & ' + cell4 + r' \\' 
            text = text + ' \n'
            text = text + r' \hline'
            text = text + ' \n'
        elif stoich < slow_percent and rate_deleted != 0:
            # Add to P_slow
            print('adding del pathway with stoich',stoich,' to the P_trash')
            is_P_trash = True
            P_trash['rate'] += rate_deleted
            P_trash['stoich'] += stoich
        
        # Is there P_trash
        if is_P_trash:
            cell1 = r'P$_{slow}$'
            cell2 = r'Cycles $<$ ' + r'{:.2f}'.format(slow_percent) + r' \%'
            cell3 = r'{:0.3e}'.format(P_trash['rate'])
            cell4 = r'{:.2f}'.format(P_trash['stoich'])
            text = text + cell1 + r' & ' + cell2 + r' & ' + cell3 + r' & ' + cell4 + r' \\' 
            text = text + ' \n'
            text = text + r' \hline'
            text = text + ' \n'

        # End of table
        text = text + r'\end{tabular}'
        text = text + ' \n'
        text = text + r'\caption{Relevant pathways for $\ce{'+s+r'}$. The rate is already evaluated for one molecule/mole of $\ce{'+s+r'}$.}'
        text = text + ' \n'
        text = text + r'\label{table:P_'+s+r'}'
        text = text + ' \n'
        text = text + r'\end{table}'
        text = text + ' \n'
        text = text + r'\end{document}'
 
        # saving file:
        with open("table_"+s+".tex",'w') as table_file:
            table_file.write(text)
        # os.system("pdflatex table_"+s+".tex")
    

def list_pathways_Tex(list_P_json:str,chem_R_json:str,filename_sav:str):
    # A simple table with the pathways written in a nice form
    header = r"""\documentclass{article}
\usepackage{makecell} % Required for multiple lines cell
\usepackage[version=4]{mhchem}
\usepackage{siunitx}
\usepackage{longtable}
"""
    # We start the text
    text = ''
    text = text + header
    # Before everything we init the text to print the pathways
    text = text + r'\begin{document}'
    text = text + ' \n'
    text = text + r'\begin{longtable}{ | c | c | }'
    text = text + ' \n'
    text = text + r' \hline'
    text = text + ' \n'
    text = text + r'Pathway & Cycle \\'
    text = text + ' \n'
    text = text + r' \hline'
    text = text + ' \n'

    # NOW we start to deal with the results!
    with open(list_P_json, 'r') as active_pathways_file:
        # Parse the JSON data and store it in a variable
        active_pathways_data = json.load(active_pathways_file)

    with open(chem_R_json, 'r') as chem_system_file:
            # Parse the JSON data and store it in a variable
            chem_system_data = json.load(chem_system_file)

    ind = 0
    for pathway in active_pathways_data:
        cell1 = r'P'+str(ind)
        cell2 = o_tools.pathway_to_latex_cell(pathway=pathway,chem_system_data=chem_system_data)
        text = text + cell1 + r' & ' + cell2 + r' \\' 
        text = text + ' \n'
        text = text + r' \hline'
        text = text + ' \n'
        ind += 1

    # End of table
    text = text + r'\end{longtable}'
    text = text + ' \n'
    text = text + r'\end{document}'

    # saving file:
    with open(filename_sav,'w') as table_file:
        table_file.write(text)
    # os.system("pdflatex table_"+s+".tex")