import json

from itertools import filterfalse
from . import branching_points as bp
from ..p__data_management import data_tools as d_tools
from ..p__data_management import data_update as up
from ..p__data_management import data_check as ch
from ..p__sub_pathways import subpathways_main as sub_main
from ..p__data_management import global_var
from ..o__cpap_output import output_tools as o_tools

def main_loop(t_min:float,f_min:float,active_p:list,deleted_p:list,chemical_species:list,BP_species:str):
    
    # main loop for connecting pathways and stuffs
    # print('Here is the list of the next species considered as branching points for a fixed minimum timescale of ',t_min)
    list_CS = bp.list_chemical_species(chemical_species=chemical_species)
    print('Here is the list of all the CS in the system: ',list_CS)

    # Test over the potential Chemical species selected as the BP of interest
    if (BP_species != 'None') and (BP_species not in list_CS):
        print('You selected BP_species ',BP_species)
        print('But ',BP_species,' is not present in the chemical system')
        print(list_CS)
        print('Please select a chemical species present or put BP_species = str(None)')
        exit()
    elif (BP_species != 'None') and (BP_species in list_CS):
        print('You selected BP_species ',BP_species)
        print('And ',BP_species,' is present in the chemical system')
        print(list_CS)
        print('Since you are interested in ',BP_species)
        print('The program will override the t_min provided.')
        print('The program will run until ',BP_species,' is used as Branching Points')
        # we overide t_min
        t_min = 1e99
        # And the list BP will be all the BP until BP_species
        # we check if BP_species is the last entry
        if list_CS.index(BP_species) != len(list_CS)-1:
            list_bp = list_CS[0:list_CS.index(BP_species)+1]
            # This is the species used as BP
            used_species = list_CS[list_CS.index(BP_species)+1:len(list_CS)]
        else:
            list_bp = list_CS
            # This is the species used as BP
            used_species = []
        print('Here is the list of all The potential BP in the system: ',list_bp)
    else:
        list_bp = bp.list_next_branching_points(t_min=t_min,chemical_species=chemical_species)
        print('Here is the list of all The potential BP in the system: ',list_bp)
        print('according to t_min: ',t_min)
        # This is the species used as BP
        used_species = []

    if global_var.chronicle_writing:
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('Here is the sorted list by lifetime of all species present in the system:')
        o_tools.write_line_chronicle(' '.join(list_CS))
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('Here is the sorted list by lifetime of the next species considered as branching points for a fixed minimum timescale of '+'{:0.3e}'.format(t_min)+':')
        o_tools.write_line_chronicle(' '.join(list_bp))
        if (BP_species != 'None') and (BP_species in list_CS):
            o_tools.write_line_chronicle('\n')
            o_tools.write_line_chronicle('You selected BP_species '+BP_species)
            o_tools.write_line_chronicle('\n')
            o_tools.write_line_chronicle('And '+BP_species+' is present in the chemical system')
            o_tools.write_line_chronicle('\n')
            o_tools.write_line_chronicle(' '.join(list_CS))
            o_tools.write_line_chronicle('\n')
            o_tools.write_line_chronicle('Since you are interested in '+BP_species)
            o_tools.write_line_chronicle('\n')
            o_tools.write_line_chronicle('The program will override the t_min provided.')
            o_tools.write_line_chronicle('\n')
            o_tools.write_line_chronicle('The program will run until '+BP_species+' is used as Branching Points')
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('##############################')
        o_tools.write_line_chronicle('Starting the pathways analysis')
        o_tools.write_line_chronicle('##############################')
        o_tools.write_line_chronicle('\n')

    # # Opening JSON file
    # ap = open('active_pathways.json')
    # dp = open('deleted_pathways.json')

    # # returns JSON object as a dictionary
    # active_p = json.load(ap)
    # deleted_p = json.load(dp)

    step = 0

    # We run the BP and Sub-BP routines until we are running out of BP
    # Connecting pathways
    # for species in list_bp:
    while list_bp:
        species = list_bp[0]
        if global_var.chronicle_writing:
            o_tools.write_line_chronicle('\n')
            o_tools.write_line_chronicle('Pathways analysis for species: '+species)
        
        # # returns JSON object as a dictionary
        # with open('active_pathways.json') as ap:
        #     active_p = json.load(ap)

        # We add the current BP to the used list to update list_bp
        used_species.append(species)
        # looking for each species from the shortest lived to the longest
        active_p,deleted_p = bp.connecting_pathways(active_pathways=active_p,species=species,list_species_done=used_species,chemical_species=chemical_species,deleted_pathways=deleted_p)

        # Printing
        if global_var.chronicle_writing:
            o_tools.write_line_chronicle('\n')
            o_tools.write_line_chronicle('##############################')
            o_tools.write_line_chronicle('Sub Pathways analysis starting')
            o_tools.write_line_chronicle('##############################')
            o_tools.write_line_chronicle('\n')
        
        # After saving, SUB-PATHWAYS analysis !!
        active_p = sub_main.main_subpathways(pathways=active_p,list_species_done=used_species,chemical_species=chemical_species)

        # Printing
        if global_var.chronicle_writing:
            o_tools.write_line_chronicle('\n')
            o_tools.write_line_chronicle('##########################')
            o_tools.write_line_chronicle('Sub Pathways analysis DONE')
            o_tools.write_line_chronicle('##########################')
            o_tools.write_line_chronicle('\n')

        # cleaning pathways that are too slow. Keeping your pathway house tight and clean.
        # The cleaning of slow Pathways MUST BE DONE AFTER the SubPathways analysis
        # Indeed, you don't know A PRIOR how new pathways might be deconstruted into EXISTING PATHWAYS
        # Printing
        if global_var.chronicle_writing:
            o_tools.write_line_chronicle('\n')
            o_tools.write_line_chronicle('#################')
            o_tools.write_line_chronicle('Cleaning Pathways')
            o_tools.write_line_chronicle('#################')
            o_tools.write_line_chronicle('\n')
        
        active_p,deleted_p = bp.cleaning_slow_pathways(active_pathways=active_p,deleted_pathways=deleted_p,f_min=f_min)

        # Maybe some checking for conservation properties?
        # Like the distribution of the chemical reaction rates onto the pathways (active & deleted)

        # Now we are at the step where we update some characteristics
        # Namely: the prod/destr rates of every species and the reactions rates
        # First, we need to update the active_p, why?
        # Because there is a possibility that when constructing new pathways
        # One of the prod/destr pathways might be related to destr/prod of the BP
        # Hence, a part of its rate is deleted and is not yet saved in active_p

        # # saving active/deleted pathways before updating the reaction/species rates
        # # d_tools.save_pathways_to_JSON(pathways=active_p,filename='active_pathways_'+species+'.json')
        # d_tools.save_pathways_to_JSON(pathways=active_p,filename='active_pathways.json')
        # # d_tools.save_pathways_to_JSON(pathways=deleted_p,filename='deleted_pathways_'+species+'.json')
        # d_tools.save_pathways_to_JSON(pathways=deleted_p,filename='deleted_pathways.json')
        
        # Updating the chemical species
        chemical_species = up.update_rates_chemical_species(active_p=active_p,deleted_p=deleted_p,chemical_species=chemical_species)

        # Updating the chemical reaction system
        up.update_rates_reaction_system(deleted_p=deleted_p)

        # Checking species flux conservation
        if global_var.chronicle_writing:
            o_tools.write_line_chronicle('\n')
            o_tools.write_line_chronicle('##########################')
            o_tools.write_line_chronicle('Checking flux conservation')
            o_tools.write_line_chronicle('##########################')
            o_tools.write_line_chronicle('\n')
        ch.check_flux_conservation_species(chemical_species=chemical_species)

        # Selecting new Branching Points
        list_bp = bp.list_next_branching_points(t_min=t_min,chemical_species=chemical_species)
        # print('This is list_bp: ',list_bp)
        list_bp = list(filterfalse(lambda x: x in used_species,list_bp))
        # print('This is used_species: ',used_species)
        # print('This is not used_species: ',[not c for c in used_species])
        # print('This is list_bp after flagged: ',list_bp)

        # Updating the step
        step += 1
        if global_var.steps_save:
            # saving active/deleted pathways before updating the reaction/species rates
            d_tools.save_pathways_to_JSON(pathways=active_p,filename='active_pathways_'+str(step)+'_'+species+'.json')
            d_tools.save_pathways_to_JSON(pathways=deleted_p,filename='deleted_pathways_'+str(step)+'_'+species+'.json')


    # Now that the main loop is over:
    # We check that the rates are conserved:
    ch.check_rates(active_pathways=active_p,deleted_pathways=deleted_p,)

    if global_var.chronicle_writing:
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('################')
        o_tools.write_line_chronicle('END OF MAIN LOOP')
        o_tools.write_line_chronicle('################')
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('################################################')
        o_tools.write_line_chronicle('list of BP used to construct the active pathways')
        o_tools.write_line_chronicle('################################################')
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle(' '.join(used_species))
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('We advise the User to check by himslef the values for flux conservation')
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('A lot a details and additional steps infos are provided in this file')
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('The User would be wise to check all the steps to check any weird occurence')
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('Especially if the results are used for a scientific publication')
        o_tools.write_line_chronicle('\n')

    print('END OF MAIN LOOP')

    return active_p,deleted_p,chemical_species

