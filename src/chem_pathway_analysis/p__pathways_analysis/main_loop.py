import json

from itertools import filterfalse
from . import branching_points as bp
from ..p__data_management import data_tools as d_tools
from ..p__data_management import data_update as up
from ..p__data_management import data_check as ch
from ..p__sub_pathways import subpathways_main as sub_main
from ..p__data_management import global_var
from ..o__cpap_output import output_tools as o_tools

def main_loop(t_min:float,f_min:float):
    
    # main loop for connecting pathways and stuffs
    print('Here is the list of the next species considered as branching points for a fixed minimum timescale of ',t_min)
    list_bp = bp.list_next_branching_points(t_min=t_min)
    print(list_bp)
    
    if global_var.chronicle_writing:
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('Here is the sorted list by lifetime of the next species considered as branching points for a fixed minimum timescale of '+'{:0.3e}'.format(t_min)+':')
        o_tools.write_line_chronicle(' '.join(list_bp))
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('Starting the pathways analysis')

    # This is the species used as BP
    used_species = []

    # We run the BP and Sub-BP routines until we are running out of BP
    while list_bp:
        # Opening JSON file
        ap = open('active_pathways.json')
        dp = open('deleted_pathways.json')

        # returns JSON object as a dictionary
        active_p = json.load(ap)
        deleted_p = json.load(dp)

        # Connecting pathways
        for species in list_bp:
            if global_var.chronicle_writing:
                o_tools.write_line_chronicle('\n')
                o_tools.write_line_chronicle('Pathways analysis for species: '+species)
            
            used_species.append(species)
            # looking for each species from the shortest lived to the longest
            active_p = bp.connecting_pathways(active_pathways=active_p,species=species,list_species_done=used_species)

            # cleaning pathways that are too slow. Keeping your pathway house tight and clean.
            active_p,deleted_p = bp.cleaning_slow_pathways(active_pathways=active_p,deleted_pathways=deleted_p,f_min=f_min)

            # Printing
            print()
            print('##############################')
            print("Sub Pathways analysis starting")
            print('##############################')
            print()

            if global_var.chronicle_writing:
                o_tools.write_line_chronicle('\n')
                o_tools.write_line_chronicle('##############################')
                o_tools.write_line_chronicle('Sub Pathways analysis starting')
                o_tools.write_line_chronicle('##############################')
                o_tools.write_line_chronicle('\n')
            
            # After saving, SUB-PATHWAYS analysis !!
            active_p = sub_main.main_subpathways(pathways=active_p,list_species_done=used_species)

            # saving
            d_tools.save_pathways_to_JSON(pathways=active_p,filename='active_pathways_'+species+'_'+'.json')
            d_tools.save_pathways_to_JSON(pathways=active_p,filename='active_pathways.json')
            d_tools.save_pathways_to_JSON(pathways=deleted_p,filename='deleted_pathways.json')

            # Printing
            print()
            print('##########################')
            print("Sub Pathways analysis done")
            print('##########################')
            print()

            if global_var.chronicle_writing:
                o_tools.write_line_chronicle('\n')
                o_tools.write_line_chronicle('##########################')
                o_tools.write_line_chronicle('Sub Pathways analysis DONE')
                o_tools.write_line_chronicle('##########################')
                o_tools.write_line_chronicle('\n')

            # Maybe some checking for conservation properties?
            # Like the distribution of the chemical reaction rates onto the pathways (active & deleted)

            # Printing
            print()
            print('!!!!!!!!!!!!!!!!!!!!!!!!')
            print("active pathways updated")
            print('!!!!!!!!!!!!!!!!!!!!!!!!')
            print()

            # Updating the chemical species
            print('Updating prod/destr rates for chemical species')
            up.update_rates_chemical_species()
            print()

            # Updating the chemical reaction system
            print('Updating rates for chemical reaction system')
            up.update_rates_reaction_system()
            print()

        # Selecting new Branching Points
        list_bp = bp.list_next_branching_points(t_min=t_min)
        # print('This is list_bp: ',list_bp)
        list_bp = list(filterfalse(lambda x: x in used_species,list_bp))
        print('This is used_species: ',used_species)
        print('This is not used_species: ',[not c for c in used_species])
        print('This is list_bp after flagged: ',list_bp)
    
    # Now that the main loop is over:
    # We check that the rates are conserved:
    ch.check_rates(active_pathways=active_p,deleted_pathways=deleted_p,)

