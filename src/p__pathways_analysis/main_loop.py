import json

from itertools import compress
from . import branching_points as bp
from p__data_management import data_tools as d_tools
from p__data_management import data_update as up
from p__sub_pathways import subpathways_main as sub_main

def main_loop(t_min:float,f_min:float):
    # main loop for connecting pathways and stuffs
    print('Here is the list of the next species considered as branching points for a fixed minimum timescale of ',t_min)
    list_bp = bp.list_next_branching_points(t_min=t_min)
    print(list_bp)

    # This is the main loop for the game!
    # The goal is to have an empty list_bp, meaning no more chemical species with a lifetime < t_min.
    # No more branching points, no more pathways to define.
    # Time to be happy and do some Science!

    while list_bp:
    # for i in range(10):
        # Opening JSON file
        ap = open('active_pathways.json')
        dp = open('deleted_pathways.json')

        # returns JSON object as a dictionary
        active_p = json.load(ap)
        deleted_p = json.load(dp)

        # Connecting pathways
        # setting up the flag list for species with no new pathways option
        flagged_species = []
        for species in list_bp:
            # looking for each species from the shortest lived to the longest
            flag,active_p = bp.connecting_pathways(active_pathways=active_p,species=species)
            flagged_species.append(flag)
            # cleaning pathways that are too slow. Keeping your pathway house tight and clean. if flagged, no cleaning necessary
            if not flag:
                active_p,deleted_p = bp.cleaning_slow_pathways(active_pathways=active_p,deleted_pathways=deleted_p,f_min=f_min)
            

            # saving
            d_tools.save_pathways_to_JSON(pathways=active_p,filename='active_pathways.json')
            d_tools.save_pathways_to_JSON(pathways=deleted_p,filename='deleted_pathways.json')

            # Printing
            print()
            print('##############################')
            print("Sub Pathways analysis starting")
            print('##############################')
            print()
            # After saving, SUB-PATHWAYS analysis !!
            sub_main.main_subpathways()

            # Printing
            print()
            print('##########################')
            print("Sub Pathways analysis done")
            print('##########################')
            print()

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

        # Selecting new Branching Points
        list_bp = bp.list_next_branching_points(t_min=t_min)
        # print('This is list_bp: ',list_bp)
        list_bp = list(compress(list_bp, flagged_species))
        # print('This is flagged_species: ',flagged_species)
        # print('This is not flagged_species: ',[not c for c in flagged_species])
        # print('This is list_bp after flagged: ',list_bp)
        

