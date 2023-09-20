import json

from . import branching_points as bp
from p__data_management import data_tools as d_tools
from p__data_management import data_update as up

def main_loop(t_min:float,f_min:float):
    # main loop for connecting pathways and stuffs
    print('Here is the list of the next species considered as branching points for a fixed minimum timescale of ',t_min)
    list_bp = bp.list_next_branching_points(t_min=t_min)
    print(list_bp)

    while list_bp:
        # Opening JSON file
        ap = open('active_pathways.json')
        dp = open('deleted_pathways.json')

        # returns JSON object as a dictionary
        active_p = json.load(ap)
        deleted_p = json.load(dp)

        # Connecting pathways
        for species in list_bp:
            active_p = bp.connecting_pathways(active_pathways=active_p,species=species)
            active_p,deleted_p = bp.cleaning_pathways(active_pathways=active_p,deleted_pathways=deleted_p,f_min=f_min)

        # saving
        d_tools.save_pathways_to_JSON(pathways=active_p,filename='active_pathways.json')
        d_tools.save_pathways_to_JSON(pathways=deleted_p,filename='deleted_pathways.json')

        
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

