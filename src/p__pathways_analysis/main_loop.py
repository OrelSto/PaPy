import json

from . import branching_points as bp

def main_loop(t_min:float):
    # main loop for connecting pathways and stuffs
    print('Here is the list of the next species considered as branching points for a fixed minimum timescale of ',t_min)
    list_bp = bp.list_next_branching_points(t_min=t_min)
    print(list_bp)

    # Opening JSON file
    ap = open('active_pathways.json')

    # returns JSON object as a dictionary
    active_p = json.load(ap)

    # Connecting pathways
    for species in list_bp:
        bp.connecting_pathways(active_pathways=active_p,species=species)