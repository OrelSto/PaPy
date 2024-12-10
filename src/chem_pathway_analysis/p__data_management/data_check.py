import json
import numpy as np
import re
from ..p__data_management import global_var
from ..o__cpap_output import output_tools as o_tools

def check_flux_conservation_species(chemical_species:list):
    # We want to check that we have the conservation of flux for each species
    for species in chemical_species:
        # We have the theoritical delta/sum of prod/destr
        theoritical_delta = abs(species["production rate"]["initial"] - species["destruction rate"]["initial"])
        delta = abs(species["production rate"]["active pathways"] + species["production rate"]["deleted pathways"] - species["destruction rate"]["active pathways"] - species["destruction rate"]["deleted pathways"])
        if theoritical_delta > 0.0:
            relative_error = abs(theoritical_delta-delta)/theoritical_delta
        else:
            relative_error = abs(theoritical_delta-delta)


        # if theoritical_delta != delta :
        if relative_error > 1e-9 :
            if global_var.chronicle_writing:
                 o_tools.write_line_chronicle("ISSUE FLUX CONSERVATION WITH "+species["name"])
                 o_tools.write_line_chronicle("theoritical_delta : "+str(theoritical_delta))
                 o_tools.write_line_chronicle("            delta : "+str(delta))
                 o_tools.write_line_chronicle("   relative error : "+str(relative_error))
        else:
            if global_var.chronicle_writing:
                 o_tools.write_line_chronicle("FLUX CONSERVED FOR "+species["name"])
                 o_tools.write_line_chronicle("theoritical_delta : "+str(theoritical_delta))
                 o_tools.write_line_chronicle("            delta : "+str(delta))
                 o_tools.write_line_chronicle("   relative error : "+str(relative_error))


def check_rates(active_pathways:dict,deleted_pathways:dict):
    pass


def elements_analysis_1st_step(chem_molecule: str) -> int:
    total_count = 0
    elements_dict = {}

    # First pattern matches every element once
    element_pattern = re.compile(r"([A][cglmrstu]|[B][aehikr]?|"
                                 r"\[3\]C?|\[1\]C?|"
                                 r"[C][adeflmnorsu]?|[D][bsy]|[E][rsu]|[F][elmr]?|"
                                 r"[G][ade]|[H][efgos]?|[I][nr]?|[K][r]?|[L][airuv]|[M][cdgnot]|[N][abdehiop]?|"
                                 r"O\[3p\]?|O\[1d\]?|"
                                 r"[O][gs]?|[P][abdmortu]?|[R][abefghnu]|[S][bcegimnr]?|[T][abcehilms]|[U]|[V]|"
                                 r"[W]|[X][e]|[Y][b]?|[Z][nr])([0-9]*)")

    matches = element_pattern.findall(chem_molecule)

    # Sum up the counts of individual elements
    for element, count in matches:
        if not count:
            count = 1
        # print('element: ',element, int(count))
        total_count += int(count)
        if element in elements_dict:
            elements_dict[element] += int(count)
        else:
            elements_dict[element] = int(count)

    return elements_dict


def elements_analysis(chem_molecule: str) -> int:
    elements_dict = {}

    elements_dict = elements_analysis_1st_step(chem_molecule)

    # Second pattern finds parentheses and multiplies elements within
    paren_pattern = re.compile(r"\((.+?)\)([0-9]+)")
    matches = paren_pattern.findall(chem_molecule)

    for subformula, multiplier in matches:
        el_tmp = {}
        print('subformula: ',subformula, int(multiplier))
        el_tmp = elements_analysis_1st_step(subformula)
        for k,v in el_tmp.items():
            elements_dict[k] += v * (int(multiplier) - 1)

    # Now we change the keys for atoms on different energy state etc
    changing_keys = {'O[3p]':'O','O[1d]':'O','O2[Dg]':'O2','N[2D]':'N','[3]C':'C','[1]C':'C'}
    elements_dict = {changing_keys[k] if k in changing_keys else k:v for k,v in elements_dict.items()}

    return elements_dict


def check_conservation(products,reactants):
    products_dict = {}
    reactants_dict = {}

    # We are gonna build the dict for all atoms in the products
    for p in products:
        p_tmp = {}
        p_tmp = elements_analysis(chem_molecule=p["compound"])
        if p_tmp:
            # then we multiply each value by stoichiometry
            p_tmp = {key:int(value * p["stoichiometry"]) for key,value in p_tmp.items()}
            # Then updating the main dict
            for element,count in p_tmp.items():
                if element in products_dict:
                    products_dict[element] += int(count)
                else:
                    products_dict[element] = int(count)
    # Sorting by value
    products_dict= {k: v for k, v in sorted(products_dict.items(), key=lambda item: item[0])}

    # We are gonna build the dict for all atoms in the reactants
    for r in reactants:
        r_tmp = {}
        r_tmp = elements_analysis(chem_molecule=r["compound"])
        if r_tmp:
            # then we multiply each value by stoichiometry
            r_tmp = {key:int(value * r["stoichiometry"]) for key,value in r_tmp.items()}
            # Then updating the main dict
            for element,count in r_tmp.items():
                if element in reactants_dict:
                    reactants_dict[element] += int(count)
                else:
                    reactants_dict[element] = int(count)
    # Sorting by value
    reactants_dict= {k: v for k, v in sorted(reactants_dict.items(), key=lambda item: item[0])}

    # print('products_dict',products_dict)
    # print('reactants_dict',reactants_dict)

    if products_dict == reactants_dict:
        flag = True
    else:
        flag = False

    return flag


def list_reaction_system_conservation(chem_reaction_system:list):
    # We check of the flag for conservation
    list_reactions_errors = []
    for reaction in chem_reaction_system:
        if not reaction["is_conserved"]:
            list_reactions_errors.append('The reaction line '+str(chem_reaction_system.index(reaction))+' does not conserve atoms')
    
    return list_reactions_errors
