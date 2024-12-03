"""
Module Name
read_reaction_system

A brief description of what this module does or provides.

Detailed Description:
- Provide additional details about the module's functionality.
- Mention any key classes, functions, or variables defined in the module.
- Explain the module's role in the larger project or system.

Usage:
- Describe how to import this module.
- Provide examples of how to use the module's features or classes.
- Mention any common use cases.

Dependencies:
- List any external libraries or modules that this module depends on.
- Include version requirements if necessary.

Author:
- Your name or the name of the module's author.

License:
- Specify the module's licensing information if applicable.

Note:
- Include any important notes, warnings, or considerations.
- This module is part of the sub-package i__user_model from CPAP

"""
import shutil

from .i__user_model import convert_reaction_system_file as i_system
from .i__user_model import convert_concentration_file as i_concentration
from .i__user_model import check_target_species as i_ts
from .p__initialization import init_pathways as p_init
from .p__data_management import data_update as up
from .p__data_management import global_var
from .p__pathways_analysis import branching_points as bp
from .p__pathways_analysis import main_loop as ml
from .o__cpap_output import output as out
from .o__cpap_output import output_tools as o_tools
from .p__data_management import data_tools as d_tools


def init_global_var(chronicle_writing:bool,steps_save:bool):
    global_var.chronicle_writing = chronicle_writing
    global_var.steps_save = steps_save


def run_cpa(timestep:float,rate_threshold:float,t_min:float,target_species:list,filename_model:str,filename_concentration:str,final_AP_file:str,final_DP_file:str,final_CS_file:str,final_SL_file:str,chronicle_writing:bool,steps_save:bool) -> None:

    # init global var
    init_global_var(chronicle_writing=chronicle_writing,steps_save=steps_save)

    # first test is to convert a given text file into a workable JSON dataset
    if global_var.chronicle_writing:
        with open('chronicles.txt', 'w') as output_file:
            output_file.write('Start of the Chemical Pathway Analysis')
            output_file.write('\n')

        o_tools.write_line_chronicle('######################')
        o_tools.write_line_chronicle('User Inputs Processing')
        o_tools.write_line_chronicle('######################')
        o_tools.write_line_chronicle('\n')
    
    i_system.convert_chemical_reaction_file(filename=filename_model)
    chemical_species = i_concentration.convert_concentration_file(filename=filename_concentration,timestep=timestep)
    i_system.adding_pseudo_reactions(chemical_species=chemical_species)

    # 2. We run the initialization
    if global_var.chronicle_writing:
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('#######################')
        o_tools.write_line_chronicle('Pathways Initialization')
        o_tools.write_line_chronicle('#######################')
        o_tools.write_line_chronicle('\n')

    active_p,deleted_p = p_init.init_pathways(json_filename="chemical_reaction_system.json")

    if global_var.steps_save:
        # saving active/deleted pathways before updating the reaction/species rates
        d_tools.save_pathways_to_JSON(pathways=active_p,filename='active_pathways_0.json')
        d_tools.save_pathways_to_JSON(pathways=deleted_p,filename='deleted_pathways_0.json')

    if global_var.chronicle_writing:
        o_tools.write_line_chronicle('Updating prod/destr rates for chemical species')
    
    chemical_species = up.update_rates_chemical_species(active_p=active_p,deleted_p=deleted_p,chemical_species=chemical_species)
    # Checking the targeted species as viable outputs
    # if global_var.chronicle_writing:
    #     o_tools.write_line_chronicle('\n')
    #     o_tools.write_line_chronicle('----------------------')
    #     o_tools.write_line_chronicle('Check targeted species')
    #     o_tools.write_line_chronicle('----------------------')
    # target_species=i_ts.check_list_target_species(target_species=target_species,t_min=t_min)

    # 3. We run the main loop
    if global_var.chronicle_writing:
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('#################')
        o_tools.write_line_chronicle('Pathways Analysis')
        o_tools.write_line_chronicle('#################')
        o_tools.write_line_chronicle('\n')
    
    active_p,deleted_p,chemical_species = ml.main_loop(t_min=t_min,f_min=rate_threshold,active_p=active_p,deleted_p=deleted_p,chemical_species=chemical_species)

    # 4. main loop done. Outputs time!!!
    # out.text_output(target_species=target_species)

    # 5 copying results files
    d_tools.save_pathways_to_JSON(pathways=active_p,filename=final_AP_file)
    d_tools.save_pathways_to_JSON(pathways=deleted_p,filename=final_DP_file)
    d_tools.save_pathways_to_JSON(pathways=chemical_species,filename=final_SL_file)
    # shutil.copy2(src='active_pathways.json',dst=final_AP_file)
    # shutil.copy2(src='deleted_pathways.json',dst=final_DP_file)
    shutil.copy2(src='chemical_reaction_system.json',dst=final_CS_file)
    # shutil.copy2(src='chemical_species.json',dst=final_SL_file)
