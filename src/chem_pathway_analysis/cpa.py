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

from .i__user_model import convert_reaction_system_file as i_system
from .i__user_model import convert_concentration_file as i_concentration
from .p__initialization import init_pathways as p_init
from .p__data_management import data_update as up
from .p__data_management import global_var
from .p__pathways_analysis import branching_points as bp
from .p__pathways_analysis import main_loop as ml
from .o__cpap_output import output_moche as out
from .o__cpap_output import output_tools as o_tools


def init_global_var(chronicle_writing:bool):
    global_var.chronicle_writing = chronicle_writing


def run_cpa(timestep:float,rate_threshold:float,t_min:float,target_specie:str,filename_model:str,filename_concentration:str,chronicle_writing=False) -> None:

    # init global var
    init_global_var(chronicle_writing=chronicle_writing)

    # first test is to convert a given text file into a workable JSON dataset

    print('######################')
    print('User Inputs Processing')
    print('######################')
    print()


    if global_var.chronicle_writing:
        with open('chronicles.txt', 'w') as output_file:
            output_file.write('Start of the Chemical Pathway Analysis')
            output_file.write('\n')

        o_tools.write_line_chronicle('######################')
        o_tools.write_line_chronicle('User Inputs Processing')
        o_tools.write_line_chronicle('######################')
        o_tools.write_line_chronicle('\n')
    
    i_system.convert_chemical_reaction_file(filename=filename_model)
    i_concentration.convert_concentration_file(filename=filename_concentration,timestep=timestep)
    i_system.adding_pseudo_reactions()

    # 2. We run the initialization
    if global_var.chronicle_writing:
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('#######################')
        o_tools.write_line_chronicle('Pathways Initialization')
        o_tools.write_line_chronicle('#######################')
        o_tools.write_line_chronicle('\n')
    print()
    print('#######################')
    print('Pathways Initialization')
    print('#######################')
    print()
    p_init.init_pathways(json_filename="chemical_reaction_system.json")
    print()
    print('Updating prod/destr rates for chemical species')
    up.update_rates_chemical_species()
    print()

    # 3. We run the main loop
    if global_var.chronicle_writing:
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('#################')
        o_tools.write_line_chronicle('Pathways Analysis')
        o_tools.write_line_chronicle('#################')
        o_tools.write_line_chronicle('\n')
    print('#################')
    print('Pathways Analysis')
    print('#################')
    print()
    ml.main_loop(t_min=t_min,f_min=rate_threshold,max_iter=3)

    # 4. main loop done. Outputs time!!!
    print('##################')
    print('Outputs formatting')
    print('##################')
    print()
    out.moche(target_specie)
