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

from i__user_model import convert_reaction_system_file as i_system
from i__user_model import convert_concentration_file as i_concentration
from p__initialization import init_pathways as p_init
from p__data_management import data_update as up
from p__pathways_analysis import branching_points as bp
from p__pathways_analysis import main_loop as ml

def cpa(timestep:float,rate_threshold:float,t_min:float) -> None:
    # first test is to convert a given text file into a workable JSON dataset
    print('######################')
    print('User Inputs Processing')
    print('######################')
    print()
    i_system.convert_chemical_reaction_file(filename='user_model_example.txt')
    i_concentration.convert_concentration_file(filename='user_concentration_example.txt',timestep=timestep)

    # 2. We run the initialization
    print()
    print('######################')
    print('Pathways Initiaziation')
    print('######################')
    print()
    p_init.init_pathways(json_filename="chemical_reaction_system.json")
    print()
    print('Updating prod/destr rates for chemical species')
    up.update_rates_chemical_species()
    print()

    # 3. We run the main loop
    print('#################')
    print('Pathways Analysis')
    print('#################')
    print()
    ml.main_loop(t_min=t_min,f_min=rate_threshold)

    # 4. main loop done. Outputs time!!!
    print('##################')
    print('Outputs formatting')
    print('##################')
    print()


if __name__=='__main__':
    # this is a stupid way to test the package and stupid values for inputs
    cpa(timestep=100.0,rate_threshold=1e-12,t_min=100.0)
