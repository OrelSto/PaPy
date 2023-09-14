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

from i__user_model import convert_reaction_system_file as i_convert
from p__initialization import init_pathways as p_init

def cpa() -> None:
    # first test is to convert a given text file into a workable JSON dataset
    i_convert.convert_chemical_reaction_file(filename='user_model_example.txt')

    # 2. We run the initialization
    p_init.init()

if __name__=='__main__':
    cpa()
