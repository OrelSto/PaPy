"""
Module Name
convert_reaction_system_file

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
- This module is part of the sub-package i__user_model from CPA

"""

import json
from ..p__pathways_analysis import branching_points as bp
from ..p__data_management import global_var
from ..o__cpap_output import output_tools as o_tools

def check_list_target_species(target_species:list,t_min:float):
    """check_list_target_specie _summary_

    _extended_summary_

    Parameters
    ----------
    reaction_equation : str
        _description_

    Returns
    -------
    _type_
        _description_
    """

    # loading the chemical species
    # Opening JSON file
    cs = open('chemical_species.json')
    # returns JSON object as a dictionary
    chemical_species = json.load(cs)

    # We get the list of BP
    list_bp = bp.list_next_branching_points(t_min=t_min)
    if global_var.chronicle_writing:
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('We have the Branching point '+str(list_bp))
        o_tools.write_line_chronicle('We have the Targeted Species '+str(target_species))

    for s in target_species:
        if global_var.chronicle_writing:
            o_tools.write_line_chronicle('\n')
            o_tools.write_line_chronicle('Looking for target specie '+s)
        if not (s in list_bp):
            if global_var.chronicle_writing:
                o_tools.write_line_chronicle(s+' is NOT in the list of Branching Points')
            for item_cs in chemical_species:
                if item_cs["name"] == s:
                    s_lifetime = item_cs["lifetime"]
                    if global_var.chronicle_writing:
                        o_tools.write_line_chronicle('The specie '+s+' cannot be processed as an output')
                        o_tools.write_line_chronicle('Because '+s+' is not a branching Point')
                        o_tools.write_line_chronicle('Lifetime of '+s+' is '+'{:0.3e}'.format(s_lifetime))
                        o_tools.write_line_chronicle('Compared to user t_min '+'{:0.3e}'.format(t_min))
                    target_species.remove(s)
        else:
            if global_var.chronicle_writing:
                o_tools.write_line_chronicle(s+' is in the list of Branching Points')

    return target_species