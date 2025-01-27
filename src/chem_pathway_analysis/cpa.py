import shutil

from .i__user_model import convert_reaction_system_file as i_system
from .i__user_model import convert_concentration_file as i_concentration
# from .i__user_model import check_target_species as i_ts
from .p__initialization import init_pathways as p_init
from .p__data_management import data_update as up
from .p__data_management import global_var
from .p__pathways_analysis import branching_points as bp
from .p__pathways_analysis import main_loop as ml
from .o__cpap_output import output as out
from .o__cpap_output import output_tools as o_tools
from .p__data_management import data_tools as d_tools


def init_global_var(chronicle_writing:bool,steps_save:bool,rate_threshold_BP_auto:float):
    global_var.chronicle_writing = chronicle_writing
    global_var.steps_save = steps_save
    global_var.rate_threshold_BP_auto = rate_threshold_BP_auto


def run_cpa(timestep:float,rate_threshold:float,t_min:float,BP_species:str,filename_model:str,filename_concentration:str,final_AP_file:str,final_DP_file:str,final_CS_file:str,final_SL_file:str,chronicle_writing:bool,steps_save:bool,rate_threshold_BP_auto:float) -> None:
    """run_cpa Running the Chemical Pathways Analysis

    _extended_summary_

    Parameters
    ----------
    timestep : float
        This is the timestep where the average chemical rates and concentrations are evaluated. Must be in the same time dimension/unit that the chemical rates.
    rate_threshold : float
        This is the minimum chemical rate considered as significant. Chemical pathways with rates lower than rate_threshold will be discarded. It can be override if the User specify a species of interest with BP_species.
    t_min : float
        This is the minimum lifetime for a chemical species to be considered. All species with a lifetime greater than t_min will be considered long-lived then not used as branching points. Must be in the same time dimension/unit that the chemical rates. This can be override using BP_species.
    BP_species : str
        This is where the User can specify one chemical species present in filename_concentration as a target. This means that the CPA will run until BP_species is used as a branching point. It overrides the value of t_min, t_min being reassign as the lifetime of BP_species.
    filename_model : str
        The name and local location of the file containing the chemical reactions and their rates.
    filename_concentration : str
        The name and local location of the file containing the chemical species and their concentration.
    final_AP_file : str
        _The name of the output JSON file containing the final Active Pathways
    final_DP_file : str
        _The name of the output JSON file containing the final Deleted Pathways
    final_CS_file : str
        _The name of the output JSON file containing the final Chemical System
    final_SL_file : str
        _The name of the output JSON file containing the final Species List
    chronicle_writing : bool
        Boolean flag that led the User save a file name chronicles.txt. This file contain each step done by the CPA algorithm and can help the User to better understand the inner working of the program or to debug it.
    steps_save : bool
        _description_
    rate_threshold_BP_auto : float
        If BP_species is specified, the User can override rate_threshold and use rate_threshold_BP_auto instead, as a percentage of BP_species total rate. Meaning that the minimal rate for a pathway considered will be expressed as a percentage of the total rate of BP_species.
    """

    # init global var
    init_global_var(chronicle_writing=chronicle_writing,steps_save=steps_save,rate_threshold_BP_auto=rate_threshold_BP_auto)

    # first test is to convert a given text file into a workable JSON dataset
    if global_var.chronicle_writing:
        with open('chronicles.txt', 'w') as output_file:
            output_file.write('Start of the Chemical Pathway Analysis')
            output_file.write('\n')

        o_tools.write_line_chronicle('######################')
        o_tools.write_line_chronicle('User Inputs Processing')
        o_tools.write_line_chronicle('######################')
        o_tools.write_line_chronicle('\n')
    
    chem_s = i_system.convert_chemical_reaction_file(filename=filename_model)
    chemical_species = i_concentration.convert_concentration_file(filename=filename_concentration,timestep=timestep,reaction_system=chem_s)
    chem_s = i_system.adding_pseudo_reactions(chemical_species=chemical_species,chemical_system=chem_s)

    # 2. We run the initialization
    if global_var.chronicle_writing:
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('#######################')
        o_tools.write_line_chronicle('Pathways Initialization')
        o_tools.write_line_chronicle('#######################')
        o_tools.write_line_chronicle('\n')

    active_p,deleted_p = p_init.init_pathways(chemical_system=chem_s)

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
    
    active_p,deleted_p,chemical_species,chem_s = ml.main_loop(t_min=t_min,f_min=rate_threshold,active_p=active_p,deleted_p=deleted_p,chemical_species=chemical_species,chemical_system=chem_s,BP_species=BP_species)

    # 4. main loop done. Outputs time!!!
    # out.text_output(target_species=target_species)

    # 5 copying results files
    d_tools.save_pathways_to_JSON(pathways=active_p,filename=final_AP_file)
    d_tools.save_pathways_to_JSON(pathways=deleted_p,filename=final_DP_file)
    d_tools.save_pathways_to_JSON(pathways=chemical_species,filename=final_SL_file)
    d_tools.save_pathways_to_JSON(pathways=chem_s,filename=final_CS_file)

def infos(timestep:float,t_min:float,filename_model:str,filename_concentration:str,chronicle_writing:bool,steps_save:bool) -> None:

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
    
    chem_s = i_system.convert_chemical_reaction_file(filename=filename_model)
    chemical_species = i_concentration.convert_concentration_file(filename=filename_concentration,timestep=timestep,reaction_system=chem_s)
    chem_s = i_system.adding_pseudo_reactions(chemical_species=chemical_species,chemical_system=chem_s)

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

    list_bp = bp.list_next_branching_points(t_min=t_min,chemical_species=chemical_species)
    print(list_bp)
    
    if global_var.chronicle_writing:
        o_tools.write_line_chronicle('\n')
        o_tools.write_line_chronicle('Here is the sorted list by lifetime of the next species considered as branching points for a fixed minimum timescale of '+'{:0.3e}'.format(t_min)+':')
        o_tools.write_line_chronicle(' '.join(list_bp))
        o_tools.write_line_chronicle('\n')
