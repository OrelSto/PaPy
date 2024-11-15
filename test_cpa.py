# run in your local directory where the github repository is copied
#$ python3 -m test_cpa

from src.chem_pathway_analysis import cpa

# this is a stupid way to test the package but it works
# Next line is the example of lehmann 2004
# cpa.run_cpa(timestep=1.0,rate_threshold=1e-12,t_min=1.0e4,
#             target_species=['O2','O3','O'],
#             filename_model='user_model_example.txt',
#             filename_concentration='user_concentration_example.txt',
#             final_AP_file='ActP_Lehmann2004.json',
#             final_DP_file='DelP_Lehmann2004.json',
#             final_CS_file='ChemS_Lehmann2004.json',
#             final_SL_file='SpecL_Lehmann2004.json',
#             chronicle_writing=True)
# cpa.out.pie_output(target_species=['O2','O3','O'],act_P_json='ActP_Lehmann2004.json',del_P_json='DelP_Lehmann2004.json',chem_R_json='ChemS_Lehmann2004.json',spec_L_json='SpecL_Lehmann2004.json')

# test article to see if P4 is deleted!
# cpa.run_cpa(timestep=1.0,rate_threshold=0.3e-9,t_min=1.0e4,
#             target_species=['O2','O3','O'],
#             filename_model='user_model_example.txt',
#             filename_concentration='user_concentration_example.txt',
#             final_AP_file='ActP_Lehmann2004_P4del.json',
#             final_DP_file='DelP_Lehmann2004_P4del.json',
#             final_CS_file='ChemS_Lehmann2004_P4del.json',
#             final_SL_file='SpecL_Lehmann2004_P4del.json',
#             chronicle_writing=True)
# cpa.out.pie_output(target_species=['O2','O3','O'],act_P_json='ActP_Lehmann2004_P4del.json',del_P_json='DelP_Lehmann2004_P4del.json',chem_R_json='ChemS_Lehmann2004_P4del.json',spec_L_json='SpecL_Lehmann2004_P4del.json')

# Test for another simple chemical model that includes Cl and ClO O3 destruction catalysis
# This model is under dev since I have no clue of the real rates and proper concentration fo Cl and ClO
# cpa.run_cpa(timestep=100.0,rate_threshold=1e-11,t_min=1.0e4,
#             filename_model='user_model_O3destruction_example.txt',
#             filename_concentration='user_concentration_O3destruction_example.txt',
#             final_AP_file='ActP_O3_Cl_example.json',
#             final_DP_file='DelP_O3_Cl_example.json',
#             final_CS_file='ChemS_O3_Cl_example.json',
#             final_SL_file='SpecL_O3_Cl_example.json',
#             chronicle_writing=True)

# test with Uranus data!
cpa.run_cpa(timestep=1.0,rate_threshold=1.0e-18,t_min=1.0e14,
            target_species=['O2','O[3p]','CH4','CH3CHO','H2O','C2H4','HCO','[3]CH2'],
            filename_model='reactions_Uranus.txt',
            filename_concentration='concentrations_Uranus.txt',
            final_AP_file='ActP_Uranus.json',
            final_DP_file='DelP_Uranus.json',
            final_CS_file='ChemS_Uranus.json',
            final_SL_file='SpecL_Uranus.json',
            chronicle_writing=True)
cpa.out.pie_output(target_species=['O2','O[3p]','CH4','CH3','CH3CHO','H2O','C2H4','HCO','[3]CH2'],act_P_json='ActP_Uranus.json',del_P_json='DelP_Uranus.json',chem_R_json='ChemS_Uranus.json',spec_L_json='SpecL_Uranus.json')
cpa.out.table_Tex(target_species=['O2','O[3p]','CH4','CH3','CH3CHO','H2O','C2H4','HCO','[3]CH2'],act_P_json='ActP_Uranus.json',del_P_json='DelP_Uranus.json',chem_R_json='ChemS_Uranus.json',spec_L_json='SpecL_Uranus.json')

# Test with the data drom ChemPath module
# cpa.run_cpa(timestep=1.0,rate_threshold=0.02e-9,t_min=10.0,
#             target_species=['HO2'],
#             filename_model='user_model_example_article.txt',
#             filename_concentration='user_concentration_example_article.txt',
#             final_AP_file='ActP_Chempath_example.json',
#             final_DP_file='DelP_Chempath_example.json',
#             final_CS_file='ChemS_Chempath_example.json',
#             final_SL_file='SpecL_Chempath_example.json',
#             chronicle_writing=True)
# cpa.out.pie_output(target_species=['HO2','H2O','OH','H2O2'],act_P_json='ActP_Chempath_example.json',del_P_json='DelP_Chempath_example.json',chem_R_json='ChemS_Chempath_example.json',spec_L_json='SpecL_Chempath_example.json')
# cpa.out.table_Tex(target_species=['HO2','H2O','OH','H2O2'],act_P_json='ActP_Chempath_example.json',del_P_json='DelP_Chempath_example.json',chem_R_json='ChemS_Chempath_example.json',spec_L_json='SpecL_Chempath_example.json')