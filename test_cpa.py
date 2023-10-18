# run in your local directory
#$ python3 -m test_cpa

from src.chem_pathway_analysis import cpa

# this is a stupid way to test the package and stupid values for inputs
cpa.run_cpa(timestep=100.0,rate_threshold=1e-12,t_min=1.0e4,filename_model='user_model_example.txt',filename_concentration='user_concentration_example.txt',chronicle_writing=True)
# cpa.run_cpa(timestep=100.0,rate_threshold=1e-10,t_min=1.0e4,filename_model='user_model_O3destruction_example.txt',filename_concentration='user_concentration_O3destruction_example.txt',chronicle_writing=True)
# test article to see if P4 is deleted!
# cpa.run_cpa(timestep=100.0,rate_threshold=0.3e-9,t_min=1.0e4)
