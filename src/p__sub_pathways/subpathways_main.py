import json
from p__sub_pathways import subpathway_analysis as sub

# This routine is set to initialize the subpathway set for each active pathway defined early on
# Steps:
# 1. reading the active_pathways.json
# 2. looping through the pathways and initialze the set of sub-pathways

def main_subpathways():
    # Open the JSON file and read its contents
    with open('active_pathways.json', 'r') as active_pathways_file:
        # Parse the JSON data and store it in a variable
        active_pathways_data = json.load(active_pathways_file)
    
    for pathway in active_pathways_data:
        # for each pathway we run the subpathway analysis
        sub.subpathway_analysis(pathway=pathway)
