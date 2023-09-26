import json

from p__initialization import init_pathways as p_init

# This is where we analyze the pathway to determine if it is a combination of subpathways

def subpathway_analysis(pathway:dict):

    # we initliaze the set of sub-pathways to the individual reactions present in the pathway
    set_SP = sub_pathway_set_init(pathway=pathway)

    # saving
    save_subpathways_to_JSON(set_SP=set_SP,filename='subpathways_tmp.json')
    return

def sub_pathway_set_init(pathway:dict):
    # we open the chemical reaction file to retrieve the actual reaction
    # Opening JSON file
    cs = open('chemical_reaction_system.json')
    # returns JSON object as a dictionary
    chemical_system = json.load(cs)

    set_SP = []
    for r in pathway["reactions"]:
        reaction = chemical_system[r["index"]]
        print('reaction index',r["index"])
        # we need to record each individual reaction into a new sub-pathway set
        sub_pathway = {}
        sub_pathway = p_init.format_first_pathway(reaction=reaction,index=r["index"])

        # Now we have the sub-pathway as an individual reaction, I need to update it a bit:
        # its rate
        # its multiplicity
        # watch out here, reaction = r and not variable reaction ! We want the multiplicity of the reaction r in pathway !!!
        update_subpathway(sub_pathway=sub_pathway,reaction=r,pathway=pathway)
        set_SP.append(sub_pathway)
    
    return set_SP

def update_subpathway(sub_pathway:dict,reaction:dict,pathway:dict):
    sub_pathway_reaction = sub_pathway["reactions"][0]
    # setting the multiplicity
    sub_pathway_reaction["multiplicity"] = reaction["multiplicity"]
    sub_pathway["rate"] = 0.0


def save_subpathways_to_JSON(set_SP:list,filename:str):
    # need a doc here?
    # Write the JSON data to an output file
    with open(filename, 'w') as output_file:
        json.dump(set_SP, output_file, indent=2)
    print("Pathways saved as",filename)