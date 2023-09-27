import json

# we make a moche output ^^'
def moche():
    with open('active_pathways.json', 'r') as active_pathways_file:
        # Parse the JSON data and store it in a variable
        active_pathways_data = json.load(active_pathways_file)

    with open('chemical_reaction_system.json', 'r') as chem_system_file:
        # Parse the JSON data and store it in a variable
        chem_system_data = json.load(chem_system_file)

    with open('output_moche.txt', 'w') as output_moche_file:
        for pathway in active_pathways_data:
            for r in pathway["reactions"]:
                list_reactant = chem_system_data[r["index"]]["reactants"]
                list_product = chem_system_data[r["index"]]["products"]
                for reactant in list_reactant:
                    if reactant["stoichiometry"] > 1:
                        output_moche_file.write(str(reactant["stoichiometry"])+ ' ' +reactant["compound"])
                    else:
                        output_moche_file.write(reactant["compound"])
                    if (list_reactant.index(reactant) + 1) == len(list_reactant):
                        output_moche_file.write(' => ')
                    else:
                        output_moche_file.write(' + ')

                for product in list_product:
                    if product["stoichiometry"] > 1:
                        output_moche_file.write(str(product["stoichiometry"])+ ' ' +product["compound"])
                    else:
                        output_moche_file.write(product["compound"])
                    if (list_product.index(product) + 1) == len(list_product):
                        output_moche_file.write(' \n')
                    else:
                        output_moche_file.write(' + ')

            output_moche_file.write('------------------')
            output_moche_file.write('\n')
            
            mask = []
            saving_prod = []
            for bp in pathway["branching points"]:
                if bp["stoichiometry"] > 1:
                    mask.append(False)
                    output_moche_file.write(str(bp["stoichiometry"])+ ' ' +bp["compound"])
                elif bp["stoichiometry"] == 0:
                    mask.append(True)
                elif bp["stoichiometry"] < 0:
                    mask.append(False)
                    if -bp["stoichiometry"] > 1:
                        saving_prod.append(str(-bp["stoichiometry"])+ ' ' +bp["compound"])
                    else:
                        saving_prod.append(bp["compound"])
                if (pathway["branching points"].index(bp) + 1) == len(pathway["branching points"]):
                    if saving_prod:
                        output_moche_file.write(' => ' + saving_prod[0])
                    output_moche_file.write(' \n')
                else:
                    if bp["stoichiometry"] != 0: 
                        output_moche_file.write(' + ')
                
            output_moche_file.write(' \n')


