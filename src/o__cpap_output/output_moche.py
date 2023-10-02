import json

# we make a moche output ^^'
def moche():
    with open('active_pathways.json', 'r') as active_pathways_file:
        # Parse the JSON data and store it in a variable
        active_pathways_data = json.load(active_pathways_file)

    with open('deleted_pathways.json', 'r') as deleted_pathways_file:
        # Parse the JSON data and store it in a variable
        deleted_pathways_data = json.load(deleted_pathways_file)

    with open('chemical_reaction_system.json', 'r') as chem_system_file:
        # Parse the JSON data and store it in a variable
        chem_system_data = json.load(chem_system_file)

    with open('output_moche.txt', 'w') as output_moche_file:
        output_moche_file.write('**************************')
        output_moche_file.write('\n')

        rate_sum = 0.0
        for pathway in active_pathways_data:
            rate_sum += pathway["rate"]
        pathway_sorted = {}
        i = 0
        for pathway in active_pathways_data:
            pathway_sorted.update({i:pathway["rate"]/rate_sum * 100})
            i += 1
        ind_pathway_sorted = sorted(pathway_sorted,key=pathway_sorted.get,reverse=True)

        rate_deleted = 0.0
        for pathway in deleted_pathways_data:
            rate_sum += pathway["rate"]
            rate_deleted += pathway["rate"]
        
        # for pathway in active_pathways_data:
        for i in ind_pathway_sorted:
            moche_writing_pathway(pathway=active_pathways_data[i],output_moche_file=output_moche_file,chem_system_data=chem_system_data,rate_sum=rate_sum)
            
        
        # Now the rate from deleted pathways
        output_moche_file.write(' RATE DELETED  : ' + '{:0.3e}'.format(rate_deleted))
        output_moche_file.write(' \n')
        output_moche_file.write(' RATE DELETED %: ' + '{:0.3f}'.format(rate_deleted/rate_sum * 100))
        output_moche_file.write(' \n')

def moche_writing_pathway(pathway:list,output_moche_file,chem_system_data,rate_sum:float):
    for r in pathway["reactions"]:
        mult = []
        if r["multiplicity"] > 1:
            mult.append(r["multiplicity"])
            output_moche_file.write('(')
        else:
            output_moche_file.write(' ')

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
                if mult:
                    output_moche_file.write(') x '+str(mult[0]))
                output_moche_file.write(' \n')
            else:
                output_moche_file.write(' + ')
            


    output_moche_file.write('------------------')
    output_moche_file.write('\n')
    
    mask = []
    saving_prod = []
    saving_react = []
    for bp in pathway["branching points"]:
        if bp["stoichiometry"] > 0:
            mask.append(False)
            if bp["stoichiometry"] > 1:
                saving_prod.append(str(bp["stoichiometry"])+ ' ' +bp["compound"])
            else:
                saving_prod.append(bp["compound"])
        elif bp["stoichiometry"] == 0:
            mask.append(True)
        elif bp["stoichiometry"] < 0:
            mask.append(False)
            if -bp["stoichiometry"] > 1:
                saving_react.append(str(-bp["stoichiometry"])+ ' ' +bp["compound"])
            else:
                saving_react.append(bp["compound"])
    if all(mask):
        output_moche_file.write(' NULL ')
        output_moche_file.write(' \n')
    else:
        output_moche_file.write(saving_react[0] + ' => ' + saving_prod[0])
        output_moche_file.write(' \n')
    
    # Now the rate
    output_moche_file.write(' \n')
    # output_moche_file.write(' RATE  : ' + str(round(pathway["rate"],3)))
    output_moche_file.write(' RATE  : ' + '{:0.3e}'.format(pathway["rate"]))
    output_moche_file.write(' \n')
    output_moche_file.write(' RATE %: ' + '{:0.3f}'.format(pathway["rate"]/rate_sum * 100))
    output_moche_file.write(' \n')
    

    output_moche_file.write(' \n')
    output_moche_file.write('**************************')
    output_moche_file.write('\n')


