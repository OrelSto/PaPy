# routines for writing stuff

import json

from io import TextIOWrapper

def writing_reaction(reaction:list,output_file:TextIOWrapper,chem_system_data:list):
    # Writing the reaction results in a given file
    mult = []
    if reaction["multiplicity"] > 1:
        mult.append(reaction["multiplicity"])
        output_file.write('(')
    else:
        output_file.write(' ')
    list_reactant = chem_system_data[reaction["index"]]["reactants"]
    list_product = chem_system_data[reaction["index"]]["products"]
    for reactant in list_reactant:
        if reactant["stoichiometry"] > 1:
            output_file.write(str(reactant["stoichiometry"])+ ' ' 
    +reactant["compound"])
        else:
            output_file.write(reactant["compound"])
        if (list_reactant.index(reactant) + 1) == len(list_reactant):
            output_file.write(' => ')
        else:
            output_file.write(' + ')
    for product in list_product:
        if product["stoichiometry"] > 1:
            output_file.write(str(product["stoichiometry"])+ ' ' 
    +product["compound"])
        else:
            output_file.write(product["compound"])
        if (list_product.index(product) + 1) == len(list_product):
            if mult:
                output_file.write(') x '+str(mult[0]))
        else:
            output_file.write(' + ')

def reaction_to_str(reaction:dict,chem_system_data:list):
    # Writing the reaction results in a given file
    result = ''
    mult = []
    if reaction["multiplicity"] > 1:
        mult.append(reaction["multiplicity"])
        result += '('
    else:
        result += ' '
    list_reactant = chem_system_data[reaction["index"]]["reactants"]
    list_product = chem_system_data[reaction["index"]]["products"]
    for reactant in list_reactant:
        if reactant["stoichiometry"] > 1:
            result += str(reactant["stoichiometry"])+ ' ' +reactant["compound"]
        else:
            result += reactant["compound"]
        if (list_reactant.index(reactant) + 1) == len(list_reactant):
            result += ' => '
        else:
            result += ' + '
    for product in list_product:
        if product["stoichiometry"] > 1:
            result += str(product["stoichiometry"])+ ' ' +product["compound"]
        else:
            result += product["compound"]
        if (list_product.index(product) + 1) == len(list_product):
            if mult:
                result += ') x '+str(mult[0])
        else:
            result += ' + '
    return result

def writing_pathway(pathway:list,output_file:TextIOWrapper,chem_system_data:list):

    for r in pathway["reactions"]:
        writing_reaction(reaction=r,output_file=output_file,chem_system_data=chem_system_data)
        output_file.write(' \n')

    output_file.write('------------------')
    output_file.write('\n')

    mask = []
    saving_prod = []
    saving_react = []
    for bp in pathway["branching points"]:
        if bp["stoichiometry"] > 0:
            mask.append(False)
            if bp["stoichiometry"] > 1:
                saving_prod += [str(bp["stoichiometry"])+ ' ' +bp["compound"]]
            else:
                saving_prod += [bp["compound"]]
        elif bp["stoichiometry"] == 0:
            mask.append(True)
        elif bp["stoichiometry"] < 0:
            mask.append(False)
            if bp["stoichiometry"] < -1:
                saving_react += [str(-bp["stoichiometry"])+ ' ' +bp["compound"]]
            else:
                saving_react += [bp["compound"]]
    if all(mask):
        output_file.write(' NULL ')
    else:
        # print('We finally write the reactants:',saving_react)
        # print('And write the products        :',saving_prod)
        output_file.write(' ' + ' + '.join(saving_react) + ' => ' + ' + '.join(saving_prod))

def pathway_to_str(pathway:list,chem_system_data:list):

    result = ''

    for r in pathway["reactions"]:
        result += reaction_to_str(reaction=r,chem_system_data=chem_system_data)
        result += ' \n'

    result += '------------------'
    result += '\n'

    mask = []
    saving_prod = []
    saving_react = []
    for bp in pathway["branching points"]:
        if bp["stoichiometry"] > 0:
            mask.append(False)
            if bp["stoichiometry"] > 1:
                saving_prod += [str(bp["stoichiometry"])+ ' ' +bp["compound"]]
            else:
                saving_prod += [bp["compound"]]
        elif bp["stoichiometry"] == 0:
            mask.append(True)
        elif bp["stoichiometry"] < 0:
            mask.append(False)
            if bp["stoichiometry"] < -1:
                saving_react += [str(-bp["stoichiometry"])+ ' ' +bp["compound"]]
            else:
                saving_react += [bp["compound"]]
    if all(mask):
        result += ' NULL '
    else:
        # print('We finally write the reactants:',saving_react)
        # print('And write the products        :',saving_prod)
        result += ' ' + ' + '.join(saving_react) + ' => ' + ' + '.join(saving_prod)
    
    return result

def write_line_chronicle(text_to_archive:str):
    # Write a a line in the chronicles
    with open('chronicles.txt', 'a') as output_file:
        if text_to_archive == '\n':
            output_file.write('\n')
        else:
            output_file.write('\n')
            output_file.write(text_to_archive)
