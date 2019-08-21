from copy import deepcopy

def check_nist(equ, nist_list, max_C = 1, max_N = 2, max_S = 1, max_Si = 1,\
    O = -1, tot = -1):

    max_H = 5

    if tot == -1:
        tot = max([max_C, max_N, max_S, max_Si])

    ele_dict = {'C':0, 'H':0, 'Ar':0, 'He':0, 'Xe':0, 'Al':0, 'Kr':0,
        'Br':0, 'Cl':0, 'F':0, 'N':0, 'O':0, 'S':0, 'Si':0}

    ref = list(deepcopy(equ))

    for equa in equ:
        rmv = False
        #check if left and right are equal
        el = list(deepcopy(equa['left']))
        er = list(deepcopy(equa['right']))
        for i in equa['left']:
            for j in equa['right']:
                if i in el and j in er and i == j:
                    el.remove(i)
                    er.remove(j)
        if el == [] and er == []:
            rmv = True
        #check if species is allowed
        for i in equa['right']:
            if i[0] != {} and 'e' not in i[0].keys() and 'heavy' not in i[0].keys():
                Sum_ele = 0

                if i[1] > 1 or i[1] < -1:

                    rmv = True
                
                species = dict(deepcopy(ele_dict))
                for ele, quant in i[0].items():
                    species[ele] += quant

                    if ele in ['C', 'N', 'Si', 'S']:
                        Sum_ele += quant
                #Calculate the max O
                if O != -1:
                    max_O = O
                else:
                    max_O = int(quant * 1.5)
                
                #limit the species scales
                if species['N'] > max_N:
                    rmv = True

                if species['C'] > max_C:
                    rmv = True

                if species['Si'] > max_Si:
                    rmv = True

                if species['S'] > max_S:
                    rmv = True
                
                if species['O'] > max_O:
                    rmv = True 
                
                if species['C'] + species['S'] + species['N'] + species['Si'] > tot:
                    rmv = True

                #convert dictionary to str, and check if that in the nist database
                spe_str = ''
                spe_copy = deepcopy(species)
                for ele, quant in species.items():
                    if quant != 0:
                        if quant == 1:
                            spe_str += ele
                        else:
                            spe_str += ele + str(quant)
                    else:
                        del spe_copy[ele]
                
                if species['H'] > max_H and len(spe_copy.keys()) == 1:
                    rmv = True
                
                if i[1] == 1:
                    spe_str += '+'
                if i[1] == -1:
                    spe_str += '-'
                
                if spe_str not in nist_list:
                    rmv = True

        if rmv == True:
            ref.remove(equa)

    return ref

def check_negative(equations, all_compounds, full_dict_lowest):
    '''
    function to check formation energy of negative ions
    '''
    species_remove = []
    equations_new = deepcopy(equations)
    species_remain = deepcopy(all_compounds)

    # loop to calculate the formation energy from neutral to ion
    for i in all_compounds:
        if i['charge'] == -1:
            comps = ''
            for comp, num in i['compound'].items():
                if num != 1:
                    comps += str(comp) + str(num)
                else:
                    comps += str(comp)
            #add charges after species
            comp_neg = comps
            if i['charge'] > 0:
                if i['charge'] == 1:
                    comp_neg += '+'
                else:
                    comp_neg += '+' + str(i['charge'])
            elif i['charge'] < 0:
                if i['charge'] == -1:
                    comp_neg += '-'
                else:
                    comp_neg += '-' + str(-i['charge'])
            #Calculate the difference between negative ions and moleculars.
            if comps in full_dict_lowest.keys() and comp_neg in full_dict_lowest.keys():
                diff = full_dict_lowest[comp_neg] - full_dict_lowest[comps]
                if diff > 0:
                    species_remove.append([i['compound'], i['charge']])
                    species_remain.remove(i)
    
    # delete those equations containing forbidden species

    for equ in equations:
        rmv = False
        for species in equ['left']:
            if species in species_remove :
                rmv = True
        for species in equ['right']:
            if species in species_remove :
                rmv = True
        if rmv == True:
            equations_new.remove(equ)
    
    return species_remain, species_remove, equations_new
