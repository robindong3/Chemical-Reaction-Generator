from initialize import initialize
from copy import deepcopy
from check import *
from tqdm import tqdm
from processes import *
from CodeToGenerateEnthapyDatabase import *
import time
import operator
from decimal import Decimal
import json

def react(inputs, only_spe = True, location = '', max_C = 1, max_N = 2, max_S = 1, max_Si = 1,\
        O = -1, tot = -1):
    # Show the time taken by this program
    a = time.time()

    # load the enthalpy data
    enth = load_enthalpy()

    # load the manully input enthalpies
    with open("EnthalpyInsertManually.txt", "r") as file:
        star_enth = []
        for line in file:
            words = line.split()
            if len(words) < 3 and len(words) > 1:
                star_enth.append(words[0])

    all_compounds = []
    new_compounds = []

    # initializing
    j = 0
    while j < len(inputs):

        # Rerugulate all elements 
        upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        lower = 'abcdefghijklmnopqrstuvwxyz'
        num = '1234567890'

        comp = inputs[j] + '/'
        gas = True
        comp_dict = {}

        elem_dict = {'C':0, 'H':0, 'Ar':0, 'He':0, 'Xe':0, 'Al':0, 'Kr':0,
                'Br':0, 'Cl':0, 'F':0, 'N':0, 'O':0, 'S':0, 'Si':0}

        for i in range(len(comp)):
        #check elements
            if comp[i] in upper and comp[i+1] not in lower:
                elem = comp[i]
                #check number of elemnts
                if comp[i+1] in num and comp[i+2] in num:
                    number = int(comp[i+1] + comp[i+2])
                elif comp[i+1] in num:
                    number = int(comp[i+1])
                else:
                    number = 1
                
                #rewrite it into a dictionary
                if elem in comp_dict.keys():
                    comp_dict[elem] += number
                else:
                    comp_dict[elem] = number

            elif comp[i] in upper and comp[i+1] in lower:
                elem = comp[i] + comp[i+1]
                #check number of elemnts
                if comp[i+2] in num and comp[i+3] in num:
                    number = int(comp[i+2] + comp[i+3])
                elif comp[i+2] in num:
                    number = int(comp[i+2])
                else:
                    number = 1
                
                #rewrite it into a dictionary
                if elem in comp_dict.keys():
                    comp_dict[elem] += number
                else:
                    comp_dict[elem] = number

            #extract out charges
            elif comp[i] == '+':
                if comp[i+1] in num:
                    charge = int(comp[i+1])
                else:
                    charge = 1

            elif comp[i] == '-':
                if comp[i+1] in num:
                    charge = - int(comp[i+1])
                else:
                    charge = -1

            if '+' not in comp and '-' not in comp:
                charge = 0
            
        for elem_new in list(comp_dict.keys()):
            elem_dict[elem_new] += int(comp_dict[elem_new])

        for elem_new in list(elem_dict.keys()):
            if elem_dict[elem_new] == 0:
                del elem_dict[elem_new]
        
        # Append all initial inputs into compound list
        comp = initialize(elem_dict, charge)
        comp_dict = comp.result()

        all_compounds.append(comp_dict)
        new_compounds.append(comp_dict)
        j += 1

    c = initialize({'e': 1}, -1)
    e = c.result()

    #load the species from nist chemistry webbook data
    with open("species.txt", "r") as file:
        array = []
        print('load species from the nist webbook database:')
        for line in tqdm(file):
            words = line.split()
            if len(words) > 2:
                array.append(words[-2])

    #empty lists contain all equations
    equations = []
    ref_equations_dict = {}
    #new list saving other equations checked without nist data.
    equations_2 = []
    n = 0

    def convert_to_str(i):
        '''
        equation convert the dictionary form of species into str form
        '''
        elem_dict = {'C':0, 'H':0, 'Ar':0, 'He':0, 'Xe':0, 'Al':0, 'Kr':0,
            'Br':0, 'Cl':0, 'F':0, 'N':0, 'O':0, 'S':0, 'Si':0}

        if 'e' not in i['compound'].keys() and 'heavy' not in i['compound'].keys():

            for elem_new in list(i['compound'].keys()):
                elem_dict[elem_new] += int(i['compound'][elem_new])

            for elem_new in list(elem_dict.keys()):
                if elem_dict[elem_new] == 0:
                    del elem_dict[elem_new]

            ref_spe = ''
            for comp, num in elem_dict.items():
                if num != 1:
                    ref_spe += str(comp) + str(num)
                else:
                    ref_spe += str(comp)
            #add charges after species
            if i['charge'] > 0:
                if i['charge'] == 1:
                    ref_spe += '+'
                else:
                    ref_spe += '+' + str(i['charge'])
            elif i['charge'] < 0:
                if i['charge'] == -1:
                    ref_spe += '-'
                else:
                    ref_spe += '-' + str(-i['charge'])
        else:
            #create empty str
            ref_spe = ''
            for comp, num in i['compound'].items():
                if num != 1:
                    ref_spe += str(comp) + str(num)
                else:
                    ref_spe += str(comp)
            #add charges after species
            if i['charge'] > 0:
                if i['charge'] == 1:
                    ref_spe += '+'
                else:
                    ref_spe += '+' + str(i['charge'])
            elif i['charge'] < 0:
                if i['charge'] == -1:
                    ref_spe += '-'
                else:
                    ref_spe += '-' + str(-i['charge'])
        return ref_spe

    while new_compounds != [] and n < 100:
        new_compounds2 = []

        for i in new_compounds:
            #add new compounds into the new_equation dict keys
            #convert compound dict to str
            ref_spe = convert_to_str(i)
            ref_equations_dict[ref_spe] = {}

        for j in all_compounds:
            for k in all_compounds:
                ref_spe_2 = convert_to_str(j)
                ref_spe_3 = convert_to_str(k)
                ref_equations_dict[ref_spe_2][ref_spe_3] = {}
                for m in ['HAS', 'EIN', 'ERC', 'EDR', 'EDS', 'EDA', 'EDI',
                    'HGN', 'HCX', 'HIR', 'HIN', 'HMM', 'HDS', 'HED', 'HDC', 'HDN']:
                    ref_equations_dict[ref_spe_2][ref_spe_3][m] = []

        #renew the dict on each step.
        equations_dict = deepcopy(ref_equations_dict)

        for i in tqdm(new_compounds):

            new_equ = []

            #add new compounds into the new_equation dict keys
            #convert compound dict to str
            ref_spe = convert_to_str(i)
            #electron collisions:

            if i['charge'] == 0:
                new_equ += EIN(i)
                if i['type'] == 'single' or i['type'] == 'multiple':
                    new_equ += EDS(i)
                    new_equ += EDA(i)
                    new_equ += EDI(i)

            if i['charge'] > 0:
                new_equ += ERC(i)
                if i['type'] == 'single' or i['type'] == 'multiple':
                    new_equ += EDR(i)

            if i['charge'] < 0:
                if i['type'] == 'single' or i['type'] == 'multiple':
                    new_equ += HIR(i, j)

            for j in all_compounds:
                if j['charge'] == 0:
                    # neutral-neutral reactions
                    if i['charge'] == 0:
                        new_equ += HAS(i, j)
                        new_equ += HIR(i, j)
                        if i['type'] == 'single' or i['type'] == 'multiple':
                            new_equ += HDS(i, j)
                            new_equ += HIN(i, j)
                    # ion-neutral reactions
                    elif i['charge'] > 0:
                        new_equ += HCX(i, j)
                        new_equ += HIR(i, j)
                    elif i['charge'] < 0:
                        new_equ += HCX(j, i)
                        new_equ += HGN(i, j)
                        if i['type'] == 'single' or i['type'] == 'multiple':
                            new_equ += HED(i, j)
                # ion-ion reactions
                elif j['charge'] > 0:
                    if i['charge'] == 0:
                        new_equ += HCX(j, i)
                        new_equ += HIR(i, j)
                        if i['type'] == 'single' or i['type'] == 'multiple':
                            new_equ += HDC(i, j)
                    elif i['charge'] > 0:
                        continue
                    elif i['charge'] < 0:
                        new_equ += HIR(i, j)
                        new_equ += HMM(j, i)
                        if i['type'] == 'single' or i['type'] == 'multiple':
                            new_equ += HDN(i, j)
                elif j['charge'] < 0:
                    if i['charge'] == 0:
                        new_equ += HCX(i, j)
                        new_equ += HGN(j, i)
                        new_equ += HIR(i, j)
                    if i['charge'] < 0:
                        continue
                    if i['charge'] > 0:
                        new_equ += HIR(i, j)
                        new_equ += HMM(i, j)

            new_equ = check_nist(new_equ, array, max_C, max_N, max_S, max_Si,\
                    O, tot)

            for equ in new_equ:

                if equ not in equations:

                    ret = True
                    s = False
                    
                    # find which sub list should the equation lies on
                    elem_list = deepcopy(equ['left'])
                    for elem in equ['left']:
                        if 'e' in elem[0].keys() or 'heavy' in elem[0].keys():
                            elem_list.remove(elem)
                        elif elem != [i['compound'], i['charge']]:
                            species_dict = {'compound': elem[0], 'charge': elem[1]}
                            ref_spe_2 = convert_to_str(species_dict)

                    if len(elem_list) == 1:
                        ref_spe_2 = ref_spe

                    # taken same equations away, forbidden double counting of N+H and H+N:
                    for equ_old in equations_dict[ref_spe][ref_spe_2][equ['type']]:
                        equ_old_ref = dict(deepcopy(equ_old))
                        equ_ref = dict(deepcopy(equ))
                        for elem in equ_old['left']:
                            if elem in equ_ref['left']:
                                equ_old_ref['left'].remove(elem)
                                equ_ref['left'].remove(elem)
                        
                        for elem in equ_old['right']:
                            if elem in equ_ref['right']:
                                equ_old_ref['right'].remove(elem)
                                equ_ref['right'].remove(elem)

                        if equ_old_ref['left'] == [] and equ_ref['right'] == []:
                            ret = False
                    spe = False
                    for ir in equ['right']:
                        if ir[0] != {} and 'e' not in ir[0].keys() and 'heavy' not in ir[0].keys():
                            # set (N)(H) = (H)(N)
                            elem_dict = {'C':0, 'H':0, 'Ar':0, 'He':0, 'Xe':0, 'Al':0, 'Kr':0,
                                'Br':0, 'Cl':0, 'F':0, 'N':0, 'O':0, 'S':0, 'Si':0}
                            
                            for elem_new in list(ir[0].keys()):
                                elem_dict[elem_new] += int(ir[0][elem_new])

                            for elem_new in list(elem_dict.keys()):
                                if elem_dict[elem_new] == 0:
                                    del elem_dict[elem_new]

                            ncomp = initialize(elem_dict, ir[1])
                            ncomp_spe = ncomp.result()

                            if ncomp_spe not in all_compounds\
                                and ncomp_spe not in new_compounds2 and 'e' not in ir[0].keys():
                                new_compounds2.append(ncomp_spe)
                                spe = True

                    if not only_spe:
                        spe = True
                    if ret == True and spe == True:
                        #add equations into the big list.
                        equations.append(equ)
                        #add equations into sub lists
                        equations_dict[ref_spe][ref_spe_2][equ['type']].append(equ)
                        if ref_spe_2 != ref_spe:
                            equations_dict[ref_spe_2][ref_spe][equ['type']].append(equ)

                    if ret == True and spe == True and s == True:
                        equations_2.append(equ)

        all_compounds += deepcopy(new_compounds2)
        new_compounds = deepcopy(new_compounds2)
        if new_compounds != []:
            print('New species discoverd in this iteration are:')
        for iii in new_compounds:
            print(iii['compound'], iii['charge'])
        n += 1

    if not only_spe:
        # write result dictionary into json file
        with open('results/' + location + 'EquationDict.json', 'w') as f:
            for i in equations:
                json.dump(i, f)
                f.write('\n')

        #write result equation into txt file
        with open('results/' + location + 'EquationStr.txt', 'w') as f:
            for i in equations:
                equation = ''
                equation += str(i['type'])
                equation += ': '
                for left in i['left']:
                    for comp, num in left[0].items():
                        equation += '(' + comp + ')' + str(num)
                    equation += '[' + str(left[1]) + ']' + ' + '
                equation = equation[:-3]
                equation += '  ==  '
                for right in i['right']:
                    for comp, num in right[0].items():
                        equation += '(' + comp + ')' + str(num)
                    equation += '[' + str(right[1]) + ']' + ' + '
                equation = equation[:-3]

                f.write(equation)
                f.write('\n')

    #write species into a .txt file
    with open('results/' + location + 'SpiciesList2.txt', 'w') as f:

        # save the enthalpies into a .txt file

        full_dict = {}
        full_dict_lowest = {}
        comp_str = []
        for i in all_compounds:
            comps = ''
            for comp, num in i['compound'].items():
                if num != 1:
                    comps += str(comp) + str(num)
                else:
                    comps += str(comp)
            #add charges after species
            if i['charge'] > 0:
                if i['charge'] == 1:
                    comps += '+'
                else:
                    comps += '+' + str(i['charge'])
            elif i['charge'] < 0:
                if i['charge'] == -1:
                    comps += '-'
                else:
                    comps += '-' + str(-i['charge'])

            enth_str = ''
            enth_list = []

            # sort the enthalpies from low to high
            for j in enth:
                if i['compound'] == j[0] and i['charge'] == j[1]:
                    enth_list.append(j[2])

            enth_list.sort(reverse = False)

            #sort the species list by their enthalpies
            full_dict[comps] = enth_list

            if enth_list != []:
                full_dict_lowest[comps] = enth_list[0]
        

        full_dict_lowest = dict(sorted(full_dict_lowest.items(),
            key=operator.itemgetter(1)))

        # move out those unexpected species
        species_remain, species_remove, equations_new = check_negative(
            equations, all_compounds, full_dict_lowest)
        
        for i in species_remove:
            comps = ''
            for comp, num in i[0].items():
                if num != 1:
                    comps += str(comp) + str(num)
                else:
                    comps += str(comp)
            #add charges after species
            if i[1] > 0:
                if i[1] == 1:
                    comps += '+'
                else:
                    comps += '+' + str(i[1])
            elif i[1] < 0:
                if i[1] == -1:
                    comps += '-'
                else:
                    comps += '-' + str(-i[1])

            # save the str form into a list
            comp_str.append(comps)
            
            enth_list = full_dict[comps]

            enth_str = ''
            for j in enth_list:
                enth_str += '     ' + str(j)

            # mark those manully inputed enthalpies
            
            if comps in star_enth:
                enth_str += '*'

            comps += enth_str
            f.write(comps)
            f.write('\n')
    
    #Write those removed species into a file
    with open('results/' + location + 'SpiciesList1.txt', 'w') as f:

        for key, elem in full_dict_lowest.items():
            if key not in comp_str:
                star = False
                if key in star_enth:
                    star = True
                key += '     ' + str(elem)
                # mark those manully inputed enthalpies
                if star:
                    key += '*'
                f.write(key)
                f.write('\n')

        # write those species have no enthalpies
        for key, value in full_dict.items():
            if value == []:
                f.write(key)
                f.write('\n')
    
    #write result equation into txt file
    if not only_spe:
        with open('results/' + location + 'EquationStr2.txt', 'w') as f:
            for i in equations_new:
                equation = ''
                equation += str(i['type'])
                equation += ': '
                for left in i['left']:
                    for comp, num in left[0].items():
                        equation += '(' + comp + ')' + str(num)
                    equation += '[' + str(left[1]) + ']' + ' + '
                equation = equation[:-3]
                equation += '  ==  '
                for right in i['right']:
                    for comp, num in right[0].items():
                        equation += '(' + comp + ')' + str(num)
                    equation += '[' + str(right[1]) + ']' + ' + '
                equation = equation[:-3]

                f.write(equation)
                f.write('\n')

    # Show the time taken by this program
    b = time.time()
    print('Total time taken:', b - a, ' seconds')

    #sort the reactions with there enthalpy change:
    if not only_spe:
        with open('results/' + location + 'EquationWithEnth.txt', 'w') as f:
            n = 0
            reaction_enth_dict = {}
            for equ in equations_new:
                left_enth = 0
                have_enth = True
                for i in equ['left']:
                    #search the enthalpy in the enthalpy list
                    comps = ''
                    if 'e' not in i[0].keys():
                        for comp, num in i[0].items():
                            if num != 1:
                                comps += str(comp) + str(num)
                            else:
                                comps += str(comp)
                        #add charges after species
                        if i[1] > 0:
                            if i[1] == 1:
                                comps += '+'
                            else:
                                comps += '+' + str(i['charge'])
                        elif i[1] < 0:
                            if i[1] == -1:
                                comps += '-'
                            else:
                                comps += '-' + str(-i['charge'])
                        # add all enthalpies on LHS        
                        if comps in full_dict_lowest.keys():
                            left_enth += full_dict_lowest[comps]
                        else:
                            have_enth = False
                
                # Calculate the RHS enthalpies
                right_enth = 0
                if have_enth == True:
                    for i in equ['right']:
                        if 'e' not in i[0].keys():
                            #search the enthalpy in the enthalpy list
                            comps = ''
                            for comp, num in i[0].items():
                                if num != 1:
                                    comps += str(comp) + str(num)
                                else:
                                    comps += str(comp)
                            #add charges after species
                            if i[1] > 0:
                                if i[1] == 1:
                                    comps += '+'
                                else:
                                    comps += '+' + str(i['charge'])
                            elif i[1] < 0:
                                if i[1] == -1:
                                    comps += '-'
                                else:
                                    comps += '-' + str(-i['charge'])
                            # add all enthalpies on LHS        
                            if comps in full_dict_lowest.keys():
                                right_enth += full_dict_lowest[comps]
                            else:
                                have_enth = False
                
                # Calculate the enthalpy of reaction
                if have_enth == True:
                    reaction_enth_dict[str(n)] = right_enth - left_enth
                if have_enth == False:
                    reaction_enth_dict[str(n)] = -99999
                n += 1

            # Sort all enthalpies
            reaction_enth_dict = dict(sorted(reaction_enth_dict.items(),
                key=operator.itemgetter(1)))

            # Write all reactions with enthalpies into the file
            for num, enthalpy in reaction_enth_dict.items():
                i = equations_new[int(num)]
                equation = ''
                equation += str(i['type'])
                equation += ': '
                # add all LHS species
                for left in i['left']:
                    for comp, num in left[0].items():
                        equation += '(' + comp + ')' + str(num)
                    equation += '[' + str(left[1]) + ']' + ' + '
                equation = equation[:-3]
                equation += '  ==  '
                #add all RHS species
                for right in i['right']:
                    for comp, num in right[0].items():
                        equation += '(' + comp + ')' + str(num)
                    equation += '[' + str(right[1]) + ']' + ' + '
                equation = equation[:-3]
                if enthalpy != -99999:
                    x = Decimal(enthalpy)
                    f.write(str(round(x,2)) + '       ')
                f.write(equation)
                f.write('\n')