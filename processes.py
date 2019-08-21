from initialize import initialize
from copy import deepcopy
from tqdm import tqdm

c = initialize({'e': 1}, -1)
e = c.result()

# define all type of processes
def HAS(a, b):
    '''
    assosication
    '''
    new = dict(a['compound'])
    for key, quant in b['compound'].items():
        if key in a['compound'].keys():
            new[key] += quant
        else:
            new[key] = quant

    equa = {'left':[[a['compound'], a['charge']], [b['compound'], b['charge']]],
           'right': [[new, a['charge'] + b['charge']]], 'type':'HAS'}
    return [equa]

def EIN(a, e = e):
    '''
    ionization
    '''
    new = dict(a['compound'])
    equa = {'left':[[a['compound'], a['charge']], [e['compound'], e['charge']]],
            'right': [[a['compound'], a['charge'] + 1],
                     [e['compound'], e['charge']], [e['compound'], e['charge']]], 
            'type':'EIN'}
    return [equa]

def ERC(a, e = e):
    '''
    recombination
    '''
    equa = {'left':[[e['compound'], e['charge']], [a['compound'], a['charge']]],
            'right':[[a['compound'], a['charge'] - 1]], 'type':'ERC'}
    return [equa]

def calculate_right_H(H_num, type_ele, a, i):
    '''
    calculate RHS of chemical equation if H is in LHS
    '''
    H_remain = H_num - i - 1
    remain = dict(a['compound'])
    remain[type_ele] = H_remain
    if H_remain == 0:
        del remain[type_ele]
    #number of double Hs
    H2 = (i + 1) // 2
    if i <= 1:
        H2 = 0
    #number of single Hs
    H1 = i + 1 - H2 * 2
    right = []
    if remain != {}:
        right.append([remain, 0])
    for j in range(H2):
        right.append([{type_ele:2}, 0])
    for j in range(H1):
        right.append([{type_ele:1}, 0])
    return right

def calculate_split(N_num, num, a, elem):
    right1 = dict(deepcopy(a['compound']))
    right2 = dict(deepcopy(a['compound']))
    right1[elem] -= N_num
    right2[elem] -= (num - N_num)
    #find all posiable N, H combinations
    for elem2, num2 in a['compound'].items():
        if elem2 != elem:
            for H_num in range(num2 + 1):
                right11 = dict(deepcopy(right1))
                right22 = dict(deepcopy(right2))
                right11[elem2] -= H_num
                right22[elem2] -= num2 - H_num
                #delete empty dictionaries
                for elem3 in right1.keys():
                    if right11[elem3] == 0:
                        del right11[elem3]
                for elem3 in right2.keys():
                    if right22[elem3] == 0:
                        del right22[elem3]
    return right11, right22

def EDR(a, e = e):
    '''
    dissociative recombination
    '''
    tot = []
    if 'H' in a['compound'].keys() or 'F' in a['compound'].keys() or 'Cl' in a['compound'].keys(): 
        for elem, num in a['compound'].items():
            #the case the compound only release Hs or single elements
            if elem == 'H' or elem == 'F' or elem == 'Cl':
                H_num = int(a['compound'][elem])
                for i in range(H_num):
                    right = calculate_right_H(H_num, elem, a, i)
                    equa = {'left':[[e['compound'], e['charge']], [a['compound'], a['charge']]],
                            'right':right, 'type':'EDR'}
                    tot.append(equa)
            #the case the compound is slpit
            else:
                if len(list(a['compound'].keys())) == 2:
                    for N_num in range(num):
                        if N_num != 0:
                            right11, right22 = calculate_split(N_num, num, a, elem)
                            equa = {'left':[[a['compound'], a['charge']]],
                                    'right':[[right11, 0],
                                            [right22, 0]], 'type':'EDR'}
                            tot.append(equa)
            
    else:
        equa = {'left':[[e['compound'], e['charge']], [a['compound'], a['charge']]],
                'right':[[a['positive'], a['charge'] - 1], [a['negative'], 0]], 'type':'EDR'}
        return [equa]
    return tot


def EDS(a, e = e):
    '''
    dissociation
    '''
    tot = []
    if 'H' in a['compound'].keys() or 'F' in a['compound'].keys() or 'Cl' in a['compound'].keys():
        # if H F Cl in the compound, the reaction most likly happen on those elements
        for elem, num in a['compound'].items():
            if elem == 'H' or elem == 'F' or elem == 'Cl':
                H_num = int(a['compound'][elem])
                for i in range(H_num):
                    right = calculate_right_H(H_num, elem, a, i)
                    right.append([e['compound'], e['charge']])
                    equa = {'left':[[e['compound'], e['charge']], [a['compound'], a['charge']]],
                            'right':right, 'type':'EDS'}
                    tot.append(equa)
            else:
                if len(list(a['compound'].keys())) == 2:
                    for N_num in range(num):
                        if N_num != 0:
                            right11, right22 = calculate_split(N_num, num, a, elem)
                            equa = {'left':[[e['compound'], e['charge']],
                                            [a['compound'], a['charge']]],
                                    'right':[[e['compound'], e['charge']],
                                             [right11, 0],
                                             [right22, 0]], 'type':'EDS'}
                            tot.append(equa)
                
    else:
        equa = {'left':[[e['compound'], e['charge']], [a['compound'], a['charge']]],
                'right':[[e['compound'], e['charge']],
                         [a['positive'], a['charge']], [a['negative'], 0]],
                'type':'EDS'}
        return [equa]
    return tot

def EDA(a, e = e):
    '''
    dissociative attachment
    '''
    tot = []
    if 'H' in a['compound'].keys() or 'F' in a['compound'].keys() or 'Cl' in a['compound'].keys()\
        or 'Br' in a['compound'].keys():
        for elem, num in a['compound'].items():
            if elem == 'H' or elem == 'F' or elem == 'Cl' or elem == 'Br':
                H_num = int(a['compound'][elem])
                for i in range(H_num):
                    right = calculate_right_H(H_num, elem, a, i)
                    right[-1][-1] -= 1 
                    equa = {'left':[[e['compound'], e['charge']], [a['compound'], a['charge']]],
                            'right':right, 'type':'EDA'}
                    tot.append(equa)
            else:
                if len(list(a['compound'].keys())) == 2:
                    for N_num in range(num):
                        if N_num != 0:
                            right11, right22 = calculate_split(N_num, num, a, elem)
                            equa = {'left':[[e['compound'], e['charge']],
                                            [a['compound'], a['charge']]],
                                    'right':[[right11, 0],
                                             [right22, -1]], 'type':'EDA'}
                            tot.append(equa)
                            equa = {'left':[[e['compound'], e['charge']],
                                            [a['compound'], a['charge']]],
                                    'right':[[right11, -1],
                                             [right22, 0]], 'type':'EDA'}
                            tot.append(equa)
    else:
        equa = {'left':[[e['compound'], e['charge']], [a['compound'], a['charge']]],
                'right':[[a['positive'], a['charge']], [a['negative'], -1]], 'type':'EDA'}
        return [equa]
    return tot

def EDI(a, e = e):
    '''
    dissociative ionization
    '''
    tot = []
    if 'H' in a['compound'].keys() or 'F' in a['compound'].keys() or 'Cl' in a['compound'].keys()\
        or 'Br' in a['compound'].keys():
        for elem, num in a['compound'].items():
            if elem == 'H' or elem == 'F' or elem == 'Cl':
                H_num = int(a['compound'][elem])
                for i in range(H_num):
                    right = calculate_right_H(H_num, elem, a, i)
                    right[0][-1] += 1
                    right.append([e['compound'], e['charge']])
                    right.append([e['compound'], e['charge']])
                    equa = {'left':[[e['compound'], e['charge']], [a['compound'], a['charge']]],
                            'right':right, 'type':'EDI'}
                    tot.append(equa)
            else:
                if len(list(a['compound'].keys())) == 2:
                    for N_num in range(num):
                        if N_num != 0:
                            right11, right22 = calculate_split(N_num, num, a, elem)
                            equa = {'left':[[e['compound'], e['charge']],
                                            [a['compound'], a['charge']]],
                                    'right':[[e['compound'], e['charge']],
                                             [e['compound'], e['charge']],
                                             [right11, 1],
                                             [right22, 0]], 'type':'EDI'}
                            tot.append(equa)
                            equa = {'left':[[e['compound'], e['charge']],
                                            [a['compound'], a['charge']]],
                                    'right':[[e['compound'], e['charge']],
                                             [e['compound'], e['charge']],
                                             [right11, 0],
                                             [right22, 1]], 'type':'EDI'}
                            tot.append(equa)
    else:
        equa = {'left':[[e['compound'], e['charge']], [a['compound'], a['charge']]],
                'right':[[a['positive'], a['charge'] + 1], [a['negative'], 0],
                         [e['compound'], e['charge']], [e['compound'], e['charge']]],
                'type':'EDI'}
        return [equa]
    return tot

def HGN(a, b):
    '''
    association and electron detachment
    '''
    new = dict(a['compound'])
    for key, quant in b['compound'].items():
        if key in a['compound'].keys():
            new[key] += quant
        else:
            new[key] = quant

    equa = {'left':[[a['compound'], a['charge']], [b['compound'], b['charge']]],
           'right': [[new, a['charge'] + b['charge'] + 1], [e['compound'], e['charge']]],
           'type':'HGN'}
    return [equa]

def HCX(a, b):
    '''
    charge transfer
    '''

    equa = {'left':[[a['compound'], a['charge']], [b['compound'], b['charge']]],
           'right': [[a['compound'], a['charge'] - 1], [b['compound'], b['charge'] + 1]],
           'type':'HCX'}
    return [equa]

def HIR(a, b):
    '''
    heavy-particle interchange
    '''
    tot = []
    for key, quant in a['compound'].items():
        if key == 'H' or key == 'Cl' or key == 'F' or key == 'O' or key == 'Br':
            for giv_num in range(quant + 1):
                new_a = dict(deepcopy(a))
                new_b = dict(deepcopy(b))
                if giv_num != 0:
                    # the compound cannot give all his elements to another (interchange)
                    if len(a['compound'].keys()) == 1:
                        if giv_num != quant:
                            if key in new_b['compound'].keys():
                                new_b['compound'][key] += giv_num
                            else:
                                new_b['compound'][key] = giv_num
                            new_a['compound'][key] -= giv_num
                            if new_a['compound'][key] == 0:
                                del new_a['compound'][key]
                    else:
                        if key in new_b['compound'].keys():
                            new_b['compound'][key] += giv_num
                        else:
                            new_b['compound'][key] = giv_num

                        new_a['compound'][key] -= giv_num
                        if new_a['compound'][key] == 0:
                                del new_a['compound'][key]
                    if a['charge'] == -b['charge']:
                        a_charge = 0
                        b_charge = 0
                    else:
                        a_charge = int(a['charge'])
                        b_charge = int(b['charge'])
                    equa = {'left':[[a['compound'], a['charge']], [b['compound'], b['charge']]],
                            'right': [[new_a['compound'], a_charge],
                                      [new_b['compound'], b_charge]],
                            'type':'HIR'}
                    tot.append(equa)

    for key, quant in b['compound'].items():
        if key == 'H' or key == 'Cl' or key == 'F' or key == 'O' or key == 'Br':
            for giv_num in range(quant + 1):
                new_a = dict(deepcopy(a))
                new_b = dict(deepcopy(b))
                if giv_num != 0:
                    # the compound cannot give all his elements to another (interchange)
                    if len(b['compound'].keys()) == 1:
                        if giv_num != quant:
                            if key in new_a['compound'].keys():
                                new_a['compound'][key] += giv_num
                            else:
                                new_a['compound'][key] = giv_num
                            new_b['compound'][key] -= giv_num
                            if new_b['compound'][key] == 0:
                                del new_b['compound'][key]
                    else:
                        if key in new_a['compound'].keys():
                            new_a['compound'][key] += giv_num
                        else:
                            new_a['compound'][key] = giv_num
                        new_b['compound'][key] -= giv_num
                        if new_b['compound'][key] == 0:
                                del new_b['compound'][key]

                    if a['charge'] == -b['charge']:
                        a_charge = 0
                        b_charge = 0
                    else:
                        a_charge = int(a['charge'])
                        b_charge = int(b['charge'])

                    equa = {'left':[[a['compound'], a['charge']], [b['compound'], b['charge']]],
                            'right': [[new_a['compound'], a_charge],
                                      [new_b['compound'], b_charge]],
                            'type':'HIR'}
                    tot.append(equa)


    return tot

def HIN(a, b, e = e):
    '''
    heavy particle collisional ionization
    '''
    equa = {'left':[[a['compound'], a['charge']], [b['compound'], b['charge']]],
            'right': [[a['compound'], a['charge'] + 1], [b['compound'], b['charge']],
                      [e['compound'], e['charge']]],
            'type':'HIN'}
    return [equa]

def HMM(a, b):
    '''
    ions recombination
    '''
    equa = {'left':[[a['compound'], a['charge']], [b['compound'], b['charge']]],
            'right':[[a['compound'], a['charge'] - 1], [b['compound'], b['charge'] + 1]],
            'type':'HMM'}
    return [equa]

def HDS(a, b):
    '''
    heavy particle collisional dissociation
    '''
    tot = []
    if 'H' in a['compound'].keys() or 'F' in a['compound'].keys() or 'Cl' in a['compound'].keys():
        for elem, num in a['compound'].items():
            if elem == 'H' or elem == 'F' or elem == 'Cl':
                H_num = int(a['compound'][elem])
                for i in range(H_num):
                    right = calculate_right_H(H_num, elem, a, i)
                    right.append([b['compound'], 0])
                    equa = {'left':[[b['compound'], 0], [a['compound'], a['charge']]],
                            'right':right, 'type':'HDS'}
                    tot.append(equa)
            else:
                if len(list(a['compound'].keys())) == 2:
                    for N_num in range(num):
                        if N_num != 0:
                            right11, right22 = calculate_split(N_num, num, a, elem)
                            equa = {'left':[[b['compound'], 0],
                                            [a['compound'], a['charge']]],
                                    'right':[[b['compound'], 0],
                                             [right11, 0],
                                             [right22, 0]], 'type':'HDS'}
                            tot.append(equa)
    else:
        equa = {'left':[[b['compound'], 0], [a['compound'], a['charge']]],
                'right':[[b['compound'], 0], [a['positive'], a['charge']], [a['negative'], 0]], 'type':'HDS'}
        return [equa]
    return tot

def HDN(a, b):
    '''
    heavy particle dissociative neutralization
    '''
    tot = []
    if 'H' in a['compound'].keys() or 'F' in a['compound'].keys() or 'Cl' in a['compound'].keys():
        for elem, num in a['compound'].items():
            if elem == 'H' or elem == 'F' or elem == 'Cl' or elem == 'O':
                H_num = int(a['compound'][elem])
                for i in range(H_num):
                    right = calculate_right_H(H_num, elem, a, i)
                    equa = {'left':[[b['compound'], b['charge']], [a['compound'], a['charge']]],
                            'right':[[b['compound'], 0]] + right, 'type':'HDN'}
                    tot.append(equa)
            else:
                if len(list(a['compound'].keys())) == 2:
                    for N_num in range(num):
                        if N_num != 0:
                            right11, right22 = calculate_split(N_num, num, a, elem)
                            equa = {'left':[[b['compound'], 0],
                                            [a['compound'], a['charge']]],
                                    'right':[[b['compound'], 0],
                                             [right11, 0],
                                             [right22, 0]], 'type':'HDN'}
                            tot.append(equa)
    else:
        equa = {'left':[[b['compound'], 0], [a['compound'], a['charge']]],
                'right':[[b['compound'], 0], [a['positive'], a['charge']], [a['negative'], 0]], 'type':'HDN'}
        return [equa]
    return tot

def HDC(a, b):
    '''
    heavy particle dissociation and charge transfer
    '''
    tot = []
    if 'H' in a['compound'].keys() or 'F' in a['compound'].keys() or 'Cl' in a['compound'].keys()\
        or 'O' in a['compound'].keys():
        if a['charge'] == 0:
#            b = {'compound':{'heavy':0}, 'charge': 1}
            for elem, num in a['compound'].items():
                if elem == 'H' or elem == 'F' or elem == 'Cl':
                    H_num = int(a['compound'][elem])
                    for i in range(H_num):
                        right = calculate_right_H(H_num, elem, a, i)
                        right[0][-1] = 1
                        equa = {'left':[[b['compound'], b['charge']], [a['compound'], a['charge']]],
                                'right':[[b['compound'], 0]] + right, 'type':'HDC'}
                        tot.append(equa)
                else:
                    if len(list(a['compound'].keys())) == 2:
                        for N_num in range(num):
                            if N_num != 0:
                                right11, right22 = calculate_split(N_num, num, a, elem)
                                equa = {'left':[[b['compound'], 0],
                                                [a['compound'], a['charge']]],
                                        'right':[[b['compound'], 0],
                                                [right11, 1],
                                                [right22, 0]], 'type':'HDC'}
                                tot.append(equa)
        elif a['charge'] > 0:
#            b = {'compound':{'heavy':0}, 'charge': 0}
            for elem, num in a['compound'].items():
                if elem == 'H' or elem == 'F' or elem == 'Cl':
                    H_num = int(a['compound'][elem])
                    for i in range(H_num):
                        right = calculate_right_H(H_num, elem, a, i)
                        right[0][-1] = 0
                        equa = {'left':[[b['compound'], b['charge']], [a['compound'], a['charge']]],
                                'right':[[b['compound'], 1]] + right, 'type':'HDC'}
                        tot.append(equa)
                else:
                    if len(list(a['compound'].keys())) == 2:
                        for N_num in range(num):
                            if N_num != 0:
                                right11, right22 = calculate_split(N_num, num, a, elem)
                                equa = {'left':[[b['compound'], 0],
                                                [a['compound'], a['charge']]],
                                        'right':[[b['compound'], 1],
                                                [right11, 0],
                                                [right22, 0]], 'type':'HDC'}
                                tot.append(equa)
   
    else:
        if a['charge'] == 0:

            equa = {'left':[[b['compound'], 1], [a['compound'], a['charge']]],
                    'right':[[b['compound'], 0], [a['positive'], a['charge'] + 1], [a['negative'], 0]], 'type':'HDC'}
            return [equa]
        elif a['charge'] > 0:

            equa = {'left':[[b['compound'], 0], [a['compound'], a['charge']]],
                    'right':[[b['compound'], 1], [a['positive'], a['charge']], [a['negative'], 0]], 'type':'HDC'}
            return [equa]
    return tot

def HED(a, b, e = e):
    '''
    heavy particle electron detachment
    '''
    equa = {'left':[[a['compound'], a['charge']], [b['compound'], 0]],
            'right':[[a['compound'], a['charge'] + 1], [b['compound'], 0],
                     [e['compound'], e['charge']]],
            'type':'HED'}
    return [equa]