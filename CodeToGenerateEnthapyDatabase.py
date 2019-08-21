import json
from tqdm import tqdm

def load_enthalpy():
    '''
    function to load enthalpy data
    both manully input data and ATcT data are load together
    produce two outputs:
        - entire_list: a list with dictionaries containing
          their enthalpy data
        - manually list: a list of species which has manually
          inpuy enthalpy
    '''
    upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    lower = 'abcdefghijklmnopqrstuvwxyz'
    num = '1234567890'
    #open the original database file
    with open("enth.txt", "r") as file:
        array = []
        
        for line in tqdm(file):
            comp_dict = {} 

            words = line.split()
            # extract only gases
            gas = False
            if len(words) > 4:
                for i in range(len(words)):
                    if words[i][0:2] == '(g':
                        comp = words[i-1] + '/'
                        gas = True
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
                comp_list = [comp_dict, charge]

                if gas:
                    #the case formation energy not equal to 0
                    if words[-5] == '±':
                        comp_list.append(float(words[-6]))
                    #the case formation energy equal to 0
                    else:
                        comp_list.append(float(words[-4]))
                    array.append(comp_list)

    # load manully input enthalpies:
    with open("EnthalpyInsertManually.txt", "r") as file:

        for line in tqdm(file):
            comp_dict = {} 

            words = line.split()
            # extract only gases

            if len(words) < 3 and len(words) > 1:
                comp = str(words[0]) + '/'
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
                comp_list = [comp_dict, charge]

                # append enthalpies into list
                comp_list.append(float(words[1]))
                array.append(comp_list)
    
    return array