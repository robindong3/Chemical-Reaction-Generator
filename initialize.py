import numpy as np
class initialize():
    def __init__(self,compound,charge):
        integers = '0123456789'
        self.compound = {}
        self.positive = {}
        self.negative = {}
        if type(compound) == dict:
            self.compound = compound
            key_list = list(compound.keys())
            if len(key_list) == 1:
                if compound[key_list[0]] == 1:
                    self.type = 'element'
                else:
                    self.type = 'multiple'
                    self.positive[key_list[0]] = compound[key_list[0]] - 1
                    self.negative[key_list[0]] = 1
            else:
                if compound[key_list[-1]] == 1:
                    self.type = 'multiple'
                    self.positive = dict(compound)
                    del self.positive[key_list[-1]]
                    self.negative[key_list[-1]] = 1
                else:
                    self.type = 'multiple'
                    self.positive = dict(compound)
                    self.positive[key_list[-1]] = compound[key_list[-1]] - 1
                    self.negative[key_list[-1]] = 1

        else:
            raise TypeError('input type should be a dictionary or a string.')

        self.charge=charge

    def result(self):
        if self.type=='element':
            lib={'compound':self.compound,'type':'element','charge':self.charge}
            return lib
        elif self.type=='multiple':
            if len(list(self.compound.keys())) == 1:
                lib={'compound':self.compound,'positive':self.positive,'negative':self.negative,
                    'type':'single','charge':self.charge}
            else:
                lib={'compound':self.compound,'positive':self.positive,'negative':self.negative,
                    'type':'multiple','charge':self.charge}
            return lib