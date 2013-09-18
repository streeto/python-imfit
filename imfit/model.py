'''
Created on Sep 18, 2013

@author: andre
'''


class Parameter(object):
    def __init__(self, name, value, limits=None, fixed=False):
        self.name = name
        self.value = value
        self.fixed = fixed
        self.limits = limits


class Function(object):
    def __init__(self, name, parameters):
        self.name = name
        self.parameters = parameters


class FunctionSet(object):
    def __init__(self, x0, y0, functions = []):
        self.x0 = x0
        self.y0 = y0
        self.functions = functions
        
    def parameterList(self):
        params = []
        params.append(self.x0)
        params.append(self.y0)
        for f in self.functions:
            params.extend(f.parameters)
        return params
        
        
class ModelDescription(object):
    def __init__(self):
        # TODO: Get rid of this test model.
        
        functions = [Function('Sersic', [Parameter('PA', 93.0217, [0, 180]),
                                         Parameter('ell', 0.37666, [0, 1]),
                                         Parameter('n', 4, fixed=True),
                                         Parameter('I_e', 1, [0, 10]),
                                         Parameter('r_e', 25, [0, 100])
                                         ]),
                     Function('Exponential', [Parameter('PA', 93.0217, [0, 180]),
                                              Parameter('ell', 0.37666, [0, 1]),
                                              Parameter('I_0', 1, [0, 10]),
                                              Parameter('h', 25, [0, 100])
                                              ]),
                     ]
        x0 = Parameter('X0', 36.0, [25, 45])
        y0 = Parameter('Y0', 32.0, [25, 45])
        self.functionSets = [FunctionSet(x0, y0, functions)]

    def load(self, fname):
        # TODO: read model from file.
        raise NotImplementedError()
    
    
    def functionSetIndices(self):
        indices = [0]
        for i in xrange(len(self.functionSets) - 1):
            indices.append(len(self.functionSets[i].functions))
        return indices
        
        
    def functionList(self):
        functions = []
        for function_set in self.functionSets:
            functions.extend([f.name for f in function_set.functions])
        return functions
    

    def parameterList(self):
        params = []
        for function_set in self.functionSets:
            params.extend(function_set.parameterList())
        return params
