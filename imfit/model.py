'''
Created on Sep 18, 2013

@author: andre
'''


__all__ = ['ParameterDescription', 'FunctionDescription', 'FunctionSetDescription']

################################################################################

class ParameterDescription(object):
    def __init__(self, name, value, limits=None, fixed=False):
        self.name = name
        self.value = value
        self.fixed = fixed
        self.limits = limits
        
        
    def __str__(self):
        if self.fixed:
            return '%s    %f     fixed' % (self.name, self.value)
        elif self.limits is not None:
            return '%s    %f     %f,%f' % (self.name, self.value, self.limits[0], self.limits[1])
        else:
            return '%s    %f' % (self.name, self.value)
            
################################################################################

class FunctionDescription(object):
    def __init__(self, name, parameters=None):
        self.name = name
        if parameters is None:
            parameters = []
        self._parameters = parameters
        
        
    def addParameter(self, p):
        if not isinstance(p, ParameterDescription):
            raise ValueError('p is not a Parameter object.')
        self._parameters.append(p)
        
    
    def parameterList(self):
        return [p for p in self._parameters]


    def __str__(self):
        lines = []
        lines.append('FUNCTION %s' % self.name)
        lines.extend(str(p) for p in self._parameters)
        return '\n'.join(lines)

################################################################################

class FunctionSetDescription(object):
    def __init__(self, x0, y0, functions=None):
        self.x0 = x0
        self.y0 = y0
        if functions is None:
            functions = []
        self._functions = functions
        
        
    def addFunction(self, f):
        if not isinstance(f, FunctionDescription):
            raise ValueError('func is not a Function object.')
        self._functions.append(f)
        
        
    def functionList(self):
        return [f.name for f in self._functions]
    
    
    def parameterList(self):
        params = []
        params.append(self.x0)
        params.append(self.y0)
        for f in self._functions:
            params.extend(f.parameterList())
        return params
    
    
    def __str__(self):
        lines = []
        lines.append(str(self.x0))
        lines.append(str(self.y0))
        lines.extend(str(f) for f in self._functions)
        return '\n'.join(lines)
        
################################################################################
        
class ModelDescription(object):
    def __init__(self, function_sets=None, options={}):
        if function_sets is None:
            function_sets = []
        self._functionSets = function_sets
        self.options = options


    @classmethod
    def example(cls):
        x0 = ParameterDescription('X0', 36.0, [25, 45])
        y0 = ParameterDescription('Y0', 32.0, [25, 45])
        fs = FunctionSetDescription(x0, y0)
        sersic = FunctionDescription('Sersic', [ParameterDescription('PA', 93.0217, [0, 180]),
                                                ParameterDescription('ell', 0.37666, [0, 1]),
                                                ParameterDescription('n', 4, fixed=True),
                                                ParameterDescription('I_e', 1, [0, 10]),
                                                ParameterDescription('r_e', 25, [0, 100]),
                                                ])
        exponential = FunctionDescription('Exponential', [ParameterDescription('PA', 93.0217, [0, 180]),
                                                          ParameterDescription('ell', 0.37666, [0, 1]),
                                                          ParameterDescription('I_0', 1, [0, 10]),
                                                          ParameterDescription('h', 25, [0, 100])
                                                          ])
        fs.addFunction(sersic)
        fs.addFunction(exponential)
        return cls([fs])
        
    
    @classmethod
    def load(cls, fname):
        from .config import parse_config_file
        return parse_config_file(fname)
    
    
    def addFunctionSet(self, fs):
        if not isinstance(fs, FunctionSetDescription):
            raise ValueError('fs is not a FunctionSet object.')
        self._functionSets.append(fs)
    
    
    def functionSetIndices(self):
        indices = [0]
        for i in xrange(len(self._functionSets) - 1):
            indices.append(len(self.functionSets[i].functions))
        return indices
        
        
    def functionList(self):
        functions = []
        for function_set in self._functionSets:
            functions.extend(function_set.functionList())
        return functions
    

    def parameterList(self):
        params = []
        for function_set in self._functionSets:
            params.extend(function_set.parameterList())
        return params


    def __str__(self):
        lines = []
        for k, v in self.options.items():
            lines.append('%s    %f' % (k, v))
        lines.extend(str(fs) for fs in self._functionSets)
        return '\n'.join(lines)
        
################################################################################
