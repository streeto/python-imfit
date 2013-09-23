'''
Created on Sep 18, 2013

@author: andre
'''


__all__ = ['ModelDescription', 'ParameterDescription', 'FunctionDescription', 'FunctionSetDescription']

################################################################################

class ParameterDescription(object):
    def __init__(self, name, value, limits=None, fixed=False):
        self._name = name
        self.setValue(value, limits, fixed)
        
    
    @property
    def name(self):
        return self._name
    
    
    def setValue(self, value, limits=None, fixed=False):
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

    
    def __getitem__(self, key):
        if not isinstance(key, str):
            raise KeyError('Parameter must be a string.')
        for p in self._parameters:
            if key == p.name:
                return p
        raise KeyError('Parameter %s not found.' % key)
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
