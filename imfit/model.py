'''
Created on Sep 18, 2013

@author: andre
'''
from copy import copy, deepcopy


__all__ = ['SimpleModelDescription', 'ModelDescription',
           'ParameterDescription', 'FunctionDescription', 'FunctionSetDescription']

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
    def __init__(self, func_type, name=None, parameters=None):
        if name is None:
            name = func_type
        self.funcType = func_type
        self.name = name
        self._parameters = []
        if parameters is not None:
            for p in parameters:
                self.addParameter(p)
        
        
    def addParameter(self, p):
        if not isinstance(p, ParameterDescription):
            raise ValueError('p is not a Parameter object.')
        self._parameters.append(p)
        
    
    def parameterList(self):
        return [p for p in self._parameters]


    def __str__(self):
        lines = []
        lines.append('FUNCTION %s # %s' % (self.funcType, self.name))
        lines.extend(str(p) for p in self._parameters)
        return '\n'.join(lines)

    
    def __getattr__(self, attr):
        return self[attr]
    
    
    def __getitem__(self, key):
        if not isinstance(key, str):
            raise KeyError('Parameter must be a string.')
        for p in self._parameters:
            if key == p.name:
                return p
        raise KeyError('Parameter %s not found.' % key)
    

    def __deepcopy__(self, memo):
        f = FunctionDescription(self.funcType, self.name)
        f._parameters = deepcopy(self._parameters)
        return f

################################################################################

class FunctionSetDescription(object):
    def __init__(self, name, functions=None):
        self.name = name
        self.x0 = ParameterDescription('X0', 0.0)
        self.y0 = ParameterDescription('Y0', 0.0)
        self._functions = []
        if functions is not None:
            for f in functions:
                self.addFunction(f)
        
        
    def addFunction(self, f):
        if not isinstance(f, FunctionDescription):
            raise ValueError('func is not a Function object.')
        if self._contains(f.name):
            raise KeyError('Function named %s already exists.' % f.name)
        self._functions.append(f)
    
    
    def _contains(self, name):
        for f in self._functions:
            if f.name == name:
                return True
        return False
    
    
    def functionList(self):
        return [f.funcType for f in self._functions]
    
    
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
    
    def __getattr__(self, attr):
        return self[attr]
    
    
    def __getitem__(self, key):
        if not isinstance(key, str):
            raise KeyError('Function must be a string.')
        for f in self._functions:
            if key == f.name:
                return f
        raise KeyError('Function %s not found.' % key)
    
    
    def __deepcopy__(self, memo):
        fs = FunctionSetDescription(self.name)
        fs.x0 = copy(self.x0)
        fs.y0 = copy(self.y0)
        fs._functions = deepcopy(self._functions)
        return fs
        
################################################################################
        
class ModelDescription(object):
    def __init__(self, function_sets=None, options={}):
        self.options = {}
        self.options.update(options)
        self._functionSets = []
        if function_sets is not None:
            for fs in function_sets:
                self.addFunctionSet(fs)


    @classmethod
    def load(cls, fname):
        '''
        Load a model description from a file. The syntax is the same
        as the imfit config file.
        '''
        from .config import parse_config_file
        return parse_config_file(fname)
    
    
    def addFunctionSet(self, fs):
        '''
        Add a function set to the model description.
        '''
        if not isinstance(fs, FunctionSetDescription):
            raise ValueError('fs is not a FunctionSet object.')
        if self._contains(fs.name):
            raise KeyError('FunctionSet named %s already exists.' % fs.name)
        self._functionSets.append(fs)
    
    
    def _contains(self, name):
        for fs in self._functionSets:
            if fs.name == name:
                return True
        return False
    
    
    def functionSetIndices(self):
        '''
        Returns the indices in the full parameters list such that
        imfit can split the parameters for in the function sets.
        '''
        indices = [0]
        for i in xrange(len(self._functionSets) - 1):
            indices.append(len(self.functionSets[i].functions))
        return indices
        
        
    def functionList(self):
        '''
        Returns the functions composing this model, as a list of strings.
        '''
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
        

    def __getattr__(self, attr):
        return self[attr]
    
    
    def __getitem__(self, key):
        if not isinstance(key, str):
            raise KeyError('FunctionSet must be a string.')
        for fs in self._functionSets:
            if key == fs.name:
                return fs
        raise KeyError('FunctionSet %s not found.' % key)
    
    
    def __deepcopy__(self, memo):
        model = ModelDescription()
        model._functionSets = deepcopy(self._functionSets)
        return model
        
################################################################################


class SimpleModelDescription(ModelDescription):
    def __init__(self):
        super(SimpleModelDescription, self).__init__()
        fs = FunctionSetDescription('fs')
        self.addFunctionSet(fs)
        
        
    @property
    def x0(self):
        return self._functionSets[0].x0
        

    @property
    def y0(self):
        return self._functionSets[0].y0
    

    def addFunction(self, f):
        self._functionSets[0].addFunction(f)
        
        
    def __getattr__(self, attr):
        return self._functionSets[0][attr]
