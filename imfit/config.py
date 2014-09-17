'''
Created on Sep 19, 2013

@author: andre
'''

from .model import ParameterDescription, FunctionDescription, FunctionSetDescription, ModelDescription

__all__ = ['parse_config_file', 'parse_config']

################################################################################

comment = '#'
x0_str = 'X0'
y0_str = 'Y0'
function_str = 'FUNCTION'
fixed_str = 'fixed'

################################################################################

def parse_config_file(fname):
    '''
    Read an Imfit model description file.
    
    Parameters
    ----------
    fname : string
        Path to the model description file.
        
    Returns
    -------
    model : :class:`~imfit.ModelDescription`
        A model description object.
        
    See also
    --------
    parse_config
    '''
    with open(fname) as fd:
        return parse_config(fd.readlines())
            
################################################################################
    
def parse_config(lines):
    '''
    Parses an Imfit model description from a list of strings.
    
    Parameters
    ----------
    fname : list of strings
        String representantion of Imfit model description.
        
    Returns
    -------
    model : :class:`~imfit.ModelDescription`
        A model description object.
        
    See also
    --------
    parse_config_file
    '''
    lines = clean_lines(lines)

    model = ModelDescription()
    
    block_start = 0
    id_fs = 0
    for i in xrange(block_start, len(lines)):
        if lines[i].startswith(x0_str):
            if block_start == 0: 
                options = read_options(lines[block_start:i])
                model.options.update(options)
            else:
                model.addFunctionSet(read_function_set('fs%2d' % id_fs, lines[block_start:i]))
                id_fs += 1
            block_start = i
    model.addFunctionSet(read_function_set('fs%d' % id_fs, lines[block_start:i+1]))
    return model

################################################################################

def clean_lines(lines):
    clean = []
    for l in lines:
        # Clean the comments.
        l = l.split(comment, 1)[0]
        # Remove leading and trailing whitespace.
        l = l.strip()
        # Skip the empty lines.
        if l == '':
            continue
        clean.append(l)
    return clean
    
################################################################################
    
def read_options(lines):
    config = {}
    for l in lines:
        # Options are key-value pairs.
        k, val = l.split(' ', 1)
        if k in [x0_str, y0_str, function_str]:
            raise ValueError('Expected option, but got %s instead.' % k)
        val = val.strip()
        config[k] = val
        
    return config
    
################################################################################

def read_function_set(name, lines):
    # A function set starts with X0 and Y0 parameters.
    x0 = read_parameter(lines[0])
    y0 = read_parameter(lines[1])
    if x0.name != x0_str or y0.name != y0_str:
        raise ValueError('A function set must begin with the parameters X0 and Y0.')
    fs = FunctionSetDescription(name)
    fs.x0 = x0
    fs.y0 = y0
    block_start = 2
    for i in xrange(block_start, len(lines)):
        # Functions for a given set start with FUNCTION.
        if i == block_start or not lines[i].startswith(function_str):
            continue
        fs.addFunction(read_function(lines[block_start:i]))
        block_start = i
    # Add the last function in the set.
    fs.addFunction(read_function(lines[block_start:i+1]))
    return fs

################################################################################
        
def read_function(lines):
    # First line contains the function name.
    test_function_str, name = lines[0].split(' ', 1)
    if test_function_str != function_str:
        raise ValueError('Function definition must begin with FUNCTION.')
    name = name.strip()
    func = FunctionDescription(name)
    # Read the function parameters.
    for i in xrange(1, len(lines)):
        func.addParameter(read_parameter(lines[i]))

    # FIXME: check function and parameters.
    return func

################################################################################

def read_parameter(line):
    ulimit = None
    llimit = None
    fixed = False
    
    # Format:
    # PAR_NAME    VALUE   ( "fixed" | LLIMIT,ULIMIT )
    name, contents = line.split(' ', 1)
    contents = contents.strip()
    value, predicate = contents.split(' ', 1)
    value = float(value)
    predicate = predicate.strip()
    
    if predicate == fixed_str:
        fixed = True

    elif ',' in predicate:
        llimit, ulimit = predicate.split(',')
        llimit = float(llimit)
        ulimit = float(ulimit)
        if llimit > ulimit:
            raise ValueError('lower limit (%f) is larger than upper limit (%f)' & (llimit, ulimit))

    return ParameterDescription(name, value, llimit, ulimit, fixed)

################################################################################
