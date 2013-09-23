'''
Created on Sep 23, 2013

@author: andre
'''
from imfit import ParameterDescription, FunctionSetDescription, ModelDescription
from imfit import getFunctionDescription


def example_model_description():
    x0 = ParameterDescription('X0', 36.0, [25, 45])
    y0 = ParameterDescription('Y0', 32.0, [25, 45])
    fs = FunctionSetDescription(x0, y0)

    sersic = getFunctionDescription('Sersic')
    sersic['PA'].setValue(93.0217, [0, 180])
    sersic['ell'].setValue(0.37666, [0, 1])
    sersic['n'].setValue(4, fixed=True)
    sersic['I_e'].setValue(1, [0, 10])
    sersic['r_e'].setValue(25, [0, 100])

    exponential = getFunctionDescription('Exponential')
    exponential['PA'].setValue(93.0217, [0, 180])
    exponential['ell'].setValue(0.37666, [0, 1])
    exponential['I_0'].setValue(1, [0, 10])
    exponential['h'].setValue(25, [0, 100])

    fs.addFunction(sersic)
    fs.addFunction(exponential)
    return ModelDescription([fs])
    

def test_model():
    desc = example_model_description()
    print desc