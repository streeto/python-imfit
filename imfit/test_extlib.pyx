import cython

cdef extern from 'imfit/add_functions.h':
    void PrintAvailableFunctions()

def print_available_functions():
    PrintAvailableFunctions()
    