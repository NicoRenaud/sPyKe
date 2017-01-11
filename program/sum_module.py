import ctypes
from ctypes.util import find_library

# find and load the library
lib_sum = ctypes.cdll.LoadLibrary(find_library('sum'))


# set the argument type
lib_sum.sum.argtypes = [ctypes.c_double,ctypes.c_double]

# set the return type
lib_sum.sum.restype = ctypes.c_double


def sum_func(arg1,arg2):
    ''' Wrapper for sum from sum.dylib '''
    return lib_sum.sum(arg1,arg2
    )

