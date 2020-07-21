import numpy as np
from ctypes import CDLL, c_int, c_double, c_long, c_float, byref




call_c_lib = CDLL("./call_function.so")

a = c_double(4)

call_c_lib.call_fortran_like_c.restype = c_double

gg = call_c_lib.call_fortran_like_c(byref(a))

print((gg))


