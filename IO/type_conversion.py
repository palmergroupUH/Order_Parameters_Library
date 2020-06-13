import numpy as np 
from ctypes import CDLL, POINTER, c_int, c_double,c_char_p,c_long,c_float,byref 

#-------------------------------------------------------------------------
#                          Python data type into ctypes                   
#-------------------------------------------------------------------------

def string_to_ctypes_string(string):

    strlength = c_int(len(string))

    string_to_bytes = string.encode() 

    string = c_char_p(string_to_bytes) 

    return string,strlength 

def int_to_ctypes_int(data): 

    try: 

        data = c_int(data) 

    except TypeError: 

        pass 

    return data 

def np_to_ctypes_array(array): 

    if ( type(array) == np.ndarray ): 
    
        return np.ctypeslib.as_ctypes(array)  

    else: 

        return None 

