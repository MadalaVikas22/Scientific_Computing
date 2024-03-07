import numpy as np
import sys
import math


def is_close(a, b):
    return abs(b-a) < sys.float_info.epsilon

def function(x):
    fn = (np.cos(x) - x)
    return fn
def bisection_method(function, a, b, tol = 0.001, max_iter = 100):
    iter = 0
    while iter < max_iter:
        c = (a + b) / 2

        if math.isclose(function(c), 0.0, abs_tol = 1.0E-6) or abs(b-a)/2 < tol:
            print(c , " is the root of the function")
        else :
            if function(a) * function(c) < 0:
                b = c
                function(b) == function(c)
            else:
                a = c
                function(a) == function(c)
            iter = iter + 1

    return None



a=0
b=2
bisection_method(function, a, b)