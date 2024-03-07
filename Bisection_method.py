# Importing necessary libraies
import numpy as np
import sys
import math

def is_close(a, b):
    return abs(b - a) < sys.float_info.epsilon

def function(x):
    fn = np.cos(x) - x
    return fn

def bisection_method(function, a, b, tol=0.001, max_iter=10):
    iter = 0
    while iter < max_iter:
        c = (a + b) / 2

        if not math.isclose(function(c), 0.0, abs_tol=1.0E-6) and abs(b - a) >= tol:
            if function(a) * function(c) < 0:
                b = c
            else:
                a = c
            iter += 1
        else:
            print(c, " is the root of the function")
            return c

    return None

# Example usage
a = 0
b = 2
tolerance = 0.0001
root = bisection_method(function, a, b, tol=tolerance)
print("Root:", root)




a=0
b=2
bisection_method(function, a, b)