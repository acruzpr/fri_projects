#!/usr/bin/env python

x = 5

# this function print x + y where y is a value you select when calling the function and x is predefined above 
def printer(y):
    print x + y

print 'function 1:'
printer(5)

# this function returns x + y where x and y is called by the user; now the predefined x is not used
def return1(x, y):
    w = x + y 
    return w

print 'function 2:'
print return1(4, 5)

