#!/usr/bin/env python

# First, we will need import the numpy module so that its tools can be used in the python script
# below are two common ways of importing modules

import numpy # ONE
# or #
import numpy as np # TWO

# The central feature of numpy is the array.  These are similar to list in python except that every element needs to be of the same type.  Creating arrays over list will make your code more efficient when doing numerical calculations. Here is how you create an array with numpy.  

a = numpy.array([1,2,3],float) # NOTE: this is how numpy is called with the import ONE 
b = np.array([1,2,3],float) # NOTE: this is how numpy is called with the import TWO 

print '1:',a
print '2:',b

# Arrays can be treated just like lists from the variables section of this tutorial 

print '3:',a[0] 
print '4:',a + b

# Note that if you add arrays of different dimensions the smaller array will be repeated as necessary to perform the operation

c = numpy.array([5],float)
print '5:',a + c

# One useful way to create new arrays is to specify arrays of a certain size and fill them with zeros or ones

print '6:',numpy.ones((2,3)) # the numbers indicate the shape of the array in this case a 2x3 matrix
print '7:',numpy.zeros(5) # The number indicates the an array of length 5 

# Additional numpy has many useful functions such as taking the dot product between two vectors, finding the mean value of an array, etc. Here are a few examples 

print '8:',numpy.vdot(a,b)
print '9:',numpy.mean(a)

# For more information on numpy's useful functions, check out there tutorial ( https://docs.scipy.org/doc/numpy-dev/user/quickstart.html#basic-operations ) or just try googling numpy and the name math operation you are interested in! 
