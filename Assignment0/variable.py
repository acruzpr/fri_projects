#!/usr/bin/env python

# Note: '#' signifies a comment in python; the comment terminates the end of the line

# Below are examples of four types of variables in python

x = 1.343  # a float or a number; Note that if this number does not have a decimal..

y = 5  # python will assume this number is an integer, like y

greeting = 'hi'  # a string 

arr = [1, 1, 2, 3, 5, 8]  # a list

matrix = [[1,2],[3,4]]  # a matrix or list of lists

# to print do the following

print 'Before variable change'
print x

# Note: a variable can change

x = 6.7

print 'After variable change'
print x
print

# Variables can also change types; below y is changed from an integrer to a string

y = 'new'

# To print a specified value from list your can do the following

print 'Below is examples of how to pull out values from a list:'
print arr[0] # prints the first value from arr; Note indexing starts at zero
print matrix[1] # prints the second row of the matrix
print matrix[0][1] # prints the second value from the first row of the matrix

