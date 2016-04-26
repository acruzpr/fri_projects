#!/usr/bin/env python

########################################################################
# The following script shows examples of various boolean operators 
# Try printing the example to see if the solution you think is correct.
########################################################################

# defining variables
x = 1
y = 'hello'

## The following operator checks for equality 
print 'The == symbol ask if two parameters are equal: 1==1 is',1==1
############Examples ##############
x == -1


# The following operator checks for inequality 
print 'The != symbol ask if two parameters are equal: 1!=1 is',1!=1
############# Examples ############
y != 'hello'


# Inequalities are written as follows
print 'You can also use inequalities: 1 > 2 is', 1 > 2
############# Examples ############
9 > 9
9 >= 9

# and, or, not
print 'Lastly, there are and statements: True and True is', True and True,',versus True and False is', True and False 
print 'or statements: True or True is', True or True,',versus True or False is', True or False, ',versus False or False is', False or False

print 'not statements: not True is', not True,',versus not False is', not False
############# Examples ############
x == 1 or y == 'hi'
x > 0 and y != 'hi'

y = 0
x or y
x and y
not y

