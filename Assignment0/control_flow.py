#!/usr/bin/env python

#########################################################################################
## Below are examples of control flow statements.  Try to figure out what this script 
## will print before running it. 
#########################################################################################

x = 0 
y = 1

# Now we can place booleans in 'if' statements; If true, then the statement will print 
print 'Example 1:'
if x:
    print 'hello1'

if not x:
    print 'hello2'

# We can add an else after the if statement so that a different set of operators can be performed in the boolean is False 
print 'Example 2:'
if x:
    print 'if'
elif y:
    print 'elif'
else:
    print 'else'


# we can also do 'for' loops where we loop through a specified set of numbers:
print 'Example 3:'
for i in range(5):
    print i

# Note: indexing in python starts at 0
# we can also loop through values in an array
print 'Example 4:'
arr = [5,'hello',6.7]
for i in arr:
    print i

# We can also use the 'while' loop;  Here the commands are iterated until the boolean is False

print 'Example 5:'
i = 0
while i < 10:
    i += 3
    print i


