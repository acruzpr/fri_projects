#!/usr/bin/env python

b = 2

f = open('write.txt', 'w')
f.write('line 1 \n')  # Note: \n tells the text file to go to the next line
f.write('line'+str(b)+ '\n') # You can convert integers or floats to strings by the function str to the left
f.write('line 3 \n')
f.close()


