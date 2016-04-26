#!/usr/bin/env python

f1 = open('write.txt', 'r')
contents = f1.read()
print contents
# or you can read in the file the following way line by line

f2 = open('write.txt', 'r')
contents2 = f2.readlines()
print contents2
for i in contents2:
    w = i.split()
    print w
