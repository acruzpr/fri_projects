#!/usr/bin/env python

# Below is a blank canvas to start your assignment
# import numpy as np
import math

print 'Hello World'

a = 23.5
b = 57.6

print a + b
print abs(a - b)
print a * b

print (a + b) <= 90
print abs(a - b) <= 90
print (a * b) <= 90

print '\nRange with for'
for i in range(1, 11):
	print i

print '\nRange with while'
i = 1
while i <= 10:
	print i
	i += 1

print '\nPrimes'
newFile = open('primes.out', 'w')
for i in range(10):
	x = 2
	prime = i >= 2
	while x < i:
		if i % x == 0:
			prime = False
			break
		x += 1
	if prime:
		print i
		newFile.write(str(i) + '\n')

def magnitude(arr):
	mag = 0
	for x in arr:
		mag += x * x
	return math.sqrt(mag)


def normalize(arr):
	mag = magnitude(arr)
	x = 0
	while x < len(arr):
		arr[x] = arr[x] / mag;
		x += 1

vector = [0.8, -0.9, 1.6]

print "\nMagnitude: "
print magnitude(vector)

normalize(vector)
print "\nNormalized Vector:"
print vector

print "\nReimann Sum:"
total = 0.0
step = .1
at = -100
while at <= 100:
	total += step * (5.0 * math.exp(-at * at / 12))
	at = at + step
print total
