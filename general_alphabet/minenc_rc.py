#!/usr/bin/env python3


from __future__ import print_function
import sys
import os
from math import log
from math import sqrt
from math import comb

def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)


#rank=int(.)

def compl(c):
	return sigma-c-1

def K(i,j):
	return int((2 * sigma - 3 - i) * i / 2 + j - 1)

def spec_l(i,j):
	return K(i,j)//sigma+1

def spec_r(i,j):
	return K(i,j)%sigma

def missing(i):
	return max(0,(sigma**(k-i)-sigma**((k+1)//2))//2)

# determine code for k-mer
# code = list of length k over {0..sigma-1}
def code(x):
	code=[]
	i=0

	# unspecific pairs
	while i<k//2:
		# specifying case found?
		if compl(x[k-i-1])>x[i]: break
		elif compl(x[k-i-1])<x[i]: return ([],i)
		# left character = 0
		code.insert(i,0)
		# right character by rank
		code.insert(i+1,compl(x[k-i-1]))
		# continue with next pair
		i+=1
		
	# specifying case
	if i==(k-1)/2: # single character left
		if compl(x[i])<x[i]: return ([],i) # not cononical
		code.insert(i,compl(x[i]))

	elif i<k//2: # specifying pair
		code.insert(i,spec_l(x[i],compl(x[k-i-1])))
		code.insert(i+1,spec_r(x[i],compl(x[k-i-1])))
		
		#remainder
		end=k-i-1
		for r in range(i+1,end):
			code.insert(r+1,int(x[r]))
	
	return (code,i)


# compose integer value (rank) from code
def intcode(code,i):
	if len(code)<=0: return -1
	c=0
	for j in range(0,k):
		# standard ranking
		c += code[k-j-1]*sigma**(j)
	
	# subtract gaps in code due to specifying pairs
	#c-=comb(sigma,2)*partial_sum((k+1)//2-1,k-3-i)
	c-=missing(i+1)
	
	# subtract gap in code due to specifying middle position
	if k%2==1: c-=(sigma**((k+1)//2))//2
	
	return c



# Usage / paramters
if len(sys.argv)<3:
    eprint("Usage: mphf_rc.py <sigma> <k>")
    eprint("<sigma> must be even")
    eprint("Output: enc(x) x for all k-mers")
    sys.exit(1)


# alphabet = {0,...,sigma-1}
sigma=int(sys.argv[1])
k=int(sys.argv[2])

# even alphabets only
if sigma%2==1:
    eprint("<sigma> must be even")
    exit(1)

# construct a list of all k-mers
l=[[i] for i in range(0,sigma)]
for i in range(1,k):
	l=[x+[i] for x in l for i in range(0,sigma)]

# encode each k-mer
for x in l:
	(cd,i)=code(x)
	if cd!=[]:
		c=intcode(cd,i)
		# output
		print(c, x)
		
		
		
		
		
		
		
