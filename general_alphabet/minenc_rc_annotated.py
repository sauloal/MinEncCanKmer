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

def compl(c, alphabet_size):
	return alphabet_size - c - 1

def K(pos_left, pos_right, alphabet_size):
	return int((2 * alphabet_size - 3 - pos_left) * pos_left / 2 + pos_right - 1)

def spec_l(pos_left, pos_right, alphabet_size):
	return K(pos_left, pos_right, alphabet_size) // alphabet_size + 1

def spec_r(pos_left, pos_right, alphabet_size):
	return K(pos_left, pos_right, alphabet_size) % alphabet_size

def missing(pos_left, alphabet_size, kmer_size):
	return max(0,(alphabet_size**(kmer_size-pos_left)-alphabet_size**((kmer_size+1)//2))//2)

# determine code for k-mer
# code = list of length k over {0..alphabet_size-1}
def code(hash_array, kmer_size, alphabet_size):
	code_array = []
	pos_left   = 0

	# unspecific pairs
	while pos_left < kmer_size//2:
		# specifying case found?
		if   compl(hash_array[kmer_size - pos_left - 1], alphabet_size) > hash_array[pos_left]: break
		elif compl(hash_array[kmer_size - pos_left - 1], alphabet_size) < hash_array[pos_left]: return ([],pos_left)
		# left character = 0
		code_array.insert(pos_left,0)
		# right character by rank
		code_array.insert(pos_left+1,compl(hash_array[kmer_size - pos_left - 1], alphabet_size))
		# continue with next pair
		pos_left += 1

	# specifying case
	if pos_left == (kmer_size-1)/2: # single character left
		if compl(hash_array[pos_left], alphabet_size) < hash_array[pos_left]: return ([],pos_left) # not cononical
		code_array.insert(pos_left, compl(hash_array[pos_left], alphabet_size))

	elif pos_left < kmer_size//2: # specifying pair
		code_array.insert(pos_left  , spec_l(hash_array[pos_left], compl(hash_array[kmer_size - pos_left - 1], alphabet_size), alphabet_size))
		code_array.insert(pos_left+1, spec_r(hash_array[pos_left], compl(hash_array[kmer_size - pos_left - 1], alphabet_size), alphabet_size))

		#remainder
		end = kmer_size - pos_left - 1
		for pos_right in range(pos_left+1, end):
			code_array.insert(pos_right+1, int(hash_array[pos_right]))

	return (code_array, pos_left)


# compose integer value (rank) from code
def intcode(code_array, pos_left, alphabet_size, kmer_size):
	if len(code_array)<=0: return -1

	code_int = 0

	for pos in range(0, kmer_size):
		# standard ranking
		code_int += code_array[kmer_size-pos-1]*alphabet_size**(pos)

	# subtract gaps in code due to specifying pairs
	#c-=comb(alphabet_size,2)*partial_sum((kmer_size+1)//2-1,kmer_size-3-pos_left)
	code_int -= missing(pos_left+1, alphabet_size, kmer_size)

	# subtract gap in code due to specifying middle position
	if kmer_size%2 == 1: code_int -= (alphabet_size**((kmer_size+1)//2))//2

	return code_int


def main():
	# Usage / paramters
	if len(sys.argv)<3:
		eprint("Usage: mphf_rc.py <alphabet_size> <kmer_size>")
		eprint("<alphabet_size> must be even")
		eprint("Output: enc(x) x for all k-mers")
		sys.exit(1)


	# alphabet = {0,...,alphabet_size-1}
	alphabet_size = int(sys.argv[1])
	kmer_size     = int(sys.argv[2])

	# even alphabets only
	if alphabet_size%2 == 1:
		eprint("<alphabet_size> must be even")
		exit(1)

	# construct a list of all k-mers
	hash_list = [[i] for i in range(0,alphabet_size)]
	for i in range(1,kmer_size):
		hash_list = [x+[i] for x in hash_list for i in range(0,alphabet_size)]

	# encode each k-mer
	for hash_array in hash_list:
		(code_array, pos_left) = code(hash_array, kmer_size, alphabet_size)
		if code_array != []:
			code_int = intcode(code_array, pos_left, alphabet_size, kmer_size)
			# output
			print(code_int, hash_array)

if __name__ == "__main__":
	main()
