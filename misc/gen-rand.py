#!/usr/bin/python3

##############################################################################
# Copyright (c) 2015 - 2016 Philipp Schubert.                                #
# All rights reserved. This program and the accompanying materials are made  #
# available under the terms of LICENSE.txt.                                  #
#                                                                            #
# Contributors:                                                              #
#     Philipp Schubert                                                       #
##############################################################################

__author__ = 'Philipp Dominik Schubert'

import sys
import random

vecs = 10
veclen = 5
min = 1
max = 100

if len(sys.argv) != 5:
	print("usage: prog <vecs n> <vec len> <min value> <max value>")
	print("for example: prog 100 5 1 20")
	exit(1)
else:
	vecs = int(sys.argv[1])
	veclen = int(sys.argv[2])
	min = int(sys.argv[3])
	max = int(sys.argv[4])

for i in range(vecs):
	if i < vecs - 1:
		print(i+1, end=",")
	else:
		print(i+1)

for i in range(veclen):
	for j in range(vecs):
		if j < vecs - 1:
			print(random.randint(min,max), end=",")
		else:
			print(random.randint(min,max))
exit(0)
