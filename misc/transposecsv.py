#!/usr/bin/python

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
import csv

if len(sys.argv) != 3:
	print("usage: prog <input> <output>")
	exit(1)

input = sys.argv[1]
output = sys.argv[2]

with open(input) as f:
    reader = csv.reader(f)
    cols = []
    for row in reader:
        cols.append(row)

with open(output, 'wb') as f:
    writer = csv.writer(f)
    for i in range(len(max(cols, key=len))):
        writer.writerow([(c[i] if i<len(c) else '') for c in cols])

exit(0)
