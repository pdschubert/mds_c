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
import csv
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) != 2:
	print("usage: prog <input>")
	exit(1)

input = sys.argv[1]

with open(input) as f:
    reader = csv.reader(f)
    x = []
    y = []
    for row in reader:
        x.append(row[1])
        y.append(row[2])

plt.scatter(x, y)
plt.show()

exit(0)
