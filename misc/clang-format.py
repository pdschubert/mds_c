#!/usr/bin/env python3

##############################################################################
# Copyright (c) 2018 Philipp Schubert.                                       #
# All rights reserved. This program and the accompanying materials are made  #
# available under the terms of LICENSE.txt.                                  #
#                                                                            #
# Contributors:                                                              #
#     Philipp Schubert                                                       #
##############################################################################

import os

SRC_DIRS = ["src/"]

cpp_extensions = (".cpp",
									".cxx",
									".c++",
									".cc",
									".cp",
									".c",
									".i",
									".ii",
									".h",
									".h++",
									".hpp",
									".hxx",
									".hh",
									".inl",
									".inc",
									".ipp",
									".ixx",
									".txx",
									".tpp",
									".tcc",
									".tpl",
									".cu",
									".cuh")

for SRC_DIR in SRC_DIRS:
	for root, dir, files in os.walk(SRC_DIR):
		for file in files:
			if file.endswith(cpp_extensions):
				os.system("clang-format -i " + root + "/" + file)
