##############################################################################
# Copyright (c) 2015 - 2016 Philipp Schubert.                                #
# All rights reserved. This program and the accompanying materials are made  #
# available under the terms of LICENSE.txt.                                  #
#                                                                            #
# Contributors:                                                              #
#     Philipp Schubert                                                       #
##############################################################################

EXE := smacofC
CUEXE := smacofCU
BUILDDIR := build
# used libraries
LIBS := -lstdc++ -lm -larmadillo
CULIBS := -lcudart
# some very important switches
# start algorithm with seeded random matrix (1) or unseeded random matrix (0)
SEED_BOOL := 1
# make use of a memory pool (pre-allocation) (1) or do not (0)
MEM_POOL := 1
# setting the compiler to use
CC := gcc
CFLAGS :=
# flag for standard C version
ifeq ($(findstring clang,$(CC)),clang)
else
	CFLAGS += -std=c11
endif
# CFLAGS += -g
CFLAGS += -Wall
CFLAGS += -Wextra
CFLAGS += -Ofast
CFLAGS += -march=native
CFLAGS += -ffunction-sections
CFLAGS += -fdata-sections
CFLAGS += -fuse-ld=gold
CFLAGS += -flto=full
CFLAGS += -fopenmp
# flags for CUDA version
CUFLAGS := -O3
CUFLAGS += -use_fast_math
# this line has to be adapted to existing hardware!
CUFLAGS += -gencode arch=compute_30,code=sm_30
CUFLAGS += -Xcompiler -std=c11
CUFLAGS += -Xcompiler -Wall
CUFLAGS += -Xcompiler -Wextra
CUFLAGS += -Xcompiler -Ofast
CUFLAGS += -Xcompiler -march=native
CUFLAGS += -Xcompiler -fopenmp

c:
	mkdir $(BUILDDIR); \
	cd src;	\
	$(CC) $(CFLAGS) -DCUDA=0 -DDEBUG=0 -DSEEDED=$(SEED_BOOL) -DMPOOL=$(MEM_POOL) *.cpp *.c -o $(EXE) $(LIBS); \
	mv $(EXE) ../$(BUILDDIR); \
	cd -; \

cuda:
	mkdir $(BUILDDIR); \
	cd src; \
	nvcc $(CUFLAGS) -DCUDA=1 -DDEBUG=0 -DSEEDED=$(SEED_BOOL) -DMPOOL=$(MEM_POOL) *.cu *.cpp *.c -o $(CUEXE) $(LIBS) $(CULIBDS); \
	mv $(CUEXE) ../$(BUILDDIR); \
	cd -; \

clang-format:
	python3 misc/clang-format.py

run-cubeC-example: c
	cd $(BUILDDIR); \
	./$(EXE) ../example-data/cube.csv cube-output.csv 0 1000 0.0001 2 2 1; \

run-cubeC-weights-example: c
	cd $(BUILDDIR); \
	./$(EXE) ../example-data/cube.csv cube-output.csv 0 1000 0.0001 2 2 1 ../example-data/cube-weights.csv; \

run-cubeCUDA-example: cuda
	cd $(BUILDDIR); \
	./$(CUEXE) ../example-data/cube.csv cube-output.csv 0 1000 0.0001 2 2 1; \

install-prerequisites:
	sudo apt-get install cmake libopenblas-dev liblapack-dev; \
	tar xvf armadillo-8.500.1.tar.xz; \
	cd armadillo-8.500.1/; \
	mkdir build; \
	cd build; \
	cmake ..; \
	make; \
	sudo make install;

doc:
	cd src; \
	doxygen doxy_config.conf

clean:
	rm -f $(BUILDDIR)/$(EXE); \
	rm -f $(BUILDDIR)/$(CUEXE); \
	rm -rf $(BUILDDIR); \
