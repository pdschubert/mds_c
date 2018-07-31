# MDS_C: The SMACOF algorithm - scaling by majorizing a convex function

## Installing the dependencies:
`$ make install-prerequisites`

## Building the programs
Build the C version:

`$ make` or `$ make c`

Build the CUDA C version:

`$ make cuda`

Switching the compiler can be done by replacing the Makefile variable CC, e.g.:

`$ make CC=clang`

Build the documentation with:

`$ make doc`

That target will clean all auto-generated sources

`$ make clean`

## Testing the build
These targets should run without any errors.

`$ make run-cubeC-example`

`$ make run-cubeC-weights-example`

This target only works when compiling for cuda of course
`$ make run-cubeCUDA-example`

The results of these example calls can be visualized using the scripts in misc/

`$ ./scatterplot <some csv>` # allows to plot your results

Your results may have to be transposed first using 

`$ ./transposecsv.py <input csv> <output csv>`

When you would like to melt down your CPU, create a random matrix using the 
`gen-rand.py` script that allows you to generate a random matrix of arbitrary
size.

## Help
Just call the compiled programs or the python helper scripts without parameters
in order to obtain a help message, e.g.:

```
$ ./smacofC
ERROR: wrong number of arguments
####################################################################
# SMACOF - Philipp D. Schubert - philipp@it-schubert.com           #
####################################################################
usage: <prog> <input> <output> <disfunc> <maxiter> <epsilon> <resultdim> <metric_p> <print>
parameter explanation:
        input - a valid path to csv file to analyse: path
        ouput - a filename for the outputfile: path
        disfunc - dissimilarity measure: ...
                0 - euclidean
                1 - cityblock
                2 - minkowski
                3 - correlation
                4 - angularseperatin
                5 - wavehedges
                6 - bahattacharyya
                7 - soergel
                8 - braycurtis
                9 - divergence
                10 - canberra
                11 - none
maxiter - maximum number of iterations to use: a positive integer
epsilon - the maximum error value which is allowed: a positive real_t
resultdim - number of dimensions for result vectors: a positive integer
minkowski_p - a value to adjust the minkowski distance: a positive integer
print - be verbose: 0 or 1
weights - a valid path to csv file containing the weights
```
