~~~~~~~~~~~~~ 
/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2010-2012                                                 |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                | 
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      Via Mesiano 77, I-38050 Trento, Italy                               |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
 |      version: 0.2                                                        |
 |                                                                          |
\*--------------------------------------------------------------------------*/
~~~~~~~~~~~~~ 

GF2 toolkit,
a library for manipulation of (full) matrices and vector with entries in F2.

to compile the library you need

- cmake
- m4ri (only for comparison with SAGE, http://m4ri.sagemath.org/)

to compile library

cd cmake
cmake ..
make

a number of test to check library is now compiled, run and enjoy.

The principal header to use is 

GF2toolkit_LU.hh

which contains the description of the class LU.
This class implement the factorrization described in

Fast matrix decomposition in F2
by Enrico bertolazzi, Anna Rimoldi
http://arxiv.org/abs/1209.5198

an example of usage is in 

testRankComputation.cc

the code shuold be easy enough for a skilled developer to understand
how to use. I am planning to add DOXYGEN documentation to the library.

Enrico Bertolazzi.
