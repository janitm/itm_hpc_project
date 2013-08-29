itm_hpc_project
===============

PROJECT ABSTRACT
================

This project aims to improve an isolated serial advection subroutine in the above mentioned climate model code, which is used at the Department of Applied Environmental Science (itm) at Stockholm University. The first step is to port the isolated subroutine from the currently existing Fortran 77 version to C. Aside from simple optimizations in advection routine, such as replacing repetitions of expensive computations and testing optimization flags of compilers, we plan to optimize the nested for loop structure of the code with OpenMP. Further optimizations we want to try are for loop-unrolling and increasing the throughput with SIMD instructions. Due to the physics of the model, there exist a lot of data-dependencies in the code, which will make it difficult to implement MPI. Moreover, the model is run on several different architectures. Nevertheless, we will keep an open mind about MPI and try to implement it if we think that it is feasible.

FOLDER STRUCTURE
================

f contains the original Fortran 77 code
c contains the C code
latex contains the LaTex code
fig contains the figures and plots generated to be included in the LaTex code
data contains the results from numerical experiments
