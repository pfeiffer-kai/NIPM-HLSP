# Primal-dual interior point method for hierarchical least-squares programs with linear constraints

This software provides an interior-point method for hierarchical least-squares programs with linear constraints.
The corresponding paper can be found here "NIPM-HLSP: An Efficient Interior-Point Method for Hierarchical Least-Squares Programs" ( https://arxiv.org/pdf/2106.13602.pdf ).

NIPM-HLSP solves hierarchical least-squares problems of the form 

min. x: ||A_l x - b_l||_2^2 , l = 1,...,p

s.t \underline{A}_{l-1} x - \underline{b}_{l-1} = w_{l-1}^*

A_l are dense matrices of any form and can be rank-deficient.

## Install required packages

Please make sure to install the following packages:

- Eigen >= 3.2.10

## Getting started

Clone the repository to the desired work folder <work_folder>
```
cd <work_folder>
mkdir build
cd build
cmake ..
ccmake .
```
set CMakeBuildType to Release
set COMPILE_TESTS to OFF
```
make && sudo make install
ccmake .
```
set COMPILE_TESTS to ON
```
make && sudo make install
```
run ./build/tests/test1

## License

Copyright (c) 2022, Kai Pfeiffer

## Authors

- Kai Pfeiffer (<pfeifferkaimartin@gmail.com>)
