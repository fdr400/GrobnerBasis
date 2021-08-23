# GrobnerBasis
Console application that implements various methods of working with polynomial systems, in particular,<br/> various optimizations of the Buchberger algoritm for finding the Grobner basis of the ideal.

### Requirements
- latest version [g++](https://gcc.gnu.org/) compiler with c++17 support.
- installed c++ libraries such as [boost](https://www.boost.org/) and [gflags](https://github.com/gflags/gflags).
- installed [gtest](https://github.com/google/googletest) library (if you have a desire to run tests).

### Installation guide
- download the repository.
- to compile the program, from the downloaded folder run `g++ main.cpp -std=c++17 -o bin/grobner -lgflags`<br/> (obligatory specify `-DFIELD` flag with number from 1 to 3 to choose the field of coefficients).
- to compile tests, from the downloaded folder run `g++ tests/main.cpp -std=c++17 -o bin/tester -lgtest -lgtest_main`.

### How to run
- to run the program, invoke `bin/grobner` with specified desired flags.
- to run tests, invoke `bin/tester`.

### Usage guide
- to look at help message, run `bin/grobner --helpon=main`
- to expand modes, run `bin/grobner --modes=expand_modes`
- to know optimizations meanings, run `bin/grobner --modes=expand_optimization`