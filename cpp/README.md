# SINATRA Pro C++ version 

## Package Details

This version of the code for implementing the SINATRA Pro pipeline was written in C++, which provides equivalent result as the Python3 version, but is more memory efficient. As part of this procedure:

1. Most athematical calculations are performed using Armadillo.
2. Inference for the Gaussian process classification (GPC) model was done using elliptical slice sampling (Murray, Prescott, and MacKay 2010).
3. Association measures are computed for the Euler characteristic curves. We use the relative centrality criterion (RATE), which is a variable selection measure for nonlinear and nonparametric statistical methods (Crawford et al. 2019; Ish-Horowicz et al. 2019).

## Dependencies

The SINATRA Pro C++ version depends on these following packages.



#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <experimental/filesystem>

## Compilation

using the follwing command

        g++ main.cpp -o main -larmadillo -fopenmp -I. -lstdc++fs

usage:

        ./main

then follow the command prompt to input the necessary parameters.


## Relevant Citations

Wai Shing Tang*, Gabriel Monteiro da Silva*, Henry Kirveslahti, Erin Skeens, Bibo Feng, Timothy Sudijono, Kevin K. Yang, Sayan Mukherjee, Brenda Rubenstein, and Lorin Crawford. Topological Data Analytic Approach for Discovering Biophysical Signatures in Protein Dynamics.

## License

GNU General Public License v3.0


