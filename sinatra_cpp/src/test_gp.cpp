#include <armadillo>
#include "gp.h"
#include <omp.h>

int main()
{
    const int number_of_threads = 4;
    omp_set_num_threads(number_of_threads);
    dmat X;
    dvec y;
    //X.load("../output/DECT_20_8_0.8_100.txt",raw_ascii);
    X.load("DECT_5_4_0.8_50.txt",raw_ascii);
    y.load("label_WT_R164S.txt",raw_ascii);
    dvec rates = find_rate_variables_with_other_sampling_methods(X,y);
    rates.save("rates_5_4_0.8_50.txt",raw_ascii);
    return 0;
}


