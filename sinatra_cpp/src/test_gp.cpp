#include <armadillo>
#include "tools.h"
#include "gp.h"
#include <omp.h>

int main()
{
    const int number_of_threads = 4;
    omp_set_num_threads(number_of_threads);
    dmat X;
    dvec y;
    X.load("DECT_10_4_0.8_60.txt",raw_ascii);
    y.load("label_WT_R164S.txt",raw_ascii);
    dvec rates = find_rate_variables_with_other_sampling_methods(X,y);
    rates.save("rates_10_4_0.8_60.txt",raw_ascii);
    return 0;
}


