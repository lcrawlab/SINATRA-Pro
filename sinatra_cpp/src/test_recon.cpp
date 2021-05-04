#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <omp.h>
#include <armadillo>

#include "tools.h"
#include "mesh.h"
#include "direction.h"
#include "euler.h"
#include "reconstruction.h"

int main()
{
    vector<vector<double> > directions = generate_equidistributed_cones(5,0.8,4,false);
    dvec rates_dvec;
    rates_dvec.load("rates_5_4_0.8_50.txt",raw_ascii);
    string mshpath_A = "../msh/WT_offset_0/";
    string mshfile = mshpath_A + "WT_frame0.msh";
    vector<double> rates_vert = reconstruct_by_sorted_threshold(mshfile,directions,rates_dvec,50,4);
    ofstream rfile;
    rfile.open("rates_vert.txt");
    for(int i=0;i<rates_vert.size();i++)
        rfile << rates_vert[i] << '\n';
    rfile.close();
    return 0;
}

