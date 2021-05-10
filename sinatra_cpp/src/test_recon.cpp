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
    /*
    vector<vector<double> > directions = generate_equidistributed_cones(10,0.8,4,false);
    dvec rates_dvec;
    rates_dvec.load("rates_10_4_0.8_60.txt",raw_ascii);
    vector<double> rates_vert_frames;
    for(int i_frame=0;i_frame<100;i_frame++)
    { 
        string mshpath_A = "../msh/WT_offset_0/";
        string mshfile = mshpath_A + "WT_frame" + to_string(i_frame) + ".msh";
        vector<double> rates_vert = reconstruct_by_sorted_threshold(mshfile,directions,rates_dvec,60,4);
        if (rates_vert_frames.empty())
            rates_vert_frames = rates_vert;
        else
            for(int i=0;i<rates_vert.size();i++)
                rates_vert_frames[i] += rates_vert[i];
    }
    ofstream rfile;
    rfile.open("rates_vert.txt");
    for(int i=0;i<rates_vert_frames.size();i++)
    {
        rates_vert_frames[i] /= 100;
        rfile << rates_vert_frames[i] << '\n';
    }
    rfile.close();
    */
    dvec rates_vert;
    rates_vert.load("rates_vert.txt");
    add_rate_pdb("../pdb/WT_offset_0/WT_frame0.pdb","rates_vert.pdb",rates_vert);
    return 0;
}

