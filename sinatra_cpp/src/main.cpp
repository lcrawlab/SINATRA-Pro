// g++ main.cpp -o main -larmadillo -fopenmp -I. -lstdc++fs

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
#include "gp.h"
#include "reconstruction.h"

using namespace std;
using namespace arma;

int main (int argc, char *argv[])
{
    
    int number_of_threads = 4;
    omp_set_num_threads(number_of_threads);
    
    string ec_type = "DECT";
    int n_filtration = 60;
    int n_cone = 1;
    int n_direction_per_cone = 1;
    double cap_radius = 0.80;
    double cutoff = 2.0;
    
    string pdbpath_A = "../../python_script/WT_R164S_65_213/pdb/WT_offset_0/";
    string mshpath_A = "../msh/WT_offset_0/";
    string pdbpath_B = "../../python_script/WT_R164S_65_213/pdb/R164S_offset_0/";
    string mshpath_B = "../msh/R164S_offset_0/";
    /** 
    cout << "Path to PDB files for protein A: ";
    cin >> pdbpath_A;
    cout << "Path to PDB files for protein B: ";
    cin >> pdbpath_B;
    cout << "Path to MSH files for protein A: ";
    cin >> mshpath_A;
    cout << "Path to MSH files for protein B: ";
    cin >> mshpath_B;
    **/
    string pdbfile_recon = pdbpath_A + "WT_frame0.pdb";
    string y_filename = "label.txt";

    string tmp;
    cout << "Type of Euler characteristic measure (DECT/ECT/SECT, default: DECT): ";
    cin >> tmp;
    if (tmp.compare("SECT") == 0) { ec_type = "SECT"; }
    else if (tmp.compare("ECT") == 0) { ec_type = "ECT"; }
    else {ec_type = "DECT"; }

    cout << "Radius cutoff for simplicial construction: ";
    cin >> cutoff;
    cout << "Number of cones: ";
    cin >> n_cone;
    cout << "Number of directions per cone: ";
    cin >> n_direction_per_cone;
    cout << "Cap radius: ";
    cin >> cap_radius;
    cout << "Number of filtration steps: ";
    cin >> n_filtration;

    string filesuffix;
    stringstream ss_filesuffix;
    ss_filesuffix << '_' << n_cone << '_' << n_direction_per_cone << '_' << setprecision(4) << cap_radius << '_' << n_filtration;
    ss_filesuffix >> filesuffix;

    pdb_to_mesh(pdbpath_A,pdbpath_B,mshpath_A,mshpath_B,cutoff);
    vector<vector<double> > directions = generate_equidistributed_cones(n_cone,cap_radius,n_direction_per_cone,false);
    cout << directions.size() << " Directions generated\n";
    compute_ec_curve_multiple_files(mshpath_A,mshpath_B,directions,n_filtration,ec_type,y_filename,filesuffix);  
    
    dmat X;
    dvec y;
    X.load(ec_type+filesuffix+".txt",raw_ascii);
    y.load(y_filename,raw_ascii);
    
    dvec rates;
    rates = find_rate_variables_with_other_sampling_methods(X,y);
    rates.save("rates_"+filesuffix+".txt",raw_ascii);    
    //rates.load("rates"+filesuffix+".txt",raw_ascii);
    vector<double> vec_rates_vert = reconstruct_multiple_mesh(mshpath_A,rates,directions,n_filtration,n_direction_per_cone);
    dvec rates_vert = conv_to<dvec>::from(vec_rates_vert);
    rates_vert.save("rates_vert"+filesuffix+".txt",raw_ascii); 
    //rates_vert.load("rates_vert"+filesuffix".txt,raw_ascii);
    add_rate_pdb(pdbfile_recon,"rates_vert"+filesuffix+".pdb",rates_vert);
    return 0;
}

