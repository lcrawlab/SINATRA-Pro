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

#include <thread>

using namespace std;
using namespace arma;

int main (int argc, char *argv[])
{
    
    string ec_type = "DECT";
    int n_filtration = 60;
    int n_cone = 4;
    int n_direction_per_cone = 4;
    double cap_radius = 0.80;
    double cutoff = 4.0;
    string y_filename = "label.txt"; 
   
    string prot_A = "protA";
    string prot_B = "protB";
    string directory = "output";
    string pdbpath_A = "data/pdb/WT_offset_0/";
    string pdbpath_B = "data/pdb/R164S_offset_0/";
    string pdbfile_recon = pdbpath_A + "WT_frame0.pdb";
    
    string input;
    cout << "Name for protein A: ";
    getline(cin,input);
    if ( !input.empty() ) {
        istringstream stream(input);
        stream >> prot_A;
    }
    cout << "Name for protein B: ";
    getline(cin,input);
    if ( !input.empty() ) {
        istringstream stream(input);
        stream >> prot_B;
    }
    cout << "Path to PDB files for protein A: ";
    getline(cin,input);
    if ( !input.empty() ) {
        istringstream stream(input);
        stream >> pdbpath_A;
    }
    cout << "Path to PDB files for protein B: ";
    getline(cin,input);
    if ( !input.empty() ) {
        istringstream stream(input);
        stream >> pdbpath_B;
    }
    cout << "PDB file for visualization protein A: ";
    getline(cin,input);
    if ( !input.empty() ) {
        istringstream stream(input);
        stream >> pdbfile_recon;
    }
    cout << "Directory for output files (default: output): ";
    getline(cin,input);
    if ( !input.empty() ) {
        istringstream stream(input);
        stream >> directory;
    }

    string mshpath_A = directory + "msh/" + prot_A + "/";
    string mshpath_B = directory + "msh/" + prot_B + "/";
   
    string tmp;
    cout << "Type of Euler characteristic measure (DECT/ECT/SECT, default: DECT): ";
    getline(cin,input);
    if ( !input.empty() ) {
        istringstream stream(input);
        stream >> tmp;
    }
    if (tmp.compare("SECT") == 0) { ec_type = "SECT"; }
    else if (tmp.compare("ECT") == 0) { ec_type = "ECT"; }
    else {ec_type = "DECT"; }

    cout << "Radius cutoff for simplicial construction (default: 4.0): ";
    cin >> cutoff;
    cout << "Number of cones: (default: 4)";
    cin >> n_cone;
    cout << "Number of directions per cone (default: 4): ";
    cin >> n_direction_per_cone;
    cout << "Cap radius: ";
    cin >> cap_radius;
    cout << "Number of filtration steps (default: 60): ";
    cin >> n_filtration;
    
    int number_of_threads = -1;
    cout << "Number of CPU cores used for calculation (default: -1, which automatically detect all cores): ";
    cin >> number_of_threads;
    if (number_of_threads == -1)
        const auto processor_count = thread::hardware_concurrency();    //may return 0 when not able to detect
    omp_set_num_threads(number_of_threads); 

    string filesuffix;
    stringstream ss_filesuffix;
    ss_filesuffix << '_' << n_cone << '_' << n_direction_per_cone << '_' << setprecision(4) << cap_radius << '_' << n_filtration;
    ss_filesuffix >> filesuffix;
    
    pdb_to_mesh(pdbpath_A,pdbpath_B,mshpath_A,mshpath_B,cutoff);
    vector<vector<double> > directions = generate_equidistributed_cones(n_cone,cap_radius,n_direction_per_cone,false);
    cout << directions.size() << " Directions generated\n";
    
    string ec_filename = directory + "/" + ec_type + filesuffix + ".txt";
    compute_ec_curve_multiple_files(mshpath_A,mshpath_B,directions,n_filtration,ec_type,y_filename,ec_filename);  
    
    dmat X;
    dvec y;
    
    X.load(directory+"/"+ec_type+filesuffix+".txt",raw_ascii);
    y.load(y_filename,raw_ascii);

    dvec rates;
    rates = find_rate_variables_with_other_sampling_methods(X,y);
    rates.save(directory+"/rates"+filesuffix+".txt",raw_ascii); 

    vector<double> vec_rates_vert = reconstruct_multiple_mesh(mshpath_A,rates,directions,n_filtration,n_direction_per_cone);
    
    dvec rates_vert = conv_to<dvec>::from(vec_rates_vert);
    rates_vert.save(directory+"/rates_vert"+filesuffix+".txt",raw_ascii); 
    add_rate_pdb(pdbfile_recon,directory+"/rates_vert"+filesuffix+".pdb",rates_vert);
    
    return 0;
}

