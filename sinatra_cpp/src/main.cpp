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

int main ()
{
    const int number_of_threads = 4;
    omp_set_num_threads(number_of_threads);

    string pdbpath_A = "../pdb/WT_offset_0/";
    string mshpath_A = "../msh/WT_offset_0/";
    string pdbpath_B = "../pdb/R164S_offset_0/";
    string mshpath_B = "../msh/R164S_offset_0/";
    
    double max_radius = 0.0;
    for(int i_frame=0;i_frame<100;i_frame++)
    {
        mesh meshA;
        string filename = pdbpath_A + "WT_frame" + to_string(i_frame) + ".pdb";
        meshA.coords = read_pdb(filename);
        double radius = meshA.calc_max_radius();
        if (radius > max_radius)
            max_radius = radius;
    }
    for(int i_frame=0;i_frame<100;i_frame++)
    {
        mesh meshA;
        string filename = pdbpath_B + "R164S_frame" + to_string(i_frame) + ".pdb";
        meshA.coords = read_pdb(filename);
        double radius = meshA.calc_max_radius();
        if (radius > max_radius)
            max_radius = radius;
    }
    cout << "Max radius = " << max_radius << '\n';

    for(int i_frame=0;i_frame<100;i_frame++)
    {
        cout << "Constructing Mesh for frame " << i_frame << "..." << '\r';
        mesh meshA;
        string filename = pdbpath_A + "WT_frame" + to_string(i_frame) + ".pdb";
        meshA.coords = read_pdb(filename);
        meshA.neighbor_grid_search(4.0);
        meshA.neighbor_to_edge();
        meshA.neighbor_to_face();
        meshA.normalize_to(max_radius);
        string mshfile = mshpath_A + "WT_frame" + to_string(i_frame) + ".msh";
        meshA.write_mesh(mshfile);
        cout << flush;
    }
    
    for(int i_frame=0;i_frame<100;i_frame++)
    {
        cout << "Constructing Mesh for frame " << i_frame << "..." << '\r';
        mesh meshA;
        string filename = pdbpath_B + "R164S_frame" + to_string(i_frame) + ".pdb";
        meshA.coords = read_pdb(filename);
        meshA.neighbor_grid_search(4.0);
        meshA.neighbor_to_edge();
        meshA.neighbor_to_face();
        meshA.normalize_to(max_radius);
        string mshfile = mshpath_B + "R164S_frame" + to_string(i_frame) + ".msh";
        meshA.write_mesh(mshfile);
        cout << flush;
    }
    cout << "Mesh Construction Complete.      \n";

    vector<vector<double> > directions = generate_equidistributed_cones(10,0.8,4,false);
    cout << directions.size() << " Directions generated\n";
    ofstream ecfile, yfile;
    ecfile.open("DECT_10_4_0.8_60.txt");
    yfile.open("label_WT_R164S.txt");
    //vector<vector<double> > ec_matrix;
    //vector<int> label;
    for(int i_frame=0;i_frame<100;i_frame++)
    {
        mesh meshA;
        string mshfile = mshpath_A + "WT_frame" + to_string(i_frame) + ".msh";
        meshA.read_mesh(mshfile);
        vector<vector<double> > ec_curves = compute_ec_curve(meshA,directions,1.0,60,"DECT");
        //vector<double> flattened(begin(ec_curves[0]), end(ec_curves[0]));
        //for(int i_dir=1;i_dir<ec_curves.size();i_dir++)
        //    flattened.insert(end(flattened), begin(ec_curves[i_dir]), end(ec_curves[i_dir]));
        //ec_matrix.push_back(flattened);
        for(int i_dir=0;i_dir<ec_curves.size();i_dir++)
        {
            for(int i_fil=0;i_fil<ec_curves[i_dir].size();i_fil++)
            {
                ecfile << setprecision(6) << ec_curves[i_dir][i_fil] << ' ';
            }
        }
        ecfile << '\n';
        yfile << "0\n";
        //label.push_back(0);
    }
    for(int i_frame=0;i_frame<100;i_frame++)
    {
        mesh meshA;
        string mshfile = mshpath_B + "R164S_frame" + to_string(i_frame) + ".msh";
        meshA.read_mesh(mshfile);
        vector<vector<double> > ec_curves = compute_ec_curve(meshA,directions,1.0,60,"DECT");
        for(int i_dir=0;i_dir<ec_curves.size();i_dir++)
        {
            for(int i_fil=0;i_fil<ec_curves[i_dir].size();i_fil++)
            {
                ecfile << setprecision(6) << ec_curves[i_dir][i_fil] << ' ';
            }
        }
        ecfile << '\n';
        yfile << "1\n";
        //vector<double> flattened(begin(ec_curves[0]), end(ec_curves[0]));
        //for(int i_dir=1;i_dir<ec_curves.size();i_dir++)
        //    flattened.insert(end(flattened), begin(ec_curves[i_dir]), end(ec_curves[i_dir]));
        //ec_matrix.push_back(flattened);
        //label.push_back(1);
    }
    cout << "EC Calculation Complete.      \n";
    /*
    dmat X;
    dvec y;
    X.load("DECT_10_4_0.8_60.txt",raw_ascii);
    y.load("label_WT_R164S.txt",raw_ascii);
    dvec rates = find_rate_variables_with_other_sampling_methods(X,y);
    rates.save("rates_10_4_0.8_60.txt",raw_ascii);
    
    vector<vector<double> > directions = generate_equidistributed_cones(10,0.8,4,false);
    dvec rates_dvec;
    rates_dvec.load("rates_10_4_0.8_60.txt",raw_ascii);
    string mshpath_A = "../msh/WT_offset_0/";
    string mshfile = mshpath_A + "WT_frame0.msh";
    vector<double> rates_vert = reconstruct_by_sorted_threshold(mshfile,directions,rates_dvec,50,4);
    ofstream rfile;
    rfile.open("rates_vert_10_4_0.8_60.txt");
    for(int i=0;i<rates_vert.size();i++)
        rfile << rates_vert[i] << '\n';
    rfile.close();
    */
    return 0;
}

