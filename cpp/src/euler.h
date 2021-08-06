#include <iostream>
#include <cmath>
#include <string>

using namespace std;

// This function computes the EC curves...
vector<vector<double> > compute_ec_curve(mesh meshA, vector<vector<double> > directions, double ball_radius, unsigned int n_filtration, string ec_type)
{
    // Initialize variables
    vector<vector<double> > ec_curves;
    for(int i_dir=0;i_dir<directions.size();i_dir++)
    {
        vector<double> empty(n_filtration,0.0);
        ec_curves.push_back(empty);
    }
    int n_vertex = meshA.coords.size();
    int n_edge = meshA.edge_list.size();
    int n_face = meshA.face_list.size();
    double binsize = 2*ball_radius/(n_filtration-1);
    
    // #pragma omp parallel
    for(int i_dir=0;i_dir<directions.size();i_dir++)  // For each direction
    {
        // Bin each vertex into a filtration bin and label the vertex with a bin number (vertex function)
        vector<unsigned int> vert_func(n_vertex,0.0);
        for(int i_vert=0;i_vert<n_vertex;i_vert++)
            vert_func[i_vert] = (unsigned int)ceil((dot(meshA.coords[i_vert],directions[i_dir])+ball_radius)/binsize);
        
        vector<double> V(n_filtration,0.0), E(n_filtration,0.0), F(n_filtration,0.0);
       
        // Count number of vertices in each filtration bin 
        for(int i_vert=0;i_vert<n_vertex;i_vert++)
        {
            V[vert_func[i_vert]] += 1.0;
        }

        // Count number of edges born in each filtration bin
        for(int i_edge=0;i_edge<n_edge;i_edge++)
        {
            int a = meshA.edge_list[i_edge][0];
            int b = meshA.edge_list[i_edge][1];
            if (vert_func[a] > vert_func[b])
                E[vert_func[a]] += 1.0;
            else
                E[vert_func[b]] += 1.0;
        }

        // Count number of faces born in each filtration bin
        for(int i_face=0;i_face<n_face;i_face++)
        {
            int max_func = 0;
            for(int j=0;j<meshA.face_list[i_face].size();j++) {
                if (vert_func[meshA.face_list[i_face][j]] > max_func)
                    max_func = vert_func[meshA.face_list[i_face][j]];
            }
            F[max_func] += 1.0;
        }

        // EC = V - E + F
        for(int i_fil=0;i_fil<n_filtration;i_fil++)
        {
            ec_curves[i_dir][i_fil] = V[i_fil] - E[i_fil] + F[i_fil];
        }
        if (not ec_type.compare("DECT")) // If ec_type == "DECT"
        {
            for(int i_fil=0;i_fil<n_filtration;i_fil++)
                ec_curves[i_dir][i_fil] /= binsize; 
        }
        else
        {
            // Integrate to compute EC curve
            double cumsum = 0.0;
            for(int i_fil=0;i_fil<n_filtration;i_fil++)
            {
                cumsum += ec_curves[i_dir][i_fil];
                ec_curves[i_dir][i_fil] = cumsum;
            }
            // If ec_type == "ECT", anything below is false
            
            if (not ec_type.compare("SECT")) // If ec_type == "SECT"
            {
                // Subtract mean from EC curve
                double mean = average(ec_curves[i_dir]);
                for(int i_fil=0;i_fil<n_filtration;i_fil++)
                    ec_curves[i_dir][i_fil] -= mean;

                // Integrate again for SECT               
                double cumsum = 0.0;
                for(int i_fil=0;i_fil<n_filtration;i_fil++)
                {
                    cumsum += ec_curves[i_dir][i_fil];
                    ec_curves[i_dir][i_fil] = cumsum;
                    ec_curves[i_dir][i_fil] *= binsize / n_filtration;
                }
            }
            else if (ec_type.compare("ECT")) // If ec_type is not ECT / DECT / SECT
            {
                cout << "Please choose from one of the EC type: ECT / DECT / SECT\n";
            }
        }
    }
    return ec_curves;
}

int compute_ec_curve_multiple_files(string mshpath_A, string mshpath_B, vector<vector<double> > directions, int n_filtration, string ec_type="DECT", string y_filename="label.txt", string ec_filename="")
{
    double ball_radius = 1.0;
    ofstream ecfile, yfile;
    ecfile.open(ec_filename);
    yfile.open(y_filename);
    for (const auto & file : experimental::filesystem::directory_iterator(mshpath_A))
    {
        string filename = file.path();
        string extension = filename.substr(filename.size()-4,4);        
        if (not extension.compare(".msh"))
        {
            string filename_short = filename.substr(filename.rfind("/") + 1);
            filename_short = filename_short.substr(0,filename_short.size()-4);
            cout << "Calculation EC for " << filename_short << "..." << '\r';
            cout << flush;
            mesh meshA;
            meshA.read_mesh(filename);            
            vector<vector<double> > ec_curves = compute_ec_curve(meshA,directions,ball_radius,n_filtration,ec_type);
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
        }
    }
    for (const auto & file : experimental::filesystem::directory_iterator(mshpath_B))
    { 
        string filename = file.path();
        string extension = filename.substr(filename.size()-4,4);        
        if (not extension.compare(".msh"))
        {
            string filename_short = filename.substr(filename.rfind("/") + 1);
            filename_short = filename_short.substr(0,filename_short.size()-4);
            cout << "Calculation EC for " << filename_short << "..." << '\r';
            cout << flush;
            mesh meshA;
            meshA.read_mesh(filename);
            vector<vector<double> > ec_curves = compute_ec_curve(meshA,directions,ball_radius,n_filtration,ec_type);
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
            yfile << "1\n";
            ecfile << '\n';
        }
    }
    ecfile.close();
    yfile.close();
    cout << "EC Calculation Complete.      \n";
    return 0;
}


