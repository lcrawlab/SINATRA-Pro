#include <iostream>
#include <cmath>
#include <string>

using namespace std;

vector<vector<double> > compute_ec_curve(mesh meshA, vector<vector<double> > directions, double ball_radius, unsigned int n_filtration, string ec_type)
{
    vector<vector<double> > ec_curves;
    for(int i_dir=0;i_dir<directions.size();i_dir++)
    {
        vector<double> empty(n_filtration,0.0);
        ec_curves.push_back(empty);
    }
    int n_vertex = meshA.coords.size();
    int n_edge = meshA.edge_list.size();
    int n_face = meshA.face_list.size();
    //vector<double> bins = linspace(-ball_radius,ball_radius,n_filtration);
    double binsize = 2*ball_radius/(n_filtration-1);
    #pragma omp parallel
    for(int i_dir=0;i_dir<directions.size();i_dir++)
    {
        vector<unsigned int> vert_func(n_vertex,0.0);
        for(int i=0;i<n_vertex;i++)
            vert_func[i] = (int)floor((dot(meshA.coords[i],directions[i_dir])+ball_radius)/binsize);
        
        vector<double> V(n_filtration,0), E(n_filtration,0), F(n_filtration,0);
        
        for(int i_vert=0;i_vert<n_vertex;i_vert++)
        {
            V[vert_func[i_vert]] += 1.0;
        }
        for(int i_edge=0;i_edge<n_edge;i_edge++)
        {
            int a = meshA.edge_list[i_edge][0];
            int b = meshA.edge_list[i_edge][1];
            if (vert_func[a] > vert_func[b])
                E[vert_func[a]] += 1.0;
            else
                E[vert_func[b]] += 1.0;
        }
        for(int i_face=0;i_face<n_face;i_face++)
        {
            int max_func = 0;
            for(int j=0;j<meshA.face_list[i_face].size();j++)
                if (vert_func[meshA.face_list[i_face][j]] > max_func)
                    max_func = vert_func[meshA.face_list[i_face][j]];
            F[max_func] += 1.0;
        }
        for(int i_fil=0;i_fil<n_filtration;i_fil++)
        {
            ec_curves[i_dir][i_fil] = V[i_fil] - E[i_fil] + F[i_fil];
        }
        if (not ec_type.compare("DECT"))
        {
            for(int i_fil=0;i_fil<n_filtration;i_fil++)
                ec_curves[i_dir][i_fil] /= binsize;
        }
        else
        {
            double cumsum = 0.0;
            for(int i_fil=0;i_fil<n_filtration;i_fil++)
            {
                cumsum += ec_curves[i_dir][i_fil];
                ec_curves[i_dir][i_fil] = cumsum;
            }
            if (not ec_type.compare("SECT"))
            {
                double mean = average(ec_curves[i_dir]);
                for(int i_fil=0;i_fil<n_filtration;i_fil++)
                    ec_curves[i_dir][i_fil] -= mean;
                double cumsum = 0.0;
                for(int i_fil=0;i_fil<n_filtration;i_fil++)
                {
                    cumsum += ec_curves[i_dir][i_fil];
                    ec_curves[i_dir][i_fil] = cumsum;
                    ec_curves[i_dir][i_fil] *= binsize / n_filtration;
                }
            }
            else if (ec_type.compare("ECT"))
            {
                cout << "Please choose from one of the EC type: ECT / DECT / SECT\n";
            }
        }
    }
    return ec_curves;
}


