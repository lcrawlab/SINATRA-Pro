#include <iostream>
#include <string>
#include <armadillo>
#include <vector>

using namespace std;
using namespace arma;

vector<double> reconstruct_by_sorted_threshold(string mshfile, vector<vector<double> > directions, dvec rates_dvec, int n_filtration, int n_direction_per_cone, double ball_radius = 1.0, bool by_rank = false, bool verbose = true)
{
    vector<double> rates = conv_to< vector<double> >::from(rates_dvec);
    mesh meshA;
    meshA.read_mesh(mshfile);
    int n_direction = directions.size();
    int n_cone = n_direction / n_direction_per_cone;
    int n_vertex = meshA.coords.size();
    double binsize = 2*ball_radius/(n_filtration-1);
    vector<double> rates_vert(n_vertex,0.0);
    for(int i_vert=0;i_vert<n_vertex;i_vert++)
    {
        double max_rates_cone = -9e99;
        for(int i=0;i<n_cone;i++)
        {
            double min_rates_direction = 9e99;
            for(int j=0;j<n_direction_per_cone;j++)
            {
                int k = i * n_direction_per_cone + j;
                unsigned int vert_func = (unsigned int)floor((dot(meshA.coords[i_vert],directions[k])+ball_radius)/binsize);
                if (rates[k*n_filtration+vert_func] < min_rates_direction)
                    min_rates_direction = rates[k*n_filtration+vert_func];                
            }
            if (min_rates_direction > max_rates_cone)
                max_rates_cone = min_rates_direction;
        }
        rates_vert[i_vert] = max_rates_cone;
    }
    return rates_vert;
}



