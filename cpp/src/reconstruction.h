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
                unsigned int vert_func = (unsigned int)ceil((dot(meshA.coords[i_vert],directions[k])+ball_radius)/binsize);
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

int add_rate_pdb (string infilename, string outfilename, dvec rates_vert)
{
    int n_vert = rates_vert.n_elem;
    uvec indices = sort_index(rates_vert);
    vector<double> score_vert(n_vert,0.0);
    for(int i=0;i<n_vert;i++)
    {
        score_vert[indices[i]] = ((double)(n_vert-i)*100.0)/n_vert;
    }
    ifstream infile;
    infile.open(infilename);
    ofstream outfile;
    outfile.open(outfilename);
    string line;
    if (infile.is_open() and outfile.is_open())
    {
        int i_atom = 0;
        while ( getline (infile,line) )
        { 
            string str1 = line.substr(0,4);
            if ( str1.compare("ATOM") == 0 )
            {
                stringstream ss;
                ss.width(6);
                ss << setprecision(3) << score_vert[i_atom];
                string str;
                ss >> str;
                for(int i=str.length();i<6;i++)
                {
                    str = str + ' ';
                }
                line.replace(61,6,str);
                outfile << line << endl;
                i_atom++;
            }            
            else
            {
                outfile << line << endl;
            }

        }
        infile.close();
        outfile.close();
    }
    return 0;
}

vector<double> reconstruct_multiple_mesh(string mshpath,dvec rates_dvec, vector<vector<double> > directions, int n_filtration, int n_direction_per_cone)
{
    vector<double> rates_vert_frames;
    int n_mesh = 0;
    for (const auto & file : experimental::filesystem::directory_iterator(mshpath))
    {
        string filename = file.path();
        string extension = filename.substr(filename.size()-4,4);
        if (not extension.compare(".msh"))
        {
            vector<double> rates_vert = reconstruct_by_sorted_threshold(filename,directions,rates_dvec,n_filtration,n_direction_per_cone);
            if (rates_vert_frames.empty())
                rates_vert_frames = rates_vert;
            else
                for(int i=0;i<rates_vert.size();i++)
                    rates_vert_frames[i] += rates_vert[i];
            n_mesh++;
        }
    }
    for(int i=0;i<rates_vert_frames.size();i++)
        rates_vert_frames[i] /= n_mesh;
    return rates_vert_frames;
}
