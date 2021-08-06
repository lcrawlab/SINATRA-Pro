#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <experimental/filesystem>

using namespace std;

bool item_in_list(auto a, vector<auto> b)
{
    for(int i=0;i<b.size();i++)
    {
        if (a == b[i])
            return 1;
    }
    return 0;
}

class mesh
{
    public:
        vector<vector<double> > coords;
        vector<vector<int> > edge_list;
        vector<vector<int> > face_list;
        vector<vector<int> > neighbor_list;
        
        double calc_max_radius()
        {
            double result = 0.0, radius;
            for(int i=0;i<coords.size();i++)
            {
                radius = norm3v(coords[i]);
                if (radius > result)
                {
                    result = radius;
                }
            }
            return result;
        }
        
        int normalize_to(double max_radius) 
        {
            for(int i=0;i<coords.size();i++)
            {
                for(int j=0;j<3;j++)
                {
                    coords[i][j] /= max_radius;
                }
            }
            return 0;
        }

        int write_mesh (string filename)
        {
            ofstream mshfile;
            mshfile.open(filename);
            mshfile << coords.size() << ' ' << edge_list.size() << ' ' << face_list.size() << '\n';
            for(int i=0;i<coords.size();i++)
                mshfile << setprecision(6) << coords[i][0] << ' ' << coords[i][1] << ' ' << coords[i][2] << '\n';
            for(int i=0;i<edge_list.size();i++)
                mshfile << edge_list[i][0] << ' ' << edge_list[i][1] << '\n';
            for(int i=0;i<face_list.size();i++)
                mshfile << face_list[i][0] << ' ' << face_list[i][1] << ' ' << face_list[i][2] << '\n';
            mshfile.close();
            return 0;
        }
        int read_mesh(string filename)
        {
            ifstream mshfile;
            mshfile.open(filename);
            string line;
            if (mshfile.is_open())
            {
                int i_line = 0;
                int n_atom, n_edge, n_face;
                while ( getline (mshfile,line) )
                { 
                    if (i_line == 0)
                    {
                        stringstream ss(line);
                        ss >> n_atom >> n_edge >> n_face;
                        coords.clear();
                        edge_list.clear();
                        face_list.clear();
                    }
                    else if (i_line <= n_atom)
                    {
                        vector<double> coord(3,0);
                        stringstream ss(line);
                        ss >> coord[0] >> coord[1] >> coord[2];
                        coords.push_back(coord);
                    }
                    else if (i_line <= n_atom + n_edge)
                    {
                        vector<int> edge(2,0);
                        stringstream ss(line);
                        ss >> edge[0] >> edge[1];
                        edge_list.push_back(edge);
                    }
                    else if (i_line <= n_atom + n_edge + n_face)
                    {
                        vector<int> face(3,0);
                        stringstream ss(line);
                        ss >> face[0] >> face[1] >> face[2];
                        face_list.push_back(face);
                    }
                    i_line++;
                }
            }
            return 0;
        }

    // Neighbor grid search to construct neighbor list for each atom
    // Note that no periodic boundary condition is considered here, it assumes the molecule is whole
    int neighbor_grid_search(double cutoff)
    {
        double sqcutoff = cutoff*cutoff;
        
        // Find boundary of the molecule i.e. min/max coordinates in x,y,z direction
        vector<double> min(3, 1e99);
        vector<double> max(3, -1e99);
        for(int i=0;i<coords.size();i++)
        {
            for(int j=0;j<3;j++)
            {
                if (coords[i][j] > max[j])
                    max[j] = coords[i][j];
                if (coords[i][j] < min[j])
                    min[j] = coords[i][j];
            }
        }

        // Divide the bounded box into cells
        double cellsize = cutoff*1.1;
        vector<int> n_block(3,0);
        for(int i=0;i<3;i++)
        {
            double n_b = ceil( (max[i] - min[i]) / cellsize);
            n_block[i] = static_cast<int>(n_b) + 1;
        }
        
        // Initialize vector for atom list in each cell
        vector<vector<vector<vector<int> > > > box_atomlist;
        for(int i=0;i<n_block[0];i++)
        {
            vector<vector<vector<int> > > vvv;
            for(int j=0;j<n_block[1];j++)
            {
                vector<vector<int> > vv;
                for(int k=0;k<n_block[2];k++)
                {
                    vector<int> empty(0,0);
                    vv.push_back(empty);
                }
                vvv.push_back(vv);
            }
            box_atomlist.push_back(vvv);  
        }
        
        // Put each atom into a cell by simple division of coordinates
        for(int i=0;i<coords.size();i++)
        {
            int ind[3];
            for(int j=0;j<3;j++)
                ind[j] = floor((coords[i][j] - min[j]) / cellsize);
            box_atomlist[ind[0]][ind[1]][ind[2]].push_back(i);
        }
        
        // Initialize neighbor list
        for(int i=0;i<coords.size();i++)
        {
            vector<int> v(0,0);
            neighbor_list.push_back(v);
        }

        // Look at each cell
        for(int i=0;i<n_block[0];i++)
            for(int j=0;j<n_block[1];j++)
                for(int k=0;k<n_block[2];k++)
                {
                    // Look at the neighboring cells (including diagonal neighbors)
                    // and look at all atoms in these neighboring cells
                    vector<int> nbbox_atom;
                    for(int a=-1;a<2;a++)
                    {
                        int p = i + a;
                        if (p < 0 or p >= n_block[0])
                            continue;
                        for(int b=-1;b<2;b++)
                        {
                            int q = j + b;
                            if (q < 0 or q >= n_block[1])
                                continue;
                            for(int c=-1;c<2;c++)
                            {
                                int r = k + c;
                                if (r < 0 or r >= n_block[2])
                                    continue;
                                for(int m=0;m<box_atomlist[p][q][r].size();m++)
                                    nbbox_atom.push_back(box_atomlist[p][q][r][m]);
                            }
                        }
                    }
                    
                    // for each atom A in the center cell
                    // if an atom B in the neighboring cells has distance < cutoff
                    // mutually put each other atom into its neighbor list
                    // i.e. put atom A into atom B neighbor list
                    //      put atom B into atom A neighbor list
                    //
                    // The method saves time on computing distances on every atom to every atom
                    // which is computationally intensive O(N^2)
                    for(int m=0;m<box_atomlist[i][j][k].size();m++)
                        for(int n=0;n<nbbox_atom.size();n++)
                        {
                            int atom_a = box_atomlist[i][j][k][m];
                            int atom_b = nbbox_atom[n];
                            if (atom_a == atom_b)
                                continue;
                            double sqr = 0.0;
                            for(int p=0;p<3;p++)
                                sqr += pow(coords[atom_a][p]-coords[atom_b][p],2.0);
                            if (sqr < sqcutoff)
                            {
                                bool exists = 0;
                                for(int m_a=0;m_a<neighbor_list[atom_a].size();m_a++)
                                {
                                    if (neighbor_list[atom_a][m_a] == atom_b)
                                    {
                                        exists = 1;
                                        break;
                                    }
                                }
                                if (not exists)
                                    neighbor_list[atom_a].push_back(atom_b);
                                exists = 0;
                                for(int m_b=0;m_b<neighbor_list[atom_b].size();m_b++)
                                {
                                    if (neighbor_list[atom_b][m_b] == atom_a)
                                    {
                                        exists = 1;
                                        break;
                                    }
                                }                           
                                if (not exists)
                                    neighbor_list[atom_b].push_back(atom_a);
                            }
                        }
                }
        return 0;
    }

    vector<vector<int> > neighbor_to_edge()
    {
        for(int n=0;n<neighbor_list.size();n++) // For each atom n
        {
            for(int i=0;i<neighbor_list[n].size();i++) // basically connect to every neighbor into edge
            {
                int a = neighbor_list[n][i];
                if (a > n)
                {
                   vector<int> edge = {n,a};
                   edge_list.push_back(edge);
                }
            }
        }
        return edge_list;
    }

    vector<vector<int> > neighbor_to_face()
    {
        for(int n=0;n<neighbor_list.size();n++) // For each atom n
        {
            for(int i=0;i<neighbor_list[n].size();i++) // look at its neighbor list
            {
                int a = neighbor_list[n][i];  // a connects n
                if (a < n)
                    continue;
                for(int j=0;j<neighbor_list[a].size();j++)  // look at a's neighbor
                {
                    int b = neighbor_list[a][j]; // b connects a
                    if (b < a)
                        continue;
                    if (item_in_list(n,neighbor_list[b]))   // if b connects back to n, construct triangular face
                    {
                        vector<int> face = {n,a,b};
                        face_list.push_back(face);
                    }
                }
            }
        }
        return face_list;
    }
};

vector<vector<double> > read_pdb (string filename)
{
    ifstream pdbfile;
    pdbfile.open(filename);
    string line;
    vector<vector<double> > coords;
    if (pdbfile.is_open())
    {
        while ( getline (pdbfile,line) )
        { 
            string str1 = line.substr(0,4);
            if ( str1.compare("ATOM") == 0 )
            {
                string substr = line.substr(30,24);
                stringstream ss(substr);
                vector<double> pos;
                for(int i=0;i<3;i++)
                {
                    double x;
                    ss >> x;
                    pos.push_back(x);
                }
                coords.push_back(pos);
            }            
        }
        pdbfile.close();
    }
    return coords;
}

int pdb_to_mesh(string pdbpath_A, string pdbpath_B, string mshpath_A, string mshpath_B, double cutoff)
{    
    cout << "Calculating radius...\n";
    double max_radius = 0.0;
    for (const auto & file : experimental::filesystem::directory_iterator(pdbpath_A))
    {
        string filename = file.path();
        string extension = filename.substr(filename.size()-4,4);
        if (not extension.compare(".pdb"))
        {
            mesh meshA;
            meshA.coords = read_pdb(filename);
            double radius = meshA.calc_max_radius();
            if (radius > max_radius)
                max_radius = radius;
        }
    }
    for (const auto & file : experimental::filesystem::directory_iterator(pdbpath_B))
    {
        string filename = file.path();
        string extension = filename.substr(filename.size()-4,4);
        if (not extension.compare(".pdb"))
        { 
            mesh meshA;
            meshA.coords = read_pdb(filename);
            double radius = meshA.calc_max_radius();
            if (radius > max_radius)
                max_radius = radius;
        }
    }
    cout << "Max radius = " << max_radius << endl;
    
    for (const auto & file : experimental::filesystem::directory_iterator(pdbpath_A))
    {
        string filename = file.path();
        string extension = filename.substr(filename.size()-4,4);
        if (not extension.compare(".pdb"))
        {
            string filename_short = filename.substr(filename.rfind("/") + 1);
            filename_short = filename_short.substr(0,filename_short.size()-4);
            cout << "Constructing mesh for " << filename_short << '\r';
            cout << flush;
            mesh meshA;
            meshA.coords = read_pdb(filename);
            meshA.neighbor_grid_search(cutoff);
            meshA.neighbor_to_edge();
            meshA.neighbor_to_face();
            meshA.normalize_to(max_radius);
            string mshfile = mshpath_A + '/' + filename_short + ".msh";
            meshA.write_mesh(mshfile);
        }
    }
    for (const auto & file : experimental::filesystem::directory_iterator(pdbpath_B))
    {
        string filename = file.path();
        string extension = filename.substr(filename.size()-4,4);
        if (not extension.compare(".pdb"))
        {
            string filename_short = filename.substr(filename.rfind("/") + 1);
            filename_short = filename_short.substr(0,filename_short.size()-4);
            cout << "Constructing mesh for " << filename_short << '\r';
            cout << flush;
            mesh meshA;
            meshA.coords = read_pdb(filename);
            meshA.neighbor_grid_search(cutoff);
            meshA.neighbor_to_edge();
            meshA.neighbor_to_face();
            meshA.normalize_to(max_radius);
            string mshfile = mshpath_B + '/' + filename_short + ".msh";
            meshA.write_mesh(mshfile);
        }
    }
    return 0;
}

