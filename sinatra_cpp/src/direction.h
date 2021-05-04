#include <cmath>

using namespace std;

vector<vector<double> > rodrigues(vector<double> z, double r, int n)
{
    z = normalize3v(z);
    vector<double> z0(3,0.0);
    bool trivial = false;
    for(int i=0;i<3;i++)
    {
       if (fabs(z[i]) < 1e-15)
       {
            trivial = true;
            z0[i] = 1.0;
       }
    }
    if (not trivial)
    {
        z0[0] = 2.0 / z[0];
        z0[1] = -1.0 / z[1];
        z0[2] = -1.0 / z[2];
    }
    z0 = normalize3v(z0);
    for(int i=0;i<3;i++)
    {
        z0[i] = z0[i] * r + z[i];
    }
    z0 = normalize3v(z0);
    vector<double> B = cross(z,z0);
    vector<double> C(3,0.0);
    double zdotz0 = dot(z,z0);
    for(int i=0;i<3;i++)
        C[i] = z[i]*zdotz0;
    vector<vector<double> > directions;
    for(int i=0;i<n;i++)
    {
        double x = 2.0 * pi * (i+1) / n;
        double cosx = cos(x);
        double sinx = sin(x);
        vector<double> direction(3,0.0);
        for(int j=0;j<3;j++)
        {
            direction[j] = z0[j] * cosx + B[j] * sinx + C[j] * (1-cosx);
        }
        direction = normalize3v(direction);
        directions.push_back(direction);
    }  
    return directions;
}


vector<vector<double> > generate_equidistributed_points(int target, int N, bool hemisphere)
{
    vector<vector<double> > points;
    double a, d, d_theta, d_phi, theta, phi;
    int M_theta, M_phi;
    if (hemisphere)
    {
        a = 2.0*pi/N;
        d = sqrt(a);
        M_theta = (int)round(pi / d * 0.5);
        d_theta = pi / M_theta * 0.5;
    }
    else
    {
        a = 4.0*pi/N;
        d = sqrt(a);
        M_theta = (int)round(pi / d);
        d_theta = pi / M_theta;
    }
    d_phi = a / d_theta;
    for(int i=0;i<M_theta;i++)
    {
        if (hemisphere)
            theta = pi * 0.5 * i / M_theta;
        else
            theta = pi * (i+0.5) / M_theta;
        M_phi = (int)round(2.0 * pi * sin(theta) / d_phi);
        for(int j=0;j<M_phi;j++)
        {
            phi = 2.0 * pi * j / M_phi;
            vector<double> point(3,0.0);
            point[0] = sin(theta) * cos(phi);
            point[1] = sin(theta) * sin(phi);
            point[2] = cos(theta);
            point = normalize3v(point);
            points.push_back(point);
        }
    }
    if (points.size() < target)
        return generate_equidistributed_points(target,N+1,hemisphere);
    else
        return points;
}

vector<vector<double> > generate_equidistributed_cones(int n_cone, double cap_radius, int n_direction_per_cone, bool hemisphere)
{
    vector<vector<double> > sphere = generate_equidistributed_points(n_cone, n_cone, hemisphere);
    vector<vector<double> > directions;
    for(int i=0;i<n_cone;i++)
    {
        directions.push_back(sphere[i]);
        if (n_direction_per_cone > 1)
        {
            vector<vector<double> > direction_in_cone = rodrigues(sphere[i],cap_radius,n_direction_per_cone-1);
            for(int j=0;j<n_direction_per_cone-1;j++)
                directions.push_back(direction_in_cone[j]);
        }
    }
    return directions;
}



