#include <cmath>
#include <vector>

using namespace std;

const double pi = M_PI;

double norm3v(vector<double> v)
{
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

vector<double> normalize3v(vector<double> v)
{
    double norm = norm3v(v);
    vector<double> w(3,0.0);
    for(int i=0;i<3;i++)
        w[i] = v[i] / norm;
    return w;
}

vector<double> cross(vector<double> u, vector<double> v)
{
    vector<double> result(3,0.0);
    result[0] = u[1] * v[2] - u[2] * v[1];
    result[1] = u[2] * v[0] - u[0] * v[2];
    result[2] = u[0] * v[1] - u[1] * v[0];
    return result;
}

double dot(vector<double> u, vector<double> v)
{
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

double average(vector<double> v)
{
    double result;
    for(int i=0;i<v.size();i++)
        result += v[i];
    return result/v.size();
}

vector<double> linspace(double min, double max, int N)
{
    vector<double> v(N,0.0);
    double step = (max-min)/(N-1);
    v[0] = min;
    for(int i=1;i<N;i++)
        v[i] = v[i-1]+step;
    return v;
}

double vmin(vector<double> v)
{
    double result=9e99;
    for(int i=0;i<v.size();i++)
        if (v[i] < result)
            result = v[i];
    return result;
}

double vmax(vector<double> v)
{
    double result=-9e99;
    for(int i=0;i<v.size();i++)
        if (v[i] > result)
            result = v[i];
    return result;
}
