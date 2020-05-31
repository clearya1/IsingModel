#include <cmath>
#include <iostream>
#include <vector>

double autocovariance(const std::vector<double>& x, int t)
{
    double r = 0;
    for (int i=0; i<x.size()-t; i++)
    {
        r += x[i]*x[i+t];
    }
    double av = 0;
    for (int i=0; i<x.size(); i++)
    {
        av += x[i];
    }
    return r/(x.size()-t)-pow(av/(x.size()),2);
}

double autocorrelation(const std::vector<double>& x, int t)
{
    return autocovariance(x,t)/autocovariance(x,0);
}

double autocorrelation_error(const std::vector<double>& x, int t)
{
    double temp = 0;
    for (int i=t; i<x.size()-t; i++)
    {
        temp += pow(autocorrelation(x, i+t)+autocorrelation(x, i-t)-autocorrelation(x, i)*autocorrelation(x, t),2);
    }
    return sqrt(temp)/(x.size()-2*t);   //only averaged over N-2t autocorrelations
}

double jackknife(const std::vector<double>& x)
{
    double av = 0;
    for (int i = 0; i<x.size(); i++)
    {
        double temp = 0;
        for (int j = 0; j<x.size(); j++)
        {
            if (i!=j)
            {
                temp += x[j];
            }
        }
        temp /= (double(x.size()-1));
        av += temp;
    }
    return av/(double(x.size()));
}

double average(const std::vector<double>& x, int start, int end)
{
    double temp = 0;
    for (int i = start; i < end; i++)
    {
        temp += x[i];
    }
    return temp/(end-start);
}

double average(const std::vector<double>& x)
{
    double temp = 0;
    for (int i = 0; i < x.size(); i++)
    {
        temp += x[i];
    }
    return temp/(x.size());
}

double error(const std::vector<double>& x, double av)
{
    double temp = 0;
    for (int i = 0; i < x.size(); i++)
    {
        temp += (x[i]-av)*(x[i]-av);
    }
    
    return sqrt((temp-av*av)/double(x.size()-1)/double(x.size()));
}

std::vector<double> abs(std::vector<double>& x)
{
    for (int i = 0; i<x.size(); i++)
    {
        x[i] = abs(x[i]);
    }
    return x;
}

std::vector<double> square(const std::vector<double>& x)
{
    
    std::vector<double> temp(x.size());
    for (int i = 0; i < x.size(); i++)
    {
        temp[i] = x[i]*x[i];
    }
    return temp;
}

std::vector<double> multiply(std::vector<double>& x, std::vector<double>& y)
{
    std::vector<double> temp(x.size());
    for (int i = 0; i < x.size(); i++) {
        temp[i] = x[i]*y[i];
    }
    return temp;
}
