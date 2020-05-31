#include <cmath>
#include <iostream>
#include<vector>
#include<fstream>
#include<string>
#include <random>

using std::cout;

#include "MetGrid.h"
#include "statistics.h"

template <class X>
void print(const MetGrid<X>& self)
{
    int n = self.size();
    for (int j=n-1; j>=0; j--)
    {
        for (int i=0; i<n; i++)
        {
            cout << self(i,j) << "  ";
        }
        cout << std::endl;
    }
}

template <class X>
std::vector<double> neighbours(const MetGrid<X>& self, int i, int j)  //BC ON!!!
{
    std::vector<double> temp(4);

    temp[0] = self(i,j+1);   //up
    temp[1] = self(i,j-1);   //down
    temp[2] = self(i-1,j);   //left
    temp[3] = self(i+1,j);   //right
    return temp;
}

template <class X>
double energy_change(const MetGrid<X>& self, int i, int j)
{
    double temp = 0;
    for (int k=0; k<neighbours(self,i,j).size(); k++)
    {
        temp += neighbours(self,i,j)[k];
    }
    return 2.0*temp*self(i,j);
}

template <class X>
double total_energy(const MetGrid<X>& self)
{
    int n = self.size();
    double tot_e = 0;
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
            tot_e -= energy_change(self,i,j)/2.0;
        }
    }
    return tot_e/2.0;
}

template <class X>
double magnetisation(const MetGrid<X>& self)
{
    int n = self.size();
    double tot_mag = 0;
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
            tot_mag += self(i,j);
        }
    }
    return tot_mag;
}

template <class X>
void metropolis(MetGrid<X>& self, char t)    //take measurements
{
    double temp_e = total_energy(self);
    double temp_m = magnetisation(self);
    int n = self.size();
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
            double e = energy_change(self,i,j);
            if (e <= 0)
            {
                self(i,j)*=(-1);
                temp_e += e;
                temp_m += 2*self(i,j);
            }
            else
            {
                double r = ((double) rand() / (RAND_MAX));
                if (exp(-1.0*e*self.beta()) > r)
                {
                    self(i,j)*=(-1);
                    temp_e += e;
                    temp_m += 2*self(i,j);
                }
            }
        }
    }
    
    self.e_list.push_back(temp_e);
    self.e_sq_list.push_back(pow(temp_e,2));
    self.mag_list.push_back(temp_m);
    self.mag_sq_list.push_back(pow(temp_m,2));
}

template <class X>
void metropolis(MetGrid<X>& self)   //for warmup
{
    int n = self.size();
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
            double e = energy_change(self,i,j);
            if (e <= 0)
            {
                self(i,j)*=(-1);
            }
            else
            {
                double r = ((double) rand() / (RAND_MAX));
                if (exp(-1.0*e*self.beta()) > r)
                {
                    self(i,j)*=(-1);
                }
            }
        }
    }
}

void four_outputs(int size, int warmup, int n)
{
    std::ofstream myenergyfile;
    myenergyfile.open("energy.dat");
    std::ofstream mymagfile;
    mymagfile.open("magnet.dat");
    std::ofstream mysuscfile;
    mysuscfile.open("susc.dat");
    std::ofstream mycapfile;
    mycapfile.open("cap.dat");
    
    std::vector<double> betalist(40);
    for (int i=0; i<9; i++) { betalist[i] = 0.1+i*0.025; }  //0.1-0.3
    for (int i=0; i<25; i++) { betalist[i+9] = 0.31+i*0.01; } //0.31-0.55
    for (int i=0; i<6; i++) { betalist[i+34] = 0.575+i*0.025; } //0.575-0.7
    int n_beta = betalist.size();
    
    for (int z = 0; z < n_beta; z++)
    {
        double beta = betalist[z];
        cout << z << std::endl;
        
        MetGrid<int> test(size, beta);
        
        cout << "Warming Up" << std::endl;
        
        for (int i=0; i<warmup; i++)
        {
            metropolis(test);
        }
        
        cout << "Taking Measurements" << std::endl;
        
        for (int i=0; i<n; i++)
        {
            metropolis(test, 't');
        }
        
        double av_energy = average(test.e_list, 0, n)/pow(size,2);
        double av_energy_sq = average(test.e_sq_list, 0, n)/pow(size,4);
        double av_mag = average(test.mag_list, 0, n)/pow(size,2);
        double av_mag_sq = average(test.mag_sq_list, 0, n)/pow(size,4);
        
        myenergyfile << beta << " " << av_energy << std::endl;
        mymagfile << beta << " " << abs(av_mag) << std::endl;
        mysuscfile << beta << " " << (av_mag_sq - av_mag*av_mag)*test.beta() << std::endl;
        mycapfile << beta << " " <<(av_energy_sq - av_energy*av_energy)*pow(test.beta(),2) << std::endl;
        
    }
    
    myenergyfile.close();
    mymagfile.close();
    mysuscfile.close();
    mycapfile.close();
    
}

template <class X>
void full_simulation(MetGrid<X>& self, int n)
{
    std::string s = std::to_string(self.size());
    std::ofstream myenergyfile;
    myenergyfile.open("full_criticalbeta_energy_L"+s+".dat");
    std::ofstream mymagfile;
    mymagfile.open("full_criticalbeta_magnet_L"+s+".dat");
    
    cout << "Warming Up" << std::endl;
    
    for (int i=0; i<100000; i++)
    {
        metropolis(self);
    }
    
    cout << "Taking Measurements" << std::endl;
    
    for (int i=0; i<n; i++)
    {
        metropolis(self, 't');
        myenergyfile << self.e_list[i] << std::endl;
        mymagfile << abs(self.mag_list[i]) << std::endl;
    }
    
    myenergyfile.close();
    mymagfile.close();
    
}

void binder_cumulants_L()
{
    int warmup = 100000;
    int n = 1000000;
    
    std::vector<double> betalist(21);
    for (int i=0; i<21; i++)    { betalist[i] = 0.30+i*0.01; } //0.3-0.5
    int n_beta = betalist.size();
    
    std::vector<int> sizelist(5);
    for (int i=0; i<5; i++)     { sizelist[i] = 8*pow(2,i); }
    
    for (int j = 4;  j < 5; j++)
    {
        std::string s = std::to_string(j);
        std::ofstream myfile;
        myfile.open("binder_cumulants_L"+s+".dat");
        
        std::ofstream myenergyfile;
        myenergyfile.open("energy_L"+s+".dat");
        std::ofstream mymagfile;
        mymagfile.open("magnet_redo_L"+s+".dat");
        std::ofstream mymagsqfile;
        mymagsqfile.open("magnet_sq_redo_L"+s+".dat");
        std::ofstream mycapfile;
        mycapfile.open("cap_L"+s+".dat");
        
        
        for (int i = 14; i < n_beta; i++)
        {
            double beta = betalist[i];
            double size = sizelist[j];
                    
            cout << j << " " << i << std::endl;
            
            MetGrid<int> test(size,beta);
            
            for (int i=0; i<warmup; i++)
            {
                metropolis(test);
            }
            cout << "Measuring" << std::endl;
            for (int i=0; i<n; i++)
            {
                metropolis(test, 't');
            }
            
            double sigma_sq = average(test.mag_sq_list);
            
            std::vector<double> sigma_four_list = square(test.mag_sq_list);
            double sigma_four = average(sigma_four_list);
            
            myfile << beta << " " << 1.0 - sigma_four/(3.0*sigma_sq*sigma_sq)  << std::endl;
            
            /* Killing two birds with one stone - plotting observs as a fn of size */
            
            double av_energy = average(test.e_list)/pow(size,2);
            double av_energy_sq = average(test.e_sq_list)/pow(size,4);
            test.mag_list = abs(test.mag_list);
            double av_mag = average(test.mag_list)/pow(size,2);
            double av_mag_sq = average(test.mag_sq_list)/pow(size,2);
            
            myenergyfile << beta << " " << av_energy << std::endl;
            mymagfile << beta << " " << av_mag << " " << error(test.mag_list, av_mag) << std::endl;
            mymagsqfile << beta << " " << av_mag_sq << " " << error(test.mag_sq_list, av_mag_sq) << std::endl;
            mysuscfile << beta << " " << (av_mag_sq - av_mag*av_mag)*test.beta() << std::endl;
            mycapfile << beta << " " <<(av_energy_sq - av_energy*av_energy)*pow(test.beta(),2) << std::endl;
            
        }
   
        myfile.close();
        
        myenergyfile.close();
        mymagfile.close();
        mymagsqfile.close();
        mycapfile.close();

        
    }
}

int main()
{
    srand(time(NULL));
    
    binder_cumulants_L();
    
    //four_outputs(32, 100000, 200000);
    
    return 0;
    
}
