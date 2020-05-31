#include <cmath>
#include <iostream>
#include<vector>
#include<fstream>
#include<string>
#include <random>

#include "WormGrid.h"
#include "statistics.h"

using std::cout;

template <class X>    //to print out the movement of the worm in a nice format
void print(const WormGrid<X>& self)
{
    for (int j=self.size()-1; j>=0; j--)
    {
        for (int i=0; i<self.size(); i++)
        {
            cout << self(i,j,1) << "   ";   //5 spaces boundary included
        }
        cout << std::endl;
        for (int i=0; i<self.size(); i++)
        {
            cout << ". " << self(i,j,0) << " ";
        }
        cout << std::endl;
    }
}

template <class X> // ensuring here that we don't double count when u and v haven't moved
void uv_separation(WormGrid<X>& self, const std::vector<int>& u, const std::vector<int>& v, int z)
{
    if (int(abs(u[0]-v[0])+abs(u[1]-v[1])) == z)
    {
        self.uv_count()++;
    }
    
}

template <class X>
void worm_step(WormGrid<X>& self, std::vector<int>& u, std::vector<int>& v)  //0,1,2,3 == left,right,down,up
{
    /* TYPE 1 */
    int direction = rand()%4;
    //cout << "direction = " << direction << std::endl;
    int n = self.size();
    int kl = 0;
    std::vector<int> old_u = u;                         //new indices: rightwards link z = 0
                                                        //             upwards link z = 1
    if (direction == 0) {                               //left
        u[0] = (u[0]-1+n)%(n);
        self(u[0], u[1], 0) = 1 - self(u[0], u[1], 0);
        kl = self(u[0], u[1], 0);

    }
    else if (direction == 1) {                          //right
        u[0] = (u[0]+1+n)%(n);
        self(old_u[0], u[1], 0) = 1 - self(old_u[0], u[1], 0);
        kl = self(old_u[0], u[1], 0);

    }
    else if (direction == 2) {                          //down
        u[1] = (u[1]-1+n)%(n);
        self(u[0], u[1], 1) = 1 - self(u[0], u[1], 1);
        kl = self(u[0], u[1], 1);

    }
    else {                                              //up
        u[1] = (u[1]+1+n)%(n);
        self(u[0], old_u[1], 1) = 1 - self(u[0], old_u[1], 1);
        kl = self(u[0], old_u[1], 1);

    }
    
    double r_num2 = double(rand())/RAND_MAX;

    double prob = pow(tanh(self.beta()),(2*kl-1));   //new kl so changed sign - def right
    
    
    if (prob < r_num2)                                                        //don't accept
    {
        
        //cout << "no change" << std::endl;
        if (direction == 0) {                               //left
            self(u[0], u[1], 0) = 1 - self(u[0], u[1], 0);
        }
        else if (direction == 1) {                          //right
            self(old_u[0], u[1], 0) = 1 - self(old_u[0], u[1], 0);
        }
        else if (direction == 2) {                          //down
            self(u[0], u[1], 1) = 1 - self(u[0], u[1], 1);
        }
        else {                                              //up
            self(u[0], old_u[1], 1) = 1 - self(u[0], old_u[1], 1);
        }
        u = old_u;
        
    }
    
    //print(self);
        
    /* TYPE 2 */
    if (u[0] == v[0] && u[1]==v[1])
    {
        uv_separation(self,u,v,0);
        //cout << "-------------------KICK------------------" << v[0] << " " << v[1] << std::endl;
        double p_kick = 0.5;
        double r_num1 = rand()/RAND_MAX;
        if (r_num1 < p_kick)
        {
            //cout << "KICK change" << std::endl;
            int r_new_x = rand()%n;
            int r_new_y = rand()%n;
            
            u[0] = r_new_x; v[0] = r_new_x;
            u[1] = r_new_y; v[1] = r_new_y;
        }
        
        //print(self);
    }
    
}

template <class X>
void worm_iteration(WormGrid<X>& self, std::vector<int>& u, std::vector<int>& v)
{
    for (int i = 0; i < self.size()*self.size(); i++)    // grid size squared to make this function implement a full "iteration"
    {
        worm_step(self,u,v);
    }

}

template <class X>
void run_worm(WormGrid<X>& self, int warmup, int n, std::vector<double>& uv_list)
{
    int new_x = rand()%self.size();
    int new_y = rand()%self.size();
    
    std::vector<int> u(2); u[0] = new_x; u[1] = new_y;
    std::vector<int> v(2); v[0] = new_x; v[1] = new_y;
    
    //cout << "Warming Up" << std::endl;
    
    for (int i = 0; i < warmup; i++)
    {
        worm_iteration(self,u,v);
    }
    
    //cout << "Measuring" << std::endl;
    
    for (int i = 0; i < n; i++)
    {
        self.uv_count() = 0;
        worm_iteration(self,u,v);

        uv_list[i] = double(self.uv_count());
    }

}


void susc_beta_L()
{
    int warmup = 100000;
    int n = 1000000;
    
    std::vector<double> betalist(21);
    for (int i=0; i<21; i++)    { betalist[i] = 0.32+i*0.01; } //0.32-5.2
    int n_beta = betalist.size();
    
    std::vector<int> sizelist(5);
    for (int i=0; i<5; i++)     { sizelist[i] = 8*pow(2,i); }
    
    for (int j = 1;  j < 5; j++)
    {
        std::string s = std::to_string(j);
        std::ofstream myfile;
        myfile.open("susc_beta_single_L"+s+".dat");
        
        for (int i = 0; i < n_beta; i++)
        {
            double beta = betalist[i];
            std::vector<double> uv_list(n);
                    
            cout << j << " " << i << std::endl;
            
            WormGrid<int> test(sizelist[j],beta);
            run_worm(test, warmup, n, uv_list);
            
            double av = average(uv_list);
            double d_error = error(uv_list, av);
            
            myfile << beta << " " << 1.0/av << " " << d_error/av/av << std::endl;
            
            if (beta == 0.44)
            {
                std::ofstream myfile2;
                myfile2.open("full_single_criticalbeta_L"+s+".dat");
                
                for (int k = 0; k < n; k++) {
                    myfile2 << uv_list[k] << std::endl;
                }
                
                myfile2.close();
            }
            
        }
   
        myfile.close();
        
    }
}

void full_sim_all_sizes(double beta)
{
    int warmup = 100000;
    int n = 1000000;

    std::vector<int> sizelist(6);
    for (int i=0; i<6; i++)     { sizelist[i] = 8*pow(2,i); }

    for (int j = 5;  j < 6; j++)
    {
     std::string s = std::to_string(j);
     std::ofstream myfile;
     myfile.open("full_single_criticalbeta_L"+s+".dat");
     
     std::vector<double> uv_list(n);
             
     cout << j << std::endl;
     
     WormGrid<int> test(sizelist[j],beta);
     run_worm(test, warmup, n, uv_list);
     
     for (int i = 0; i < n; i++) {
         myfile << uv_list[i] << std::endl;
     }

     myfile.close();
     
    }
}


int main()
{
    srand(time(NULL));
    
    //binder_cumulants_L();
    
    full_sim_all_sizes(0.4407);
    
    //susc_beta_L();
    
    return 0;
    
}

