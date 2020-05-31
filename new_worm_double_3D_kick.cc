#include <cmath>
#include <iostream>
#include<vector>
#include<fstream>
#include<string>
#include <random>

#include "WormGrid_2_3D.h"
#include "statistics.h"

using std::cout;

bool separation(const std::vector<int>& u, const std::vector<int>& v)
{
    if (int(abs(u[0]-v[0])+abs(u[1]-v[1])+abs(u[2]-v[2])) == 0)
    {
        return true;
    }
    else {return false;}
}

int overlaps(const std::vector<std::vector <int> >& ends)
{
    int temp = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = i+1; j<4; j++) {
            if (separation(ends[i],ends[j])) {
                //test = true;
                temp++;
            }
        }
    }
    
    return temp;
    
}

template <class X> // ensuring here that we don't double count when u and v haven't moved
int uv_separation(WormGrid_2_3D<X>& self, const std::vector<std::vector <int> >& ends)
{
    int temp = overlaps(ends);
    
    self.any_2_count()+=temp;

    if(temp==2)
      self.uvwz_count()++;
    if(temp==6)
      self.uvwz_count()+=temp;
    
    return temp;
    
}

template <class X>
void worm_step(WormGrid_2_3D<X>& self, std::vector<std::vector <int> >& ends)  //0,1,2,3,4,5 == left,right,down,up,front - z--,back z++
{
    
    //cout << u[0] << " " << u[1] << std::endl;
    /* TYPE 1 */
    int end_point = rand()%4;    //which end point to move
    int direction = rand()%6;    //which direction to move chosen endpoint
    
    std::vector<int> old = ends[end_point];    //identifying which end to move
    std::vector<int>& trial = ends[end_point];

    int n = self.size();
    int kl = 0;
                                                        //new indices: rightwards link k = 0
                                                        //             upwards link k = 1
                                                        //             lateral link k = 2
    if (direction == 0) {                               //left
        trial[0] = (trial[0]-1+n)%(n);
        self(trial[0], trial[1], trial[2], 0) = 1 - self(trial[0], trial[1], trial[2], 0);
        kl = self(trial[0], trial[1], trial[2], 0);

    }
    else if (direction == 1) {                          //right
        trial[0] = (trial[0]+1+n)%(n);
        self(old[0], trial[1], trial[2], 0) = 1 - self(old[0], trial[1], trial[2], 0);
        kl = self(old[0], trial[1], trial[2], 0);

    }
    else if (direction == 2) {                          //down
        trial[1] = (trial[1]-1+n)%(n);
        self(trial[0], trial[1], trial[2], 1) = 1 - self(trial[0], trial[1], trial[2], 1);
        kl = self(trial[0], trial[1], trial[2], 1);

    }
    else if (direction == 3) {                          //up
        trial[1] = (trial[1]+1+n)%(n);
        self(trial[0], old[1], trial[2], 1) = 1 - self(trial[0], old[1], trial[2], 1);
        kl = self(trial[0], old[1], trial[2], 1);

    }
    else if (direction == 4) {                          //front
        trial[2] = (trial[2]-1+n)%(n);
        self(trial[0], trial[1], trial[2], 2) = 1 - self(trial[0], trial[1], trial[2], 2);
        kl = self(trial[0], trial[1], trial[2], 2);
    }
    else {                                              //back
        trial[2] = (trial[2]+1+n)%(n);
        self(trial[0], trial[1], old[2], 2) = 1 - self(trial[0], trial[1], old[2], 2);
        kl = self(trial[0], trial[1], old[2], 2);
    }
    
    double r_num2 = double(rand())/RAND_MAX;
    double prob = pow(tanh(self.beta()),(2*kl-1));   //new kl so changed sign - def right

    if (prob < r_num2)                                                        //don't accept
    {
        
        //cout << "no change" << std::endl;
        if (direction == 0) {                               //left
            self(trial[0], trial[1], trial[2], 0) = 1 - self(trial[0], trial[1], trial[2], 0);
        }
        else if (direction == 1) {                          //right
            self(old[0], trial[1], trial[2], 0) = 1 - self(old[0], trial[1], trial[2], 0);
        }
        else if (direction == 2) {                          //down
            self(trial[0], trial[1], trial[2], 1) = 1 - self(trial[0], trial[1], trial[2], 1);
        }
        else if (direction == 3) {                          //up
            self(trial[0], old[1], trial[2], 1) = 1 - self(trial[0], old[1], trial[2], 1);
        }
        else if (direction == 4) {
            self(trial[0], trial[1], trial[2], 2) = 1 - self(trial[0], trial[1], trial[2], 2);
        }
        else {
            self(trial[0], trial[1], old[2], 2) = 1 - self(trial[0], trial[1], old[2], 2);
        }
        
        ends[end_point] = old;
        
    }

    int combs = uv_separation(self,ends);
    
    /* TYPE 2 - Kick substep - A lot of if and else statements, incorporate 4C2, 4C3, 4C4 into probabilities*/
    
    if (combs > 1) {
        double p_kick = 0.5*combs;
        double r_num1 = rand()/RAND_MAX;
        
        if (r_num1 < p_kick)
        {
            std::vector<std::vector <int> > pairs;
            for (int i = 0; i < 3; i++) {
                for (int j = i+1; j<4; j++) {
                    if (separation(ends[i],ends[j])) {
                        std::vector<int> temp(2); temp[0] = i; temp[1] = j;
                        pairs.push_back(temp);
                    }
                }
            }
            
            int rand_1 = rand()%combs;   //pick a random pair to kick
                        
            //cout << "KICK change" << std::endl;
            int r_new_x = rand()%n;
            int r_new_y = rand()%n;
            int r_new_z = rand()%n;
            
            ends[pairs[rand_1][0]][0] = r_new_x; ends[pairs[rand_1][1]][0] = r_new_x;
            ends[pairs[rand_1][0]][1] = r_new_y; ends[pairs[rand_1][1]][1] = r_new_y;
            ends[pairs[rand_1][0]][2] = r_new_z; ends[pairs[rand_1][1]][2] = r_new_z;
        }
    }
}

template <class X>
void worm_iteration(WormGrid_2_3D<X>& self, std::vector<std::vector<int> >& ends)
{
    for (int i = 0; i < self.size()*self.size()*self.size(); i++)    // grid size squared to make this function implement a full "iteration"
    {
        worm_step(self,ends);
    }

}

template <class X>
void run_worm(WormGrid_2_3D<X>& self, int warmup, int n, std::vector<double>& uvwz_list, std::vector<double>& any_2_list)
{
    int new_x1 = rand()%self.size();
    int new_y1 = rand()%self.size();
    int new_z1 = rand()%self.size();
    
    int new_x2 = rand()%self.size();
    int new_y2 = rand()%self.size();
    int new_z2 = rand()%self.size();
    
    std::vector<int> u(3); u[0] = new_x1; u[1] = new_y1; u[2] = new_z1;
    std::vector<int> v(3); v[0] = new_x1; v[1] = new_y1; v[2] = new_z1;
    
    std::vector<int> w(3); w[0] = new_x2; w[1] = new_y2; w[2] = new_z2;
    std::vector<int> z(3); z[0] = new_x2; z[1] = new_y2; w[2] = new_z2;
    
    std::vector<std::vector<int> > ends; ends.push_back(u); ends.push_back(v); ends.push_back(w); ends.push_back(z);
    
    //cout << "Warming Up" << std::endl;
    
    for (int i = 0; i < warmup; i++)
    {
        worm_iteration(self,ends);
    }
    
    //cout << "Measuring" << std::endl;
    
    //std::ofstream myfile;
    //myfile.open("binder_cumulant_worm2_no_kick.dat");
    
    for (int i = 0; i < n; i++)
    {
        self.uvwz_count() = 0;
        self.any_2_count() = 0;
        
        worm_iteration(self,ends);
        
        uvwz_list[i] = double(self.uvwz_count());
        any_2_list[i] = double(self.any_2_count());
        
        //myfile << i << " " << self.size()*self.size()*average(uvwz_list, 0, i) << " " << average(any_2_list,0,i)  << std::endl;
    }
    
    //myfile.close();

}


void susc_beta_L()
{
    int warmup = 200000;
    int n = 2000000;
    
    std::vector<double> betalist(11);
    for (int i=0; i<11; i++)    { betalist[i] = 0.31+i*0.02; } //0.3-0.5
    int n_beta = betalist.size();
    
    std::vector<int> sizelist(5);
    for (int i=0; i<5; i++)     { sizelist[i] = 8*pow(2,i); }
    
    for (int j = 0;  j < 5; j++)
    {
        std::string s = std::to_string(j);
        std::ofstream myfile;
        myfile.open("susc_beta_zoom_L"+s+".dat");
        
        
        for (int i = 0; i < n_beta; i++)
        {
            double beta = betalist[i];
            std::vector<double> uvwz_list(n);
            std::vector<double> any_2_list(n);
                    
            cout << j << " " << i << std::endl;
            
            WormGrid_2_3D<int> test(sizelist[j],beta);
            run_worm(test, warmup, n, uvwz_list, any_2_list);
            
            double av = average(uvwz_list);
            myfile << beta << " " << av << " " << error(uvwz_list, av) << std::endl;
            
        }
   
        myfile.close();
        
    }
}

void full_sim_all_sizes(double beta)
{
    int warmup = 100000;
    int n = 1000000;

    std::vector<int> sizelist(6);
    for (int i=0; i<6; i++)     { sizelist[i] = 4*pow(2,i); }

    for (int j = 0;  j < 1; j++)
    {
     std::string s = std::to_string(4);
     std::ofstream myfile;
     myfile.open("full_sim_criticalbeta_L"+s+".dat");
     
     std::vector<double> uvwz_list(n);
     std::vector<double> any_2_list(n);
             
     cout << j << std::endl;
     
     WormGrid_2_3D<int> test(sizelist[j],beta);
     run_worm(test, warmup, n, uvwz_list, any_2_list);
     
     for (int i = 0; i < n; i++) {
         myfile << uvwz_list[i] << std::endl;
     }

     myfile.close();
     
    }
}

void binder_cumulants_L()
{
    int warmup = 100000;
    int n = 1000000;
    
    std::vector<double> betalist(5);
    for (int i=0; i<5; i++)    { betalist[i] = 0.17+i*0.02; }
    int n_beta = betalist.size();
    
    std::vector<int> sizelist(3);
    for (int i=0; i<3; i++)     { sizelist[i] = 8*pow(2,i); }
    
    for (int j = 2;  j < 3; j++)
    {
        std::string s = std::to_string(j);
        std::ofstream myfile;
        myfile.open("3D_binder_cumulant_worm2_kick_extra_L"+s+".dat");
        
        std::ofstream myfile2;
        myfile2.open("3D_susc_beta_worm2_kick_extra_L"+s+".dat");
        
        for (int i = 0; i < 1; i++)
        {
            //double beta = betalist[i];
            double beta = 0.26;
            int size = sizelist[j];
                    
            cout << j << " " << i << std::endl;
            
            WormGrid_2_3D<int> test(size,beta);
            
            std::vector<double> uvwz_list(n);
            std::vector<double> any_2_list(n);
            
            run_worm(test, warmup, n, uvwz_list, any_2_list);
            
            double uv = average(any_2_list);
            double uvwz = average(uvwz_list);
            
            double uv_err = error(any_2_list,uv);
            double uvwz_err = error(uvwz_list,uvwz);
            
            double divsq = uvwz/uv/uv;
            double b_error = divsq*sqrt(pow(uvwz_err/uvwz,2)+pow(2.0*uv_err/uv,2))*4.0;
            
            double div = uv/uvwz;
            double x_error = div*sqrt(pow(uvwz_err/uvwz,2)+pow(uv_err/uv,2));
            
            myfile << beta << " " << 1.0 - 4.0*divsq*size*size*size << " " << b_error*size*size*size << std::endl;
            
            myfile2 << beta << " " << div/size/size/size << " " << x_error/size/size/size << std::endl;
            
        }
   
        myfile.close();
        myfile2.close();
        
    }
}

int main()
{
    srand(time(NULL));
    
    binder_cumulants_L();
    
    return 0;
    
}
