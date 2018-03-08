#include <string.h>
#include <cstdio>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <cstdlib>
#include <random>
#include <ctime>


using namespace std;
void Print(double *mat, int nlin, int ncol );
void Print(int *mat, int nlin, int ncol );
void copy(double *in, double *out, int nelem);

// generates a random double number d, uniformly distributed beween 0 and 1 (0 <= d < 1) 
inline double drand(){
  return rand()/(double)RAND_MAX;
}
// generates a random integer number i, uniformly distributed beween 0 and n (0 <= i < n)
inline int irand(int n){
  return rand()%n;
}

struct Gene_network{
        int nspins;
        int ngrad;
        int nparam;
        double *J;
        double *C;
        Gene_network(): nspins(0), ngrad(0), nparam(0), J(0), C(0){}
        Gene_network(int ng, int nm);
        Gene_network(int ng, int nm, double *M);
        ~Gene_network();
        void Init(int ng, int nm);
        void Init_rand(int ng, int nm, double abs_min, double abs_max);
        void Fill(double *M);
        void Print();
        int Check_bounds(double minim, double maxim);
};

struct Spatial_grid{
        int nsites;
        int length;
        int thickness;
        int nlinks;
        int *sites;
        int *index;
        int *number;
        Spatial_grid(): nsites(0), thickness(0), nlinks(0), number(0), index(0), sites(0){}
        ~Spatial_grid();   
        void Init(int L, int l);
        int construction_1D(int ns); // 1D : NN=i, i+1, i-1
        int construction_2D_rectangle_square(int L, int l); // NN= square centered in i,j
        int construction_2D_rectangle_cross(int L, int l); // NN=cross centered in i,j
        void Print(); // Display for each site who are the neighbors
        
};

struct Gradient{
        int ngrad;
        int nsites;
        double *values;    
        Gradient(): ngrad(0), nsites(0), values(0){}
        Gradient(int ngr, int nsi);    
        Gradient(int ngr, int nsi, double *data);    
        ~Gradient();
        void Init(int ngr, int nsi);
        void Fill(double *data);
        void Print();
        void Construct_simple_gradient(int nsites);
        void Construct_opp_gradient(int nsites);      
        
};

struct Parameters{
        int nspins;
        int nsites;
        int ngrad;
        double temperature;
        Gradient gradient;
        Gene_network network;
        Spatial_grid neighbors;    
        Parameters(): nspins(0), nsites(0), ngrad(0), temperature(0){}    
        Parameters(int nsp, int nsi, int ngr, double T);
        ~Parameters();
        void Init(int nsp, int nsi, int ngr, double T);
        void Fill_gradient( double *H);
        void Fill_param(double *par);
        void Print_conditions();
        
};

struct Spins{
        int nspins;
        int nsites;
        double *state;       
        Spins(): nspins(0), nsites(0), state(0){};    
        Spins(int nsp, int nsi);
        Spins(int nsp, int nsites, double *data);
        ~Spins();
        void Init(int nsp, int nsi);
        void Fill(double *data);
        void Fill_rand();    
        int Add(Spins &spin_2);
        void Print();
        void Switch_one(double *heff, double T, int nuc);
        void Switch_mean(double *heff, double T);
        void Calculate_heff_1nuc(double *heff, Parameters &system, int nuc);
        void Calculate_heff(double *heff, Parameters &system);    
        int Test_stability(Spins &spin_mean, double critere);
        double Quality_max();
};

void MFSAexp_asym(Spins &spin, Parameters &system, int n_iterations);

void MFSAexp_asym_rec(Spins & spin, Parameters &system, int n_iterations, double *spin_rec);

void MCMC(Spins &spin, Parameters &system, int n_equilibrium, int n_recording, int step_recording);
