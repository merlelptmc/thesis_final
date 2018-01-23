#include"library_droso.h"


int main(){
        
        int nsites=50;
        int ngenes=1;
        int ngrad=1;
        int n_eq=1000000;
        int n_rec = 100;
        int step_rec = 50;
        
        double temperature=1;
        double network[(ngenes+ngrad)*ngenes] = {5, -15};
        
        
        //Initialization of the system
        Parameters system(ngenes, nsites, ngrad, temperature);
        system.gradient.Construct_simple_gradient(nsites);
        system.network.Fill(network);
        system.network.Print();
        system.neighbors.construction_1D(nsites);    
        Spins spin(ngenes, nsites);
        spin.Fill_rand();
        
        spin.Print();
        
        MCMC(spin, system, n_eq, n_rec, step_rec);
        
        spin.Print();
        
               
        
}
