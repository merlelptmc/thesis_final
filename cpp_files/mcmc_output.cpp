#include"library_droso.h"


int main(){
        
        int nsites=40;
        int ngenes=2;
        int ngrad=1;
        int n_eq=1000000;
        int n_rec = 100;
        int step_rec = 50;
        
        double temperature=1;
        double network[(ngenes+ngrad)*ngenes] = {5, -2.1, 0, 6, -15, -12};
        
        
        //Initialization of the system
        Parameters system(ngenes, nsites, ngrad, temperature);
        system.gradient.Construct_simple_gradient(nsites);
        system.network.Fill(network);
        system.neighbors.construction_1D(nsites);    
        Spins spin(ngenes, nsites);
        spin.Fill_rand();
        
        system.Print_conditions();
        cout << endl << "N equilibration = " << n_eq << endl;
        cout << "Recording " << n_rec << " every " << step_rec << " calculations points after equilibration." << endl;
        
        cout << endl << "Initial Spin" << endl;
        spin.Print();
        
        MCMC(spin, system, n_eq, n_rec, step_rec);
        
        double quality = spin.Quality_max();
        
        cout << endl << "Final Spin" << endl;
        spin.Print();        
        cout << "Quality of the pattern (0<Q<1) = " << quality << endl << endl;
               
        
}
