#include"library_droso.h"


int main(){
        
        int nspins = 1;
        int ngrad=1;
        int nsites = 100;
        double temperature = 1;
        //     double dauto = 1;
        //     double dneig = 1;
        int n_iterations = 1000;
        
        Parameters system(nspins, nsites, ngrad, temperature);
        double network[system.network.nparam] = {-3.14, -1.87};
        Print(network, 2,1);
        system.gradient.Construct_simple_gradient(nsites);
        system.network.Fill(network);
        system.neighbors.construction_1D(nsites);    
        Spins spin(nspins, nsites);
        
        double *spin_out;
        double *spin_record;
        
        system.Print_conditions();
        
        MFSAexp_asym_rec(spin, system, n_iterations, spin_record);
        double quality = spin.Quality_max();
        
        cout << endl << "Final Spin" << endl;
        spin.Print();        
        cout << "Quality of the pattern (0<Q<1) = " << quality << endl << endl;
        
}


