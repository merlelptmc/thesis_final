#include"library_droso.h"

int main(){
        
    int nspins = 2;
    int nsites = 100;
    double temperature = 1;
    int n_iterations = 1000; //number of iterations in the MFSA
    //parameters of autocrine-paracrine if needed
//     double dauto = pr_in[0];
//     double dneig = pr_in[1];
    int ngrad=1;
    double basis_network[nspins*(nspins+ngrad)] = {5, -2, 0,5 , -15, -10};
    
    int n_directions = 2; //number of dimensions we explore
    double pr_exploration[n_directions*4] = {1, 2, -10, -10, 1, 1, 10, 10}; // how we explore: have to change the code in C++ to make this matrix more readable, the first numbers are the index of the parameters we variate
            
    int n_points_tot =1;
    int n_points[n_directions];
    for(int i=0; i<n_directions; i++){
        n_points[i] = ((pr_exploration[3*n_directions+i] - pr_exploration[1*n_directions+i]) / pr_exploration[2*n_directions+i]) +1;
        n_points_tot *= n_points[i];
    }

    //Initialization of the system
    Parameters system(nspins, nsites, ngrad, temperature);
    system.gradient.Construct_simple_gradient(nsites);
    system.network.Fill(basis_network);
    system.neighbors.construction_1D(nsites);    
    Spins spin(nspins, nsites );
    
    //Creation of output "matrices"
    double *pr_network = (double*)calloc((nspins+ngrad)*nspins*n_points_tot,sizeof(double));
    double *pr_quality = (double*)calloc((nspins+ngrad)*nspins*n_points_tot,sizeof(double));
    
    //Display informations on the landscape exploration
    cout << "Begining to explore" << endl << "Number of points to calculate = " << n_points_tot << endl << "Number of directions = " << n_directions << endl;
    for(int i=0; i< n_directions; i++){
        cout << "Direction " << pr_exploration[0*n_directions+i] << " de " << pr_exploration[1*n_directions +i] << " à " << pr_exploration[3*n_directions+i] << " par pas de " << pr_exploration[2*n_directions+i]<< endl;
    }
    cout << endl;
    
    //Actual exploration
    for(int i=0; i<n_points_tot; i++){
            
        //Update of the network
        for(int j=0; j<n_directions;j++){
            int modulo = n_points[j];
            int divid=1;
            for(int k=j+1; k<n_directions; k++){
                divid *= n_points[k];
            }
            
            system.network.J[ (int)pr_exploration[j]-1 ] = pr_exploration[1*n_directions+j] + i/divid%(modulo)* pr_exploration[2*n_directions+j];
        }
        
        if(i%(100) ==0 ){ //regular check out
            cout << i << "ème point " << endl;
        }
        
        //Calculation of the equilibrium state by MFSA and the associated quality for this network
        MFSAexp_asym(spin, system, n_iterations); 
        double Q = spin.Quality_max();
        
        //Saving the data
        for(int m=0; m < nspins+ngrad; m++){
            for(int n=0; n<nspins; n++){
                pr_network[i*nspins*(nspins+ngrad) + n*(nspins+ngrad) + m] = system.network.J[ n*(nspins+ngrad)+m];
            }
        }
        pr_quality[i] = Q;
        
    }
    
}


