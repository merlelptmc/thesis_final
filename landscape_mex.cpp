#include"library_droso.h"
#include "mex.h"
#include "matrix.h"

// int main(){
//     
//  
//     int n=5;
//     int m=10;
//     cout << n << "   " << m << endl;
//     
//     Gene_network net=Gene_network(n,m);
//     net.Print();
//     
//     
//     
// }


void mexFunction( int nlhs, mxArray *plhs[],     
                  int nrhs, const mxArray *prhs[] ) { 
    
//     prhs[0] = nspins
//     prhs[1] = nsites
//     prhs[2] = temperature
//     prhs[3] = genetic network
//     prhs[4] = [diffusion_auto diffusion_neighbors]
//     prhs[5] = n_iterations dans mfsa 
//     prhs[6] = [variables à modifier début pas fin]
    
    double *pr_in = mxGetPr(prhs[0]);
    int nspins = pr_in[0];
    pr_in =mxGetPr(prhs[1]);
    int nsites = pr_in[0];
    pr_in =mxGetPr(prhs[2]);
    double temperature = pr_in[0];
    pr_in =mxGetPr(prhs[5]);
    int n_iterations = pr_in[0];
    pr_in =mxGetPr(prhs[4]);
    double dauto = pr_in[0];
    double dneig = pr_in[1];
    int ngrad=1;
    
    int n_directions = mxGetM(prhs[6]);
    if(mxGetN(prhs[6]) != 4){
        mexErrMsgTxt("Error in the matrix describing the exploration to do");
    }
    double *pr_exploration = mxGetPr(prhs[6]);
    int n_points_tot =1;
    int n_points[n_directions];
    for(int i=0; i<n_directions; i++){
        n_points[i] = ((pr_exploration[3*n_directions+i] - pr_exploration[1*n_directions+i]) / pr_exploration[2*n_directions+i]) +1;
        n_points_tot *= n_points[i];
    }
            
    Parameters system(nspins, nsites, 1, temperature);
    system.gradient.Construct_simple_gradient(nsites);
    system.network.Fill(mxGetPr(prhs[3]));
    system.neighbors.construction_1D(nsites);    
    Spins spin(nspins, nsites, dauto, dneig);
        
    nlhs = 2;
    int nb_dim = 3;
    mwSize dim[nb_dim];
    dim[0]= nspins+ngrad;
    dim[1] = nspins;
    dim[2] = n_points_tot;
    plhs[0] = mxCreateNumericArray(nb_dim, dim, mxDOUBLE_CLASS,mxREAL);
    double *pr_network = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(1,n_points_tot, mxREAL);
    double *pr_quality = mxGetPr(plhs[1]);
    
    cout << "Begining to explore" << endl << "Number of points to calculate = " << n_points_tot << endl << "Number of directions = " << n_directions << endl;
    for(int i=0; i< n_directions; i++){
        cout << "Direction " << pr_exploration[0*n_directions+i] << " de " << pr_exploration[1*n_directions +i] << " à " << pr_exploration[3*n_directions+i] << " par pas de " << pr_exploration[2*n_directions+i]<< endl;
    }
    cout << endl;
    
    for(int i=0; i<n_points_tot; i++){
        for(int j=0; j<n_directions;j++){
            int modulo = n_points[j];
            int divid=1;
            for(int k=j+1; k<n_directions; k++){
                divid *= n_points[k];
            }
            
            system.network.J[ (int)pr_exploration[j]-1 ] = pr_exploration[1*n_directions+j] + i/divid%(modulo)* pr_exploration[2*n_directions+j];
        }
        if(i%(10000) ==0 ){
            cout << i << "ème point " << endl;
//             system.network.Print;
        }
        MFSAexp_asym(spin, system, n_iterations);
        double Q = spin.Quality_max();
        
        for(int m=0; m < nspins+ngrad; m++){
            for(int n=0; n<nspins; n++){
                pr_network[i*nspins*(nspins+ngrad) + n*(nspins+ngrad) + m] = system.network.J[ n*(nspins+ngrad)+m];
            }
        }
        pr_quality[i] = Q;
        
    }
    
//     system.network.J[0]= 5;
//     system.network.J[1] = -15;
//     spin.Fill_rand();
// //     spin.Print();
//     MFSAexp_asym(spin, system, n_iterations);    
// //     system.network.Print();
//     spin.Print();
//     double Q = spin.Quality_max();
//     cout << Q << endl;
    
    
}


