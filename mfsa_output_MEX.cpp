#include "mex.h"
#include "matrix.h"
#include "library_droso.h"

// to use as [spin_eq, spin_rec, quality] = mfsa_out_MEX(nspins,nsites,temperature, network, [da dn], n_iterations);

void mexFunction( int nlhs, mxArray *plhs[],     
                  int nrhs, const mxArray *prhs[] ) { 
    
//     prhs[0] = nspins
//     prhs[1] = nsites
//     prhs[2] = temperature
//     prhs[3] = genetic network
//     prhs[4] = [diffusion_auto diffusion_neighbors]
//     prhs[5] = n_iterations dans mfsa 
    
//     plhs[0] = equilibrium state
//     plhs[1] = recording
//     plhs[2] = quality
    
    double *pr_in = mxGetPr(prhs[0]);
    int nspins = pr_in[0];
    pr_in =mxGetPr(prhs[1]);
    int nsites = pr_in[0];
    pr_in =mxGetPr(prhs[2]);
    double temperature = pr_in[0];
    pr_in =mxGetPr(prhs[4]);
    double dauto = pr_in[0];
    double dneig = pr_in[1];
    pr_in =mxGetPr(prhs[5]);
    int n_iterations = pr_in[0];
//     cout << n_iterations << endl;
    int ngrad=1;
    
    Parameters system(nspins, nsites, 1, temperature);
    system.gradient.Construct_simple_gradient(nsites);
    system.network.Fill(mxGetPr(prhs[3]));
    system.neighbors.construction_1D(nsites);    
    Spins spin(nspins, nsites, dauto, dneig);
        
    system.Print_conditions();
    nlhs = 2;
    
    plhs[0] = mxCreateDoubleMatrix(nspins,nsites, mxREAL);
    double *pr_spin_out = mxGetPr(plhs[0]);

    int nb_dim = 3;
    mwSize dim[nb_dim];
    dim[0]= nspins;
    dim[1] = nsites;
    dim[2] = n_iterations;
    plhs[1] = mxCreateNumericArray(nb_dim, dim, mxDOUBLE_CLASS,mxREAL);
    double *pr_record = mxGetPr(plhs[1]);

    plhs[2] = mxCreateDoubleMatrix(1,1, mxREAL);
    double *pr_quality = mxGetPr(plhs[2]);
    
    spin.Fill_rand();
    MFSAexp_asym_rec(spin, system, n_iterations, pr_record);
    copy(spin.state, pr_spin_out, nspins*nsites);
    pr_quality[0] = spin.Quality_max();
    
}


