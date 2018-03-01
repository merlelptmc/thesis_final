#include "mex.h"
#include "matrix.h"
#include "/users/invites/merle/thesis_final_folder/cpp_files/libraries/library_droso.h"

// to use as [spin_eq, spin_rec, quality] = 
//                     mcmc_out_MEX(nspins,nsites,temperature, network, [n_eq n_rec step_rec ],  spin_ini);

void mexFunction( int nlhs, mxArray *plhs[],     
                  int nrhs, const mxArray *prhs[] ) { 
        
        //     prhs[0] = nspins
        //     prhs[1] = nsites
        //     prhs[2] = temperature
        //     prhs[3] = genetic network
        //     prhs[4] = [n_equilibration n_recording step_rec]
        //     prhs[5] = spin_initial
        
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
        int n_eq = pr_in[0];
        int n_rec = pr_in[1];
        int step_rec = pr_in[2];
        int ngrad=1;
        
        Parameters system(nspins, nsites, 1, temperature);
        system.gradient.Construct_simple_gradient(nsites);
        system.network.Fill(mxGetPr(prhs[3]));
        system.neighbors.construction_1D(nsites);    
        Spins spin(nspins, nsites);
        spin.Fill(mxGetPr(prhs[5]));
        
        system.Print_conditions();
        spin.Print();        
        
        nlhs = 2;        
        plhs[0] = mxCreateDoubleMatrix(nspins,nsites, mxREAL);
        double *pr_spin_out = mxGetPr(plhs[0]);
        
        int nb_dim = 3;
        mwSize dim[nb_dim];
        dim[0]= nspins;
        dim[1] = nsites;
        dim[2] = n_rec;
        plhs[1] = mxCreateNumericArray(nb_dim, dim, mxDOUBLE_CLASS,mxREAL);
        double *pr_record = mxGetPr(plhs[1]);
        
        plhs[2] = mxCreateDoubleMatrix(1,1, mxREAL);
        double *pr_quality = mxGetPr(plhs[2]);
        

        MCMC(spin, system, n_eq, n_rec, step_rec );
        
        double quality = spin.Quality_max();
        
        cout << endl << "Final Spin" << endl;
        spin.Print();        
        cout << "Quality of the pattern (0<Q<1) = " << quality << endl << endl;
        copy(spin.state, pr_spin_out, nspins*nsites);
        pr_quality[0] = spin.Quality_max();
        
                          
        }
                  
