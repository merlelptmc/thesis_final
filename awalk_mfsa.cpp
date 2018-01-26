#include"library_droso.h"

int main(){
    
    int nspins = 2;
    int ngrad = 1;
    double criterion = 0.45;
    double sigma = 0.1;
    double temperature = 0.5;
    double dauto =1;
    double dneigh =1;
    int nsites = 100;
    int n_step = 100;
    double bnd =15;
    double min_bd = -bnd;
    double max_bd = bnd;
    
    string title = "mfaw_" + to_string(nspins)+"s" + to_string(ngrad)+"g_T" to_string(temperature) ;
    
    ofstream parameters;    
    parameters.open( title + ".txt");
    parameters << "nspins = " << 2 << cd 
    
    Parameters system(nspins, nsites, ngrad, temperature);
    if( ngrad == 1){
        system.gradient.Construct_simple_gradient(nsites);
    }
    if( ngrad ==2){
        system.gradient.Construct_opp_gradient(nsites);        
    }
    system.network.Fill(mxGetPr(prhs[0]));
    system.neighbors.construction_1D(nsites);
    Spins spin(nspins, nsites, dauto, dneig);
    
    cout << "System created" << endl << "Creating the output .txt files" << endl;
    ofstream netall, netgood, qualall, qualgood; 
    
    
    cout<< "Initiating the normal distribution generator" << endl;
    unsigned int mersenneSeed = rand()%100000;
    mt19937_64 generator; // it doesn't matter if I use the 64 or the std::mt19937
    generator.seed(mersenneSeed);
    normal_distribution<double> normal; // default is 0 mean and 1.0 std0;
    
    double over_range;
    sigma /= sqrt(system.network.nparam);
    
    cout << "Starting point :" << endl;
    system.network.Print();
    MFSAexp_asym(spin, system, n_iterations);
    double Q = spin.Quality_max();
    over_range = system.network.Check_bounds(min_bd, max_bd);
    if(over_range){
        cout << endl << " Starting point out of range ! stopping now" << endl;
    }
    else{
        cout << "Quality starting point = " << Q ; 
        if(Q<criterion){
            cout << endl << "Quality not sufficient !! Criterion = " << criterion << "stopping now" << endl;
        }
        else
        {
            cout << endl <<"Beneath range and quality sufficient" << endl;
            
            Gene_network net_save(nspins,ngrad, system.network.J);
            int k_kept = 0;
            int index[n_step];
            
            cout << "Starting the loop " << endl ;
            for(int i=0; i<n_step; i++){
                if( i %(n_step/10) == 0){
                    cout << "step n " << i << endl;
                }
                
                over_range =1;
                while( over_range ==1){
                    for(int j=0; j< system.network.nparam; j++){
                        system.network.J[j] = net_save.J[j] + sigma*normal(generator);
                    }
                    over_range = check_bounds(system.network, min_bd, max_bd);
                }
                MFSAexp_asym(spin, system, n_iterations);
                double Q = spin.Quality_max();
                for(int j=0; j< system.network.nparam; j++){
                    pr_out[j*n_step+i] = system.network.J[j];
                }
                pr_quality[i] = Q;
                //                 cout << "i = " << i << "    Q = " << Q << endl;
                
                if( Q > criterion){
                    index[k_kept] = i;
                    net_save.Fill(system.network.J);
                    k_kept++;
                }
                else{
                    system.network.Fill(net_save.J);
                }
            }
            
            if( k_kept > n_step){
                mexErrMsgTxt("Error in the count of kept parameters sets ! ");
            }
            else{
                
                plhs[2] = mxCreateDoubleMatrix(k_kept, system.network.nparam, mxREAL);
                double *pr_kept= mxGetPr(plhs[2]);
                
                plhs[3] = mxCreateDoubleMatrix(k_kept,1, mxREAL);
                double *pr_quality_kept = mxGetPr(plhs[3]);
                
                for(int i=0; i<k_kept; i++){
                    for(int j=0; j<system.network.nparam; j++){
                        //                                             cout << pr_kept[j*k_kept+i] << "  " << index[i] << "  " << j*n_step+index[i] <<  "   " << pr_out[ j*n_step + index[i]] << endl;
                        pr_kept[j*k_kept + i] = pr_out[j*n_step + index[i]];
                    }
                }
                
                for(int i=0; i<k_kept; i++){
                    for(int j=0; j<system.network.nparam; j++){
                        //                         cout << pr_quality_kept[i] << "  " << index[i] <<  "   " << pr_quality[index[i]] << endl;
                        pr_quality_kept[i] = pr_quality[index[i]];
                    }
                }
                
            }
            
            
        }
    }
    
    
    return 0;
}











