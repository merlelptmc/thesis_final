#include"library_droso.h"

/// Returns in score the ratio of bins of the histogram where the number is > average
void update_score(int &score, int *histogram, int nbins, int average){
        score= 0;
        double average2 = ((double)average)/nbins;
        for(int i = 0 ; i < nbins; i++){
                if( histogram[i] > average2){
                        score ++;
                }
        }
};

void update_sigma(double &sigma, double sigma_min, double sigma_max, double Qmax,  double Qnew){
        sigma = sigma_max * exp( log( sigma_min /sigma_max) * Qnew / Qmax);
}

int main(int argc, char *argv[]){
        
        int nbins= 100;
        double edges[nbins];
        for(int i=0; i< nbins; i++){
                edges[i] = (double)i/nbins;
        }
        
        //         Print(edges,nbins);
        int    *histogram   = (int*)   calloc( nbins, sizeof(int));    // set to zero
        double *log_density = (double*)calloc( nbins, sizeof(double)); // set to zero
        double density_new = 0;
        double density_old = 0;
        double lnf = 1;
        double f = exp(lnf);
        int average = 0;
        int score   = 0;
        
        //Parameters
        int nspins = 1;
        int nsites = 100;
        double temperature = 1;
        int n_mfsa = 10; //number of iterations in the MFSA
        int ngrad=1;
        //         double basis_network[nspins*(nspins+ngrad)] = {0,0,0,0,0,0};
        double basis_network[nspins*(nspins+ngrad)] = {0, 0};
        
        //Initialization of the system
        Parameters system(nspins, nsites, ngrad, temperature);
        system.gradient.Construct_simple_gradient(nsites);
        system.network.Fill(basis_network);
        system.neighbors.construction_1D(nsites);    
        Spins spin(nspins, nsites );
        
        Parameters system_tested(nspins,nsites,ngrad,temperature);
        system_tested.gradient.Construct_simple_gradient(nsites);
        system_tested.neighbors.construction_1D(nsites);    
        
        //Initialization random walk parameters
        cout<< "Initiating the normal distribution generator" << endl << endl;
        unsigned int mersenneSeed = rand()%100000;
        mt19937_64 generator; // it doesn't matter if I use the 64 or the std::mt19937
        generator.seed(mersenneSeed);
        normal_distribution<double> normal; // default is 0 mean and 1.0 std;
        double sigma = 0;
        double sigma_min = 0.1;
        double sigma_max = 0.3;
        double proba = 0; 
        double a = 0;
        
        //Calculation of the first point
        MFSAexp_asym(spin, system, n_mfsa); 
        double Qold = spin.Quality_max();
        
        int bin_old = floor(Qold*nbins*2);
        histogram[bin_old]++;
        log_density[bin_old] += lnf;
        int bin_new;
        
        average++;
        update_score(score, histogram, nbins, average);
        
        double min_bd = -15;
        double max_bd = 15;
        double over_range=0;
        
        //Opening the results file
        ostringstream prefix, smin, smax, T;
        ofstream result;
        smin << setprecision(2) << sigma_min;
        smax << setprecision(2) << sigma_max;
        T << setprecision(2) << temperature;
        string title = "histogram_" + to_string(nspins) +"s" + to_string(ngrad) + "g_T" + T.str() + "smin" + smin.str() +"smax" + smax.str() + ".txt";
        if( argc > 1){
                prefix << argv[1];                
                title = prefix.str() + "_" +  title;
        }
        result.open(title);
        
        cout << title << endl << endl ; 
        
        int n_loop=0;
        
        while((f-1) > 1e-2){
                
                n_loop++;
                cout << "Iteration " << n_loop << " : " << setprecision(8) << "f =" << f << endl;
                int k=0;
                
                while(score < 0.6*nbins){
                        
                        k++;
                        if( k%10000 == 0){
                                cout << "k = " << k << ", average = " << average << ", score = " << score << ", sigma = "<< sigma << endl;
                                cout << "Histogram";
                                Print(histogram,1,nbins);
                                cout << "Density";
                                Print(log_density,1,nbins);
                                cout << "Average = " << average << " Score = " << score << endl << endl;
                        }
                        
                        //Add random values to the actual network parameters
                        system_tested.network.Fill(system.network.J);
                        over_range = 1;
                        while( over_range == 1){
                                for(int j = 0; j < system_tested.network.nparam; j++){
                                        system_tested.network.J[j] = system.network.J[j] + sigma*normal(generator);
                                }
                                over_range = system_tested.network.Check_bounds(min_bd, max_bd);
                                //                         if (over_range==1){
                                //                                 system_tested.network.Print();
                                //                         }
                        }
                        
                        // Calculate Q
                        MFSAexp_asym(spin, system_tested, n_mfsa); 
                        double Qnew = spin.Quality_max();
                        
                        bin_new = floor(Qnew*nbins*2);
                        bin_old = floor(Qold*nbins*2);
                        proba   = exp(log_density[bin_old]-log_density[bin_new]);
                        if (proba <1){
                          a = drand();
                          if( a < proba){
                                  system.network.Fill(system_tested.network.J);
                                  Qold    = Qnew;
                                  bin_old = bin_new;
                                  histogram[bin_new] ++;
                                  log_density[bin_new] += lnf;
                                  average ++;
                                  update_score(score, histogram, nbins, average);
                                  update_sigma(sigma, sigma_min, sigma_max, 0.5, Qnew);
                          }
                          else{
                                  histogram[bin_old] ++;
                                  log_density[bin_old] += lnf;
                          }
                        }
                        else{
                                system.network.Fill(system_tested.network.J);
                                Qold = Qnew;
                                histogram[bin_old] ++;
                                log_density[bin_old] += lnf;
                                average ++;
                                update_score(score, histogram, nbins, average);
                                update_sigma(sigma, sigma_min, sigma_max, 0.5, Qnew);
                                
                        }
                        
                }
                
                cout << "Histogram";
                Print(histogram,1,nbins);
                cout << "Density";
                Print(log_density,1,nbins);
                
                f = sqrt(f);
                lnf = log(f);
                score   = 0;
                average = 0;
                double min = log_density[0];
                memset(histogram, 0, nbins*sizeof(int));
//                 for(int i=0; i< nbins; i++){
//                         histogram[i] = 0;
//                         if (log_density[i] < min){
//                                 min = log_density[i];
//                         }
//                 }
//                 for(int i=0; i<nbins; i++){
//                         log_density[i] -= min;
//                 }
                
                result << n_loop << "\t" << setprecision(10) << f << "\t" ;
                for(int i=0; i<nbins; i++){
                        result << log_density[i] << "\t";
                }
                result << endl;                
                
        }
        
        
        
        
        return 0;
}
