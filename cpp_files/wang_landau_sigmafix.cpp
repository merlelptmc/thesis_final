#include"library_droso.h"


void update_score(double &score, double *histogram, int nbins, double average){
        score= 0;
        for(int i=0; i< nbins; i++){
                if( histogram[i] > average){
                        score ++;
                }
        }
        score /= nbins;        
};

void update_sigma(double &sigma, double Qnew){
        double sig_max = 0.5;
        sigma = sig_max*exp(log( 0.1 /sig_max)* Qnew / 0.5);
        
}

int main(){
        
        int nbins= 50;
        double edges[nbins];
        for(int i=0; i< nbins; i++){
                edges[i] = (double)i/nbins;
        }
        
        //         Print(edges,nbins);
        double *histogram=(double*)calloc( nbins, sizeof(double));
        double *log_density=(double*)calloc( nbins, sizeof(double));
        for(int i=0; i < nbins; i++){
                log_density[i]=1;
        }
        double density_new=0;
        double density_old=0;
        double f = exp(1);
        double average =0;
        double score=0;
        
        //Parameters
        int nspins = 1;
        int nsites = 100;
        double temperature = 1;
        int n_mfsa = 1000; //number of iterations in the MFSA
        int ngrad=1;
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
        
        //Initialization random
        cout<< "Initiating the normal distribution generator" << endl << endl;
        unsigned int mersenneSeed = rand()%100000;
        mt19937_64 generator; // it doesn't matter if I use the 64 or the std::mt19937
        generator.seed(mersenneSeed);
        normal_distribution<double> normal; // default is 0 mean and 1.0 std;
        double sigma=0.1;
        double proba=0; double a=0;
        
        //Calculation of the first point
        MFSAexp_asym(spin, system, n_mfsa); 
        double Qold = spin.Quality_max();
        
        int bin = floor(Qold*nbins*2);
        histogram[bin] ++;
        log_density[bin] += log(f);
        
        average += (double)1/nbins;
        update_score(score,histogram,nbins, average);
        
        double min_bd = -15;
        double max_bd = 15;
        double over_range=0;
        
        ofstream result;
        result.open("histogram_sigmamax0.5.txt");
        
        int n_loop =0;
        while((f-1) > 10^(-3)){
                
                cout << setprecision(8) << "f =" << f << endl;
                int k=0;
                n_loop++;                
                
                while(score < 0.6){
                        
                        k++;
                        if( k%100 == 0){
                                cout << "k = " << k << ", average = " << average << ", score = " << score << ", sigma = "<< sigma << endl;
                        }
                        
                        //Add random values to the actual network parameters
                        system_tested.network.Fill(system.network.J);
                        over_range =1;
                        while( over_range ==1){
                                for(int j=0; j< system_tested.network.nparam; j++){
                                        system_tested.network.J[j] = system.network.J[j] + sigma*normal(generator);
                                }
                                over_range = system_tested.network.Check_bounds(min_bd, max_bd);
                       }
                        
                        // Calculate Q
                        MFSAexp_asym(spin, system_tested, n_mfsa); 
                        double Qnew = spin.Quality_max();
                        
                        bin = floor(Qnew*nbins*2);
                        density_new = exp(log_density[bin]);
                        bin = floor(Qold*nbins*2);
                        density_old = exp(log_density[bin]);
                        if (density_old/density_new <1){
                                proba = density_old/density_new;
                                a =(double)(rand()%(10^8))/(double)(10^8);
                                if( a < proba){
                                        system.network.Fill(system_tested.network.J);
                                        Qold=Qnew;
                                        bin = floor(Qold*nbins*2);
                                        histogram[bin] ++;
                                        log_density[bin] += log(f);
                                        average += (double)1/nbins;
                                        update_score(score,histogram,nbins, average);
                                        update_sigma(sigma, Qnew);
                                }
                                else{
                                        bin = floor(Qold*nbins*2);
                                        histogram[bin] ++;
                                        log_density[bin] += log(f);
                                        
                                }
                        }
                        else{
                                system.network.Fill(system_tested.network.J);
                                Qold=Qnew;
                                bin = floor(Qold*nbins*2);
                                histogram[bin] ++;
                                log_density[bin] += log(f);
                                average += (double)1/nbins;
                                update_score(score,histogram,nbins, average);
                                
                                
                        }
                        
                }
                
                Print(histogram,1,nbins);
                
                f = sqrt(f);
                score =0;
                average=0;
                double min = log_density[0];
                for(int i=0; i< nbins; i++){
                        histogram[i] = 0;
                        if (log_density[i] < min){
                                min = log_density[i];
                        }
                }
                for(int i=0; i<nbins; i++){
                        log_density[i] -= min;
                }
                
                result << n_loop << "\t" << f << "\t" ;
                for(int i=0; i<nbins; i++){
                        result << log_density[i] << "\t";
                }
                result << endl;                
                
        }
        
        
        
        
        return 0;
}
