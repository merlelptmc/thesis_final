#include"library_droso.h"

bool fexists(const std::string& filename) {
        ifstream ifile(filename.c_str());
        return (bool)ifile;
}

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
        
        // Parameters of the Wang Landau
        int nbins= 100;        
        double crit_fmin = 1e-3;
        double crit_score_flat =0.8;
        
        // Quantities followed during WL
        int    *histogram   = (int*)   calloc( nbins, sizeof(int));    // set to zero
        double *log_density = (double*)calloc( nbins, sizeof(double)); // set to zero
        
        // Quantities used during WL 
        double lnf = 1;
        double f = exp(lnf);
        double density_new = 0;
        double density_old = 0;
        double Qnew=0;
        double Qold=0;
        int bin_old=0;
        int bin_new=0;
        int average = 0;
        int score   = 0;
        double proba = 0; 
        double a = 0;
        
        //Parameters of th studied system
        int nspins = 1; // number of genes/spins
        int nsites = 100; // number of nuclei/sites
        double temperature = 1;
        int n_mfsa = 1000; //number of iterations in the MFSA algorithm
        int ngrad=1; // number of maternel gradients if ngrad=2, need to change the graident filling function
        
        // Boundaries of the parameters
        double min_bd = -15;
        double max_bd = 15;
        double over_range=0;
        
        //Initialization of the system
        Parameters system(nspins, nsites, ngrad, temperature);
        system.gradient.Construct_simple_gradient(nsites);
        system.network.Init_rand(nspins,nsites,max_bd);
        system.neighbors.construction_1D(nsites);    
        Spins spin(nspins, nsites );
        
        // Initialization of the tested system during a WL step
        Parameters system_tested(nspins,nsites,ngrad,temperature);
        system_tested.gradient.Construct_simple_gradient(nsites);
        system_tested.neighbors.construction_1D(nsites);    
        
        //Initialization of random walk parameters
        cout<< "Initiating the normal distribution generator" << endl << endl;
        unsigned int mersenneSeed = rand()%100000;
        mt19937_64 generator; // it doesn't matter if I use the 64 or the std::mt19937
        generator.seed(mersenneSeed);
        normal_distribution<double> normal; // default is 0 mean and 1.0 std;
        double sigma = 5;
        
        //Calculation of the first point and Initializationof histogram and log_density
        MFSAexp_asym(spin, system, n_mfsa); 
        Qold = spin.Quality_max();
        bin_old = floor(Qold*nbins*2);
        histogram[bin_old]++;
        log_density[bin_old] += lnf;
        
        average++;
        update_score(score, histogram, nbins, average);
        
        //Filling parameters file and opening the results file
        ostringstream prefix, sig, T, nb;
        ofstream result, parameters;
        sig << fixed << setprecision(1) << sigma;
        T << fixed <<  setprecision(1) << temperature;
        nb << setprecision(0) << nbins;
        string title = "WL_" + to_string(nspins) +"s" + to_string(ngrad) + "g_T" + T.str()  + "_nbins" + nb.str() + "_sigmafix" + sig.str();
        if( argc > 1){
                prefix << argv[1];                
                title = prefix.str() + "_" +  title;
        }
        bool test=fexists(title + "_parameters.txt");
        if(test){
                cout << "Please change the name of the output files, a paramters file with the same name already exists." << endl << endl;
                return 0;
        }
        
        cout << "Prefix is " << title << endl << endl;
        
        parameters.open(title + "_parameters.txt");        
        time_t t = time(0);   // get time now
        struct tm * now = localtime( & t );
        stringstream ss;
        ss << (now->tm_year + 1900) << '-'<< (now->tm_mon + 1) << '-'
        <<  now->tm_mday  << endl;
        
        parameters << "date : " << ss.str() << endl << "system:" << endl;
        parameters << "nspins = " << nspins << endl << "ngradients = " << ngrad << endl; 
        parameters << "nsites = " << nsites << endl << "temperature = " << temperature << endl << endl;
        parameters << "wang-landau : " << endl  << "sigma = " << sigma << endl << "nbins = " << nbins << endl << "flatening critera = "<< crit_score_flat  << endl << "ending critera : f-1 < " << crit_fmin << endl << "boundaries = " << min_bd << " "<< max_bd << endl << endl;
        parameters << "calculation : MFSA " << endl << "niterations = " << n_mfsa << endl << endl;
        parameters << "Starting network : " << endl;
        for(int i=0; i<system.network.nparam; i++){
                parameters << system.network.J[i] << "\t";
                if((i+1)%(nspins)==0){ 
                        parameters << endl;
                }
                
        }
        cout << "Parameters registered." << endl;
        parameters.close();
        result.open(title + "_histograms.txt");
        
        int n_loop=0;        
        while((f-1) > crit_fmin){
                
                n_loop++;
                cout << "Iteration " << n_loop << " : " << setprecision(8) << "f =" << f << endl;
                int k=0;
                
                while(score < crit_score_flat*nbins){
                        
                        k++;
                        if( k%100 == 0){
                                cout << "(" << n_loop << ',' << k << ")" << ", average = " << average << ", score = " << score << ", sigma = "<< sigma << endl;
                                
                                                                cout << "Histogram";
                                                                Print(histogram,1,nbins);
                                                                cout << "Density";
                                                                Print(log_density,1,nbins);
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
                        Qnew = spin.Quality_max();
                        
                        bin_new = floor(Qnew*nbins*2);
                        proba   = exp(log_density[bin_old]-log_density[bin_new]);
                        
                        if( drand() < proba){
                                system.network.Fill(system_tested.network.J);
                                bin_old=bin_new;
                        }
                        
                        histogram[bin_old] ++;
                        log_density[bin_old] += lnf; 
                        average ++;
                        update_score(score, histogram, nbins, average);
                        
                        
                }
                
                
                cout << "Histogram";
                Print(histogram,1,nbins);
                cout << "Density";
                Print(log_density,1,nbins);
                
                //                 f = pow(f,0.8);
                f=sqrt(f);
                cout << f << endl;
                lnf = log(f);
                score   = 0;
                average = 0;
                double min = log_density[0];
                memset(histogram, 0, nbins*sizeof(int));
                
                // Removing the minimum value of log(density) to keep registered values low
                for(int i=0; i< nbins; i++){
                        if (log_density[i] < min){
                                min = log_density[i];
                        }
                }
                for(int i=0; i<nbins; i++){
                        log_density[i] -= min;
                }
                
                cout << "Density";
                Print(log_density,1,nbins);
                
                result << n_loop << "\t" << setprecision(10) << f << "\t" ;
                for(int i=0; i<nbins; i++){
                        result << log_density[i] << "\t";
                }
                result << endl;                
                
                
        }                             
        
        return 0;
}
