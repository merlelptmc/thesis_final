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

void update_sigma(double &sigma, double sig_min, double sig_max, double Qnew, double Qmax){
        sigma = sig_max*exp(log( sig_min /sig_max)* Qnew / Qmax);
        
}

int main(int argc, char *argv[]){
        
        srand(time(NULL));
        
        // Parameters of the Wang Landau
        int nbins= 150;        
        double crit_fmin = 1e-4;
        double crit_score_flat =0.4;
        
        // Quantities followed during WL
        int    *histogram   = (int*)   calloc( nbins, sizeof(int));    // set to zero
        double *log_density = (double*)calloc( nbins, sizeof(double)); // set to zero
        
        // Quantities used during WL 
        double lnf = 1;
        double f = exp(lnf);
        double pow_f= 0.9;
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
        int nspins = 2; // number of genes/spins
        int nsites = 100; // number of nuclei/sites
        double temperature = 1;
        int n_mfsa = 1000; //number of iterations in the MFSA algorithm
        int ngrad=1; // number of maternel gradients if ngrad=2, need to change the graident filling function
        double Qmax=0.75;
        
        // Boundaries of the parameters
        double min_bd = -10;
        double max_bd = 10;
        double over_range=0;
        
        //Initialization of the system
        Parameters system(nspins, nsites, ngrad, temperature);
        system.gradient.Construct_simple_gradient(nsites);
        system.network.Init_rand(nspins,ngrad,min_bd, max_bd);
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
        double sigma_min = 0.1;
        double sigma_max = 4;
        
        //Calculation of the first point and Initializationof histogram and log_density
        MFSAexp_asym(spin, system, n_mfsa); 
        Qold = spin.Quality_max();
        bin_old = floor(Qold*nbins/Qmax);
        histogram[bin_old]++;
        log_density[bin_old] += lnf;
        
        average++;
        update_score(score, histogram, nbins, average);
        double sigma =0;
        update_sigma(sigma, sigma_min, sigma_max, Qold, Qmax);
        
        //Filling parameters file and opening the results file
        ostringstream prefix, sig_min, sig_max, T, nb;
        ofstream result, parameters, walk;
        sig_min << fixed << setprecision(1) << sigma_min;
        sig_max << fixed << setprecision(1) << sigma_max;
        T << fixed <<  setprecision(1) << temperature;
        nb << setprecision(0) << nbins;
        string title = "WL_" + to_string(nspins) +"s" + to_string(ngrad) + "g_T" + T.str()  + "_nbins" + nb.str() + "_sigma" + sig_min.str() + "to" + sig_max.str() ;
        if( argc > 1){
                prefix << argv[1];                
                title = prefix.str() + "_" +  title;
        }
        bool test=fexists(title + "_parameters.txt");
        cout << "Prefix is " << title << endl << endl;
        if(test){
                cout << "Please change the name of the output files, a paramters file with the same name already exists." << endl << endl;
                return 0;
        }
        
        
        parameters.open(title + "_parameters.txt");        
        time_t t = time(0);   // get time now
        struct tm * now = localtime( & t );
        stringstream ss;
        ss << (now->tm_year + 1900) << '-'<< (now->tm_mon + 1) << '-'
        <<  now->tm_mday  << endl;
        
        parameters << "date : " << ss.str() << endl << "system:" << endl;
        parameters << "nspins = " << nspins << endl << "ngradients = " << ngrad << endl << "nsites = " << nsites << endl << "temperature = " << temperature << endl << "Qmax = " << Qmax << endl;
        parameters << "boundaries = " << min_bd << " "<< max_bd << endl << endl;
        parameters << "wang-landau : " << endl  << "sigma min = " << sigma_min << endl << "sigma max = " << sigma_max << endl << "nbins = " << nbins << endl << "flatening critera = "<< crit_score_flat  << endl << "ending critera : f-1 < " << crit_fmin << endl <<  "first value of lnf = " << lnf << endl << "pow_f = " << pow_f << endl << endl;
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
        walk.open(title + "_walk.txt");
        
        int n_loop=0;        
        clock_t timer0 = clock();
        clock_t timer = timer0;
        double elapsed_secs;
        
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
                        
                        bin_new = floor(Qnew*nbins/Qmax);
                        proba   = exp(log_density[bin_old]-log_density[bin_new]);
                        
                        if( drand() < proba){
                                system.network.Fill(system_tested.network.J);
                                bin_old=bin_new;
                                update_sigma(sigma, sigma_min, sigma_max, Qnew, Qmax);
                        }
                        
                        histogram[bin_old] ++;
                        log_density[bin_old] += lnf; 
                        average ++;
                        update_score(score, histogram, nbins, average);

                        
                        for(int i=0; i< system.network.nparam; i++){
                                walk << system.network.J[i] << "\t" ;
                        }
                        walk << endl;
                        
                        
                }
                
                
                //                 cout << "Histogram";
                //                 Print(histogram,1,nbins);
                //                 cout << "Density";
                //                 Print(log_density,1,nbins);
                
                //                 f = pow(f,0.8);
                f=pow(f, pow_f);
                //                 cout << f << endl;
                lnf = log(f);
                score   = 0;
                average = 0;
                memset(histogram, 0, nbins*sizeof(int));
                
                // Removing the minimum value of log(density) to keep registered values low
                double min = log_density[0];
                for(int i=0; i< nbins; i++){
                        if (log_density[i] < min){
                                min = log_density[i];
                        }
                }
                for(int i=0; i<nbins; i++){
                        log_density[i] -= min;
                }
                
                timer = clock() -timer;
                elapsed_secs = double(timer) / CLOCKS_PER_SEC ;
                cout <<  "Elapsed time for this loop is " << elapsed_secs/60 << "min." << endl;
                
                //                 cout << "Density";
                //                 Print(log_density,1,nbins);
                
                result << n_loop << "\t" << setprecision(10) << f << "\t" << elapsed_secs << "\t" << k << "\t"   ;
                for(int i=0; i<nbins; i++){
                        result << log_density[i] << "\t";
                }
                result << endl;    
                
                
        }                             
        
        // Measure execution time
        clock_t timer_tot = clock() -timer0;
        elapsed_secs = double(timer_tot) / CLOCKS_PER_SEC ;
        
        //         cout << endl << endl <<  "Elapsed time total is " << elapsed_secs/60 << "min." << endl << "End of the program." << endl << endl;
        
        return 0;
}
