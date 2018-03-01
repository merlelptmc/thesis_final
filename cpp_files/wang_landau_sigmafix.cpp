#include"library_droso.h"

bool fexists(const std::string& filename) {
        ifstream ifile(filename.c_str());
        return (bool)ifile;
}

void update_score(int &score, int *histogram, int nbins, int average){
        score= 0;
        double average2 = ((double)average)/nbins;
        for(int i = 0 ; i < nbins; i++){
                if( histogram[i] > average2){
                        score ++;
                }
        }
};

void update_sigma(double &sigma, double Qnew){
        double sig_max = 0.5;
        sigma = sig_max*exp(log( 0.1 /sig_max)* Qnew / 0.5);
        
}

int main(int argc, char *argv[]){
        
        int nbins= 100;
        double edges[nbins];
        for(int i=0; i< nbins; i++){
                edges[i] = (double)i/nbins;
        }
        
        //         Print(edges,nbins);
        int *histogram=(int*)calloc( nbins, sizeof(int));
        double *log_density=(double*)calloc( nbins, sizeof(double));
        double density_new=0;
        double density_old=0;
        double f = exp(1);
        double crit_fmin = pow(10,-5);
        int average =0;
        int score=0;
        double crit_score_flat = 0.8;
        
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
        
        //Initialization random generator
        cout<< "Initiating the normal distribution generator" << endl << endl;
        unsigned int mersenneSeed = rand()%100000;
        mt19937_64 generator; // it doesn't matter if I use the 64 or the std::mt19937
        generator.seed(mersenneSeed);
        normal_distribution<double> normal; // default is 0 mean and 1.0 std;
        double sigma=0.5;
        double proba=0; double a=0;
        
        //Calculation of the first point
        MFSAexp_asym(spin, system, n_mfsa); 
        double Qold = spin.Quality_max();
        double Qnew=0;
        
        int bin = floor(Qold*nbins*2);
        histogram[bin] ++;
        log_density[bin] += log(f);
        
        average ++;
        update_score(score,histogram,nbins, average);
        
        double min_bd = -15;
        double max_bd = 15;
        double over_range=0;
        
        //Opening the results file
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
        result.open(title + "_histograms.txt");
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
        
        cout << title << endl << endl ; 
        
        int n_loop =0;
        while((f-1) > crit_fmin){
                
                cout << setprecision(8) << "f =" << f << endl;
                int k=0;
                n_loop++;                
                
                while(score < crit_score_flat*nbins){
                        
                        k++;
                        if( k%10000 == 0){
                                cout << "(" << n_loop << ',' << k << ")" << ", average = " << average << ", score = " << score << ", sigma = "<< sigma << endl;
                                
//                                 cout << "Histogram";
//                                 Print(histogram,1,nbins);
//                                 cout << "Density";
//                                 Print(log_density,1,nbins);
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
                        
                        system.network.Print();
                        system_tested.network.Print();
                        cout << "Q " << Qnew << endl;
                        
                        bin = floor(Qnew*nbins*2);      
                        density_new = exp(log_density[bin]);
                        bin = floor(Qold*nbins*2);
                        density_old = exp(log_density[bin]);
                        cout << "old " << density_old << " new " << density_new ;
                        if (density_old/density_new <1){
                                proba = density_old/density_new;
                                a =(double)(rand()%(10^8))/(double)(10^8);
                                cout << endl << " proba " << proba << " a " << a ;
                                if( a < proba){
                                        system.network.Fill(system_tested.network.J);
                                        Qold=Qnew;
                                        bin = floor(Qold*nbins*2);
                                        histogram[bin] ++;
                                        log_density[bin] += log(f);
                                        average ++;
                                        update_score(score,histogram,nbins, average);
                                        cout << " taken !" << endl;
                                }
                                else{
                                        bin = floor(Qold*nbins*2);
                                        histogram[bin] ++;
                                        log_density[bin] += log(f);
                                        cout << " left out ! " << endl;
                                        
                                }
                        }
                        else{
                                system.network.Fill(system_tested.network.J);
                                Qold=Qnew;
                                bin = floor(Qold*nbins*2);
                                histogram[bin] ++;
                                log_density[bin] += log(f);
                                average ++;
                                update_score(score,histogram,nbins, average);
                                cout << " taken !" << endl;
                                
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
