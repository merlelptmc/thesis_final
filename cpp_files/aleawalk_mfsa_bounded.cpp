#include"library_droso.h"

bool fexists(const std::string& filename) {
        ifstream ifile(filename.c_str());
        return (bool)ifile;
}


int main(int argc, char *argv[] ){
        
        // first argument is the program name
        // second argument will be a prefix for the recording
        // for example : ./myexp or ./mydir/myexp 
        
        // parameters 
        int nspins = 2;
        int ngrad = 1;
        double criterion_mn = 0.68;
        double criterion_mx = 0.75;
        double sigma = 0.1;
        double temperature = 0.5;
        int nsites = 100;
        int n_step = 5000000;
        double n_mfsa=1000;
        double bnd =15;
        double min_bd = -bnd;
        double max_bd = bnd;
        
        // creation of the structures
        Parameters system(nspins, nsites, ngrad, temperature);
        double network[system.network.nparam] = {5, -2.1, 0 , 6, -15, -12};        
        if( ngrad == 1){
                system.gradient.Construct_simple_gradient(nsites);
        }
        if( ngrad ==2){
                system.gradient.Construct_opp_gradient(nsites);        
        }
        system.network.Fill(network);
        system.neighbors.construction_1D(nsites);
        Spins spin(nspins, nsites);
        cout << endl << "System initialized."<< endl << endl;
        
        // creation of the string giving the output files name
        ostringstream T, mx, mn, crit_mn, crit_mx, prefix;
        T << fixed << setprecision(1) << temperature;
        mx << fixed << setprecision(0) << max_bd;
        mn << fixed << setprecision(0) << min_bd; 
        crit_mn << fixed << setprecision(3) << criterion_mn;
        crit_mx << fixed << setprecision(3) << criterion_mx;
        prefix << argv[1];
        string title = "aw_mfsa_" + to_string(nspins)+"s" + to_string(ngrad) + "g_T" + T.str() + "_bnd" + mn.str() + "to" + mx.str() + "_crit" + crit_mn.str() + "to" + crit_mx.str();
        if( argc > 1){
                title = prefix.str() + "_" +  title;
        }
        cout << "file prefix is : " << title << endl;
        
        time_t t = time(0);   // get time now
        struct tm * now = localtime( & t );
        stringstream ss;
        ss << (now->tm_year + 1900) << '-'<< (now->tm_mon + 1) << '-'
        <<  now->tm_mday  << endl;
        
        // recording the parameters of the walk in an annex file
        ofstream parameters;
        bool test=fexists(title + "_parameters.txt");
        if(test){
                cout << "Please change the name of the output files, a paramters file with the same name already exists." << endl << endl;
                return 0;
        }
        parameters.open( title + "_parameters.txt");
        parameters << "date : " << ss.str() << "system:" << endl << endl;
        parameters << "nspins = " << nspins << endl << "ngradients = " << ngrad << endl; 
        parameters << "nsites = " << nsites << endl << "temperature = " << temperature << endl << endl;
        parameters << "walk : " << endl << "step size sigma = " << sigma << endl << "quality criterion min = " << criterion_mn << endl << "quality criterion max = "<< criterion_mx  << endl << "boundaries = " << min_bd << " "<< max_bd << endl << endl;
        parameters << "calculation : MFSA " << endl << "niterations = " << n_mfsa << endl << endl;
        parameters << "Starting network : " << endl;
        for(int i=0; i<system.network.nparam; i++){
                parameters << system.network.J[i] << "\t";
                if((i+1)%(nspins)==0){ 
                        parameters << endl;
                }
                
        }
        cout << "Parameters registered." << endl;
        
        // opening the file containing the results
        ofstream walk;
        test=fexists(title+"results.txt");
        if(test){
                cout << "Please change the name of the output files, a results file with the same name already exists." << endl << endl;
                return 0;
        }
        walk.open( title + "_results.txt");
        if(!walk.is_open()){
                cout << "Could not open the .txt file to register data !" << endl;
                return 0;
        }
        walk << fixed << setprecision(3);
        
        // initializing the normal distribution random number generator
        cout<< "Initiating the normal distribution generator" << endl << endl;
        unsigned int mersenneSeed = rand()%100000;
        mt19937_64 generator; // it doesn't matter if I use the 64 or the std::mt19937
        generator.seed(mersenneSeed);
        normal_distribution<double> normal; // default is 0 mean and 1.0 std;
        
        double over_range;
        sigma /= sqrt(system.network.nparam);
        
        // calculating the values for the starting point and checking that it is in the boundaries and below the choosen criterion 
        cout << "Starting point :";
        system.network.Print();
        cout << endl;
        MFSAexp_asym(spin, system, n_mfsa);
        double Q = spin.Quality_max();
        over_range = system.network.Check_bounds(min_bd, max_bd);
        clock_t timer = clock();
        if(over_range){
                cout << endl << " Starting point out of range ! stopping now" << endl;
        }
        else{
                cout << "Quality starting point = " << Q ; 
                if(Q<criterion_mn || Q > criterion_mx){
                        cout << endl << "Quality not in range !! Criterion must be between " << criterion_mn << " and " << criterion_mx  << ", stopping now." << endl;
                }
                else
                {
                        cout << endl <<"Beneath range and quality sufficient" << endl<< endl;
                        
                        Gene_network net_save(nspins,ngrad,system.network.J);
                        
                        cout << "Starting the loop " << endl << "nstep = " << n_step << endl << endl;
                        for(int i=0; i<n_step; i++){
                                if( i %(n_step/10) == 0){
                                        cout << "step n " << i << endl;
                                }
                                
                                over_range =1;
                                while( over_range ==1){
                                        for(int j=0; j< system.network.nparam; j++){
                                                system.network.J[j] = net_save.J[j] + sigma*normal(generator);
                                        }
                                        over_range = system.network.Check_bounds(min_bd, max_bd);
                                }
                                MFSAexp_asym(spin, system, n_mfsa);
                                double Q = spin.Quality_max();
                                for(int j=0; j< system.network.nparam; j++){
                                        walk << system.network.J[j] << "\t\t ";
                                }
                                walk <<  Q << endl;
                                
                                if( Q > criterion_mn && Q < criterion_mx ){
                                        net_save.Fill(system.network.J);
                                }
                                else{
                                        system.network.Fill(net_save.J);
                                }
                        }
                        
                }
        }
        
        timer = clock() -timer;
        double elapsed_secs = double(timer) / CLOCKS_PER_SEC ;
        
        cout << endl << endl <<  "Elapsed time for " << n_step << " points is " << elapsed_secs/60 << "min." << endl << endl<< "End of the program." << endl << endl;
        
        return 0;
}











