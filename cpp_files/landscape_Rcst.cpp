#include"library_droso.h"

bool fexists(const std::string& filename) {
        ifstream ifile(filename.c_str());
        return (bool)ifile;
}

int main(int argc, char *argv[]){
        
        int nspins = 1;
        int nsites = 100;
        double temperature = 0.5;
        int n_mfsa = 1000; //number of iterations in the MFSA
        int ngrad=1;
        double *basis_network = (double*)calloc(nspins*(nspins+ngrad),sizeof(double));
//         double basis_network[nspins*(nspins+ngrad)] = {0, 0, 0}; //, 5 , -15, -10};

        double R=2.5;
        double norm=0;
        
        int n_directions = 2; //number of dimensions we explore
//         double exploration[n_directions*4] = {1, 2, 3,-10, -10, -10,0.1, 0.1, 0.1, 10,10, 10}; // how we explore: have to change the code in C++ to make this matrix more readable, the first numbers are the index of the parameters we variate, then the begining value for each, the step for each and the ending value for each
        double exploration[n_directions*4] = {1, 2, -10, -10, 0.1, 0.1, 10, 10}; // how we explore: have to change the code in C++ to make this matrix more readable, the first numbers are the index of the parameters we variate, then the begining value for each, the step for each and the ending value for each

        //Calculation of the total number of points to calculate
        int n_points_tot =1;
        int n_points[n_directions];
        for(int i=0; i<n_directions; i++){
                n_points[i] = ((exploration[3*n_directions+i] - exploration[1*n_directions+i]) / exploration[2*n_directions+i]) +1;
                n_points_tot *= n_points[i];
        }
        
        //Initialization of the system
        Parameters system(nspins, nsites, ngrad, temperature);
        if(ngrad ==1){
        system.gradient.Construct_simple_gradient(nsites);
        }
        if(ngrad==2){
                system.gradient.Construct_opp_gradient(nsites);
        }
        system.network.Fill(basis_network);
        system.neighbors.construction_1D(nsites);    
        Spins spin(nspins, nsites );
        
        
        // creation of the string giving the output files name
        ostringstream T, par, prefix;
        T << fixed << setprecision(1) << temperature;
        string title = "landscape_mfsa_" + to_string(nspins)+"s" + to_string(ngrad) + "g_T" + T.str() + "varparam";
        for(int i=0; i< n_directions; i++){
                par << setprecision(1) << exploration[i];
        }
        title = title + par.str();
        if( argc > 1){
                prefix << argv[1];                
                title = prefix.str() + "_" +  title;
        }
        cout << endl << "file prefix is : " << title << endl;
        
        
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
        parameters << "landscape : " << endl; /*<< "step size sigma = " << sigma << endl << "quality criterion min = " << criterion_mn << endl << "quality criterion max = "<< criterion_mx  << endl << "boundaries = " << min_bd << " "<< max_bd << endl << endl;*/
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
        ofstream landscape;
        test=fexists(title+"results.txt");
        if(test){
                cout << "Please change the name of the output files, a results file with the same name already exists." << endl << endl;
                return 0;
        }
        landscape.open( title + "_results.txt");
        if(!landscape.is_open()){
                cout << "Could not open the .txt file to register data !" << endl;
                return 0;
        }
        landscape << fixed << setprecision(3);
        
        
        //Display informations on the landscape exploration
        cout << endl << "Begining to explore" << endl << "Number of points to calculate = " << n_points_tot << endl << "Number of directions = " << n_directions << endl;
        for(int i=0; i< n_directions; i++){
                cout << "Direction " << exploration[0*n_directions+i] << " de " << exploration[1*n_directions +i] << " à " << exploration[3*n_directions+i] << " par pas de " << exploration[2*n_directions+i]<< endl;
        }
        cout << endl;
        
        cout << "Starting parameters set: " << endl;
        system.network.Print();
        cout << endl;
        
        //Actual exploration
        clock_t timer = clock();
        for(int i=0; i<n_points_tot; i++){
                
                //Update of the network
                for(int j=0; j<n_directions;j++){
                        int modulo = n_points[j];
                        int divid=1;
                        for(int k=j+1; k<n_directions; k++){
                                divid *= n_points[k];
                        }
                        
                        system.network.J[ (int)exploration[j]-1 ] = exploration[1*n_directions+j] + i/divid%(modulo)* exploration[2*n_directions+j];
                }
                
                norm=0;
                for(int j=0; j<system.network.nparam; j++){
                        norm += abs(system.network.J[j]);
                }
                for(int j=0; j<system.network.nparam; j++){
                        system.network.J[j] *= R/norm;
                }
                
                if(i%(n_points_tot/50) ==0 ){ //regular check out
                        cout << i << "ème point ";
                        system.network.Print();
                }
                
                //Calculation of the equilibrium state by MFSA and the associated quality for this network
                MFSAexp_asym(spin, system, n_mfsa); 
                double Q = spin.Quality_max();
                
                //Saving the data
                for(int m=0; m < system.network.nparam; m ++){
                        landscape << system.network.J[m] << "\t\t";
                }
                landscape << setprecision(5) << Q << endl;
                
        }
        
        timer = clock() -timer;
        double elapsed_secs = double(timer) / CLOCKS_PER_SEC ;
        
        cout << endl << endl <<  "Elapsed time for " << n_points_tot << " points is " << elapsed_secs/60 << "min." << endl << endl<< "End of the program." << endl << endl;
        
}


