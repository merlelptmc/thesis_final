#include"library_droso.h"
using namespace std;


// Function to display a matrix of size nlin*ncol
void Print(double *mat, int nlin, int ncol ){
        cout << endl;
        for(int i=0; i<nlin; i++){
                for (int j=0; j<ncol; j++){
                        cout <<  setprecision(4) << mat[j*nlin+i] << "\t" ;
                }
                cout << endl;
        }
        cout << endl;
}


//Functions associated to the gene network 
Gene_network::Gene_network(int ng, int nm){
        Init(ng,nm);
}

Gene_network::Gene_network(int ng, int nm, double *M){
        Init(ng,nm);
        Fill(M);
}

Gene_network::~Gene_network(){
        free(J);
}

void Gene_network::Init(int ng, int nm){
        nspins=ng;
        ngrad=nm;
        nparam=nspins*(nspins+ngrad);
        J = (double*)calloc(nparam,sizeof(double));
        C = J +nspins*nspins;
}

void Gene_network::Init_rand(int ng, int nm, int abs_max){
        Init(ng,nm);
        srand(time(0));
        for(int i=0; i<nparam; i++){
                J[i] = -abs_max + rand() %(2*abs_max);
        }
}

void Gene_network::Fill(double *M){
        for(int i=0; i<nparam; i++){
                J[i] = M[i];
        }
}

void Gene_network::Print(){
        cout << endl <<  "Gene network :" << endl << "J = " << endl;
        for(int i=0; i< nspins; i++){
                for(int j=0; j<nspins;j++){
                        cout << setprecision(3) << J[i*nspins+j] << "\t";
                }
                cout << endl;
        }
        cout <<  endl << "C = " << endl;
        for(int i=0; i< ngrad; i++){
                for(int j=0; j< nspins; j++){
                        cout << setprecision(3) << C[i*nspins + j] << "\t" ;
                }
                cout << endl;
        }
}


//Functions associated to the structure Spatial_grid
Spatial_grid::~Spatial_grid(){
        free(sites);
        free(index);
        free(number);
        nsites=0;
        length=0;
        thickness=0;
        nlinks=0;
        index=0;
        number=0;
        sites=0;
}

void Spatial_grid::Init(int L, int l){
        length=L;
        thickness=l;
        nsites=l*L;
        nlinks=0;
        index  = (int*)malloc(nsites * sizeof(int));
        number = (int*)malloc(nsites * sizeof(int));
}

int Spatial_grid::construction_1D(int ns){
        
        if(ns <2){
                cout << endl << "Error, neighbors matrix must be constructed on more than one site." << endl;
        }
        
        Init(ns,1);
        nlinks = 3*nsites - 2;
        sites  = (int*)malloc(nlinks* sizeof(int));
        int iL=0;
        int iS=0;
        
        index[iS] = iL;
        number[iS]=2;
        sites[iL]=iS;
        sites[iL+1]=iS+1;
        iL += number[iS];
        
        for(int iS=1; iS<nsites-1; iS++){
                index[iS] = iL;
                number[iS] = 3;
                sites[iL]=iS-1;
                sites[iL+1]=iS;
                sites[iL+2]=iS+1;
                iL = iL+ 3;
        }
        
        iS=nsites-1;       
        index[iS] = iL;
        number[iS]=2;
        sites[iL]=iS-1; 
        sites[iL+1]=iS;
        iL += number[iS];
        iS++;
        
        if(iS != nsites){
                cout << "Nsites=" << nsites << " iS=" << iS << endl;
                return 1;
        }
        if(iL!=nlinks){
                cout << "Nlinks=" << nlinks << " iL=" << iL << endl;
                return 1;
        }
        return 0;
}  

int Spatial_grid::construction_2D_rectangle_square(int L, int l){
        
        Init(L, l);
        
        if(l<2 || L<2){
                cout << endl <<  "Error, to construct a 2D neighbors matrix, both dimensions must be greater than 1" << endl;
                return 1;
        }
        
        nlinks =  9* nsites - 6*(l+L) +4;
        sites = (int*)malloc(nlinks*sizeof(int));
        int iL=0;
        int xS=0;
        int yS=0;
        int iS = xS*l + yS;
        
        index[iS] = iL;
        number[iS] = 4;
        sites[iL] = iS;
        sites[iL+1] = iS+1;
        sites[iL+2] = iS+l;
        sites[iL+3] = iS+l+1;
        iL += number[iS];
        
        for(int xS=1; xS< L-1; xS++){
                iS = xS*l + yS;
                //             cout << "iS = " << iS << endl;
                index[iS] = iL;
                number[iS] = 6;
                sites[iL] = iS;
                sites[iL+1] = iS+1;
                sites[iL+2] = iS+l;
                sites[iL+3] = iS-l;
                sites[iL+4] = iS+l+1;
                sites[iL+5] = iS-l+1;
                iL += number[iS];
        }
        
        xS=L-1;
        iS = xS*l + yS;
        index[iS] = iL;
        number[iS] = 4;
        sites[iL] = iS;
        sites[iL+1] = iS+1;
        sites[iL+2] = iS-l;
        sites[iL+3] = iS-l+1;
        iL += number[iS];
        
        for(int i=1; i<l-1; i++){            
                for(int j=0; j<L; j++){
                        
                        xS = j;
                        yS = i;
                        iS = xS*l + yS;
                        //                 cout << "xS = " << xS << " yS = " << yS << " iS = " << iS << endl;
                        
                        if(j == 0){
                                index[iS] = iL;
                                number[iS] = 6;
                                sites[iL] = iS;
                                sites[iL+1] = iS+1;
                                sites[iL+2] = iS-1;
                                sites[iL+3] = iS+l;
                                sites[iL+4] = iS+l+1;
                                sites[iL+5] = iS+l-1;
                                iL += number[iS];
                        }  
                        else if(j == L-1){
                                index[iS] = iL;
                                number[iS] = 6;
                                sites[iL] = iS;
                                sites[iL+1] = iS+1;
                                sites[iL+2] = iS-1;
                                sites[iL+3] = iS-l;
                                sites[iL+4] = iS-l+1;
                                sites[iL+5] = iS-l-1;
                                iL += number[iS];
                        }
                        else{
                                index[iS]=iL;
                                number[iS]=9;
                                sites[iL] = iS;
                                sites[iL] = iS;
                                sites[iL+1] = iS+1;
                                sites[iL+2] = iS-1;
                                sites[iL+3] = iS+l;
                                sites[iL+4] = iS+l+1;
                                sites[iL+5] = iS+l-1;
                                sites[iL+6] = iS-l;
                                sites[iL+7] = iS-l+1;
                                sites[iL+8] = iS-l-1;
                                iL += number[iS];
                                
                        }               
                        
                }
        }
        
        xS=0;
        yS=l-1;
        iS=xS*l+yS;
        
        index[iS] = iL;
        number[iS] = 4;
        sites[iL] = iS;
        sites[iL+1] = iS-1;
        sites[iL+2] = iS+l;
        sites[iL+3] = iS+l-1;
        iL += number[iS];
        
        for(int xS=1; xS< L-1; xS++){
                iS = xS*l + yS;
                index[iS] = iL;
                number[iS] = 6;
                sites[iL] = iS;
                sites[iL+1] = iS-1;
                sites[iL+2] = iS+l;
                sites[iL+3] = iS+l-1;
                sites[iL+4] = iS-l;
                sites[iL+5] = iS-l-1;
                iL += number[iS];
        }
        
        xS = L-1;
        iS = xS*l + yS;
        index[iS] = iL;
        number[iS] = 4;
        sites[iL] = iS;
        sites[iL+1] = iS-1;
        sites[iL+2] = iS-l;
        sites[iL+3] = iS-l-1;
        iL += number[iS];
        
        if(iL!=nlinks){
                cout << "Nlinks=" << nlinks << " iL=" << iL << endl;
                return 1;
        }
        
        return 0;
}

int Spatial_grid::construction_2D_rectangle_cross(int L, int l){
        
        Init(L, l);
        
        if(l<2 || L<2){
                cout << endl <<  "Error, to construct a 2D neighbors matrix, both dimensions must be greater than 1" << endl;
                return 1;
        }
        
        nlinks =  5*nsites - 2*(l+L);
        sites = (int*)malloc(nlinks*sizeof(int));
        int iL=0;
        int xS=0;
        int yS=0;
        int iS = xS*l + yS;
        
        index[iS] = iL;
        number[iS] = 3;
        sites[iL] = iS;
        sites[iL+1] = iS+1;
        sites[iL+2] = iS+l;
        iL += number[iS];
        
        for(int xS=1; xS< L-1; xS++){
                iS = xS*l + yS;
                index[iS] = iL;
                number[iS] = 4;
                sites[iL] = iS;
                sites[iL+1] = iS+1;
                sites[iL+2] = iS+l;
                sites[iL+3] = iS-l;
                iL += number[iS];
        }
        
        xS=L-1;
        iS = xS*l + yS;
        index[iS] = iL;
        number[iS] = 3;
        sites[iL] = iS;
        sites[iL+1] = iS+1;
        sites[iL+2] = iS-l;
        iL += number[iS];
        
        for(int yS=1; yS<l-1; yS++){            
                for(int xS=0; xS<L; xS++){
                        
                        iS = xS*l + yS;
                        
                        if(xS == 0){
                                index[iS] = iL;
                                number[iS] = 4;
                                sites[iL] = iS;
                                sites[iL+1] = iS+1;
                                sites[iL+2] = iS-1;
                                sites[iL+3] = iS+l;
                                iL += number[iS];
                        }  
                        else if(xS == L-1){
                                index[iS] = iL;
                                number[iS] = 4;
                                sites[iL] = iS;
                                sites[iL+1] = iS+1;
                                sites[iL+2] = iS-1;
                                sites[iL+3] = iS-l;
                                iL += number[iS];
                        }
                        else{
                                index[iS]=iL;
                                number[iS]=5;
                                sites[iL] = iS;
                                sites[iL] = iS;
                                sites[iL+1] = iS+1;
                                sites[iL+2] = iS-1;
                                sites[iL+3] = iS+l;
                                sites[iL+4] = iS-l;
                                iL += number[iS];
                                
                        }               
                        
                }
        }
        
        xS=0;
        yS=l-1;
        iS=xS*l+yS;
        
        index[iS] = iL;
        number[iS] = 3;
        sites[iL] = iS;
        sites[iL+1] = iS-1;
        sites[iL+2] = iS+l;
        iL += number[iS];
        
        for(int xS=1; xS< L-1; xS++){
                iS = xS*l + yS;
                index[iS] = iL;
                number[iS] = 4;
                sites[iL] = iS;
                sites[iL+1] = iS-1;
                sites[iL+2] = iS+l;
                sites[iL+3] = iS-l;
                iL += number[iS];
        }
        
        xS = L-1;
        iS = xS*l + yS;
        index[iS] = iL;
        number[iS] = 3;
        sites[iL] = iS;
        sites[iL+1] = iS-1;
        sites[iL+2] = iS-l;
        iL += number[iS];
        
        if(iL!=nlinks){
                cout << "Nlinks=" << nlinks << " iL=" << iL << endl;
                return 1;
        }
        
        return 0;
}

void Spatial_grid::Print(){
        cout << endl << "Site \t\t Spatial_grid " << endl;
        for(int i=0; i<thickness; i++){
                for(int j=0; j<length; j++){
                        cout << "i = " << i << "  j = " << j  << "   site " << j*thickness+i << "  number = " << number[j*thickness+i] << endl ;
                        for(int k=0; k < number[j*thickness+i]; k++){
                                cout << sites[index[j*thickness +i] + k] << " \t" ;
                        }
                        cout << endl;
                }
                cout << endl;
        }
}


//Functions associated to the structure Gradient
Gradient::Gradient(int ngr, int nsi){
        Init(ngr, nsi);
}

Gradient::Gradient(int ngr, int nsi, double *data){
        Init(ngr, nsi);
        Fill(data);
}

Gradient::~Gradient(){
        free(values);
}

void Gradient::Init(int ngr, int nsi){
        ngrad = ngr;
        nsites= nsi;
        values = (double*)calloc(nsites*ngrad, sizeof(double));
}

void Gradient::Fill(double *data){
        for(int i=0; i<nsites*ngrad; i++){
                values[i] = data[i];
        }
}

void Gradient::Print(){
        cout << "Sites\t" ;
        for(int i=0; i<nsites; i++){
                cout << i << "\t" ;
        }
        cout << endl;
        for(int i=0; i<ngrad;i++){
                cout << "Gene " << i+1 << "\t";
                for(int j=0; j<nsites; j++){
                        cout << setprecision(3) << values[j*ngrad+i] << "\t";
                }
                cout << endl;
        }
        cout << endl;
} 

void Gradient::Construct_simple_gradient(int nsites){
        Init(1,nsites);
        for(int i=0; i<nsites; i++){
                values[i] = 1- (double)i/((double)nsites-1);
        }
}

//Functions associated to the structure Parameters
Parameters::Parameters(int nsp, int nsi, int ngr, double T){
        Init(nsp, nsi, ngr,T);
}

Parameters::~Parameters(){}

void Parameters::Init(int nsp, int nsi, int ngr, double T){
        nspins= nsp;
        nsites= nsi;
        ngrad= ngr;
        temperature = T;
        neighbors.Init(1, nsites);
        network.Init(nspins, ngrad);
        gradient.Init(1,nsites);
}

void Parameters::Fill_gradient( double *H){
        gradient.Fill(H);
}

void Parameters::Fill_param(double *par){
        network.Fill(par);
}

void Parameters::Print_conditions(){
        cout << endl << "Parameters : " << endl << "Temperature = " << temperature  << endl << "N genes = " << nspins << endl << "N sites = " << nsites << endl << "N maternals = " << ngrad << endl;
        network.Print();
        cout << endl;
        
}

// Functions associated to the structure Spins
Spins::Spins(int nsp, int nsi){
        Init(nsp, nsi);
}

Spins::Spins(int nsp, int nsi, double dauto, double dneigh){
        Init(nsp,nsi);
        Fill_diffusion(dauto, dneigh);
}

Spins::Spins(int nsp, int nsites, double *data){
        Init(nspins, nsites);
        Fill(data);
}

Spins::~Spins(){
        free(state);
}

void Spins::Init(int nsp, int nsi){
        nspins = nsp;
        nsites = nsi;
        state = (double*)calloc(nspins*nsites, sizeof(double));
        diff_auto = 1;
        diff_neigh = 1;
}

void Spins::Fill_diffusion(double da, double dn){
        diff_auto = da;
        diff_neigh = dn;
}

void Spins::Fill(double *data){
        for(int i =0; i< nsites*nspins; i++){
                state[i] = data[i];
        }
}

void Spins::Fill_rand(){
        srand(time(0));
        for(int i=0; i < nsites*nspins; i++){
                state[i] = rand()%2;
        }
}

int Spins::Add(Spins &spin_2){
        if(spin_2.nspins != nspins || spin_2.nsites != nsites){
                //             mexErrMsgTxt("Error, adding spin matrix of different dimensions");
                return -1;
        }
        for(int i=0; i< nsites*nspins; i++){
                state[i] += spin_2.state[i];
        }            
        return 0;
}

void Spins::Print(){
        cout << "Sites\t" ;
        for(int i=0; i<nsites; i++){
                cout << i << "\t" ;
        }
        cout << endl;
        for(int i=0; i<nspins;i++){
                cout << "Gene " << i+1 << "\t";
                for(int j=0; j<nsites; j++){
                        cout << setprecision(2) << state[j*nspins+i] << "\t";
                }
                cout << endl;
        }
        cout << endl;
}   

void Spins::Switch_one(double *heff, double T, int nuc){
        double proba[nspins];
        for(int i=0; i< nspins;i++){
                double ex = exp(-heff[i]/T);
                proba[i] = 1/(1+ex);
                cout << "proba = " << proba[i] << "   " ;
                double a =(double)(rand()%1000)/(double)1000;
                cout << "a= " << a << endl;
                if(a< proba[i]){
                        state[i*nsites+nuc]=1;
                }
                else{
                        state[i*nsites+nuc]=0;
                }
        }
}            

void Spins::Switch_mean(double *heff, double T){
        for(int i=0; i<nsites*nspins;i++){
                double ex = exp(-heff[i]/T);
                state[i] = 1/(1+ex);
        }
}  

void Spins::Calculate_heff_1nuc(double *heff, Parameters &system, int nuc){
        double CH[nsites*nspins]= {0};
        for(int j = 0; j<nspins; j++){
                heff[j]=0;
                for(int k = 0; k< system.ngrad; k++){
                        heff[j] += system.gradient.values[nuc*system.ngrad+k] * system.network.C[k*nspins+j];
                }
        }
        
        double *spin_neigh;
        spin_neigh=(double*)calloc(nspins, sizeof(double));
        for(int j=0; j< system.neighbors.number[nuc]; j++){
                if(system.neighbors.sites[ system.neighbors.index[nuc] +j] == nuc){
                        for(int k=0; k< nspins; k++){
                                spin_neigh[k] += diff_auto*state[nuc*nspins+k];
                        }
                }
                else{
                        for(int k=0; k< nspins; k++){
                                spin_neigh[k] += diff_neigh*state[system.neighbors.sites[ system.neighbors.index[nuc] +j]*nspins+k];
                        }
                }
        }
        
        for(int j=0; j<nspins;j++){
                for(int k=0; k<nspins; k++){
                        heff[j] += system.network.J[k*nspins+j]*spin_neigh[k];
                }
        }
        free(spin_neigh);   
}

void Spins::Calculate_heff(double *heff, Parameters &system){
        double CH[nsites*nspins]= {0};
        for(int i = 0; i<nsites; i++){
                for(int j = 0; j<nspins; j++){
                        heff[i*nspins+j]=0;
                        for(int k = 0; k< system.ngrad; k++){
                                heff[i*nspins+j] += system.gradient.values[i*system.ngrad+k] * system.network.C[k*nspins+j];
                        }
                }
        }
        
        double *spin_neigh;
        spin_neigh=(double*)calloc(nsites*nspins, sizeof(double));
        for(int i=0;i<nsites; i++){
                for(int j=0; j< system.neighbors.number[i]; j++){
                        if(system.neighbors.sites[ system.neighbors.index[i] +j] == i){
                                for(int k=0; k< nspins; k++){
                                        spin_neigh[i*nspins +k] += diff_auto*state[i*nspins+k];
                                }
                        }
                        else{
                                for(int k=0; k< nspins; k++){
                                        spin_neigh[i*nspins +k] += diff_neigh*state[system.neighbors.sites[ system.neighbors.index[i] +j]*nspins+k];
                                }
                        }
                }
                
                for(int j=0; j<nspins;j++){
                        for(int k=0; k<nspins; k++){
                                heff[i*nspins+j] += system.network.J[k*nspins+j]*spin_neigh[i*nspins+k];
                        }
                }
        }     
        free(spin_neigh);
}

int Spins::Test_stability(Spins &spin_mean, double critere){
        if(spin_mean.nspins != nspins || spin_mean.nsites != nsites){
                //             /*mex*/ErrMsgTxt("Error, adding spin matrix of different dimensions");
        }
        else{
                double test_osc=0;
                for(int i=0; i< nsites*nspins; i++){
                        test_osc += fabs(spin_mean.state[i]- state[i]);
                }
                test_osc /= nspins*nsites;
                if(test_osc < critere)
                        return 1;
                else
                        return 0;
        }
}

double Spins::Quality_max(){
        double Q=0;
        for(int i=0; i<nsites ; i++)
        {
                for(int j=0; j<nsites; j++)
                {
                        double diff = fabs(state[i*nspins+0]-state[j*nspins+0]);
                        for(int k=1; k<nspins ; k++)
                        {
                                double diff2 = fabs(state[i*nspins+k]-state[j*nspins+k]);
                                if (diff2>diff)
                                        diff = diff2;
                        }
                        Q += diff;
                }
        }
        Q /= pow(nsites,2);
        return Q;        
}

void MFSAexp_asym(Spins &spin, Parameters &system, int n_iterations){
        
        double Tmin = system.temperature;
        double T=1000*Tmin;
        double *heff;
        heff = (double*)calloc(system.nspins*system.nsites, sizeof(double));
        double a = pow((1e-3),1./n_iterations); 
        
        Spins spin_mean(spin.nspins, spin.nsites);    
        for(int i=0; i<n_iterations; i++){
                T *= a;
                system.temperature =T;
                spin.Calculate_heff(heff, system);
                spin.Switch_mean(heff, T);
                
                if(i > 19*n_iterations/20){
                        spin_mean.Add(spin);
                }
        }
        
        for(int i=0; i<spin.nspins*spin.nsites; i++){
                spin_mean.state[i] *= 20/(double)n_iterations;
        }
        
        int test = spin.Test_stability(spin_mean, 0.02);
        if(test ==0){
                spin.Fill(spin_mean.state);
        }
        system.temperature=Tmin;
        
        free(heff);
}

void MCMC(Spins &spin, Parameters &system, int n_equilibrium, int n_recording, int step_recording){
        double T=system.temperature;
        double *heff;
        heff = (double*)calloc(system.nspins*system.nsites, sizeof(double));
        Spins spin_mean(spin.nspins, spin.nsites);
        int nuc=0;
        
        for(int i=0; i< n_equilibrium; i++){
                nuc = rand()%(spin.nsites);
                spin.Calculate_heff_1nuc(heff, system, nuc);
//                 Print(heff,1,spin.nspins);
                spin.Switch_one(heff, T, nuc);
        }
        
        for(int i=0; i<n_recording*step_recording; i++){
                nuc = rand()%(spin.nsites); 
                cout << nuc << endl;
                spin.Calculate_heff_1nuc(heff, system, nuc);
                spin.Switch_one(heff, T, nuc);
                if(i%(step_recording) ==0){
                        cout << "add spin mean " << endl;
                        spin_mean.Add(spin);
                }                
        }
        
        spin_mean.Print();
        
        
        for(int i=0; i<spin.nspins*spin.nsites; i++){
                spin_mean.state[i] = spin_mean.state[i]/(double)n_recording;
        }
        
        spin.Fill(spin_mean.state);
        
}

