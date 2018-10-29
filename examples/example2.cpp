#include <iostream>
#include <iomanip>
#include <fstream>
#include "nuFATE.h"


int main(){

/*
 * flavor: Select a neutrino flavor (1, 2, 3, for electron, muon, and tau, respectively. the number is negative for anti-neutrinos)
 *  gamma: Spectral index of the incoming flux
 *  file: path to file containing the cross sections.
*/
    int flavor = 2;
    double gamma = 2;
    std::string fluxfile = "../resources/1PeV_flux.txt";
    bool include_secondaries = false;
    std::string file = "../resources/NuFATECrossSections.h5";
    //Initialize an instance of the nuFATE class with these parameters.
    nufate::nuFATE object(flavor, gamma, file, include_secondaries);

    std::ifstream fin;
    fin.open(fluxfile);
    std::vector<double> phi_in;
    double flux;
    while (fin >> flux) { 
        phi_in.push_back(flux);
    }
    object.SetInitialFlux(phi_in);

    //Result is a struct that stores the solution to the cascade equation.
    nufate::Result result;
    //get_eigs() solves the cascade equation
    result = object.getEigensystem();

    int NumNodes;
    NumNodes = object.getNumNodes();
    std::vector<double> eval = result.eval;
    std::shared_ptr<double> evec = result.evec;
    std::vector<double> ci = result.ci;
    std::vector<double> energy_nodes = result.energy_nodes_;
    std::vector<double> phi_0 = result.phi_0_;
    //Calculate earth column density for a given zenith
    double Na = 6.0221415e23;
    //double zenith = 2.2689280276;
    double zenith = 3.1415;
    double t;
    t = object.getEarthColumnDensity(zenith) * Na;
    //t = 6.59272246299e33 * Na;
    std::cerr << "t is " << t << std::endl;

    std::vector<double> abs;
    std::vector<double> phi_sol;
    //Get Attenuation!
    if(not include_secondaries){
      //Without Secondaries
      for(int i=0; i<NumNodes; i++){
        double sum = 0.;
        for (int j=0; j<NumNodes;j++){
          abs.push_back(ci[j] * exp(-t*eval[j]));
          sum+= abs[j] *  *(evec.get()+i*NumNodes+j) * (std::pow(energy_nodes[i],-2)); // arrival flux
        }
        phi_sol.push_back(sum);
      }
      //Print Solution
      std::cout << "Solution = " << std::endl;
      std::ofstream fout;
      fout.open("example_out2_cxx.txt", std::ios::out);
      for(int i =0; i<NumNodes; i++){
        std::cout << std::fixed << std::setprecision(6) << log10(energy_nodes[i]) << ", " << std::scientific << std::setprecision(6) << phi_sol[i] << std::endl;
        fout << std::fixed << std::setprecision(6) << log10(energy_nodes[i]) << ", " << std::scientific << std::setprecision(6) << phi_sol[i] << std::endl;
      }
      fout.close();

    } else{
      //With Secondaries
        int rsize = 2*NumNodes;
        for(int i=0; i<rsize; i++){
          double sum = 0.;
          abs.clear();
          for (int j=0; j<rsize;j++){
            abs.push_back(ci[j] * exp(-t*eval[j]));
            sum+= (abs[j] * *(result.evec.get()+i*rsize+j)) / phi_0[i];
          }
        phi_sol.push_back(sum);
        }

        //Print Solution
        std::cout << "Solution Including secondaries= " << std::endl;
        for(int i =0; i<NumNodes; i++){
          std::cout << phi_sol[i] << std::endl;
        }
    }

    return 0;
}
