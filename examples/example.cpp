#include <iostream>
#include "nuFATE.h"


int main(){

/*
 * flavor: Select a neutrino flavor (1, 2, 3, for electron, muon, and tau, respectively. the number is negative for anti-neutrinos)
 *  gamma: Spectral index of the incoming flux
 *  file: path to file containing the cross sections.
*/
    int flavor = -2;
    double gamma = 1.2;
    double pedestal_index = 2.0;
    bool include_secondaries = false;
    std::string file = "../resources/NuFATECrossSections.h5";

    //Initialize an instance of the nuFATE class with these parameters.
    nufate::nuFATE object(flavor, gamma, file, include_secondaries);

    int NumNodes;
    NumNodes = object.getNumNodes();

    //
    // test 1 : replace input flux
    //
    gamma = 3.65;
    pedestal_index = 3.7;
    // prepare input flux
    std::vector<double> initial_flux;
    const std::vector<double> &energies = object.getEnergyNodes();
    for (int i=0; i<NumNodes; ++i) {
        initial_flux.push_back(std::pow(energies[i], -gamma));
    }
    object.setInitialFlux(initial_flux, pedestal_index);
    //

    //Result is a struct that stores the solution to the cascade equation.
    nufate::Result result;
    //get_eigs() solves the cascade equation
    result = object.getEigensystem();

    std::vector<double> eval = result.eval;
    std::shared_ptr<double> evec = result.evec;
    std::vector<double> ci = result.ci;
    std::vector<double> energy_nodes = result.energy_nodes_;
    std::vector<double> phi_0 = result.phi_0_;
    //Calculate earth column density for a given zenith
    double Na = 6.0221415e23;
    double zenith = 120.*3.1415 / 180.;
    double t;
    t = object.getEarthColumnDensity(zenith) * Na;

    std::vector<double> abs;
    std::vector<double> phi_sol;
    //Get Attenuation!
    if(not include_secondaries){
      //Without Secondaries
      for(int i=0; i<NumNodes; i++){
        double sum = 0.;
        for (int j=0; j<NumNodes;j++){
          abs.push_back(ci[j] * exp(-t*eval[j]));
          sum+= abs[j] *  *(evec.get()+i*NumNodes+j) * (std::pow(energy_nodes[i],-pedestal_index) / std::pow(energy_nodes[i],-gamma));
        }
        phi_sol.push_back(sum);
      }
      //Print Solution
      std::cout << "Solution = " << std::endl;
      for(int i =0; i<NumNodes; i++){
        std::cout << log10(energy_nodes[i]) << ", " <<  phi_sol[i] << std::endl;
      }

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
          std::cout << log10(energy_nodes[i]) << ", " <<  phi_sol[i] << std::endl;
        }
    }

    return 0;
}
