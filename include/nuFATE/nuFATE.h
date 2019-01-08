#include <vector>
#include <string>
#include <math.h>
#include <cmath>
#include <memory>
#include <iostream>

#include "hdf5.h"
#include "hdf5_hl.h"
#include "H5Gpublic.h"
#include "H5Fpublic.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_integration.h>

#ifndef NUFATE_H
#define NUFATE_H

namespace nufate{

///\class Result
///\brief Simple class to hold the results
class Result {
  public :
   std::vector<double> eval;
   std::shared_ptr<double> evec;
   std::vector<double> ci;
   std::vector<double> energy_nodes_;
   std::vector<double> phi_0_;

   // for pybinding
   std::vector<double> get_eval() { return eval; }
   double              get_evec() { return *evec; }
   std::vector<double> get_ci()   { return ci; }
   std::vector<double> get_energy_nodes() { return energy_nodes_; }
   std::vector<double> get_phi_0() { return phi_0_; }

   // for debug
   void Print() {
      std::cout << "***** Result print *****"  << std::endl;
      std::cout << "ci ---"  << std::endl;
      for (unsigned int i=0; i<ci.size(); ++i) {
          std::cout << ci[i] << " " ;
      }
      std::cout << std::endl;
      std::cout << "energy_nodes ---"  << std::endl;
      for (unsigned int i=0; i<energy_nodes_.size(); ++i) {
          std::cout << energy_nodes_[i] << " ";
      }
      std::cout << std::endl;
   }
};

///\class nuFATE
///\brief nuFATE main class
class nuFATE {
  private:
    const double GF = 1.16e-5;
    const double hbarc = 1.97e-14;
    const double GW = 2.085;
    const double MW = 80.385e0;
    const double mmu = 0.106e0;
    const double me = 511.e-6;
    const double pi = M_PI;
    const double MZ = 91.18;
    const double s2t = 0.23;
    const double gL =  s2t-0.5;
    const double REarth = 6371.;
    const double gR = s2t;
  private:
    /// flavor index for input neutrino
    int newflavor_;
    /// gamma index for input neutrino (for power law flux)
    double newgamma_;
    /// pedestal power law index to make the input flux flat, default is 2.0
    double pedestal_index_;
    /// array of input neutrino flux in energy bins
    /// size of the energy bins must be same as the size of
    /// energy bins in cross section 
    std::vector<double> input_flux_;
    ///  name of h5 file for cross section
    std::string newh5_filename_;
  public:
    /// \brief Constructor
    /// @param flavor Position of the system.
    /// @param gamma spectral index of input flux.
    /// @param h5_filename name of hdf5 file containing the cross sections.
    /// @param include_secondaries if true secondaries are added to the flux propagation.
    nuFATE(int flv, double gamma, std::string h5_filename, bool include_secondaries);
    /// \brief Constructor
    /// @param flavor Position of the system.
    /// @param gamma spectral index of input flux.
    /// @param energy_nodes energy nodes in GeV.
    /// @param sigma_array total cross section at each energy node in cm^2.
    /// @param dsigma_dE square array of differential cross section in cm^2/GeV.
    /// @param include_secondaries if true secondaries are added to the flux propagation.
    nuFATE(int flv, double gamma, std::vector<double> energy_nodes, std::vector<double> sigma_array, std::vector<std::vector<double>> dsigma_dE, bool include_secondaries);
    /// \brief Eigensystem calculator
    Result getEigensystem();
    /// \brief Function to get Earth column density
    /// @param theta Zenith angle in radians.
    double getEarthColumnDensity(double theta);
    /// \brief Function to get flavor
    int getFlavor() const;
    /// \brief Function to get input spectral index
    double getGamma() const;
    /// \brief Function to get pedestal power law index
    double getPedestalIndex() const;
    /// \brief Function to get input filename
    std::string getFilename() const;
    /// \brief Function to get number of energy nodes
    double getNumNodes() const;
    /// \brief Function to toggle secondaries
    const std::vector<double> &getEnergyNodes() const { return energy_nodes_; }
    /// \brief Function to toggle secondaries
    void setAddSecondaries(bool opt) { add_secondary_term_ = opt;}
    /// \brief Function to set initial power law flux
    /// @param gamma gamma index of power law
    /// @param pedestal_index power law index that makes power law spectrum flat. For example, if gamma index is close to 2.0, set 2.0 for pedestal_index.
    void setInitialPowerLawFlux(double gamma, double pedestal_index);
    /// \brief Function to set pedestal power law for input flux
    /// @param pedestal_index power law index that makes input flux flat. If input flux is atmospheric flux, for example, set pedestal_index 3.7.
    void setPedestalIndex(double pedestal_index);
    /// \brief Function to set initial user-defined flux. 
    /// @param v input energy flux at each energy node, must be same size as NumNodes_.
    /// @param pedestal_index power law index that makes input flux flat. If input flux is atmospheric flux, for example, set pedestal_index 3.7.
    void setInitialFlux(const std::vector<double> &v, double pedestal_index);

  protected:
    void AddAdditionalTerms();
    void LoadCrossSectionFromHDF5();
    void SetCrossSectionsFromInput(std::vector<std::vector<double>> dsigma_dE);
    void set_glashow_total();
    void set_glashow_partial();
    void set_RHS_matrices(std::shared_ptr<double> RMatrix_, std::shared_ptr<double> dxs_array);
    static double rho_earth(double theta, void * p); //returns density for a given angle in radians along x distance
    double readDoubleAttribute(hid_t, std::string) const;
    unsigned int readUIntAttribute(hid_t, std::string) const;
    std::vector<double> logspace(double min,double max,unsigned int samples) const;
  private:
    void AllocateMemoryForMembers(unsigned int num_nodes);
    void SetEnergyBinWidths();
    std::string grpdiff_;
    std::string grptot_;
    hid_t file_id_;
    hid_t root_id_;
    unsigned int NumNodes_;
    unsigned int rsize_;
    double Emax_;
    double Emin_;
    int dxsdim_[2];
    std::vector<double> energy_nodes_;
    std::vector<double> sigma_array_;
    std::vector<double> DeltaE_;
    std::vector<double> phi_0_;
    std::vector<double> phi_0_pedestal_;
    std::vector<double> glashow_total_;
    std::vector<double> sig3_array_;
    std::shared_ptr<double> glashow_partial_;
    std::shared_ptr<double> dxs_array_;
    std::shared_ptr<double> tau_array_;
    std::shared_ptr<double> sec_array_;
    std::shared_ptr<double> regen_array_;
    std::shared_ptr<double> RHSMatrix_;
    std::shared_ptr<double> RHSMatrix1_;
    std::shared_ptr<double> RHSMatrix2_;
    std::shared_ptr<double> RHSMatrix3_;
    std::shared_ptr<double> RHSMatrix4_;
    std::shared_ptr<double> RHregen_;
    std::shared_ptr<double> Enuin_;
    std::shared_ptr<double> Enu_;
    std::shared_ptr<double> selectron_;
    std::shared_ptr<double> den_;
    std::shared_ptr<double> t1_;
    std::shared_ptr<double> t2_;
    std::shared_ptr<double> t3_;
  private:
    bool include_secondaries_ = false;
    bool memory_allocated_ = false;
    bool initial_flux_set_ = false;
    bool total_cross_section_set_ = false;
    bool differential_cross_section_set_ = false;
    bool RHS_set_ = false;
    bool add_tau_regeneration_ = true;
    bool add_glashow_term_= true;
    bool add_secondary_term_ = true;
};

}

#endif
