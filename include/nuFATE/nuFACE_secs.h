#include <vector>
#include <string>
#include <math.h>
#include <cmath>
#include <memory>

#include "hdf5.h"
#include "hdf5_hl.h"
#include "H5Gpublic.h"
#include "H5Fpublic.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_integration.h>

#ifndef NUFACE_SECS_H
#define NUFACE_SECS_H

namespace nufate {

struct Result {
  double* eval;
  double* evec;
  double* ci;
  std::vector<double> energy_nodes_;
  std::vector<double> phi_0_;
};

class nuFACE_secs {
  private:
    const double GF = 1.16e-5;
    const double hbarc = 1.97e-14;
    const double GW = 2.085;
    const double MW = 80.385e0;
    const double mmu = 0.106e0;
    const double me = 511.e-6;
    const double pi = 3.14159265358979323846;
    const double MZ = 91.18;
    const double REarth = 6371.;
    const double s2t = 0.23;
    const double gL =  s2t-0.5;
    const double gR = s2t;
  private:
    int newflavor_;
    double newgamma_;
    std::string newh5_filename_;
  public:
    nuFACE_secs(int, double, std::string);
    Result get_eigs();
    double get_t_earth(double theta); //calculates column density for a given angle in radians
    int getFlavor() const;
    double getGamma() const;
    std::string getFilename() const;
  protected:
    void set_glashow_total();
    void set_glashow_partial();
    void set_RHS_matrices();
    unsigned int readUIntAttribute(hid_t, std::string) const;
    std::vector<double> logspace(double min,double max,unsigned int samples) const;
    static double rho_earth(double theta, void * p); //returns density for a given angle in radians along x distance
    double readDoubleAttribute(hid_t, std::string) const;
  private:
    hid_t file_id_;
    hid_t root_id_;
    std::string grpdiff_;
    std::string grptot_;
    unsigned int NumNodes_;
    unsigned int rsize_;
    double Emax_;
    double Emin_;
    std::vector<double> energy_nodes_;
    std::vector<double> sigma_array_;
    std::vector<double> phi_0_;
    std::vector<double> DeltaE_;
    std::vector<double> glashow_total_;
    std::vector<double> sig3_array_;
    std::shared_ptr<double> dxs_array_;
    std::shared_ptr<double> sec_array_;
    std::shared_ptr<double> regen_array_;
    std::shared_ptr<double> RHSMatrix_;
    std::shared_ptr<double> RHSMatrix1_;
    std::shared_ptr<double> RHSMatrix2_;
    std::shared_ptr<double> RHSMatrix3_;
    std::shared_ptr<double> RHSMatrix4_;
    std::shared_ptr<double> glashow_partial_;
    std::shared_ptr<double> Enuin_;
    std::shared_ptr<double> Enu_;
    std::shared_ptr<double> selectron_;
    std::shared_ptr<double> den_;
    std::shared_ptr<double> t1_;
    std::shared_ptr<double> t2_;
    std::shared_ptr<double> t3_;
    int dxsdim_[2];
};

} // close namespace

#endif