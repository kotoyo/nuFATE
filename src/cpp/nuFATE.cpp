#include "nuFATE.h"
#include <iostream>

namespace nufate{

static double UNPHYSICAL_NUMBER = 1e-30;

nuFATE::nuFATE(int flavor, double gamma, std::string h5_filename, bool include_secondaries) : newflavor_(flavor), newgamma_(gamma), newh5_filename_(h5_filename), include_secondaries_(include_secondaries) {
    //A few sanity checks
    if(include_secondaries_ and (newflavor_ == 3 or newflavor_== -3))
      throw std::runtime_error("nuFATE::nuFATE Cannot Include secondaries for tau's.");
    if(not (newflavor_ == 3 or newflavor_ == -3 or newflavor_ == 2 or newflavor_ == -2 or newflavor_ == 1 or newflavor_ == -1))
      throw std::runtime_error("nuFATE::nuFATE flavor has to be plus or minus 1,2,3.");
    //open h5file containing cross sections (xsh5)
    file_id_ = H5Fopen(h5_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    root_id_ = H5Gopen(file_id_, "/", H5P_DEFAULT);
    grptot_ = "/total_cross_sections";
    grpdiff_ = "/differential_cross_sections";
    hid_t group_id = H5Gopen(root_id_, grptot_.c_str(), H5P_DEFAULT);
    //assign some important variables
    Emax_ = readDoubleAttribute(group_id, "max_energy");
    Emin_ = readDoubleAttribute(group_id, "min_energy");
    NumNodes_ = readUIntAttribute(group_id, "number_energy_nodes");
    // allocate memory
    AllocateMemoryForMembers(NumNodes_);
    //set energy nodes and deltaE
    energy_nodes_ = logspace(Emin_, Emax_, NumNodes_);
    // calculate and set energy bin widths
    SetEnergyBinWidths();
    // load cross sections from file
    LoadCrossSectionFromHDF5();
    // set the initial flux
    setInitialPowerLawFlux(newgamma_);
    // allocate gsl buffers for evaluating eigenvalues
    allocate_gsl_buffers();
}

nuFATE::nuFATE(int flavor, double gamma, std::vector<double> energy_nodes, std::vector<double> sigma_array, std::vector<std::vector<double>> dsigma_dE, bool include_secondaries):
  newflavor_(flavor), newgamma_(gamma), energy_nodes_(energy_nodes), sigma_array_(sigma_array), include_secondaries_(include_secondaries)
{
  NumNodes_ = energy_nodes_.size();
  Emax_ = energy_nodes_.back();
  Emin_ = energy_nodes_.front();
  if(sigma_array.size() != NumNodes_)
    throw std::runtime_error("nuFATE::nuFATE Total cross section array does not match energy nodes size.");
  if(dsigma_dE.size() != NumNodes_ or dsigma_dE.front().size() != NumNodes_)
    throw std::runtime_error("nuFATE::nuFATE Differential cross section array does not match energy nodes size.");
  AllocateMemoryForMembers(NumNodes_);
  SetEnergyBinWidths();
  setInitialPowerLawFlux(newgamma_);
  SetCrossSectionsFromInput(dsigma_dE);
  // allocate gsl buffers for evaluating eigenvalues
  allocate_gsl_buffers();
}

void nuFATE::SetCrossSectionsFromInput(std::vector<std::vector<double>> dsigma_dE){
    for(unsigned int i = 0; i<NumNodes_; i++){
        for(unsigned int j=0; j<NumNodes_; j++)
        *(dxs_array_.get()+i*NumNodes_+j) = dsigma_dE[i][j];
    }
    total_cross_section_set_ = true;
    differential_cross_section_set_ = true;
    // RHS matrices depends on cross sections and now
    // need to be recalculated.
    RHS_set_ = false;
}

void nuFATE::AllocateMemoryForMembers(unsigned int NumNodes){
  // if the energy nodes are different or if memory has not been allocated. Allocate it.
  if(NumNodes_ != NumNodes or (not memory_allocated_)){
    NumNodes_ = NumNodes;
    //allocate memory that will be used in functions below
    if(include_secondaries_){
       RHSMatrix_ = std::shared_ptr<double>((double *)malloc((2*NumNodes_)*(2*NumNodes_)*sizeof(double)),free);
       RHSMatrix1_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
       RHSMatrix2_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
       RHSMatrix3_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
       RHSMatrix4_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
       regen_array_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
       sec_array_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
       sig3_array_ = std::vector<double>(NumNodes_);
    }
    else{
       RHSMatrix_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    }

  glashow_total_ = std::vector<double>(NumNodes_);
  sigma_array_ = std::vector<double>(NumNodes_);
  glashow_partial_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
  Enuin_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
  Enu_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
  den_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
  selectron_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
  t1_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
  t2_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
  t3_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
  dxs_array_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);

  }
  memory_allocated_ = true;
}

void nuFATE::SetEnergyBinWidths(){
  if(energy_nodes_.size() == 0)
    throw std::runtime_error("nuFATE::SetEnergyBinWidths Energy nodes need to be set before the widths can be set.");
  DeltaE_ = std::vector<double>(NumNodes_);
  for(unsigned int i = 0; i < NumNodes_-1;i++){
      DeltaE_[i] = log(energy_nodes_[i+1]) - log(energy_nodes_[i]);
  }
}

//reads an attribute of type double from h5 object
double nuFATE::readDoubleAttribute(hid_t object, std::string name) const{
    double target;
    hid_t attribute_id = H5Aopen(object,name.c_str(),H5P_DEFAULT);
    herr_t status = H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &target);
    if(status<0)
        throw std::runtime_error("Failed to read attribute '"+name+"'");
    H5Aclose(attribute_id);
    return target;
}

//reads an attribute of type unsigned int from h5 object
unsigned int nuFATE::readUIntAttribute(hid_t object, std::string name) const{
    unsigned int target;
    hid_t attribute_id = H5Aopen(object,name.c_str(),H5P_DEFAULT);
    herr_t status = H5Aread(attribute_id, H5T_NATIVE_UINT, &target);
    if(status<0)
        throw std::runtime_error("Failed to read attribute '"+name+"'");
    H5Aclose(attribute_id);
    return target;
}

//given a minimum, maximum, and number of points, this function returns equally spaced array in logspace
std::vector<double> nuFATE::logspace(double Emin,double Emax,unsigned int div) const {
    if(div==0)
        throw std::length_error("number of samples requested from logspace must be nonzero");
    std::vector<double> logpoints(div);
    double Emin_log,Emax_log;
    Emin_log = log10(Emin);
    Emax_log = log10(Emax);
    double step_log = (Emax_log - Emin_log)/double(div-1);
    logpoints[0]=Emin;
    double EE = Emin_log+step_log;
    for(unsigned int i=1; i<div-1; i++, EE+=step_log)
        logpoints[i] = std::pow(10,EE);
    logpoints[div-1]=Emax;
    return logpoints;
}

//sets the contribution from glashow
void nuFATE::set_glashow_total(){
    for(unsigned int i=0; i<NumNodes_; i++){
        glashow_total_[i] = 2.*me*energy_nodes_[i];
        double x = glashow_total_[i];
        glashow_total_[i] = 1. /3.*std::pow(GF,2)*x/pi*std::pow((1.-(std::pow(mmu,2)-std::pow(me,2))/x),2)/(std::pow((1.-x/std::pow(MW,2)),2)+std::pow(GW,2)/std::pow(MW,2))*0.676/0.1057*std::pow(hbarc,2);
    }
    return;
}

void nuFATE::set_glashow_partial(){
    for (unsigned int i =0; i<NumNodes_;i++){
        for(unsigned int j = 0; j<NumNodes_;j++){
            *(Enuin_.get()+i*NumNodes_+j) = energy_nodes_[j];
            *(Enu_.get()+i*NumNodes_+j) = energy_nodes_[i];
            *(Enu_.get()+i*NumNodes_+j) = 1 - *(Enu_.get()+i*NumNodes_+j)/ *(Enuin_.get()+i*NumNodes_+j);
            *(selectron_.get()+i*NumNodes_+j) = 2.*me* *(Enuin_.get()+i*NumNodes_+j);
            *(den_.get()+i*NumNodes_+j) = std::pow(1. - *(selectron_.get()+i*NumNodes_+j)/std::pow(MW,2),2) + std::pow(GW,2)/std::pow(MW,2);
            *(t1_.get()+i*NumNodes_+j) = std::pow(gR,2)/std::pow((1.+ *(Enu_.get()+i*NumNodes_+j)* *(selectron_.get()+i*NumNodes_+j)/std::pow(MZ,2)),2);
            *(t2_.get()+i*NumNodes_+j) = gL/(1.+ *(Enu_.get()+i*NumNodes_+j)* *(selectron_.get()+i*NumNodes_+j)/std::pow(MZ,2)) + (1. - *(selectron_.get()+i*NumNodes_+j)/std::pow(MW,2))/ *(den_.get()+i*NumNodes_+j);
            *(t3_.get()+i*NumNodes_+j) = GW/MW/ *(den_.get()+i*NumNodes_+j);
            if (*(Enu_.get()+i*NumNodes_+j) >= 0.){
                *(glashow_partial_.get()+i*NumNodes_+j) = (std::pow(GF,2)* *(selectron_.get()+i*NumNodes_+j)/pi*(*(t1_.get()+i*NumNodes_+j)+ (std::pow(*(t2_.get()+i*NumNodes_+j) , 2) + std::pow(*(t3_.get()+i*NumNodes_+j),2)) * std::pow((1.-*(Enu_.get()+i*NumNodes_+j)),2))*std::pow(hbarc,2))/ *(Enuin_.get()+i*NumNodes_+j);
            } else {
                *(glashow_partial_.get()+i*NumNodes_+j) = 0.;
            }
        }
    }
    return;
}

void nuFATE::set_RHS_matrices(std::shared_ptr<double> RMatrix, std::shared_ptr<double> dxsarray) {

    // make sure scaling_flux is set
    setScalingFlux(scaling_index_);

    if(include_secondaries_){
      for(unsigned int i=0; i<NumNodes_; i++){
        for(unsigned int j=0; j<NumNodes_; j++){
          *(RHSMatrix1_.get()+i*NumNodes_+j) = 0.;
          *(RHSMatrix2_.get()+i*NumNodes_+j) = 0.;
          *(RHSMatrix3_.get()+i*NumNodes_+j) = 0.;
          *(RHSMatrix4_.get()+i*NumNodes_+j) = 0.;
        }
      }

      for (unsigned int i=0; i<NumNodes_; i++){
        for(unsigned int j=i+1; j<NumNodes_; j++){
          // nue or numu NC
          *(RHSMatrix1_.get()+i*NumNodes_+j) = DeltaE_[j-1] * *(dxs_array_.get()+j*NumNodes_+i) * std::pow(energy_nodes_[j],-1) * scaling_flux_[i];
          // tau regen + tau NC
          *(RHSMatrix4_.get()+i*NumNodes_+j) = DeltaE_[j-1] * (*(dxs_array_.get()+j*NumNodes_+i) + *(regen_array_.get()+j*NumNodes_+i)) * std::pow(energy_nodes_[j],-1) * scaling_flux_[i];
        }
      }

      for (unsigned int i=0; i<NumNodes_; i++){
        *(RHSMatrix1_.get()+i*NumNodes_+i) = *(RHSMatrix1_.get()+i*NumNodes_+i) - sigma_array_[i];
        *(RHSMatrix4_.get()+i*NumNodes_+i) = *(RHSMatrix4_.get()+i*NumNodes_+i) - sig3_array_[i];
        for(unsigned int j=i+1; j<NumNodes_; j++){
          // nue/mu production
          *(RHSMatrix2_.get()+i*NumNodes_+j) = DeltaE_[j-1] * *(sec_array_.get()+j*NumNodes_+i) * std::pow(energy_nodes_[j],-1) * scaling_flux_[i];
        }
      }

      rsize_ = 2*NumNodes_;
      for (unsigned int i =0; i<NumNodes_;i++){
        unsigned int x = i+NumNodes_;
        for(unsigned int j =0; j<NumNodes_;j++){
          unsigned int y = j+NumNodes_;
          *(RHSMatrix_.get()+i*rsize_+j) = *(RHSMatrix1_.get() + i *NumNodes_+j);
          *(RHSMatrix_.get()+x*rsize_+j) = *(RHSMatrix3_.get() + i *NumNodes_+j);
          *(RHSMatrix_.get() +i*rsize_+y) = *(RHSMatrix2_.get() + i *NumNodes_+j);
          *(RHSMatrix_.get() +x*rsize_+y) =  *(RHSMatrix4_.get() + i *NumNodes_+j);
        }
      }

    } else{
      for(unsigned int i = 0; i < NumNodes_; i++){
        for(unsigned int j= i+1; j < NumNodes_; j++){
          double e1 = 1./ energy_nodes_[j];
          double e2 = scaling_flux_[i];
          *(RMatrix.get()+i*NumNodes_+j) = DeltaE_[j - 1] * *(dxsarray.get()+j * dxsdim_[1]+i) * e1 * e2;
        }
      }
    }
    return;
}

void nuFATE::LoadCrossSectionFromHDF5(){
    hid_t group_id;
    group_id = H5Gopen(root_id_, grptot_.c_str(), H5P_DEFAULT);

    if (newflavor_ == -1) {
        hsize_t sarraysize[1];
        H5LTget_dataset_info(group_id,"nuebarxs", sarraysize,NULL,NULL);
        if(sarraysize[0] != NumNodes_)
          throw std::runtime_error("nuFATE::LoadCrossSectionFromHDF5 Total cross section array does not match number of energy nodes.");
        H5LTread_dataset_double(group_id, "nuebarxs", sigma_array_.data());
    }  else if (newflavor_ == -2){
        hsize_t sarraysize[1];
        H5LTget_dataset_info(group_id,"numubarxs", sarraysize,NULL,NULL);
        if(sarraysize[0] != NumNodes_)
          throw std::runtime_error("nuFATE::LoadCrossSectionFromHDF5 Total cross section array does not match number of energy nodes.");
        H5LTread_dataset_double(group_id, "numubarxs", sigma_array_.data());
    }  else if (newflavor_ == -3){
        hsize_t sarraysize[1];
        H5LTget_dataset_info(group_id,"nutaubarxs", sarraysize,NULL,NULL);
        if(sarraysize[0] != NumNodes_)
          throw std::runtime_error("nuFATE::LoadCrossSectionFromHDF5 Total cross section array does not match number of energy nodes.");
        H5LTread_dataset_double(group_id, "nutaubarxs", sigma_array_.data());
    }  else if (newflavor_ == 1){
        hsize_t sarraysize[1];
        H5LTget_dataset_info(group_id,"nuexs", sarraysize,NULL,NULL);
        if(sarraysize[0] != NumNodes_)
          throw std::runtime_error("nuFATE::LoadCrossSectionFromHDF5 Total cross section array does not match number of energy nodes.");
        H5LTread_dataset_double(group_id, "nuexs", sigma_array_.data());
    }  else if (newflavor_ == 2){
        hsize_t sarraysize[1];
        H5LTget_dataset_info(group_id,"numuxs", sarraysize,NULL,NULL);
        if(sarraysize[0] != NumNodes_)
          throw std::runtime_error("nuFATE::LoadCrossSectionFromHDF5 Total cross section array does not match number of energy nodes.");
        H5LTread_dataset_double(group_id, "numuxs", sigma_array_.data());
    }  else if (newflavor_ == 3){
        hsize_t sarraysize[1];
        H5LTget_dataset_info(group_id,"nutauxs", sarraysize,NULL,NULL);
        if(sarraysize[0] != NumNodes_)
          throw std::runtime_error("nuFATE::LoadCrossSectionFromHDF5 Total cross section array does not match number of energy nodes.");
        H5LTread_dataset_double(group_id, "nutauxs", sigma_array_.data());
    }

    hsize_t dxarraysize[2];
    group_id = H5Gopen(root_id_, grpdiff_.c_str(), H5P_DEFAULT);
   if (newflavor_ > 0){
        H5LTget_dataset_info(group_id,"dxsnu", dxarraysize,NULL,NULL);
        dxsdim_[0] = dxarraysize[0];
        dxsdim_[1] = dxarraysize[1];
        if((unsigned int)(dxsdim_[0]) != NumNodes_ or (unsigned int)(dxsdim_[1]) != NumNodes_)
          throw std::runtime_error("nuFATE::LoadCrossSectionFromHDF5 Differential arrays do not match number of energy nodes.");
        H5LTread_dataset_double(group_id, "dxsnu", dxs_array_.get());
    } else {
        H5LTget_dataset_info(group_id,"dxsnubar", dxarraysize,NULL,NULL);
        dxsdim_[0] = dxarraysize[0];
        dxsdim_[1] = dxarraysize[1];
        if((unsigned int)(dxsdim_[0]) != NumNodes_ or (unsigned int)(dxsdim_[1]) != NumNodes_)
          throw std::runtime_error("nuFATE::LoadCrossSectionFromHDF5 Differential arrays do not match number of energy nodes.");
        H5LTread_dataset_double(group_id, "dxsnubar", dxs_array_.get());
    }

   if(include_secondaries_){
       if (newflavor_ > 0){
          group_id = H5Gopen(root_id_, grptot_.c_str(), H5P_DEFAULT);
          hsize_t sarraysize[1];
          H5LTget_dataset_info(group_id,"nutauxs", sarraysize,NULL,NULL);
          if((unsigned int)(sarraysize[0]) != NumNodes_)
            throw std::runtime_error("nuFATE::LoadCrossSectionFromHDF5 Tau arrays do not match number of energy nodes.");
          H5LTread_dataset_double(group_id, "nutauxs", sig3_array_.data());
          std::string grptau = "/tau_decay_spectrum";
          group_id = H5Gopen(root_id_, grptau.c_str(), H5P_DEFAULT);
          hsize_t tauarraysize[2];
          H5LTget_dataset_info(group_id,"tfull", tauarraysize,NULL,NULL);
          if((unsigned int)(tauarraysize[0]) != NumNodes_ or (unsigned int)(tauarraysize[1]) != NumNodes_)
            throw std::runtime_error("nuFATE::LoadCrossSectionFromHDF5 Tau arrays do not match number of energy nodes.");
          H5LTread_dataset_double(group_id, "tfull", regen_array_.get());
          hsize_t secarraysize[2];
          H5LTget_dataset_info(group_id,"secfull", secarraysize,NULL,NULL);
          if((unsigned int)(secarraysize[0]) != NumNodes_ or (unsigned int)(secarraysize[1]) != NumNodes_)
            throw std::runtime_error("nuFATE::LoadCrossSectionFromHDF5 Secondary arrays do not match number of energy nodes.");
          H5LTread_dataset_double(group_id, "secfull", sec_array_.get());

      } else {
          group_id = H5Gopen(root_id_, grptot_.c_str(), H5P_DEFAULT);
          hsize_t sarraysize[1];
          H5LTget_dataset_info(group_id,"nutaubarxs", sarraysize,NULL,NULL);
          if((unsigned int)(sarraysize[0]) != NumNodes_)
            throw std::runtime_error("nuFATE::LoadCrossSectionFromHDF5 Tau arrays do not match number of energy nodes.");
          H5LTread_dataset_double(group_id, "nutaubarxs", sig3_array_.data());
          std::string grptau = "/tau_decay_spectrum";
          group_id = H5Gopen(root_id_, grptau.c_str(), H5P_DEFAULT);
          hsize_t tauarraysize[2];
          H5LTget_dataset_info(group_id,"tbarfull", tauarraysize,NULL,NULL);
          if((unsigned int)(tauarraysize[0]) != NumNodes_ or (unsigned int)(tauarraysize[1]) != NumNodes_)
            throw std::runtime_error("nuFATE::LoadCrossSectionFromHDF5 Tau arrays do not match number of energy nodes.");
          H5LTread_dataset_double(group_id, "tbarfull", regen_array_.get());
          hsize_t secarraysize[2];
          H5LTget_dataset_info(group_id,"secbarfull", secarraysize,NULL,NULL);
          if((unsigned int)(secarraysize[0]) != NumNodes_ or (unsigned int)(secarraysize[1]) != NumNodes_)
            throw std::runtime_error("nuFATE::LoadCrossSectionFromHDF5 Secondary arrays do not match number of energy nodes.");
          H5LTread_dataset_double(group_id, "secbarfull", sec_array_.get());
      }
   }

   total_cross_section_set_ = true;
   differential_cross_section_set_ = true;

   // RHS matrices depends on cross sections and now
   // need to be recalculated.
   RHS_set_ = false;
}

std::vector<double> nuFATE::getRelativeAttenuation(double number_of_targets) 
{
   // calculate eigensystem first
   getEigensystem();
   double t = number_of_targets;
   std::vector<double> phi_sol;

   unsigned int rsize = NumNodes_; 
   if (include_secondaries_) {
      rsize = 2*NumNodes_;
   }
   phi_sol.resize(rsize);

   double abs;
   for (unsigned int i=0; i<rsize; i++){
      double sum = 0.;
      for (unsigned int j=0; j<rsize; j++){
         abs = r1_.ci[j] * exp(-t * r1_.eval[j]);
         // phi_0_ = initial_flux * scaling_flux
         // arrval_flux = abs (dot) eigenvec / scaling_flux 
         // attenuation = arrival_flux / initial_flux = abs (dot) eigenvec / (scaling_flux*initial_flux) 
         sum += abs *  *((r1_.evec).get()+i*rsize+j);
      }
      phi_sol[i] = sum / r1_.phi_0_[i];
   }
   return phi_sol;
}

void nuFATE::setInitialPowerLawFlux(double gamma)
{
    if (newgamma_ == gamma && initial_flux_set_ == true) {
        // if fluxes are already calculated, just return.
        return;
    }
    
    initial_flux_set_ = false;
    newgamma_ = gamma;

    // calculate scaling flux first.
    setScalingFlux(scaling_index_);

    if(include_secondaries_){

        phi_0_ = std::vector<double>(2*NumNodes_);
        for (unsigned int i = 0; i < NumNodes_; i++){
            phi_0_[i] = std::pow(energy_nodes_[i],- newgamma_) * scaling_flux_[i] ;
            phi_0_[i+NumNodes_] = phi_0_[i];
        }

    } else {
        phi_0_ = std::vector<double>(NumNodes_);
        for (unsigned int i = 0; i < NumNodes_; i++){
          phi_0_[i] = std::pow(energy_nodes_[i],- newgamma_) * scaling_flux_[i];
        }
    }
    initial_flux_set_ = true;
}

void nuFATE::setScalingFlux(double scaling_index) {

    if (scaling_flux_set_ == true && scaling_index_ == scaling_index) {
        return;
    }

    scaling_flux_set_ = false;
    scaling_index_ = scaling_index;

    // Calculate scaling flux. 
    // That removes std::pow operation in setInitialFlux 
    // function, which may speed up calculation a little
    // when the setInitialFlux function is called in 
    // for_loop...
    if (include_secondaries_) {
        scaling_flux_ = std::vector<double>(2*NumNodes_);
        for (unsigned int i = 0; i < NumNodes_; i++){
            scaling_flux_[i] = std::pow(energy_nodes_[i],scaling_index_);
            scaling_flux_[i+NumNodes_] = scaling_flux_[i];
        }
    } else {
        scaling_flux_ = std::vector<double>(NumNodes_);
        for (unsigned int i = 0; i < NumNodes_; i++){
            scaling_flux_[i] = std::pow(energy_nodes_[i],scaling_index_);
        }
    }

    // RHS matrices depends on scaling flux. Reset flag.
    RHS_set_ = false;

    scaling_flux_set_ = true;
}

void nuFATE::setInitialFlux(const std::vector<double> &flux)
{
    if (flux.size() != NumNodes_) {
        throw std::runtime_error("nuFATE::nuFATE number of energy nodes of input flux doesn't match with energy nodes of cross section.");
    }

    initial_flux_set_ = false;

    // size is OK. set phi_0_ vector.
    // prepare scaling flux.
    setScalingFlux(scaling_index_);

    if(include_secondaries_){
        phi_0_ = std::vector<double>(2*NumNodes_);
        for (unsigned int i = 0; i < NumNodes_; i++){
            phi_0_[i] = flux[i]*scaling_flux_[i];
            phi_0_[i+NumNodes_] = phi_0_[i];
        }

    } else {
        phi_0_ = std::vector<double>(NumNodes_);
        for (unsigned int i = 0; i < NumNodes_; i++){
            phi_0_[i] = flux[i]*scaling_flux_[i];
        }
    }

    // set unphysical number to newgamma_ to avoid to be used for analysis
    newgamma_ = UNPHYSICAL_NUMBER;

    // ci depends on initial flux. reset flag.
    ci_set_ = false;

    initial_flux_set_ = true;

}

void nuFATE::AddAdditionalTerms(){
    hid_t group_id;
    if (newflavor_ == -3 and add_tau_regeneration_){
        std::string grptau = "/tau_decay_spectrum";
        group_id = H5Gopen(root_id_, grptau.c_str(), H5P_DEFAULT);
        hsize_t tauarraysize[2];
        H5LTget_dataset_info(group_id,"tbarfull", tauarraysize,NULL,NULL);
        size_t dim1 = tauarraysize[0];
        size_t dim2 = tauarraysize[1];
        tau_array_ = std::shared_ptr<double>((double *)malloc(dim1*dim2*sizeof(double)),free);
        H5LTread_dataset_double(group_id, "tbarfull", tau_array_.get());
        RHregen_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
        set_RHS_matrices(RHregen_, tau_array_);
        for (unsigned int i = 0; i<NumNodes_; i++){
            for(unsigned int j=0; j<NumNodes_;j++)
            *(RHSMatrix_.get()+i*NumNodes_+j) = *(RHSMatrix_.get()+i*NumNodes_+j) + *(RHregen_.get()+i*NumNodes_+j);
        }
    } else if (newflavor_ == 3 and add_tau_regeneration_){
        std::string grptau = "/tau_decay_spectrum";
        group_id = H5Gopen(root_id_, grptau.c_str(), H5P_DEFAULT);
        hsize_t tauarraysize[2];
        H5LTget_dataset_info(group_id,"tfull", tauarraysize,NULL,NULL);
        size_t dim1 = tauarraysize[0];
        size_t dim2 = tauarraysize[1];
        tau_array_ = std::shared_ptr<double>((double *)malloc(dim1*dim2*sizeof(double)),free);
        H5LTread_dataset_double(group_id, "tfull", tau_array_.get());
        RHregen_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
        set_RHS_matrices(RHregen_, tau_array_);
        for (unsigned int i = 0; i<NumNodes_; i++){
            for(unsigned int j=0; j<NumNodes_;j++){
               *(RHSMatrix_.get()+i*NumNodes_+j) = *(RHSMatrix_.get()+i*NumNodes_+j) + *(RHregen_.get()+i*NumNodes_+j);
            }
        }
    } else if (newflavor_ == -1 and add_glashow_term_){
        set_glashow_total();
        set_glashow_partial();

        if(include_secondaries_){
          for (unsigned int i = 0; i<NumNodes_;i++){
            *(glashow_partial_.get()+i*NumNodes_+i) = (*(glashow_partial_.get()+i*NumNodes_+i) - glashow_total_[i])/2.;
            for(unsigned int j = 0; j<NumNodes_;j++){
              if(i!=j)
                *(glashow_partial_.get()+i*NumNodes_+j) = (*(glashow_partial_.get()+i*NumNodes_+j)/2.);
            }
          }
          for(unsigned int i=0; i<NumNodes_;i++){
            for(unsigned int j =0; j<NumNodes_;j++){
              *(RHSMatrix_.get()+i*rsize_+j) = *(RHSMatrix_.get()+i*rsize_+j) + *(glashow_partial_.get()+i*NumNodes_+j);
            }
          }

        } else{
            for (unsigned int i = 0; i < NumNodes_; i++){
              sigma_array_[i] = sigma_array_[i] + glashow_total_[i]/2.;
              for(unsigned int j=0; j<NumNodes_;j++){
                *(RHSMatrix_.get() +i*NumNodes_+j) = *(RHSMatrix_.get() +i*NumNodes_+j) + *(glashow_partial_.get() + i*NumNodes_+j)/2.;
              }
            }
        }
    }
}

void nuFATE::allocate_gsl_buffers() 
{
   unsigned int msize = NumNodes_;
   if(include_secondaries_){
       msize = 2*NumNodes_;
   }
   eval_ = gsl_vector_complex_alloc (msize);
   evec_ = gsl_matrix_complex_alloc (msize, msize);
   w_    = gsl_eigen_nonsymmv_alloc (msize);
   ci_   = gsl_vector_alloc(msize);
   p_    = gsl_permutation_alloc(msize);
   V_    = gsl_matrix_alloc (msize,msize);
}

void nuFATE::free_gsl_buffers() 
{
   gsl_vector_complex_free(eval_);
   gsl_matrix_complex_free(evec_);
   gsl_eigen_nonsymmv_free(w_);
   gsl_vector_free(ci_);
   gsl_permutation_free (p_);
   gsl_matrix_free(V_);
}

Result nuFATE::getEigensystem(){
    if(not initial_flux_set_)
      throw std::runtime_error("nuFATE::getEigensystem initial flux not set.");
    if(not total_cross_section_set_)
      throw std::runtime_error("nuFATE::getEigensystem total cross section not set.");
    if(not differential_cross_section_set_)
      throw std::runtime_error("nuFATE::getEigensystem differential cross section not set.");

    // Set RHS matrices.
    // RHS matrices depends on:
    // 1. Cross Sections
    // 2. Energy Nodes
    // 3. Scaling Flux
    // and does not depend on the initial flux.
    // Because AddAdditionalTerms function calls file
    // open, calculate RHS matrices only when any of
    // above is touched.
    if (RHS_set_ == false) {
       set_RHS_matrices(RHSMatrix_, dxs_array_);
       if (add_secondary_term_) {
          AddAdditionalTerms();
       }

       if(not include_secondaries_){
          for (unsigned int i = 0; i < NumNodes_; i++){
             *(RHSMatrix_.get()+i*NumNodes_+i) = *(RHSMatrix_.get()+i*NumNodes_+i) - sigma_array_[i];
          }
       }
       RHS_set_ = true;

       // once RHSMatrix is modified, eigenvalues and
       // eigenvectors need to be recalculated.
       eigenvalue_and_vector_set_ = false;
    }

    unsigned int msize;
    if(include_secondaries_){
       msize = 2*NumNodes_;
    } else{
       msize = NumNodes_;
    }

    // Among 5 parameters of Result class, 
    // eigenvector and eigenvalues depends on RHSMatrix
    // and energynodes only. Separate the calculation 
    // from solving ci.
    if (eigenvalue_and_vector_set_ == false) {

       // copy energy_nodes to result
       r1_.energy_nodes_ = energy_nodes_;

       // copy RHSMatrix to gsl matrix
       gsl_matrix_view m = gsl_matrix_view_array(RHSMatrix_.get(), msize, msize);

       // calculate eigenvalues and right eivenvectors of 
       // matrix m
       gsl_eigen_nonsymmv (&m.matrix, eval_, evec_, w_);

       // sort eigenvalues and corresponding eigenvectors
       // in ascending order of magnitude
       gsl_eigen_nonsymmv_sort(eval_, evec_, GSL_EIGEN_SORT_ABS_ASC);

       // prepare buffer for eigenvectors
       r1_.evec = std::shared_ptr<double>((double *)malloc(msize*msize*sizeof(double)),free);

       r1_.eval.resize(msize);
       for(unsigned int i = 0; i<msize;i++){
          // copy eigenvalues to result
          gsl_complex eval_i
              = gsl_vector_complex_get(eval_, i);
          r1_.eval[i] = gsl_complex_abs(eval_i);
 
          // copy eigenvectors to result
          for(unsigned int j=0; j<msize;j++){
              gsl_vector_complex_view evec_i = gsl_matrix_complex_column(evec_, j);
              gsl_complex z = gsl_vector_complex_get(&evec_i.vector, i);
              double value = GSL_REAL(z);
              gsl_matrix_set(V_ , i, j, value);
              *(r1_.evec.get() + i * msize + j) = value;
          }
       }

       // factorize matrix V into the LU decomposition
       // PA = LU, V will be updated
       int s;
       gsl_linalg_LU_decomp(V_, p_, &s);

       eigenvalue_and_vector_set_ = true;
    }

    if (ci_set_ == false) {

       // copy phi_0_ to gsl vector
       gsl_vector_view b = gsl_vector_view_array(&phi_0_.front(), msize);

       // solve the square system Ax = b using the LU
       // decomposition of V into (LU, p) given by 
       // gsl_linalg_LU_decomp
       gsl_linalg_LU_solve(V_, p_, &b.vector, ci_);

       // copy phi_0_ and ci to result
       r1_.phi_0_ = phi_0_;
       r1_.ci.resize(msize);
       for (unsigned int i = 0; i<msize;i++) {
           double newval = gsl_vector_get(ci_, i);
           r1_.ci[i] = newval;
       }

       // done!
       ci_set_ = true;
    } 

    return r1_;
}

struct rho_earth_params{double theta;};

double nuFATE::rho_earth(double x, void * p){
    double RE = 6371.;
    struct rho_earth_params * params = (struct rho_earth_params *)p;
    double theta = (params->theta);
    double xmax = 2.*std::abs(RE*cos(theta));
    double R = std::pow(RE,2) + std::pow((xmax-x),2) + 2.*RE*(xmax-x)*cos(theta);
    double r = std::pow(R,0.5);
    double p1;
    double p2;
    double p3;

    if (r<1221.){
        p1 = -0.0002177;
        p2 = -4.265e-06;
        p3 = 1.309e+04;
    } else if (r<3480.){
        p1 = -0.0002409;
        p2 = 0.1416;
        p3 = 1.234e+04;
    } else if (r<5721.){
        p1 = -3.764e-05;
        p2 = -0.1876;
        p3 = 6664.;
    } else if (r<5961.){
        p1 = 0.;
        p2 = -1.269;
        p3 = 1.131e+04;
    } else if (r<6347.){
        p1 = 0.;
        p2 = -0.725;
        p3 = 7887.;
    } else if (r<6356.){
        p1 = 0.;
        p2 = 0.;
        p3 = 2900.;
    } else if (r<6368.){
        p1 = 0.;
        p2 = 0.;
        p3 = 2600.;
    } else {
        p1 = 0.;
        p2 = 0.;
        p3 = 1020.;
    }

    double rho = p1*std::pow(r,2)+p2*r+p3;
    return rho*1.0e-3;
}

double nuFATE::getEarthColumnDensity(double theta){
    double t;
    if (theta < pi/2.){
       t = 0;
    } else {
      double kmtocm = 1.0e5;
      double result, error;
      gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
      gsl_function F;
      struct rho_earth_params params = {theta};
      F.function = &nuFATE::rho_earth;
      F.params = &params;
      double xmax = 2.*std::abs(REarth*cos(theta));

      gsl_integration_qags(&F, 0, xmax, 1.0e-18, 1.0e-3, 1000, w, &result, &error);
      t = result*kmtocm;
    }

   return t;
}


int nuFATE::getFlavor() const {
    return newflavor_;
}

double nuFATE::getGamma() const {
    if (newgamma_ == UNPHYSICAL_NUMBER) {
        std::cout << "You set user-defined input function. Gamma " << newgamma_ << " is unphysical." << std::endl;
    }
    return newgamma_;
}

std::string nuFATE::getFilename() const {
    return newh5_filename_;
}

double nuFATE::getNumNodes() const {
    return NumNodes_;
}

} //close namespace
