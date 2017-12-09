#include "nuFATE.h"
#include <iostream>

namespace nufate{

nuFACE::nuFACE(int flavor, double gamma, std::string h5_filename) : newflavor_(flavor), newgamma_(gamma), newh5_filename_(h5_filename) {

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
    //allocate memory that will be used in functions below
    energy_nodes_.resize(NumNodes_);
    DeltaE_ = std::vector<double>(NumNodes_);
    glashow_total_ = std::vector<double>(NumNodes_);
    glashow_partial_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    RHSMatrix_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    Enuin_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    Enu_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    den_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    selectron_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    t1_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    t2_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    t3_ = std::shared_ptr<double>((double *)malloc(NumNodes_*NumNodes_*sizeof(double)),free);
    //set energy nodes and deltaE
    energy_nodes_ = logspace(Emin_, Emax_, NumNodes_);
//    for(unsigned int i = 0; i < NumNodes_;i++){
 //       std::cout << energy_nodes_[i] << std::endl; 
//    }
    for(unsigned int i = 0; i < NumNodes_-1;i++){
        DeltaE_[i] = log(energy_nodes_[i+1]) - log(energy_nodes_[i]);
//        std::cout <<  DeltaE_[i] << std::endl;
    }

}

//reads an attribute of type double from h5 object
double nuFACE::readDoubleAttribute(hid_t object, std::string name) const{
    double target;
    hid_t attribute_id = H5Aopen(object,name.c_str(),H5P_DEFAULT);
    herr_t status = H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &target);
    if(status<0)
        throw std::runtime_error("Failed to read attribute '"+name+"'");
    H5Aclose(attribute_id);
    return target;
}
//reads an attribute of type unsigned int from h5 object
unsigned int nuFACE::readUIntAttribute(hid_t object, std::string name) const{
    unsigned int target;
    hid_t attribute_id = H5Aopen(object,name.c_str(),H5P_DEFAULT);
    herr_t status = H5Aread(attribute_id, H5T_NATIVE_UINT, &target);
    if(status<0)
        throw std::runtime_error("Failed to read attribute '"+name+"'");
    H5Aclose(attribute_id);
    return target;
}
//given a minimum, maximum, and number of points, this function returns equally spaced array in logspace
std::vector<double> nuFACE::logspace(double Emin,double Emax,unsigned int div) const {
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
void nuFACE::set_glashow_total(){
    for(unsigned int i=0; i<NumNodes_; i++){
        glashow_total_[i] = 2.*me*energy_nodes_[i];
        double x = glashow_total_[i];
        glashow_total_[i] = 1. /3.*std::pow(GF,2)*x/pi*std::pow((1.-(std::pow(mmu,2)-std::pow(me,2))/x),2)/(std::pow((1.-x/std::pow(MW,2)),2)+std::pow(GW,2)/std::pow(MW,2))*0.676/0.1057*std::pow(hbarc,2);
        std::cout << glashow_total_[i] << std::endl;
    }
    return;
}

void nuFACE::set_glashow_partial(){
    for (unsigned int i =0; i<NumNodes_;i++){
        for(unsigned int j = 0; j<NumNodes_;j++){
            *(Enuin_.get()+i*NumNodes_+j) = energy_nodes_[i];
            *(Enu_.get()+i*NumNodes_+j) = energy_nodes_[j];
            *(Enu_.get()+i*NumNodes_+j) = 1 - *(Enu_.get()+i*NumNodes_+j)/ *(Enuin_.get()+i*NumNodes_+j);
            *(selectron_.get()+i*NumNodes_+j) = 2.*me* *(Enuin_.get()+i*NumNodes_+j);
            *(den_.get()+i*NumNodes_+j) = std::pow(1. - *(selectron_.get()+i*NumNodes_+j)/std::pow(MW,2),2) + std::pow(GW,2)/std::pow(MW,2);
            *(t1_.get()+i*NumNodes_+j) = std::pow(gR,2)/std::pow((1.+ *(Enu_.get()+i*NumNodes_+j)* *(selectron_.get()+i*NumNodes_+j)/std::pow(MZ,2)),2);
            *(t2_.get()+i*NumNodes_+j) = gL/(1.+ *(Enu_.get()+i*NumNodes_+j)* *(selectron_.get()+i*NumNodes_+j)/std::pow(MZ,2)) + (1. - *(selectron_.get()+i*NumNodes_+j)/std::pow(MW,2))/ *(den_.get()+i*NumNodes_+j);
            *(t3_.get()+i*NumNodes_+j) = GW/MW/ *(den_.get()+i*NumNodes_+j);
            if (*(Enu_.get()+i*NumNodes_+j) > 0.){
                *(glashow_partial_.get()+i*NumNodes_+j) = (std::pow(GF,2)* *(selectron_.get()+i*NumNodes_+j)/pi*(*(t1_.get()+i*NumNodes_+j)+ (std::pow(*(t2_.get()+i*NumNodes_+j) , 2) + std::pow(*(t3_.get()+i*NumNodes_+j),2)) * std::pow((1.-*(Enu_.get()+i*NumNodes_+j)),2))*std::pow(hbarc,2))/ *(Enuin_.get()+i*NumNodes_+j);
            } else {
                *(glashow_partial_.get()+i*NumNodes_+j) = 0.;
            }
        }
    }
    return;
}

void nuFACE::set_RHS_matrices(std::shared_ptr<double> RMatrix, std::shared_ptr<double> dxsarray) const{
    for(unsigned int i = 0; i < NumNodes_; i++) 
    {
        for(unsigned int j= i+1; j < NumNodes_; j++){
            double e1 = 1./ energy_nodes_[j];
            double e2 = energy_nodes_[i] * energy_nodes_[i];
            *(RMatrix.get()+i*NumNodes_+j) = DeltaE_[j - 1] * *(dxsarray.get()+j * dxsdim_[1]+i) * e1 * e2;
        }
    }
    return;
}

Result nuFACE::get_eigs() {

    hid_t group_id;
    //file_id = H5Fopen(newh5_filename_.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    //root_id = H5Gopen(file_id, "/", H5P_DEFAULT);
    //std::string grptot = "/total_cross_sections";
    group_id = H5Gopen(root_id_, grptot_.c_str(), H5P_DEFAULT);

    if (newflavor_ == -1) {
        hsize_t sarraysize[1];
        H5LTget_dataset_info(group_id,"nuebarxs", sarraysize,NULL,NULL);
        sigma_array_ = std::vector<double>(sarraysize[0]);
        H5LTread_dataset_double(group_id, "nuebarxs", sigma_array_.data());
    }  else if (newflavor_ == -2){
        hsize_t sarraysize[1];
        H5LTget_dataset_info(group_id,"numubarxs", sarraysize,NULL,NULL);
        sigma_array_ = std::vector<double>(sarraysize[0]);
        H5LTread_dataset_double(group_id, "numubarxs", sigma_array_.data());
//       for(unsigned long i = 0; i<sigma_array_.size();i++){
//           std::cout << sigma_array_[i] << std::endl;
//        }
    }  else if (newflavor_ == -3){
        hsize_t sarraysize[1];
        H5LTget_dataset_info(group_id,"nutaubarxs", sarraysize,NULL,NULL);
        sigma_array_ = std::vector<double>(sarraysize[0]);
        H5LTread_dataset_double(group_id, "nutaubarxs", sigma_array_.data());
    }  else if (newflavor_ == 1){
        hsize_t sarraysize[1];
        H5LTget_dataset_info(group_id,"nuexs", sarraysize,NULL,NULL);
        sigma_array_ = std::vector<double>(sarraysize[0]);
        H5LTread_dataset_double(group_id, "nuexs", sigma_array_.data());
    }  else if (newflavor_ == 2){
        hsize_t sarraysize[1];
        H5LTget_dataset_info(group_id,"numuxs", sarraysize,NULL,NULL);
	      sigma_array_ = std::vector<double>(sarraysize[0]);
        H5LTread_dataset_double(group_id, "numuxs", sigma_array_.data());
    }  else if (newflavor_ == 3){
        hsize_t sarraysize[1];
        H5LTget_dataset_info(group_id,"nutauxs", sarraysize,NULL,NULL);
 	      sigma_array_ = std::vector<double>(sarraysize[0]);
        H5LTread_dataset_double(group_id, "nutauxs", sigma_array_.data());
    }

    hsize_t dxarraysize[2];
    group_id = H5Gopen(root_id_, grpdiff_.c_str(), H5P_DEFAULT);
   if (newflavor_ > 0){
        H5LTget_dataset_info(group_id,"dxsnu", dxarraysize,NULL,NULL);    
        size_t dim1 = dxarraysize[0];
        size_t dim2 = dxarraysize[1];
        dxsdim_[0] = dxarraysize[0];
        dxsdim_[1] = dxarraysize[1];
        dxs_array_ = std::shared_ptr<double>((double *)malloc(dim1*dim2*sizeof(double)),free);
        H5LTread_dataset_double(group_id, "dxsnu", dxs_array_.get());
    } else {
        H5LTget_dataset_info(group_id,"dxsnubar", dxarraysize,NULL,NULL);
        size_t dim1 = dxarraysize[0];
        size_t dim2 = dxarraysize[1];
        dxsdim_[0] = dxarraysize[0];
        dxsdim_[1] = dxarraysize[1];
        dxs_array_ = std::shared_ptr<double>((double *)malloc(dim1*dim2*sizeof(double)),free);
        H5LTread_dataset_double(group_id, "dxsnu", dxs_array_.get());
    }

    set_RHS_matrices(RHSMatrix_, dxs_array_);
    if (newflavor_ == -3){
        std::string grptau = "/tau_decay_spectrum";
        group_id = H5Gopen(root_id_, grptau.c_str(), H5P_DEFAULT);
        hsize_t tauarraysize[2];
        H5LTget_dataset_info(group_id,"tbarfull", tauarraysize,NULL,NULL);
        size_t dim1 = tauarraysize[0];
        size_t dim2 = tauarraysize[1];
        tau_array_ = std::shared_ptr<double>((double *)malloc(dim1*dim2*sizeof(double)),free);
        H5LTread_dataset_double(group_id, "tbarfull", tau_array_.get());
        set_RHS_matrices(RHregen_, tau_array_); 
        for (unsigned int i = 0; i<NumNodes_; i++){
            for(unsigned int j=0; j<NumNodes_;j++)
            *(RHSMatrix_.get()+i*NumNodes_+j) = *(RHSMatrix_.get()+i*NumNodes_+j) + *(RHregen_.get()+i*NumNodes_+j);
        }
    } else if(newflavor_ == 3){
        std::string grptau = "/tau_decay_spectrum";
        group_id = H5Gopen(root_id_, grptau.c_str(), H5P_DEFAULT);
        hsize_t tauarraysize[2];
        H5LTget_dataset_info(group_id,"tfull", tauarraysize,NULL,NULL);
        size_t dim1 = tauarraysize[0];
        size_t dim2 = tauarraysize[1];
        tau_array_ = std::shared_ptr<double>((double *)malloc(dim1*dim2*sizeof(double)),free);
        H5LTread_dataset_double(group_id, "tbarfull", tau_array_.get());
        set_RHS_matrices(RHregen_, tau_array_);
        for (unsigned int i = 0; i<NumNodes_; i++){
            for(unsigned int j=0; j<NumNodes_;j++)
            *(RHSMatrix_.get()+i*NumNodes_+j) = *(RHSMatrix_.get()+i*NumNodes_+j) + *(RHregen_.get()+i*NumNodes_+j);
        }
    } else if(newflavor_ == -1){
        set_glashow_total();
        for (unsigned int i = 0; i < NumNodes_; i++){
            sigma_array_[i] = sigma_array_[i] + glashow_total_[i]/2.; 
            *(RHSMatrix_.get() +i) = *(RHSMatrix_.get() +i) + *(glashow_partial_.get() + i)/2.;

        }
    }
    phi_0_ = std::vector<double>(NumNodes_);
    for (unsigned int i = 0; i < NumNodes_; i++){
        phi_0_[i] = std::pow(energy_nodes_[i],(2.-newgamma_));
    }
    for (unsigned int i = 0; i < NumNodes_; i++){
        *(RHSMatrix_.get()+i*NumNodes_+i) = *(RHSMatrix_.get()+i*NumNodes_+i) + sigma_array_[i];    
    }
//    
//    for(unsigned int i=0;i<NumNodes_;i++){
//      for(unsigned int j=0;j<NumNodes_;j++){
//        std::cout << *(RHSMatrix_.get()+i * NumNodes_+j) << std::endl;
//     }
//    }

    //compute eigenvalues and eigenvectors

    gsl_matrix_view m = gsl_matrix_view_array(RHSMatrix_.get(), NumNodes_, NumNodes_);
    gsl_vector_complex *eval = gsl_vector_complex_alloc (NumNodes_);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc (NumNodes_, NumNodes_);
    gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (NumNodes_);

    gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);

    int s;
    gsl_vector *ci = gsl_vector_alloc(NumNodes_);
    gsl_permutation *p = gsl_permutation_alloc(NumNodes_);
    gsl_vector_view b = gsl_vector_view_array (&phi_0_.front(), NumNodes_);

    gsl_linalg_LU_decomp (&m.matrix, p, &s);
    gsl_linalg_LU_solve (&m.matrix, p, &b.vector, ci);

    //free unneeded memory
    gsl_permutation_free (p);
    gsl_eigen_nonsymmv_free (w);

    struct Result r1;
    r1.eval = eval->data;
    r1.evec = evec->data;
    r1.ci = ci->data;
    r1.energy_nodes_ = energy_nodes_;
    r1.phi_0_ = phi_0_;

    return r1;
}

struct rho_earth_params {double theta;};

double nuFACE::rho_earth(double x, void * p){
    double RE = 6371.;
    struct rho_earth_params * params = (struct my_f_params *)p;
    double theta = (params->theta);
    double xmax = 2.*abs(RE*cos(theta));
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
        p2 = -4.265e-06;
        p3 = 1.309e+04;
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

double nuFACE::get_t_earth(double theta){
    double t;
    if (theta < pi/2.){
       t = 0;
    } else {
      double kmtocm = 1.0e5;
      double result, error;
      gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
      gsl_function F;
      struct my_f_params params = {theta};
      F.function = &nuFACE::rho_earth;
      F.params = &params;
      double xmax = 2.*abs(REarth*cos(theta));

      gsl_integration_qags(&F, 0, xmax, 1.0e-18, 1.0e-3, 1000, w, &result, &error);
      t = result*kmtocm;
    }

   return t;
}


int nuFACE::getFlavor() const {
    return newflavor_;
}

double nuFACE::getGamma() const {
    return newgamma_;
}

std::string nuFACE::getFilename() const {
    return newh5_filename_;
}

double nuFACE::getNumNodes() const {
    return NumNodes_;
}

} //close namespace
