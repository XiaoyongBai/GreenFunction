//
//  GuzinaMultLayerT.h
//  GreenFunction
//
//  Created by Xiaoyong Bai on 5/14/16.
//
//

#ifndef GuzinaMultLayerT_hpp
#define GuzinaMultLayerT_hpp

#include <iostream>
#include <complex>
#include "GuzinaBiMatT.h"

namespace GreenFunction{
    
class GuzinaMultLayerT
{
public:
    
    GuzinaMultLayerT();
    ~GuzinaMultLayerT();
    
    
    //set configuration of layers
    void SetLayers(int num_layer, double* depth);
    
    //set material constants
    //
    //input:
    //		E=Young's modulus
    //      nu=Possion's ratio
    // 		density=density
    void SetMaterial(double *E, double *nu, double *density);

    
    //set positions of source and observation
    //
    //input:
    //      src=coordinate of source point
    //      dest=coordinate of observation point
    //      normal=normal vector of the observation surface
    void SetGeometry(double* src, double* dest, double* normal);

    //drive the computation after all the settings
    //
    //input:
    //      w=circular frequency
    void Compute(double w);
    
    void ComputeTraction();
    
    
    //Get pointer to the Green's functions.
    //For displacement and tractions, the first index stands for component, and 2nd is loading direction.
    //SX is stress corresponding to loading in X-direction. Similar for Y and Z.
    //
    const std::complex<double>* GetTraction();
    const std::complex<double>* GetDisplacement();
    void GetStress(std::complex<double>** SX, std::complex<double>** SY, std::complex<double>** SZ);

private:
    //find which layer contains the point
    int WhichLayer(double z);
    
    //compute parameters of a, b, c, d, e, and f
    void Parameters_P_SV(void);
    void Parameters_SH(void);
    
    //input:
    //  w=normalized frequency
    //  p_start and p_end are two points specifying the straigt line of integration
    void ContourIntegral(double r, double z, double s, double w,
                         std::complex<double> p_start, std::complex<double> p_end, int num_int);
    
    //sign function
    double Sign(double x);
        
private:

    //fNumLayers=number of layers.
    //NOTE: this number doesn't include the half-space.
    //      It corresponds to n in Figure 5.4 of Guzina's thesis.
    //depth of interfaces, including the free surface and the artificial interface of half-space.
    int fNumLayers;
    double* fDepth;
    
    // fE=Young's modulus
    // fNu=Possion's ratio
    // fMu=Shear modulus
    // fRho=density
    double *fE, *fNu, *fMu, *fRho;
    double *C1, *C2; //wave speeds
    double cscd2; //square of ratio of Cs to Cd (the two wave speeds)
    double L_Mu; //Lambda/Mu
    
    //static solution, or singular part
    GuzinaBiMatT fBiMat;
    
    //fSrc, fDest=coordinate of source and observation point
    // fN=normal of the surface on which the traction is computing
    double fSrc[3], fDest[3];
    double fN[3];
    int fL_source, fL_dest; //number of layer containning the source point and ovservation point
        
    
    
    //fDis=part 2 of Pak displacement
    //fTrac=part 2 of Pak traction
    //fStress_*=part 2 ofPak stress due to loading in direction *
    std::complex<double> fDis[9]; //the 1st index is component, and 2nd is loading direction
    std::complex<double> fTrac[9]; // the 1st index is component, and 2nd is loading direction
    std::complex<double> fStress_X[9];
    std::complex<double> fStress_Y[9];
    std::complex<double> fStress_Z[9];
    
    ///////////////////////
    //temporary parameters
    //////////////////////
    std::complex<double> *a_n1, *a_0, *a_p1;
    std::complex<double> *b_n1, *b_0, *b_p1;
    std::complex<double> *c_n1, *c_0, *c_p1;
    std::complex<double> *d_n1, *d_0, *d_p1;
    std::complex<double> *e_n1, *e_0, *e_p1;
    std::complex<double> *f_n1, *f_0, *f_p1;
    
};//end of definition of class PakT
    


    
}//end of namespace



#endif /* PakT_hpp */
