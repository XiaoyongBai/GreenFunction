//
//  PakT.h
//  GreenFunction
//
//  Created by Xiaoyong Bai on 5/14/16.
//
//

#ifndef PakT_hpp
#define PakT_hpp

#include <iostream>
#include <complex>
#include "MindlinT.h"

namespace GreenFunction{
    
class PakT
{
public:
    
    PakT();
    ~PakT();
        
    //set material constants
    //
    //input:
    //		E=Young's modulus
    //      nu=Possion's ratio
    // 		density=density
    void SetMaterial(double E, double nu, double density);
    
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
    //compute  numerically-evaluated part of Pak parameters
    //Note!!!! asymptotic part is extracted from the parameters
    //
    //input:
    //  w=normalized frequency
    void PakParameters(std::complex<double> xi, double z, double s, double w);
    
    //input:
    //  w=normalized frequency
    //  p_start and p_end are two points specifying the straigt line of integration
    void ContourIntegral(double r, double z, double s, double w,
                         std::complex<double> p_start, std::complex<double> p_end, int num_int);
    
    //sign function
    double Sign(double x);
        
private:

    // fE=Young's modulus
    // fNu=Possion's ratio
    // fMu=Shear modulus
    // fRho=density
    double fE, fNu, fMu, fRho;
    double C1, C2; //wave speeds
    double cscd2; //square of ratio of Cs to Cd (the two wave speeds)
    double L_Mu; //Lambda/Mu
    
    //static solution, or singular part
    MindlinT fMindlin;
    
    //fSrc, fDest=coordinate of source and observation point
    // fN=normal of the surface on which the traction is computing
    double fSrc[3], fDest[3];
    double fN[3];
        
    
    
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
    
    //part 2 of Pak parameters and theire derivatives with respect to z.
    //Pak parameters are computed by funciton PakParameters with normalized inputs.
    std::complex<double> fGamma1, fGamma2, fGamma3, fOmega1, fOmega2;
    std::complex<double> fDG1, fDG2, fDG3, fDO1, fDO2;
    
    
    std::complex<double> UX[3];
    std::complex<double> UY[3];
    std::complex<double> UZ[3];
    std::complex<double> SX[9];
    std::complex<double> SY[9];
    std::complex<double> SZ[9];
    
};//end of definition of class PakT
    


    
}//end of namespace



#endif /* PakT_hpp */
