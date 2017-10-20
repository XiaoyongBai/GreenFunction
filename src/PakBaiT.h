//
//  PakBaiT.h
//  GreenFunction
//
//  Created by Xiaoyong Bai on 5/14/16.
//
//

#ifndef PakBaiT_hpp
#define PakBaiT_hpp

#include <iostream>
#include <complex>

#include "PakT.h"

namespace GreenFunction{
    
class PakBaiT
{
public:
    
    PakBaiT();
    ~PakBaiT();
        
    //set material constants
    //
    //input:
    //		E=Young's modulus
    //      nu=Possion's ratio
    // 		density=density
    void SetMaterial(double E, double nu, double density);
    
    //num=length of t
    //t=time when Green's function is evaluted
    void SetTime(int num, double* t);
    
    //set width of cubic B-Spline functions
    void SetBWidth(double T);
    
    //set positions of source and observation
    //
    //input:
    //      src=coordinate of source point
    //      dest=coordinate of observation point
    //      normal=normal vector of the observation surface
    void SetGeometry(double* src, double* dest, double* normal);

    //drive the computation after all the settings
    //
    void Compute();
    
    //Fourier transfor of B-Spline function
    std::complex<double> FTBSpline(double w);
    
    //Get pointer to the Green's functions.
    //For displacement and tractions, the first index stands for component, and 2nd is loading direction.
    //
    void GetTraction(int& num, double** traction);
    void GetDisplacement(int& num, double** displacement);
    void GetStress(std::complex<double>** SX, std::complex<double>** SY, std::complex<double>** SZ);

        
private:
    
    double fDT; //width of Cubic B-Spline basis function
    
    int fNumT; //number of time points
    double* ft; //times
    
    double* fDis;
    double* fTrac;
    
    ///////////////////////
    //temporary parameters
    //////////////////////
    int fNumW; //number of frequencies
    
    double* fW; //frequencies to be computed
    
    PakT* fPak;
    
    std::complex<double>* fFTB;//fourier transform of cubic B-Spline funcitons
    
};//end of definition of class PakT
    


    
}//end of namespace



#endif /* PakT_hpp */
