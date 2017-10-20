//
//  GuzinaBaiT.h
//  GreenFunction
//
//  Created by Xiaoyong Bai on 7/15/16.
//
//

#ifndef GuzinaBaiT_hpp
#define GuzinaBaiT_hpp

#include <iostream>
#include <complex>

namespace GreenFunction{
    
class GuzinaBaiT
{
public:
    
    GuzinaBaiT();
    ~GuzinaBaiT();
    
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
    
    //integration of B(w)*exp(i*w*t) in subinterval
    std::complex<double> Integral_Sub(double w1, double w2, double t);
    
    //Get pointer to the Green's functions.
    //For displacement and tractions, the first index stands for component, and 2nd is loading direction.
    //
    void GetTraction(int& num, double** traction);
    void GetDisplacement(int& num, double** displacement);
    double* GetStaticTraction(void);
    void GetStress(int& num, double ** SX, double ** SY, double ** SZ);

        
private:
    
    double fSrc[3], fDest[3], fNormal[3];
    
    double fDT; //width of Cubic B-Spline basis function
    
    int fNumT; //number of time points
    double* ft; //times
    
    double* fDis;
    double* fTrac;
    double fTrac_static[9];
    
    double* fSX;//stress corresonding load in X-direction
    double* fSY;//stress corresonding load in Y-direction
    double* fSZ;//stress corresonding load in Z-direction
    
    ///////////////////////
    //temporary parameters
    //////////////////////
    int fNumW; //number of frequencies
    
    double* fW; //frequencies to be computed
    
    std::complex<double> *fDis_FD, *fTrac_FD;
    std::complex<double> *fSX_FD, *fSY_FD, *fSZ_FD; //Stresses in frequency domain
    
};//end of definition of class PakT
    


    
}//end of namespace



#endif /* GuzinaBaiT_hpp */
