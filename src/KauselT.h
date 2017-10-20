//
//  KauselT.h
//  GreenFunction
//
//  Created by Xiaoyong Bai on 06/15/17.
//
//

#ifndef KauselT_hpp
#define KauselT_hpp

#include <iostream>
#include <complex>

namespace GreenFunction{
    
class KauselT
{
public:
        
    KauselT();
    ~KauselT();
        
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
        
    //set width of B-spline functions
    void SetBWidth(double T);
        
    //set positions of source and observation
    //
    //input:
    //      src=coordinate of source point
    //      dest=coordinate of observation point
    //Note!!! both source and receivr must be located on free surface
    void SetGeometry(double* src, double* dest);
    
    
    //drive the computation after all the settings
    //
    void Compute();
    
        
    //Get pointer to the Green's functions.
    //For displacement and tractions, the first index stands for component, and 2nd is loading direction.
    //
    void GetTraction(int& num, double** traction);
    void GetDisplacement(int& num, double** displacement);
    void GetDisplacement_Heavi(int& num, double** time, double** displacement);

    void RayleighWave(void);
    
    std::complex<double> Ellip_K( std::complex<double> n);
    std::complex<double> Ellip_PI( std::complex<double> m,std::complex<double> n);

    double Heaviside(double t);
    void CubicB1(double& B1, double& DB1, double t, double T);
    
private:
    
    int fNumT_BSpline; //number of time points
    double* ft_BSpline; //times
    
    int fNumT_Heavi;
    double* ft_Heavi;
    
    double fDT; //width of Cubic B-Spline basis function 
    
    double* fDis_Heavi; //Displacement results related to Heaviside load;
    double* fDis_BSpline;//Displacement results related to BSpline load;
    double* fTrac; //Traction are zeros, regardless of loading.
    
    double fX;
    double fY;
    double fr;
    double fsTheta;
    double fcTheta;
    
    
    
    // fE=Young's modulus
    // fNu=Possion's ratio
    // fMu=Shear modulus
    // fRho=density
    double fE, fNu, fMu, fRho;
    double Cp, Cs; //wave speeds
    double cscd2; //square of ratio of Cs to Cd (the two wave speeds)
    double fa;
    double L_Mu; //Lambda/Mu
    
    //Intermittent variables for the computation
    double K1_s; //squares of Rayleigh function roots;
    std::complex<double> K2_s;
    std::complex<double> K3_s;
    
    double gamma; //Rayleigh wave speed;

    
    
    };//end of definition of class KauselT
    
    
    
    
}//end of namespace



#endif /* KauselT_hpp */
