//
//  GuzinaBiMatT.h
//  GreenFunction
//
//  Created by Xiaoyong Bai on 5/14/16.
//
//

#ifndef BuzinaBiMatT_hpp
#define BuzinaBiMatT_hpp

#include <iostream>

namespace GreenFunction{
    
class GuzinaBiMatT
{
public:
    
    GuzinaBiMatT();
    ~GuzinaBiMatT();
        
    //set material constants
    //
    //input:
    //		E=Young's modulus
    //      nu=Possion's ratio
    // 		density=density
    void SetMaterial(double E_1, double nu_1, double E_2, double nu_2);
    
    //set positions of source and observation
    //
    //input:
    //      src=coordinate of source point
    //      dest=coordinate of observation point
    //      normal=normal vector of the observation surface
    void SetGeometry(double* src, double* dest, double* normal);

    //drive the computation after all the settings
    void Compute();
    void ComputeTraction();
    
    //compute Guzina parameters
    void GuzinaParameters(double& Omega1, double& Omega2, double& Gamma1, double& Gamma2, double& Gamma3,
                          double& DO1,    double& DO2,    double& DG1,    double& DG2,    double& DG3,
                          int n, int l, double r, double z, double s);
    
    //Definite integral of combination of Bessel, exponential and polynomial functions.
    //See equation (4.55) in Guzina's thesis
    double BesselInt(int n, int l, double r, double d);
    
    //Get pointer to the Green's functions.
    //For displacement and tractions, the first index stands for component, and 2nd is loading direction.
    //SX is stress corresponding to loading in X-direction. Similar for Y and Z.
    //
    const double* GetTraction();
    const double* GetDisplacement();
    void GetStress(double ** SX, double ** SY, double ** SZ);

    double Sign(double x);
    
private:

    // fE1, fNu1, fMu1=Material parameters for the half-space containing the source point
    // fE1, fNu1, fMu1=Material parameters for the other half-space
    double fE1, fE2, fNu1, fNu2, fMu1, fMu2, fLambda1, fLambda2;
    
    //auxilliary constants
    double o1v1, o1v2, o2v1, o2v2, t4v1, t4v2, t2v2, o4v2, f4v2, mx1, mx2;
    double c[10];
    
    // fSrc, fDest=coordinate of source and observation point
    // fN=normal of the surface on which the traction is computing
    double fSrc[3];
    double fDest[3];
    double fN[3];
    
        
    //fDis=displacement
    //fTrac=traction
    //fStress_X=stress due to loading in x-direction
    double fDis[9];
    double fTrac[9];
    double fStress_X[9];
    double fStress_Y[9];
    double fStress_Z[9];
    
};//end of definition of class MindlinT
    


    
}//end of namespace



#endif /* BuzinaBiMatT_hpp */
