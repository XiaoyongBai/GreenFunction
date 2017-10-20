//
//  MindlinT.h
//  GreenFunction
//
//  Created by Xiaoyong Bai on 5/14/16.
//
//

#ifndef MindlinT_hpp
#define MindlinT_hpp

#include <iostream>

namespace GreenFunction{
    
class MindlinT
{
public:
    
    MindlinT();
    ~MindlinT();
        
    //set material constants
    //
    //input:
    //		E=Young's modulus
    //      nu=Possion's ratio
    // 		density=density
    void SetMaterial(double E, double nu);
    
    //set positions of source and observation
    //
    //input:
    //      src=coordinate of source point
    //      dest=coordinate of observation point
    //      normal=normal vector of the observation surface
    void SetGeometry(double* src, double* dest, double* normal);

    //drive the computation after all the settings
    void Compute();
    void ComputeDisplacement();
    void ComputeStress();
    void ComputeTraction();
    
    
    //Get pointer to the Green's functions.
    //For displacement and tractions, the first index stands for component, and 2nd is loading direction.
    //SX is stress corresponding to loading in X-direction. Similar for Y and Z.
    //
    const double* GetTraction();
    const double* GetDisplacement();
    void GetStress(double ** SX, double ** SY, double ** SZ);

    
private:

    // fE=Young's modulus
    // fNu=Possion's ratio
    // fMu=Shear modulus
    // fSrc, fDest=coordinate of source and observation point
    // fR1=distance from source point to observation point
    // fR2=distance from mirror of source point to observation point
    // fN=normal of the surface on which the traction is computing
    
    double fE, fNu, fMu;
    double fSrc[3];
    double fDest[3];
    double fR1, fR2;
    double fN[3];
    
    //Rotation matrix for loading in Y-directions
    double fRotation[9];
    double fRotationT[9];
        
    /**
    * fDT_Kelvin_1, ..=Time integrations of Kelvin's tractions
    * fDU_Stokes_1, ..=Time integrations of Stokes' displacements
    * fDT_Stokes_1, ..=Time integrations of Stokes' tractions */
    double fDis[9];
    double fTrac[9];
    double fStress_X[9];
    double fStress_Y[9];
    double fStress_Z[9];
    
};//end of definition of class MindlinT
    


    
}//end of namespace



#endif /* MindlinT_hpp */
