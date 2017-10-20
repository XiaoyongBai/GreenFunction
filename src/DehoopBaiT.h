//
//  DeHoopBaiT.h
//  GreenFunction
//
//  Created by Xiaoyong Bai on 06/19/17.
//
//

#ifndef DehoopBaiT_hpp
#define DehoopBaiT_hpp

#include <iostream>
#include <complex>

#include "GuzinaBiMatT.h"

namespace GreenFunction{
    
class DehoopBaiT
{
public:
        
    DehoopBaiT();
    ~DehoopBaiT();
    
    //set layer configurations
    //input:
    //      z=depth of the layer interface
    //		E=Young's modulus
    //      nu=Possion's ratio
    // 		density=density
    void SetLayers(int nl, double* z_interface, double* E, double* nu, double* density);
    
        
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
    //      Normal=surface normal
    void SetGeometry(double* src, double* dest, double* normal);
    
    
    //drive the computation after all the settings
    //
    void Compute();
    
    void ComputeStaticPart();
    
    //Compute BSpline resuslts by synthesis of the Heaviside results
    void Heavi2BSpline();
        
    //Get pointer to the Green's functions.
    //For displacement and tractions, the first index stands for component, and 2nd is loading direction.
    //
    void GetTraction(int& num, double** traction);
    void GetDisplacement(int& num, double** displacement);
    void GetDisplacement_Heavi(int& num, double** time, double** displacement);
    
    void GetStress_BSpline(int& num, double** time, double** SX, double** SY, double** SZ);
    void GetStress_Heavi(int& num, double** time, double** SX_1, double** SX_2, double** SY_1, double** SY_2,
                         double** SZ_1, double** SZ_2);
    
    const double* GetStaticTraction();
    
    double Heaviside(double t);
    
    
    //B1=Bspline base function
    //DB1=The first derivative of B1
    //DDB1=The second derivative of B1
    void CubicB1(double& B1, double& DB1, double& DDB1, double t, double T);
    
    //Generate both transmission and relfection waves
    void GenerateWaves(void);
    void DeleteWaves(void);
    
    //Expand waves (WL1) in level L1 to waves (WL2) in level L2=L1+1
    //N1=number of waves in WL1
    //N2=number of waves in WL2
    void WaveExpand(int Level1, int N1, int *ID1, int* Layer1, int& Level, int& N2, int *&ID2, int*& Layer2);
    
    //Find z is in which layer
    int WhichLayer(double zcoord);
    
    //Compute arrival time of a specified wave
    //ns=number of segments;
    //r=horizontal distance;
    //d=vertical distances;
    //c=wave speeds;
    //id=wave ids;
    //layer=layers the wave propagates in;
    double WaveArrival(int ns, double r, double* d, double* c, int* id, int* layer);
    
    //!!note that the input should be normalized
    //t=time of interested
    //r=horiztonal distance
    //ns=number of segement, i.e., the length of d and c.
    //d=vertical distances
    //c=wave speeds
    //cPhi=cos(Phi)
    //xi, dxi_dt=output
    void Find_xi_Newton(double t, double r, int ns, double* d, double* c, std::complex<double> cPhi,
                        std::complex<double>& xi, std::complex<double>& dxi_dt);
    
    
    void WaveIntegrand(std::complex<double> xi, int ns, int* ID, int* Type, int* Direction, int* Layer,
                       std::complex<double> phi, std::complex<double>* M, std::complex<double>* N);
    
    //Sd, Su=source term for p-sv wave
    //su, sd=source term for sh-wave
    void Source_X(std::complex<double> xi, double mu, double rho, double cd, double cs, int ID,
                  std::complex<double>& Source_sv, std::complex<double>& Source_sh);
    void Source_Z(std::complex<double> xi, double mu, double rho, double cd, double cs, int ID,
                  std::complex<double>& Source_sv);
    
    
    //L1, L2 the two layer near the interface
    //D1, D2 = direction of two segement
    //Type1, Type2=type of the two segments
    void Propagator(std::complex<double>xi, int L1, int L2, int D1, int D2, int Type1, int Type2,
               std::complex<double>& SV_factor, std::complex<double>& SH_factor);
    
    //Formulate CC matrix
    void CC(std::complex<double> xi, std::complex<double>* cc_sv, std::complex<double>* cc_sh);
    
private:
    
    int fNumT_BSpline; //number of time points
    double* ft_BSpline; //times
    double ft_max;//the maximum time of interest
    
    int fNumT_Heavi;
    double* ft_Heavi;
    
    double fDT; //width of Cubic B-Spline basis function 
    
    double* fDis_Heavi; //Displacement results related to Heaviside load;
    double* fDis_BSpline;//Displacement results related to BSpline load;
    
    double *fSX_Heavi_p1, *fSY_Heavi_p1, *fSZ_Heavi_p1; //Stresses for Heaviside loading
    double *fSX_Heavi_p2, *fSY_Heavi_p2, *fSZ_Heavi_p2;
    
    double *fSX_BSpline, *fSY_BSpline, *fSZ_BSpline;
    
    double *fTrac_BSpline;
    
    double fX;
    double fY;
    double fZ; //depth of receiver point
    double fS; //depth of source point
    double fr; //Horizontal distance
    double fsTheta;
    double fcTheta;
    
    double fNormal[3];
    
    int fS_layer;//Layer contains the source point
    int fZ_layer;//Layer contains the receiver point
    
    // fE=Young's modulus
    // fNu=Possion's ratio
    // fMu=Shear modulus
    // fRho=density
    int fNumLayer;
    double* fZ_Interface;
    
    
    double *fE, *fNu, *fMu, *fRho;
    double *fCp, *fCs; //wave speeds

    //Information related to the waves
    int fMinLevel; //Level = number of segments
    int fMaxLevel;
    
    int* fLevel;
    int* fN;
    int** fID; //ID of the wave group. Each level corresponds to a pointer.
    int** fType; //Type of the wave group. Type 1=Compressional, Type 2=shear;
    int** fDirection; //Direction of the wave. 1=upgoing, 2=downgoing;
    int** fLayer;//Layers of the waves
    double** fC;//wave speeds related to the waves
    double** fd;//vertical distances of the waves
    
    double** fArrival; //Arrival time of each wave
    
    
    //normalized model
    double fa; //reference legnth
    double fMu_ref; //reference shear modulus;
    double fRho_ref; //reference density;
    double fCs_ref; //reference wave speed;
    double ft_ref; //reference time;
    
    double *fMuBar, *fLamdaBar, *fRhoBar, *fCpBar, *fCsBar; //Normalized material parameters;
    
    double frBar; //Normalized horizontal distance;
    double** fCBar; //Normalized wave speed
    double** fdBar; //Normalized vrtical distance;
    double** fArrivalBar; //NOrmalized arrival time;
    
    //Static part of the stresses and displacements
    GuzinaBiMatT fGuzinaBiMat;
    
    };//end of definition of class DehoopBaiT
    
    
    
    
}//end of namespace



#endif /* DehoopBaiT */
