//
//  PakT.cpp
//  GreenFunction
//
//  Created by Xiaoyong Bai on 5/13/16.
//
//

#include "GuzinaMultLayerT.h"
#include "MathOperationT.h"
#include "BesselT.h"

#include "cmath"
#include <algorithm>


using namespace GreenFunction;
using namespace std;

GuzinaMultLayerT::GuzinaMultLayerT()
{
    // do nothing
}


GuzinaMultLayerT::~GuzinaMultLayerT()
{
    // do nothing
    delete [] fDepth;
    delete [] fE;
    delete [] fNu;
    delete [] fMu;
    delete [] C1;
    delete [] C2;
    
    delete [] a_n1;
    delete [] a_0;
    delete [] a_p1;
    
    delete [] b_n1;
    delete [] b_0;
    delete [] b_p1;
    
    delete [] c_n1;
    delete [] c_0;
    delete [] c_p1;
    
    delete [] d_n1;
    delete [] d_0;
    delete [] d_p1;
    
    delete [] e_n1;
    delete [] e_0;
    delete [] e_p1;
    
    delete [] f_n1;
    delete [] f_0;
    delete [] f_p1;
}


void GuzinaMultLayerT::SetLayers(int num_layer, double *depth)
{
    fNumLayers=num_layer;
    
    if (!fDepth)  delete [] fDepth;
    fDepth=new double[fNumLayers+2];
    
    MathOperationT::MemCopy(fNumLayers+2, depth, fDepth);
    
    if (!fE)  delete [] fE;
    fE=new double[fNumLayers+1];
    
    if (!fNu)  delete [] fNu;
    fNu=new double[fNumLayers+1];
    
    if (!fMu)  delete [] fMu;
    fMu=new double[fNumLayers+1];
    
    if (!C1)  delete [] C1;
    C1=new double[fNumLayers+1];
    
    if (!C2)  delete [] C2;
    C2=new double[fNumLayers+1];
    
    if (!a_n1)  delete [] a_n1;
    a_n1=new complex<double>[fNumLayers+1];
    if (!a_0)  delete [] a_0;
    a_0=new complex<double>[fNumLayers+1];
    if (!a_p1)  delete [] a_p1;
    a_p1=new complex<double>[fNumLayers+1];
    
    if (!b_n1)  delete [] b_n1;
    b_n1=new complex<double>[fNumLayers+1];
    if (!b_0)  delete [] b_0;
    b_0=new complex<double>[fNumLayers+1];
    if (!b_p1)  delete [] b_p1;
    b_p1=new complex<double>[fNumLayers+1];
    
    if (!c_n1)  delete [] c_n1;
    c_n1=new complex<double>[fNumLayers+1];
    if (!c_0)  delete [] c_0;
    c_0=new complex<double>[fNumLayers+1];
    if (!c_p1)  delete [] c_p1;
    c_p1=new complex<double>[fNumLayers+1];
    
    if (!d_n1)  delete [] d_n1;
    d_n1=new complex<double>[fNumLayers+1];
    if (!d_0)  delete [] d_0;
    d_0=new complex<double>[fNumLayers+1];
    if (!d_p1)  delete [] d_p1;
    d_p1=new complex<double>[fNumLayers+1];
    
    if (!e_n1)  delete [] e_n1;
    e_n1=new complex<double>[fNumLayers+1];
    if (!e_0)  delete [] e_0;
    e_0=new complex<double>[fNumLayers+1];
    if (!e_p1)  delete [] e_p1;
    e_p1=new complex<double>[fNumLayers+1];
    
    if (!f_n1)  delete [] f_n1;
    f_n1=new complex<double>[fNumLayers+1];
    if (!f_0)  delete [] f_0;
    f_0=new complex<double>[fNumLayers+1];
    if (!f_p1)  delete [] f_p1;
    f_p1=new complex<double>[fNumLayers+1];
}

void GuzinaMultLayerT::SetMaterial(double* E, double* nu, double* rho)
{
    MathOperationT::MemCopy(fNumLayers+1, E, fE);
    MathOperationT::MemCopy(fNumLayers+1, nu, fNu);
    MathOperationT::MemCopy(fNumLayers+1, rho, fRho);
    
    for (int li=0; li<fNumLayers+1; li++) {
        fMu[li]=fE[li]/(2.0*(1.0+fNu[li]));
        
        double lamda=fMu[li]*2.0*fNu[li]/(1.0-2.0*fNu[li]);
        C1[li]=sqrt((lamda+2.0*fMu[li])/fRho[li]);
        C2[li]=sqrt(fMu[li]/fRho[li]);
    }
}


void GuzinaMultLayerT::SetGeometry(double* src, double* dest, double* normal)
{
    
    MathOperationT::MemCopy(3, src, fSrc);
    MathOperationT::MemCopy(3, dest, fDest);
    MathOperationT::MemCopy(3, normal, fN);
    
    fL_source=WhichLayer(fSrc[2]);
    fL_dest=WhichLayer(fDest[2]);
    
    double nn=MathOperationT::VecNormal(3, fN);
    MathOperationT::VecScale(3, fN, 1/nn);
    
}



void GuzinaMultLayerT::Parameters_P_SV()
{

    if (fL_source==fL_dest) { //source and observation points in the same layer
        
        
    }else{
        
    }
}




//integration with Trapzoidal rule
void GuzinaMultLayerT::ContourIntegral(double r, double z, double s, double w,
                           std::complex<double> p_start, std::complex<double> p_end, int num_int)
{

}



void GuzinaMultLayerT::Compute(double w)
{

    //****************
    ComputeTraction();
}


void GuzinaMultLayerT::ComputeTraction()
{
    complex<double> traction[3];
    
    MathOperationT::multAB(3, 3, fStress_X, 3, 1, fN, traction);
    for (int i=0; i<3; i++) {
        fTrac[0*3+i]=traction[i];
    }
    
    MathOperationT::multAB(3, 3, fStress_Y, 3, 1, fN, traction);
    for (int i=0; i<3; i++) {
        fTrac[1*3+i]=traction[i];
    }
    
    MathOperationT::multAB(3, 3, fStress_Z, 3, 1, fN, traction);
    for (int i=0; i<3; i++) {
        fTrac[2*3+i]=traction[i];
    }
}


const complex<double>* GuzinaMultLayerT::GetDisplacement()
{
    return fDis;
}


const complex<double>* GuzinaMultLayerT::GetTraction()
{
    return fTrac;
}


void GuzinaMultLayerT::GetStress(std::complex<double>** SX, std::complex<double>** SY, std::complex<double>** SZ)
{
    *SX=fStress_X;
    *SY=fStress_Y;
    *SZ=fStress_Z;
}


int GuzinaMultLayerT::WhichLayer(double z)
{
    for (int i=0; i<fNumLayers+1; i++) {
        if (z>=fDepth[i] && z<fDepth[i+1]) {
            return i+1;
        }
    }
    
    if (z>=fDepth[fNumLayers+2]) {
        cout<<"!!!!!!!!!!!!!!!Warning: point is out of all layers" << endl;
        cout<<"z is " << z << " depth is " << fDepth[fNumLayers+2] <<endl;
        return fNumLayers+1;
    }
    
    cout << "!!!!!error: configurations of layers incorrect \n";
    return 0;
}


double GuzinaMultLayerT::Sign(double x)
{
    if (x>0) return 1.0;
    if (x<0) return -1.0;
    return 0.0;
}











