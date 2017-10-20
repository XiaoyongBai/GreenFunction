//
//  PakT.cpp
//  GreenFunction
//
//  Created by Xiaoyong Bai on 5/13/16.
//
//

#include "PakT.h"
#include "MindlinT.h"
#include "MathOperationT.h"
#include "BesselT.h"

#include "cmath"
#include <algorithm>


using namespace GreenFunction;
using namespace std;

PakT::PakT()
{
    // do nothing
}


PakT::~PakT()
{
    // do nothing
}


void PakT::SetMaterial(double E, double nu, double rho)
{
    fE=E;
    fNu=nu;
    fRho=rho;
    
    fMu=fE/(2*(1+nu));
    
    double lamda=fMu*2*nu/(1-2*nu);
    C1=sqrt((lamda+2*fMu)/rho);
    C2=sqrt(fMu/rho);
    
    L_Mu= lamda/fMu;
    cscd2=(1-2*nu)/(2*(1-nu)); //(C2/C1)^2
    
    //set Mateiral for fMindlin
    fMindlin.SetMaterial(E, nu);
}


void PakT::SetGeometry(double* src, double* dest, double* normal)
{
    
    MathOperationT::MemCopy(3, src, fSrc);
    MathOperationT::MemCopy(3, dest, fDest);
    MathOperationT::MemCopy(3, normal, fN);
    
    double nn=MathOperationT::VecNormal(3, fN);
    MathOperationT::VecScale(3, fN, 1/nn);
    
    fMindlin.SetGeometry(src, dest, normal);
}




void PakT::PakParameters(std::complex<double> xi, double z, double s, double w)
{
    z=w*z;
    s=w*s;
    
    double d1=abs(z-s);
    double d2=(z+s);
    double d3=z-s;
    
    complex<double> xi_2=xi*xi;
    
    complex<double> alpha=sqrt(xi_2-cscd2);
    complex<double> beta=sqrt(xi_2-1.0);
    complex<double> alpha_2=xi_2-cscd2;
    complex<double > beta_2=xi_2-1.0;
    
    complex<double> exp_alpha_d1=exp(-alpha*d1);
    complex<double> exp_beta_d1=exp(-beta*d1);
    complex<double> exp_alpha_d2=exp(-alpha*d2);
    complex<double> exp_beta_d2=exp(-beta*d2);
    complex<double> exp_alpha_z=exp(-alpha*z);
    complex<double> exp_alpha_s=exp(-alpha*s);
    complex<double> exp_beta_z=exp(-beta*z);
    complex<double> exp_beta_s=exp(-beta*s);
    
    complex<double> Rp=pow(2.0*xi_2-1.0,2)+4.0*xi_2*alpha*beta;
    complex<double> Rn=pow(2.0*xi_2-1.0,2)-4.0*xi_2*alpha*beta;
    
    
    ///////////////////////////////////
    //compute pak parameters
    fGamma1=       (xi_2/alpha)*exp_alpha_d1-beta*exp_beta_d1;
    fGamma1=fGamma1-(Rp/Rn)*((xi_2/alpha)*exp_alpha_d2+beta*exp_beta_d2);
    fGamma1=fGamma1+(4.0*xi_2*beta*(2.0*xi_2-1.0)/Rn)*(exp_beta_z*exp_alpha_s+exp_beta_s*exp_alpha_z);
    fGamma1=fGamma1/2.0;
    fGamma1=fGamma1*xi;
    
    fGamma2=(exp_beta_d1+exp_beta_d2)/(2.0*beta);
    fGamma2=fGamma2*xi;
    
    fGamma3=        Sign(z-s)*(exp_alpha_d1-exp_beta_d1)*xi/2.0;
    fGamma3=fGamma3+(xi/2.0)*(Rp/Rn)*(exp_alpha_d2+exp_beta_d2);
    fGamma3=fGamma3-(2.0*xi*(2.0*xi_2-1.0)/Rn)*(alpha*beta*exp_beta_z*exp_alpha_s+xi_2*exp_beta_s*exp_alpha_z);
    fGamma3=fGamma3*xi;
    
    fOmega1=       -1.0*Sign(z-s)*(exp_alpha_d1-exp_beta_d1)*xi/2.0;
    fOmega1=fOmega1+(xi/2.0)*(Rp/Rn)*(exp_alpha_d2+exp_beta_d2);
    fOmega1=fOmega1-(2.0*xi*(2.0*xi_2-1.0)/Rn)*(xi_2*exp_beta_z*exp_alpha_s+alpha*beta*exp_beta_s*exp_alpha_z);
    fOmega1=fOmega1*xi;
    
    fOmega2=      (-alpha*exp_alpha_d1+(xi_2/beta)*exp_beta_d1)/2.0;
    fOmega2=fOmega2-(1/2.0)*(Rp/Rn)*(alpha*exp_alpha_d2+(xi_2/beta)*exp_beta_d2);
    fOmega2=fOmega2+(2.0*xi_2*alpha*(2.0*xi_2-1.0)/Rn)*(exp_beta_z*exp_alpha_s+exp_beta_s*exp_alpha_z);
    fOmega2=fOmega2*xi;
    
    ///////////////////////////////////
    //Asymptotic part
    complex<double> exp_d1=exp(-xi*d1);
    complex<double> exp_d2=exp(-xi*d2);
    
    complex<double> Gamma1_asym, Gamma2_asym, Gamma3_asym, Omega1_asym, Omega2_asym;
    
    Gamma1_asym=(-d1*xi+(3.0-4.0*fNu))*exp_d1 + (2.0*z*s*xi_2-(3.0-4.0*fNu)*d2*xi+(5.0-12.0*fNu+8.0*fNu*fNu))*exp_d2;
    Gamma1_asym=Gamma1_asym/(8.0*(1.0-fNu));
    
    Gamma2_asym= (exp_d1+exp_d2)/2.0;
    
    Gamma3_asym= -d3*xi*exp_d1+(-2.0*z*s*xi_2-(3.0-4.0*fNu)*d3*xi+(4.0-4.0*fNu)*(1.0-2.0*fNu))*exp_d2;
    Gamma3_asym=Gamma3_asym/(8.0*(1.0-fNu));
    
    Omega1_asym= d3*xi*exp_d1+(-2.0*z*s*xi_2+(3.0-4.0*fNu)*d3*xi+(4.0-4.0*fNu)*(1.0-2.0*fNu))*exp_d2;
    Omega1_asym=Omega1_asym/(8.0*(1.0-fNu));
    
    Omega2_asym=(d1*xi+(3-4*fNu))*exp_d1 + (2.0*z*s*xi_2+(3.0-4.0*fNu)*d2*xi+(5.0-12.0*fNu+8.0*fNu*fNu))*exp_d2;
    Omega2_asym=Omega2_asym/(8.0*(1.0-fNu));
    
    //subtraction
    fGamma1=fGamma1-Gamma1_asym;
    fGamma2=fGamma2-Gamma2_asym;
    fGamma3=fGamma3-Gamma3_asym;
    fOmega1=fOmega1-Omega1_asym;
    fOmega2=fOmega2-Omega2_asym;
    
    
    ///////////////////////////////////////
    //compute derivatives
    fDG1=    -Sign(z-s)*xi_2*exp_alpha_d1+Sign(z-s)*beta_2*exp_beta_d1;
    fDG1=fDG1+(Rp/Rn)*(xi_2*exp_alpha_d2+beta_2*exp_beta_d2);
    fDG1=fDG1+((4.0*xi_2*beta*(2.0*xi_2-1.0))/Rn)*(-beta*exp_beta_z*exp_alpha_s-alpha*exp_beta_s*exp_alpha_z);
    fDG1=fDG1*xi/2.0;
    
    fDG2=-(Sign(z-s)*exp_beta_d1+exp_beta_d2)*xi/2.0;
    
    fDG3=     (xi/2.0)*(-alpha*exp_alpha_d1+beta*exp_beta_d1);
    fDG3=fDG3+(xi/2.0)*(Rp/Rn)*(-alpha*exp_alpha_d2-beta*exp_beta_d2);
    fDG3=fDG3+((2.0*xi*(2.0*xi_2-1.0))/Rn)*(alpha*beta_2*exp_beta_z*exp_alpha_s+alpha*xi_2*exp_beta_s*exp_alpha_z);
    fDG3=fDG3*xi;
    
    fDO1=     (xi/2.0)*(alpha*exp_alpha_d1-beta*exp_beta_d1);
    fDO1=fDO1-(xi/2.0)*(Rp/Rn)*(alpha*exp_alpha_d2+beta*exp_beta_d2);
    fDO1=fDO1+((2.0*xi*(2.0*xi_2-1.0))/Rn)*(beta*xi_2*exp_beta_z*exp_alpha_s+alpha_2*beta*exp_beta_s*exp_alpha_z);
    fDO1=fDO1*xi;
    
    fDO2=     Sign(z-s)*(alpha_2/2.0)*exp_alpha_d1-Sign(z-s)*(xi_2/2.0)*exp_beta_d1;
    fDO2=fDO2+(1.0/2.0)*(Rp/Rn)*(alpha_2*exp_alpha_d2+xi_2*exp_beta_d2);
    fDO2=fDO2-((2.0*xi_2*alpha*(2.0*xi_2-1.0))/Rn)*(beta*exp_beta_z*exp_alpha_s+alpha*exp_beta_s*exp_alpha_z);
    fDO2=fDO2*xi;
    
    //asymptotic part of DG1,...
    complex<double> DG1_asym, DG2_asym, DG3_asym, DO1_asym, DO2_asym;
    
    DG1_asym=         exp_d1*Sign(z-s)*(d1*xi-(4.0-4.0*fNu));
    DG1_asym=DG1_asym-exp_d2*(2*z*s*xi_2-((3.0-4.0*fNu)*d2+2.0*s)*xi+8.0*pow(1-fNu,2));
    DG1_asym=DG1_asym*xi/(8.0*(1.0-fNu));
    
    DG2_asym=-(exp_d1*Sign(z-s)+exp_d2)*xi/2.0;
    
    DG3_asym=         exp_d1*(xi*d1-1.0);
    DG3_asym=DG3_asym+exp_d2*(2.0*z*s*xi_2+xi*((3.0-4.0*fNu)*d3-2.0*s)-(8.0*fNu*fNu-16.0*fNu+7.0));
    DG3_asym=DG3_asym*xi/(8.0*(1.0-fNu));
    
    DO1_asym=         exp_d1*(1.0-xi*d1);
    DO1_asym=DO1_asym+exp_d2*(2.0*z*s*xi_2-xi*(z*(3.0-4.0*fNu)-s*(1.0-4.0*fNu))-(1.0-8.0*fNu+8.0*fNu*fNu));
    DO1_asym=DO1_asym*xi/(8.0*(1.0-fNu));
    
    DO2_asym=        -exp_d1*Sign(z-s)*(xi*d1+2.0-4.0*fNu);
    DO2_asym=DO2_asym-exp_d2*(2.0*z*s*xi_2+xi*(d2*(3.0-4.0*fNu)-2.0*s)+2.0*pow(1.0-2.0*fNu,2));
    DO2_asym=DO2_asym*xi/(8.0*(1.0-fNu));
    
    fDG1=fDG1-DG1_asym;
    fDG2=fDG2-DG2_asym;
    fDG3=fDG3-DG3_asym;
    fDO1=fDO1-DO1_asym;
    fDO2=fDO2-DO2_asym;
}


//integration with Trapzoidal rule
void PakT::ContourIntegral(double r, double z, double s, double w,
                           std::complex<double> p_start, std::complex<double> p_end, int num_int)
{
    //initialization of the integral
    for (int i=0; i<3; i++) {
        UX[i]=0;
        UZ[i]=0;
    }
    for (int i=0; i<9; i++) {
        SX[i]=0;
        SZ[i]=0;
    }
    
    
    complex<double> d_xi=(p_end-p_start)/(num_int+0.0);
    
    for (int ii=1; ii<=num_int; ii++){
    
        complex<double> xi=p_start+(ii-0.5)*d_xi;
    
        PakParameters(xi, z, s, w);
    
        complex<double>  J0=BesselT::BesselJ0(xi*w*r);
        complex<double>  J1=BesselT::BesselJ1(xi*w*r);
        complex<double>  J2=BesselT::BesselJ(2,xi*w*r);
    
        UX[0]=UX[0]+(fGamma1*(J0-J2)+fGamma2*(J0+J2))*d_xi;
        UX[1]=UX[1]+(fGamma1*(J0+J2)+fGamma2*(J0-J2))*d_xi;
        UX[2]=UX[2]+fOmega1*J1*d_xi;
    
        UZ[0]=UZ[0]+fGamma3*J1*d_xi;
        UZ[2]=UZ[2]+fOmega2*J0*d_xi;
    
        SX[0*3+0]=SX[0*3+0]+(L_Mu*fDO1-(L_Mu+2.0)*xi*fGamma1)*J1*d_xi;
        SX[0*3+1]=SX[0*3+1]+xi*fGamma2*J1*d_xi;
        SX[0*3+2]=SX[0*3+2]+(fDG2+fDG1+xi*fOmega1)*J0*d_xi+(fDG2-fDG1-xi*fOmega1)*J2*d_xi;
        SX[1*3+1]=SX[1*3+1]+L_Mu*(fDO1-xi*fGamma1)*J1*d_xi;
        SX[1*3+2]=SX[1*3+2]+(fDG2+fDG1+xi*fOmega1)*J0*d_xi-(fDG2-fDG1-xi*fOmega1)*J2*d_xi;
        SX[2*3+2]=SX[2*3+2]+((L_Mu+2.0)*fDO1-L_Mu*xi*fGamma1)*J1*d_xi;
    
        SZ[0*3+0]=SZ[0*3+0]+(L_Mu*fDO2-(L_Mu+2.0)*xi*fGamma3)*J0*d_xi;
        SZ[0*3+2]=SZ[0*3+2]-(fDG3+xi*fOmega2)*J1*d_xi;
        SZ[1*3+1]=SZ[1*3+1]+L_Mu*(fDO2-xi*fGamma3)*J0*d_xi;
        SZ[2*3+2]=SZ[2*3+2]+ ((L_Mu+2.0)*fDO2-L_Mu*xi*fGamma3)*J0*d_xi;
    }
    
    for (int i=0; i<3; i++) {
        UX[i]=UX[i]*w;
        UY[i]=UX[i];
        UZ[i]=UZ[i]*w;
    }

    
    SX[1*3+0]=SX[0*3+1];
    SX[2*3+0]=SX[0*3+2];
    SX[2*3+1]=SX[1*3+2];
    
    double w2=w*w;
    
    for (int i=0; i<9; i++) {
        SX[i]=SX[i]*w2;
        SY[i]=SX[i];
        SZ[i]=SZ[i]*w2;
    }
    
    SZ[2*3+0]=SZ[0*3+2];
}



void PakT::Compute(double w)
{
    //cout<<"frequency =" << w<<endl;
    
    fMindlin.Compute();
    
    const double* U_static=fMindlin.GetDisplacement();
    double* SX_static;
    double* SY_static;
    double* SZ_static;
    
    fMindlin.GetStress(&SX_static, &SY_static, &SZ_static);
    
    //For very small frequency, static solution is used.
    if (w<1e-3) {
        for (int i=0; i<9; i++) {
            fDis[i]=U_static[i];
            fStress_X[i]=SX_static[i];
            fStress_Y[i]=SY_static[i];
            fStress_Z[i]=SZ_static[i];
        }
        
        //****************
        ComputeTraction();
        return;
    }
    
    
    double s=fSrc[2]; //depth of the source point
    double x=fDest[0]-fSrc[0];
    double y=fDest[1]-fSrc[1];
    double z=fDest[2]; //depth of the field point
    
    double r=sqrt(x*x+y*y);
    double stheta=y/r;
    double ctheta=x/r;
    if (r==0) {
        stheta=0;
        ctheta=1;
    }
    
    //specify the reference length a
    double a;
    
    if (r>=abs(z-s)){
        a=r;
    }
    else{
        a=abs(z-s);
    }
    
    
    //compute Dimensionless parameters
    double rBar=r/a;
    double zBar=z/a;
    double sBar=s/a;
    double wBar=w*a/C2;
    
    //specify contour path
    double HMAX=10.0;
    double h;
    if (wBar*rBar<=HMAX){
        h=1.0;
    }
    else{
        h=HMAX/(wBar*rBar);
    }
    
    
    complex<double> I(0.0, 1.0);
    complex<double> tail=5.0;
    complex<double> p0, p1, p2, p3; //Four points on the triangular path
    p0=0.0;
    p1=1.0+ h* I;
    p2=2.0;
    p3=tail;
    
    
    //number of integration points
    int NN=80;
    
    //initializtion
    complex<double> UX_temp[3];
    complex<double> UY_temp[3];
    complex<double> UZ_temp[3];
    complex<double> SX_temp[9];
    complex<double> SY_temp[9];
    complex<double> SZ_temp[9];
    
    for (int i=0; i<3; i++) {
        UX_temp[i]=0.0;
        UY_temp[i]=0.0;
        UZ_temp[i]=0.0;
    }
    for (int i=0; i<9; i++) {
        SX_temp[i]=0.0;
        SY_temp[i]=0.0;
        SZ_temp[i]=0.0;
    }
    
    
    
    //integration over the first segment
    int num_int=ceil(NN*real(p1-p0)*max(1.0, wBar));
    ContourIntegral(rBar, zBar, sBar, wBar, p0, p1, num_int);
    
    for (int i=0; i<3; i++) {
        UX_temp[i]+=UX[i];
        UY_temp[i]+=UY[i];
        UZ_temp[i]+=UZ[i];
    }
    for (int i=0; i<9; i++) {
        SX_temp[i]+=SX[i];
        SY_temp[i]+=SY[i];
        SZ_temp[i]+=SZ[i];
    }
    
    //integration over the second segment
    num_int=ceil(NN*real(p2-p1)*max(1.0, wBar));
    ContourIntegral(rBar, zBar, sBar, wBar, p1, p2, num_int);
    
    for (int i=0; i<3; i++) {
        UX_temp[i]+=UX[i];
        UY_temp[i]+=UY[i];
        UZ_temp[i]+=UZ[i];
    }
    for (int i=0; i<9; i++) {
        SX_temp[i]+=SX[i];
        SY_temp[i]+=SY[i];
        SZ_temp[i]+=SZ[i];
    }
    
    //integration over the third segment
    num_int=ceil(NN*real(p3-p2)*max(1.0, wBar));
    ContourIntegral(rBar, zBar, sBar, wBar, p2, p3, num_int);
    
    for (int i=0; i<3; i++) {
        UX_temp[i]+=UX[i];
        UY_temp[i]+=UY[i];
        UZ_temp[i]+=UZ[i];
    }
    for (int i=0; i<9; i++) {
        SX_temp[i]+=SX[i];
        SY_temp[i]+=SY[i];
        SZ_temp[i]+=SZ[i];
    }
    
    
    //compute regular part of displacement
    complex<double> UX_part2[3], UY_part2[3], UZ_part2[3];
    
    UX_part2[0]=UX_temp[0]*ctheta/(4.0*M_PI);
    UX_part2[1]=UX_temp[1]*(-stheta)/(4.0*M_PI);
    UX_part2[2]=UX_temp[2]*ctheta/(2.0*M_PI);
    
    UY_part2[0]=UY_temp[0]*stheta/(4.0*M_PI);
    UY_part2[1]=UY_temp[1]*ctheta/(4.0*M_PI);
    UY_part2[2]=UY_temp[2]*stheta/(2.0*M_PI);
    
    UZ_part2[0]=UZ_temp[0]/(-2.0*M_PI);
    UZ_part2[1]=0;
    UZ_part2[2]=UZ_temp[2]/(2.0*M_PI);
    
    //compute regular part of SX
    complex<double> SX_part2[9], SY_part2[9], SZ_part2[9];
    
    if (rBar==0.0) {
        rBar+=1e-20;
    }
    
    SX_part2[0*3+0]=SX_temp[0*3+0]*ctheta/(2.0*M_PI)-(2.0/rBar)*(UX_part2[0]-UX_temp[1]*ctheta/(4.0*M_PI));
    SX_part2[0*3+1]=SX_temp[0*3+1]*stheta/(2.0*M_PI)-(2.0/rBar)*(UX_part2[1]+UX_temp[0]*stheta/(4.0*M_PI));
    SX_part2[0*3+2]=SX_temp[0*3+2]*ctheta/(4.0*M_PI);
    SX_part2[1*3+1]=SX_temp[1*3+1]*ctheta/(2.0*M_PI)+(2.0/rBar)*(UX_part2[0]-UX_temp[1]*ctheta/(4.0*M_PI));
    SX_part2[1*3+2]=SX_temp[1*3+2]*(-stheta)/(4.0*M_PI);
    SX_part2[2*3+2]=SX_temp[2*3+2]*ctheta/(2.0*M_PI);
    SX_part2[1*3+0]=SX_part2[0*3+1];
    SX_part2[2*3+0]=SX_part2[0*3+2];
    SX_part2[2*3+1]=SX_part2[1*3+2];
    
    //compute regular part of SY
    SY_part2[0*3+0]=SY_temp[0*3+0]*stheta/(2.0*M_PI)-(2.0/rBar)*(UY_part2[0]-UY_temp[1]*stheta/(4.0*M_PI));
    SY_part2[0*3+1]=SY_temp[0*3+1]*(-ctheta)/(2.0*M_PI)-(2.0/rBar)*(UY_part2[1]+UY_temp[0]*(-ctheta)/(4.0*M_PI));
    SY_part2[0*3+2]=SY_temp[0*3+2]*stheta/(4.0*M_PI);
    SY_part2[1*3+1]=SY_temp[1*3+1]*stheta/(2.0*M_PI)+(2.0/rBar)*(UY_part2[0]-UY_temp[1]*stheta/(4.0*M_PI));
    SY_part2[1*3+2]=SY_temp[1*3+2]*ctheta/(4.0*M_PI);
    SY_part2[2*3+2]=SY_temp[2*3+2]*stheta/(2.0*M_PI);
    SY_part2[1*3+0]=SY_part2[0*3+1];
    SY_part2[2*3+0]=SY_part2[0*3+2];
    SY_part2[2*3+1]=SY_part2[1*3+2];
    
    
    
    //compute regular part of SZ
    for (int i=0; i<9; i++) {
        SZ_part2[i]=SZ_temp[i]/(2*M_PI);
    }
    
    if (rBar!=0) {
        SZ_part2[0*3+0]+=-2.0*UZ_part2[0]/rBar;
        SZ_part2[1*3+1]+= 2.0*UZ_part2[0]/rBar;
    }

    
    //tranfer form cylindrical to Cartesian system
    fDis[0*3+0] = UX_part2[0]*ctheta-UX_part2[1]*stheta;
    fDis[1*3+0] = UX_part2[0]*stheta+UX_part2[1]*ctheta;
    fDis[2*3+0] = UX_part2[2];
    
    fDis[0*3+1] = UY_part2[0]*ctheta-UY_part2[1]*stheta;
    fDis[1*3+1] = UY_part2[0]*stheta+UY_part2[1]*ctheta;
    fDis[2*3+1] = UY_part2[2];
    
    fDis[0*3+2] = UZ_part2[0]*ctheta-UZ_part2[1]*stheta;
    fDis[1*3+2] = UZ_part2[0]*stheta+UZ_part2[1]*ctheta;
    fDis[2*3+2] = UZ_part2[2];
    
    
    
    double R[]={ctheta, -stheta, 0, stheta,  ctheta, 0,  0, 0, 1};
    double RT[]={ctheta, stheta, 0,  -stheta,  ctheta, 0,  0, 0, 1};
    
    
    complex<double> MM[9]; //temporary matrix for multiplications
    
    MathOperationT::multAB(3, 3, R, 3, 3, SX_part2, MM);
    MathOperationT::multAB(3, 3, MM, 3, 3, RT, fStress_X);
    
    MathOperationT::multAB(3, 3, R, 3, 3, SY_part2, MM);
    MathOperationT::multAB(3, 3, MM, 3, 3, RT, fStress_Y);
    
    MathOperationT::multAB(3, 3, R, 3, 3, SZ_part2, MM);
    MathOperationT::multAB(3, 3, MM, 3, 3, RT, fStress_Z);

    
    //back to non-normalized
    double aa=a*a;
    
    for ( int i=0; i<9; i++) {
        fDis[i]=fDis[i]/(a*fMu)+U_static[i];
        fStress_X[i]=fStress_X[i]/aa+SX_static[i];
        fStress_Y[i]=fStress_Y[i]/aa+SY_static[i];
        fStress_Z[i]=fStress_Z[i]/aa+SZ_static[i];
    }
    
    //****************
    ComputeTraction();
    
}


void PakT::ComputeTraction()
{
    complex<double> traction[3];
    
    MathOperationT::multAB(3, 3, fStress_X, 3, 1, fN, traction);
    for (int i=0; i<3; i++) {
        fTrac[i*3+0]=traction[i];
    }
    
    MathOperationT::multAB(3, 3, fStress_Y, 3, 1, fN, traction);
    for (int i=0; i<3; i++) {
        fTrac[i*3+1]=traction[i];
    }
    
    MathOperationT::multAB(3, 3, fStress_Z, 3, 1, fN, traction);
    for (int i=0; i<3; i++) {
        fTrac[i*3+2]=traction[i];
    }
    
    
    
    //MathOperationT::PrintVector(3, fSrc, "source for pak");
    //MathOperationT::PrintVector(3, fDest, "dest for pak");
    //MathOperationT::PrintVector(3, fN, "Normal for pak");
    //MathOperationT::PrintMatrix(3, 3, fDis, "displacement");
    //MathOperationT::PrintMatrix(3, 3, fTrac, "traction");
}


const complex<double>* PakT::GetDisplacement()
{
    return fDis;
}


const complex<double>* PakT::GetTraction()
{
    return fTrac;
}


void PakT::GetStress(std::complex<double>** SX, std::complex<double>** SY, std::complex<double>** SZ)
{
    *SX=fStress_X;
    *SY=fStress_Y;
    *SZ=fStress_Z;
}


double PakT::Sign(double x)
{
    if (x>0) return 1.0;
    if (x<0) return -1.0;
    return 0.0;
}











