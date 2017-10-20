//
//  GuzinaBiMatT.cpp
//  GreenFunction
//
//  Created by Xiaoyong Bai on 5/13/16.
//
//

#include "GuzinaBiMatT.h"
#include "MathOperationT.h"
#include "cmath"

using namespace GreenFunction;
using namespace std;

GuzinaBiMatT::GuzinaBiMatT()
{
    //do nothing

    
}

GuzinaBiMatT::~GuzinaBiMatT()
{
    //do nothing
}

//!!!Note: E_1 and nu_1 must be of the layer containing the source;
void GuzinaBiMatT::SetMaterial(double E_1, double nu_1, double E_2, double nu_2)
{
    fE1=E_1;
    fNu1=nu_1;
    fMu1=fE1/(2.0*(1.0+nu_1));
    fLambda1=fE1*fNu1/((1.0+fNu1)*(1.0-2.0*fNu1));
    
    fE2=E_2;
    fNu2=nu_2;
    fMu2=fE2/(2.0*(1.0+nu_2));
    fLambda2=fE2*fNu2/((1.0+fNu2)*(1.0-2.0*fNu2));
    
    //compute auxilliary constants
    o1v1 = 1.0-1.0*fNu1;
    o1v2 = 1.0-1.0*fNu2;
    o2v1 = 1.0-2.0*fNu1;
    o2v2 = 1.0-2.0*fNu2;
    t4v1 = 3.0-4.0*fNu1;
    t4v2 = 3.0-4.0*fNu2;
    t2v2 = 3.0-2.0*fNu2;
    o4v2 = 1.0-4.0*fNu2;
    f4v2 = 5.0-4.0*fNu2;
    mx1  = fMu1+t4v1*fMu2;
    mx2  = fMu2+t4v2*fMu1;
    
}


void GuzinaBiMatT::SetGeometry(double* src, double* dest, double* normal)
{
    MathOperationT::MemCopy(3, src, fSrc);
    MathOperationT::MemCopy(3, dest, fDest);
    MathOperationT::MemCopy(3, normal, fN);
    
    double nn=MathOperationT::VecNormal(3, fN);
    MathOperationT::VecScale(3, fN, 1/nn);
    
    //initialize the array of auxiliary coefficients
    c[0] = fMu1*o2v1*t4v2-fMu2*o2v2*t4v1;
    
    double s=fSrc[2];
    double z=fDest[2];
    if (s<0) {
        s=-s;
        z=-z;
    }
    
    if (z<0.0) {
        c[1] = fMu2/(2.0*mx1*mx2);
        c[2] = fMu2/(fMu1+fMu2);
        c[3] = fMu1*2.0*o1v1*t4v2+fMu2*2.0*o1v2*t4v1;
    }else{
        c[1] = 1.0/(8.0*o1v2*mx1*mx2);
        c[2] = 1.0/(2.0*(fMu1+fMu2));
        c[3] = 2.0*mx1*(fMu1-fMu2)*z*s;
        c[4] = mx1*(fMu1-fMu2);
        c[5] = fMu1*fMu1*t4v2*t4v2-fMu2*fMu2*t4v1*(5.0-12.0*fNu2+8.0*fNu2*fNu2)+2.0*fMu1*fMu2*o2v1*t4v2*o2v2;
        c[6] = fMu1*fMu1*t4v2+fMu2*fMu2*t4v1*(1.0-8.0*fNu2+8.0*fNu2*fNu2)-2.0*fMu1*fMu2*o2v1*t4v2*o2v2;
        c[7] = 2.0*fMu1*fMu1*t4v2*o2v2-2.0*fMu2*fMu2*t4v1*o2v2*o2v2-4.0*fMu1*fMu2*fNu2*t4v2*o2v1;
        c[8] = 4.0*o1v2*(fMu1*fMu1*t4v2-2.0*fMu2*fMu2*t4v1*o1v2+fMu1*fMu2*t4v2*o2v1);
        c[9] = fMu1*fMu1*t4v2-fMu2*fMu2*t4v1*(7.0-16.0*fNu2+8.0*fNu2*fNu2)+2.0*fMu1*fMu2*o2v1*t4v2*t2v2;
    }

}


void GuzinaBiMatT::Compute()
{
    double s=fSrc[2]; //depth of the source point
    double x=fDest[0]-fSrc[0];
    double y=fDest[1]-fSrc[1];
    double z=fDest[2]; //depth of the field point
    double r=sqrt(x*x+y*y);

    //Special treatment for source in the upper medium
    if (s<0) {
        s=-s;
        z=-z;
        y=-y;
    }
    
    double stheta, ctheta;
    if (r==0) {
        stheta=0;
        ctheta=1;
    }else{
        stheta=y/r;
        ctheta=x/r;
    }
    
    double O1_01, O2_01, G1_01, G2_01, G3_01;
    double O1_02, O2_02, G1_02, G2_02, G3_02;
    double O1_11, O2_11, G1_11, G2_11, G3_11;
    double O1_12, O2_12, G1_12, G2_12, G3_12;
    double O1_21, O2_21, G1_21, G2_21, G3_21;
    double O1_22, O2_22, G1_22, G2_22, G3_22;
    
    double DO1_0, DO2_0, DG1_0, DG2_0, DG3_0;
    double DO1_1, DO2_1, DG1_1, DG2_1, DG3_1;
    double DO1_2, DO2_2, DG1_2, DG2_2, DG3_2;
    
    
    GuzinaParameters(O1_01, O2_01, G1_01, G2_01, G3_01, DO1_0, DO2_0, DG1_0, DG2_0, DG3_0, 0, 1, r, z, s);
    GuzinaParameters(O1_02, O2_02, G1_02, G2_02, G3_02, DO1_0, DO2_0, DG1_0, DG2_0, DG3_0, 0, 2, r, z, s);
    GuzinaParameters(O1_11, O2_11, G1_11, G2_11, G3_11, DO1_1, DO2_1, DG1_1, DG2_1, DG3_1, 1, 1, r, z, s);
    GuzinaParameters(O1_12, O2_12, G1_12, G2_12, G3_12, DO1_1, DO2_1, DG1_1, DG2_1, DG3_1, 1, 2, r, z, s);
    GuzinaParameters(O1_21, O2_21, G1_21, G2_21, G3_21, DO1_2, DO2_2, DG1_2, DG2_2, DG3_2, 2, 1, r, z, s);
    GuzinaParameters(O1_22, O2_22, G1_22, G2_22, G3_22, DO1_2, DO2_2, DG1_2, DG2_2, DG3_2, 2, 2, r, z, s);
    
    double fac=1/(4.0*M_PI*fMu2);
    
    //compute displacements in cylindrical system
    double UX_r=fac*ctheta*(G2_01+G1_01+r*G2_21-r*G1_21);
    double UX_t=fac*(-stheta)*(G2_01+G1_01-r*G2_21+r*G1_21);
    double UX_z=2.0*fac*ctheta*r*O1_11;
    
    double UY_r=fac*stheta*(G2_01+G1_01+r*G2_21-r*G1_21);
    double UY_t=fac*(ctheta)*(G2_01+G1_01-r*G2_21+r*G1_21);
    double UY_z=2.0*fac*stheta*r*O1_11;
    
    double UZ_r=-2.0*fac*r*G3_11;
    double UZ_t=0.0;
    double UZ_z=2.0*fac*O2_01;
    
    
    //compute stresses in cylindrical system
    double Lb, Mb;
    if (z<0) {
        Lb=fLambda1;
        Mb=fMu1;
    }else{
        Lb=fLambda2;
        Mb=fMu2;
    }
    
    double SX[9], SY[9], SZ[9];
    
    SX[0*3+0]=2.0*fac*ctheta*(Lb*r*DO1_1-(Lb+2*Mb)*r*G1_12)-2.0*Mb*2.0*fac*ctheta*(G2_21-G1_21);
    SX[0*3+1]=2.0*fac*stheta*Mb*r*G2_12-2.0*Mb*2.0*fac*stheta*(G2_21-G1_21);
    SX[0*3+2]=fac*Mb*ctheta*(DG2_0+DG1_0+O1_02+r*DG2_2-r*DG1_2-r*O1_22);
    SX[1*3+0]=SX[0*3+1];
    SX[1*3+1]=2.0*fac*Lb*ctheta*(r*DO1_1-r*G1_12)+2.0*Mb*2.0*fac*ctheta*(G2_21-G1_21);
    SX[1*3+2]=-fac*Mb*stheta*(DG2_0+DG1_0+O1_02-r*DG2_2+r*DG1_2+r*O1_22);
    SX[2*3+0]=SX[0*3+2];
    SX[2*3+1]=SX[1*3+2];
    SX[2*3+2]=2.0*fac*ctheta*((Lb+2*Mb)*r*DO1_1-Lb*r*G1_12);
 
    SY[0*3+0]=2.0*fac*stheta*(Lb*r*DO1_1-(Lb+2*Mb)*r*G1_12)-2.0*Mb*2.0*fac*stheta*(G2_21-G1_21);
    SY[0*3+1]=2.0*fac*(-ctheta)*Mb*r*G2_12-2.0*Mb*2.0*fac*(-ctheta)*(G2_21-G1_21);
    SY[0*3+2]=fac*Mb*stheta*(DG2_0+DG1_0+O1_02+r*DG2_2-r*DG1_2-r*O1_22);
    SY[1*3+0]=SY[0*3+1];
    SY[1*3+1]=2.0*fac*Lb*stheta*(r*DO1_1-r*G1_12)+2.0*Mb*2.0*fac*stheta*(G2_21-G1_21);
    SY[1*3+2]=-fac*Mb*(-ctheta)*(DG2_0+DG1_0+O1_02-r*DG2_2+r*DG1_2+r*O1_22);
    SY[2*3+0]=SY[0*3+2];
    SY[2*3+1]=SY[1*3+2];
    SY[2*3+2]=2.0*fac*stheta*((Lb+2*Mb)*r*DO1_1-Lb*r*G1_12);
    
    SZ[0*3+0]=2.0*fac*(Lb*DO2_0-(Lb+2*Mb)*G3_02)+2.0*Mb*2.0*fac*G3_11;
    SZ[0*3+1]=0.0;
    SZ[0*3+2]=-2.0*fac*Mb*(r*DG3_1+r*O2_12);
    SZ[1*3+0]=SZ[0*3+1];
    SZ[1*3+1]=2.0*fac*Lb*(DO2_0-G3_02)-2.0*Mb*2.0*fac*G3_11;;
    SZ[1*3+2]=0.0;
    SZ[2*3+0]=SZ[0*3+2];
    SZ[2*3+1]=SZ[1*3+2];
    SZ[2*3+2]=2.0*fac*((Lb+2*Mb)*DO2_0-Lb*G3_02);
    
    
    //transfter to Cartesian system
    fDis[0*3+0] = UX_r*ctheta-UX_t*stheta;
    fDis[1*3+0] = UX_r*stheta+UX_t*ctheta;
    fDis[2*3+0] = UX_z;
    
    fDis[0*3+1] = UY_r*ctheta-UY_t*stheta;
    fDis[1*3+1] = UY_r*stheta+UY_t*ctheta;
    fDis[2*3+1] = UY_z;
    
    fDis[0*3+2] = UZ_r*ctheta-UZ_t*stheta;
    fDis[1*3+2] = UZ_r*stheta+UZ_t*ctheta;
    fDis[2*3+2] = UZ_z;
    
    double R[]={ctheta, -stheta, 0, stheta,  ctheta, 0,  0, 0, 1};
    double RT[]={ctheta, stheta, 0,  -stheta,  ctheta, 0,  0, 0, 1};
    
    double MM[9]; //temporary matrix for multiplications
    
    MathOperationT::multAB(3, 3, R, 3, 3, SX, MM);
    MathOperationT::multAB(3, 3, MM, 3, 3, RT, fStress_X);
    
    MathOperationT::multAB(3, 3, R, 3, 3, SY, MM);
    MathOperationT::multAB(3, 3, MM, 3, 3, RT, fStress_Y);
    
    MathOperationT::multAB(3, 3, R, 3, 3, SZ, MM);
    MathOperationT::multAB(3, 3, MM, 3, 3, RT, fStress_Z);
    
    
    if (fSrc[2]<0) {
        fDis[1*3+0] = -fDis[1*3+0];
        fDis[2*3+0] = -fDis[2*3+0];
        fDis[0*3+1] = -fDis[0*3+1];
        fDis[0*3+2] = -fDis[0*3+2];
        
        double R2[]={1, 0, 0, 0,  -1, 0,  0, 0, -1};
        double R2_negative[]={-1, 0, 0, 0,  1, 0,  0, 0, 1};
    
        MathOperationT::multAB(3, 3, R2, 3, 3, fStress_X, MM);
        MathOperationT::multAB(3, 3, MM, 3, 3, R2, fStress_X);
        
        MathOperationT::multAB(3, 3, R2_negative, 3, 3, fStress_Y, MM);
        MathOperationT::multAB(3, 3, MM, 3, 3, R2, fStress_Y);
        
        MathOperationT::multAB(3, 3, R2_negative, 3, 3, fStress_Z, MM);
        MathOperationT::multAB(3, 3, MM, 3, 3, R2, fStress_Z);
    }
    
    
    ComputeTraction();
}



void GuzinaBiMatT::GuzinaParameters(double &Omega1, double &Omega2, double &Gamma1, double &Gamma2, double &Gamma3,
                                    double &DO1,    double &DO2,    double &DG1,    double &DG2,    double &DG3,
                                    int n, int l, double r, double z, double s)
{
    double d1=std::abs(z-s);
    double d2=z+s;
    double d3=z-s;
    
    double szs=Sign(d3);
    
    double J1_nl=BesselInt(n, l, r, d1); //J1(n,l);
    double J1_nl_n1=BesselInt(n, l-1, r, d1); //J1(n, l-1);
    
    double J2_nl=BesselInt(n, l, r, d2); //J2(n,l);
    double J2_nl_n1=BesselInt(n, l-1, r, d2); //J2(n, l-1);
    double J2_nl_p1=BesselInt(n, l+1, r, d2); //J2(n, l+1);
    
    double J1_n1=BesselInt(n, 1, r, d1);
    double J1_n2=BesselInt(n, 2, r, d1);
    
    double J2_n1=BesselInt(n, 1, r, d2);
    double J2_n2=BesselInt(n, 2, r, d2);
    double J2_n3=BesselInt(n, 3, r, d2);
    
    if (z<0.0) { //observation point and the source are in different media
        Omega1 = c[1] * ((z*mx2-s*mx1)*J1_nl - c[0]*J1_nl_n1);
        Omega2 =-c[1] * ((z*mx2-s*mx1)*J1_nl - c[3]*J1_nl_n1);
        Gamma1 = c[1] * ((z*mx2-s*mx1)*J1_nl + c[3]*J1_nl_n1);
        Gamma2 = c[2] * J1_nl_n1;
        Gamma3 =-c[1] * ((z*mx2-s*mx1)*J1_nl + c[0]*J1_nl_n1);
        
        DO1 =  c[1] * ((z*mx2-s*mx1)*J1_n2 - (c[0]-mx2)*J1_n1);
        DO2 = -c[1] * ((z*mx2-s*mx1)*J1_n2 - (c[3]-mx2)*J1_n1);
        DG1 =  c[1] * ((z*mx2-s*mx1)*J1_n2 + (c[3]+mx2)*J1_n1);
        DG2 =  c[2] * J1_n1;
        DG3 = -c[1] * ((z*mx2-s*mx1)*J1_n2 + (c[0]+mx2)*J1_n1);
    }else{ //observation point and the source are in the same media
        Omega1 =  c[1] * ( mx1*mx2*d3*J1_nl + c[3]*J2_nl_p1 - c[4]*t4v2*d3*J2_nl - 4.0*fMu2*o1v2*c[0]*J2_nl_n1);
        Omega2 =  c[1] * ( mx1*mx2*(d1*J1_nl+t4v2*J1_nl_n1) - c[3]*J2_nl_p1 - c[4]*t4v2*d2*J2_nl - c[5]*J2_nl_n1);
        Gamma1 =  c[1] * (-mx1*mx2*(d1*J1_nl-t4v2*J1_nl_n1) - c[3]*J2_nl_p1 + c[4]*t4v2*d2*J2_nl - c[5]*J2_nl_n1);
        Gamma2 =  c[2] * ((fMu1+fMu2)*J1_nl_n1-(fMu1-fMu2)*J2_nl_n1);
        Gamma3 =  c[1] * (-mx1*mx2*d3*J1_nl + c[3]*J2_nl_p1 + c[4]*t4v2*d3*J2_nl - 4.0*fMu2*o1v2*c[0]*J2_nl_n1);
        
        DO1 =  c[1] * (-mx1*mx2*(d1*J1_n2-J1_n1) - c[3]*J2_n3 + c[4]*(z*t4v2-s*o4v2)*J2_n2 - c[6]*J2_n1);
        DO2 =  c[1] * (-mx1*mx2*szs*(d1*J1_n2+2.0*o2v2*J1_n1) + c[3]*J2_n3 + c[4]*(z*t4v2+s*o4v2)*J2_n2 + c[7]*J2_n1);
        DG1 =  c[1] * ( mx1*mx2*szs*(d1*J1_n2-4.0*o1v2*J1_n1) + c[3]*J2_n3 - c[4]*(z*t4v2+s*f4v2)*J2_n2 + c[8]*J2_n1);
        DG2 =  c[2] * (-(fMu1+fMu2)*szs*J1_n1+(fMu1-fMu2)*J2_n1);
        DG3 =  c[1] * ( mx1*mx2*(d1*J1_n2-J1_n1) - c[3]*J2_n3 - c[4]*(z*t4v2-s*f4v2)*J2_n2 + c[9]*J2_n1);
    }
    
}



double GuzinaBiMatT::BesselInt(int n, int l, double r, double d)
{
    double rs   = r*r;
    double ds   = d*d;
    double dist = sqrt(ds+rs);
    double coef = 1.0/(dist+d);
    double coef_n=pow(coef, n+0.0);
    
    double jx;
    
    //evaluation of jx(0,l), jx(1,l)/r, jx(2,l)/r^2

    switch (l) {
        case 0:
            jx = coef_n/dist;
            break;
            
        case 1:
            jx = coef_n/pow(dist, 3.0) * (n*dist+d);
            break;
        
        case 2:
            jx = coef_n/pow(dist, 5.0) * (n*dist*(n*dist+3.0*d) + 2.0*ds - rs);
            break;
            
        case 3:
            jx = coef_n/pow(dist, 7.0) * (n*dist*(pow(n*dist, 2)+11.0*ds-4.0*rs)
                                          + 3.0*d*(2.0*pow(n*dist, 2.0) + 2.0*ds - 3.0*rs));
            break;
            
        default:
            break;
    }
    

    //conversion jx(2,k)/r = (jx(2,k)/r^2) * r
    if (n==2)  jx = r * jx;
    
    return jx;
}


void GuzinaBiMatT::ComputeTraction()
{
    double traction[3];
    
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
}

const double* GuzinaBiMatT::GetDisplacement()
{
    return fDis;
}

const double* GuzinaBiMatT::GetTraction()
{
    return fTrac;
}

void GuzinaBiMatT::GetStress(double ** SX, double **SY, double ** SZ)
{
    *SX=fStress_X;
    *SY=fStress_Y;
    *SZ=fStress_Z;
}

double GuzinaBiMatT::Sign(double x)
{
    if (x>0) return 1.0;
    if (x<0) return -1.0;
    return 0.0;
}










