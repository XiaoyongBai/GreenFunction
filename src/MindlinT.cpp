//
//  MindlinT.cpp
//  GreenFunction
//
//  Created by Xiaoyong Bai on 5/13/16.
//
//

#include "MindlinT.h"
#include "MathOperationT.h"
#include "cmath"

using namespace GreenFunction;


MindlinT::MindlinT()
{
    //do nothing
    fRotation[0]= 0; fRotation[1]= 1; fRotation[2]=0;
    fRotation[3]=-1; fRotation[4]= 0; fRotation[5]=0;
    fRotation[6]= 0; fRotation[7]= 0; fRotation[8]=1;
    
    fRotationT[0]= 0; fRotationT[1]=-1; fRotationT[2]=0;
    fRotationT[3]= 1; fRotationT[4]= 0; fRotationT[5]=0;
    fRotationT[6]= 0; fRotationT[7]= 0; fRotationT[8]=1;
}

MindlinT::~MindlinT()
{
    //do nothing
}

void MindlinT::SetMaterial(double E, double nu)
{
    fE=E;
    fNu=nu;
    
    fMu=fE/(2*(1+nu));
}

void MindlinT::SetGeometry(double* src, double* dest, double* normal)
{
    MathOperationT::MemCopy(3, src, fSrc);
    MathOperationT::MemCopy(3, dest, fDest);
    MathOperationT::MemCopy(3, normal, fN);
    
    double RV[3]={fDest[0]-fSrc[0], fDest[1]-fSrc[1], fDest[2]-fSrc[2]};
    fR1=MathOperationT::VecNormal(3, RV);
    
    RV[2]=fDest[2]+fSrc[2];
    fR2=MathOperationT::VecNormal(3, RV);

    double nn=MathOperationT::VecNormal(3, fN);
    MathOperationT::VecScale(3, fN, 1.0/nn);
}


void MindlinT::Compute()
{
    ComputeDisplacement();
    ComputeStress();
    ComputeTraction();
}


void MindlinT::ComputeDisplacement()
{
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
    
    double factor=(16*M_PI*fMu*(1-fNu));
    
    ////////////////////////////////////////////
    //Displacement when load in X-direction
    //Mindlin_1936
    double ux_x, uy_x, uz_x;
    
    ux_x=(3-4*fNu)/fR1 + pow(x,2)/pow(fR1,3); //full-space part
    ux_x=ux_x + 1/fR2 + (3-4*fNu)*pow(x,2)/pow(fR2,3) + (2*s*z/pow(fR2,3))*(1-3*pow(x,2)/pow(fR2,2))
                      +(4*(1-fNu)*(1-2*fNu)/(fR2+z+s))*(1-pow(x, 2)/(fR2*(fR2+z+s))); //complementary part
    ux_x=ux_x/factor;
    
    uy_x=     1/pow(fR1,3); //full-space part
    uy_x=uy_x+(3-4*fNu)/pow(fR2,3)-6*s*z/pow(fR2,5)-4*(1-fNu)*(1-2*fNu)/(fR2*pow(fR2+z+s,2)); //complementary part
    uy_x=uy_x*x*y/factor;
    
    uz_x=     (z-s)/pow(fR1,3); //full-space part
    uz_x=uz_x+(3-4*fNu)*(z-s)/pow(fR2,3)-6*s*z*(z+s)/pow(fR2,5)+4*(1-fNu)*(1-2*fNu)/(fR2*(fR2+z+s));
    uz_x=uz_x*x/factor;
    
    ////////////////////////////////////////////
    //Displacement when load in Y-direction
    double ux_y, uy_y, uz_y;
    
    double y_p=y; //coordinate at rotated coordinate system
    y=-x;
    x=y_p;
    
    ux_y=       (3-4*fNu)/fR1 + pow(x,2)/pow(fR1,3); //full-space part
    ux_y=ux_y + 1/fR2 + (3-4*fNu)*pow(x,2)/pow(fR2,3) + (2*s*z/pow(fR2,3))*(1-3*pow(x,2)/pow(fR2,2))
              +(4*(1-fNu)*(1-2*fNu)/(fR2+z+s))*(1-pow(x, 2)/(fR2*(fR2+z+s))); //complementary part
    ux_y=ux_y/factor;
    
    uy_y=     1/pow(fR1,3); //full-space part
    uy_y=uy_y+(3-4*fNu)/pow(fR2,3)-6*s*z/pow(fR2,5)-4*(1-fNu)*(1-2*fNu)/(fR2*pow(fR2+z+s,2)); //complementary part
    uy_y=uy_y*x*y/factor;
    
    uz_y=     (z-s)/pow(fR1,3); //full-space part
    uz_y=uz_y+(3-4*fNu)*(z-s)/pow(fR2,3)-6*s*z*(z+s)/pow(fR2,5)+4*(1-fNu)*(1-2*fNu)/(fR2*(fR2+z+s));
    uz_y=uz_y*x/factor;
    
    ////////////////////////////////////////////
    //Displacement when load in Z-direction
    double ur_z, uz_z;
    
    ur_z=   (z-s)/pow(fR1,3);     //full-space part
    ur_z=ur_z+(3-4*fNu)*(z-s)/pow(fR2,3)-4*(1-fNu)*(1-2*fNu)/(fR2*(fR2+z+s))+6*z*s*(z+s)/pow(fR2,5); //complementary part
    ur_z=ur_z*r/factor;
    
    uz_z=      (3-4*fNu)/fR1+pow((z-s),2)/pow(fR1,3); //full-space part
    uz_z=uz_z+ (8*pow(1-fNu,2)-(3-4*fNu))/fR2 + ((3-4*fNu)*pow(z+s,2)-2*z*s)/pow(fR2,3)
             + 6*z*s*pow(z+s,2)/pow(fR2,5); //complementary part
    uz_z=uz_z/factor;
    
    
    
    fDis[0*3+0]=ux_x;
    fDis[1*3+0]=uy_x;
    fDis[2*3+0]=uz_x;
    
    fDis[0*3+1]=-uy_y;
    fDis[1*3+1]=ux_y;
    fDis[2*3+1]=uz_y;
    
    fDis[0*3+2]=ur_z*ctheta;
    fDis[1*3+2]=ur_z*stheta;
    fDis[2*3+2]=uz_z;
}


void MindlinT::ComputeStress()
{
    double s=fSrc[2]; //depth of the source point
    double x=fDest[0]-fSrc[0];
    double y=fDest[1]-fSrc[1];
    double z=fDest[2]; //depth of the field point
    
    double factor=1/(8*M_PI*(1-fNu));
    
    double x2=pow(x,2);
    double y2=pow(y,2);
    double R1_3=pow(fR1,3);
    double R1_5=pow(fR1,5);
    double R2_2=pow(fR2,2);
    double R2_3=pow(fR2,3);
    double R2_5=pow(fR2,5);
    double R2_7=pow(fR2,7);
    
    ////////////////////////////////////////////
    //Stress when load in Z-direction
    fStress_Z[0*3+0]=                (1-2*fNu)*(z-s)/R1_3-3*x2*(z-s)/R1_5;
    fStress_Z[0*3+0]=fStress_Z[0*3+0]+(1-2*fNu)*(3*(z-s)-4*fNu*(z+s))/R2_3;
    fStress_Z[0*3+0]=fStress_Z[0*3+0]-(3*(3-4*fNu)*x2*(z-s)-6*s*(z+s)*((1-2*fNu)*z-2*fNu*s))/R2_5;
    fStress_Z[0*3+0]=fStress_Z[0*3+0]-30*s*x2*z*(z+s)/R2_7;
    fStress_Z[0*3+0]=fStress_Z[0*3+0]-(4*(1-fNu)*(1-2*fNu)/(fR2*(fR2+z+s)))*(1-x2/(fR2*(fR2+z+s))-x2/R2_2);
    
    fStress_Z[1*3+1]=                 (1-2*fNu)*(z-s)/R1_3-3*y2*(z-s)/R1_5;
    fStress_Z[1*3+1]=fStress_Z[1*3+1]+(1-2*fNu)*(3*(z-s)-4*fNu*(z+s))/R2_3;
    fStress_Z[1*3+1]=fStress_Z[1*3+1]-(3*(3-4*fNu)*y2*(z-s)-6*s*(z+s)*((1-2*fNu)*z-2*fNu*s))/R2_5;
    fStress_Z[1*3+1]=fStress_Z[1*3+1]-30*s*y2*z*(z+s)/R2_7;
    fStress_Z[1*3+1]=fStress_Z[1*3+1]-(4*(1-fNu)*(1-2*fNu)/(fR2*(fR2+z+s)))*(1-y2/(fR2*(fR2+z+s))-y2/R2_2);
    
    fStress_Z[2*3+2]=                -(1-2*fNu)*(z-s)/R1_3-3*pow(z-s,3)/R1_5;
    fStress_Z[2*3+2]=fStress_Z[2*3+2]+(1-2*fNu)*(z-s)/R2_3;
    fStress_Z[2*3+2]=fStress_Z[2*3+2]-(3*(3-4*fNu)*z*pow(z+s,2)-3*s*(z+s)*(5*z-s))/R2_5;
    fStress_Z[2*3+2]=fStress_Z[2*3+2]-30*s*z*pow(s+z,3)/R2_7;
    
    fStress_Z[0*3+1]=                -3*(z-s)/R1_5-3*(3-4*fNu)*(z-s)/R2_5;
    fStress_Z[0*3+1]=fStress_Z[0*3+1]+(4*(1-fNu)*(1-2*fNu)/(R2_2*(fR2+s+z)))*(1/(fR2+s+z)+1/fR2)-30*s*z*(s+z)/R2_7;
    fStress_Z[0*3+1]=fStress_Z[0*3+1]*x*y;
    
    fStress_Z[0*3+2]=                -(1-2*fNu)/R1_3-3*pow(z-s,2)/R1_5;
    fStress_Z[0*3+2]=fStress_Z[0*3+2]+(1-2*fNu)/R2_3-(3*(3-4*fNu)*z*(z+s)-3*s*(3*z+s))/R2_5-30*s*z*pow(s+z,2)/R2_7;
    fStress_Z[0*3+2]=fStress_Z[0*3+2]*x;
    
    fStress_Z[1*3+2]=                -(1-2*fNu)/R1_3-3*pow(z-s,2)/R1_5;
    fStress_Z[1*3+2]=fStress_Z[1*3+2]+(1-2*fNu)/R2_3-(3*(3-4*fNu)*z*(z+s)-3*s*(3*z+s))/R2_5-30*s*z*pow(z+s,2)/R2_7;
    fStress_Z[1*3+2]=fStress_Z[1*3+2]*y;
    
    fStress_Z[1*3+0]=fStress_Z[0*3+1];
    fStress_Z[2*3+0]=fStress_Z[0*3+2];
    fStress_Z[2*3+1]=fStress_Z[1*3+2];

    MathOperationT::VecScale(9, fStress_Z, factor);
    
    ////////////////////////////////////////////
    //Stress when load in X-direction
    fStress_X[0*3+0]=                -(1-2*fNu)/R1_3-3*x2/R1_5;
    fStress_X[0*3+0]=fStress_X[0*3+0]+(1-2*fNu)*(5-4*fNu)/R2_3-3*(3-4*fNu)*x2/R2_5;
    fStress_X[0*3+0]=fStress_X[0*3+0]-(4*(1-fNu)*(1-2*fNu)/(fR2*pow(fR2+z+s,2))*(3-x2*(3*fR2+z+s)/(R2_2*(fR2+z+s))));
    fStress_X[0*3+0]=fStress_X[0*3+0]+(6*s/R2_5)*(3*s-(3-2*fNu)*(z+s)+5*x2*z/R2_2);
    fStress_X[0*3+0]=fStress_X[0*3+0]*x;

    fStress_X[1*3+1]=                 (1-2*fNu)/R1_3-3*y2/R1_5+(1-2*fNu)*(3-4*fNu)/R2_3-3*(3-4*fNu)*y2/R2_5;
    fStress_X[1*3+1]=fStress_X[1*3+1]-4*(1-fNu)*(1-2*fNu)/(fR2*pow(fR2+s+z,2))*(1-y2*(3*fR2+z+s)/(R2_2*(fR2+z+s)));
    fStress_X[1*3+1]=fStress_X[1*3+1]+(6*s/R2_5)*(s-(1-2*fNu)*(z+s)+5*y2*z/R2_2);
    fStress_X[1*3+1]=fStress_X[1*3+1]*x;
    
    fStress_X[2*3+2]=                 (1-2*fNu)/R1_3-3*pow(z-s,2)/R1_5-(1-2*fNu)/R2_3-3*(3-4*fNu)*pow(z+s,2)/R2_5;
    fStress_X[2*3+2]=fStress_X[2*3+2]+(6*s/R2_5)*(s+(1-2*fNu)*(z+s)+5*z*pow(z+s,2)/R2_2);
    fStress_X[2*3+2]=fStress_X[2*3+2]*x;
    
    fStress_X[0*3+1]=                -(1-2*fNu)/R1_3-3*x2/R1_5+(1-2*fNu)/R2_3-3*(3-4*fNu)*x2/R2_5;
    fStress_X[0*3+1]=fStress_X[0*3+1]-4*(1-fNu)*(1-2*fNu)/(fR2*pow(fR2+z+s,2))*(1-x2*(3*fR2+s+z)/(R2_2*(fR2+z+s)));
    fStress_X[0*3+1]=fStress_X[0*3+1]-(6*z*s/R2_5)*(1-5*x2/R2_2);
    fStress_X[0*3+1]=fStress_X[0*3+1]*y;
    
    fStress_X[0*3+2]=                -(1-2*fNu)*(z-s)/R1_3-3*x2*(z-s)/R1_5;
    fStress_X[0*3+2]=fStress_X[0*3+2]+(1-2*fNu)*(z-s)/R2_3-3*(3-4*fNu)*x2*(z+s)/R2_5;
    fStress_X[0*3+2]=fStress_X[0*3+2]-(6*s/R2_5)*(z*(z+s)-(1-2*fNu)*x2-5*x2*z*(z+s)/R2_2);
    
    fStress_X[1*3+2]=               -3*(z-s)/R1_5-3*(3-4*fNu)*(z+s)/R2_5+6*s/R2_5*(1-2*fNu+5*z*(z+s)/R2_2);
    fStress_X[1*3+2]=fStress_X[1*3+2]*x*y;
    
    fStress_X[1*3+0]=fStress_X[0*3+1];
    fStress_X[2*3+0]=fStress_X[0*3+2];
    fStress_X[2*3+1]=fStress_X[1*3+2];
    
    MathOperationT::VecScale(9, fStress_X, factor);
    
    ////////////////////////////////////////////
    //Stress when load in Y-direction
    double y_p=y; //coordinate at rotated coordinate system
    y=-x;
    x=y_p;
    
    x2=x*x;
    y2=y*y;
    
    double Stress_Y[9];
    Stress_Y[0*3+0]=                -(1-2*fNu)/R1_3-3*x2/R1_5;
    Stress_Y[0*3+0]=Stress_Y[0*3+0]+(1-2*fNu)*(5-4*fNu)/R2_3-3*(3-4*fNu)*x2/R2_5;
    Stress_Y[0*3+0]=Stress_Y[0*3+0]-(4*(1-fNu)*(1-2*fNu)/(fR2*pow(fR2+z+s,2))*(3-x2*(3*fR2+z+s)/(R2_2*(fR2+z+s))));
    Stress_Y[0*3+0]=Stress_Y[0*3+0]+(6*s/R2_5)*(3*s-(3-2*fNu)*(z+s)+5*x2*z/R2_2);
    Stress_Y[0*3+0]=Stress_Y[0*3+0]*x;
    
    Stress_Y[1*3+1]=                 (1-2*fNu)/R1_3-3*y2/R1_5+(1-2*fNu)*(3-4*fNu)/R2_3-3*(3-4*fNu)*y2/R2_5;
    Stress_Y[1*3+1]=Stress_Y[1*3+1]-4*(1-fNu)*(1-2*fNu)/(fR2*pow(fR2+s+z,2))*(1-y2*(3*fR2+z+s)/(R2_2*(fR2+z+s)));
    Stress_Y[1*3+1]=Stress_Y[1*3+1]+(6*s/R2_5)*(s-(1-2*fNu)*(z+s)+5*y2*z/R2_2);
    Stress_Y[1*3+1]=Stress_Y[1*3+1]*x;
    
    Stress_Y[2*3+2]=                 (1-2*fNu)/R1_3-3*pow(z-s,2)/R1_5-(1-2*fNu)/R2_3-3*(3-4*fNu)*pow(z+s,2)/R2_5;
    Stress_Y[2*3+2]=Stress_Y[2*3+2]+(6*s/R2_5)*(s+(1-2*fNu)*(z+s)+5*z*pow(z+s,2)/R2_2);
    Stress_Y[2*3+2]=Stress_Y[2*3+2]*x;
    
    Stress_Y[0*3+1]=                -(1-2*fNu)/R1_3-3*x2/R1_5+(1-2*fNu)/R2_3-3*(3-4*fNu)*x2/R2_5;
    Stress_Y[0*3+1]=Stress_Y[0*3+1]-4*(1-fNu)*(1-2*fNu)/(fR2*pow(fR2+z+s,2))*(1-x2*(3*fR2+s+z)/(R2_2*(fR2+z+s)));
    Stress_Y[0*3+1]=Stress_Y[0*3+1]-(6*z*s/R2_5)*(1-5*x2/R2_2);
    Stress_Y[0*3+1]=Stress_Y[0*3+1]*y;
    
    Stress_Y[0*3+2]=                -(1-2*fNu)*(z-s)/R1_3-3*x2*(z-s)/R1_5;
    Stress_Y[0*3+2]=Stress_Y[0*3+2]+(1-2*fNu)*(z-s)/R2_3-3*(3-4*fNu)*x2*(z+s)/R2_5;
    Stress_Y[0*3+2]=Stress_Y[0*3+2]-(6*s/R2_5)*(z*(z+s)-(1-2*fNu)*x2-5*x2*z*(z+s)/R2_2);
    
    Stress_Y[1*3+2]=               -3*(z-s)/R1_5-3*(3-4*fNu)*(z+s)/R2_5+6*s/R2_5*(1-2*fNu+5*z*(z+s)/R2_2);
    Stress_Y[1*3+2]=Stress_Y[1*3+2]*x*y;
    
    Stress_Y[1*3+0]=Stress_Y[0*3+1];
    Stress_Y[2*3+0]=Stress_Y[0*3+2];
    Stress_Y[2*3+1]=Stress_Y[1*3+2];
    
    double MM[9];
    
    MathOperationT::multAB(3, 3, fRotationT, 3, 3, Stress_Y, MM);
    MathOperationT::multAB(3, 3, MM, 3, 3, fRotation, fStress_Y);
    MathOperationT::VecScale(9, fStress_Y, factor);
}


void MindlinT::ComputeTraction()
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

const double* MindlinT::GetDisplacement()
{
    return fDis;
}

const double* MindlinT::GetTraction()
{
    return fTrac;
}

void MindlinT::GetStress(double ** SX, double **SY, double ** SZ)
{
    *SX=fStress_X;
    *SY=fStress_Y;
    *SZ=fStress_Z;
}















