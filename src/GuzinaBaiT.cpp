//
//  GuzinaBaiT.cpp
//  GreenFunction
//
//  Created by Xiaoyong Bai on 5/13/16.
//
//

#include "GuzinaBaiT.h"
#include "MathOperationT.h"

#include "cmath"
#include "complex"
#include <algorithm>


using namespace GreenFunction;
using namespace std;

extern "C" void haldgreen4_(long* dyngreen, long* kdgl, double* omega, double* xyzs, double* xyzo,
                            complex<double>* ugreen, complex<double>* tdgreen, complex<double>* tsgreen);

GuzinaBaiT::GuzinaBaiT()
{
    //Allocate memberss
    fNumT=1;
    ft=new double[fNumT];
    fDis=new double[fNumT*9];
    fTrac=new double[fNumT*9];
    fSX=new double[fNumT*9];
    fSY=new double[fNumT*9];
    fSZ=new double[fNumT*9];
    
    
    fNumW=3000;
    fW=new double[fNumW];
    
    fDis_FD=new complex<double>[fNumW*9];
    fTrac_FD=new complex<double>[fNumW*9];
    
    fSX_FD=new complex<double>[fNumW*9];
    fSY_FD=new complex<double>[fNumW*9];
    fSZ_FD=new complex<double>[fNumW*9];
}


GuzinaBaiT::~GuzinaBaiT()
{
    //Deallocate
    if (ft) delete [] ft;
    if (fW) delete [] fW;
    
    if (fDis) delete [] fDis;
    if (fTrac) delete [] fTrac;
    
    if (fSX) delete [] fSX;
    if (fSY) delete [] fSY;
    if (fSZ) delete [] fSZ;
    
    if (fDis_FD) delete [] fDis_FD;
    if (fTrac_FD) delete [] fTrac_FD;
    
    if (fSX_FD) delete [] fSX_FD;
    if (fSY_FD) delete [] fSY_FD;
    if (fSZ_FD) delete [] fSZ_FD;

}

void GuzinaBaiT::SetGeometry(double* src, double* dest, double* normal)
{
    MathOperationT::MemCopy(3, src, fSrc);
    MathOperationT::MemCopy(3, dest, fDest);
    MathOperationT::MemCopy(3, normal, fNormal);
}


void GuzinaBaiT::SetBWidth(double T)
{
    fDT=T;
    
    double w_max=25.0/T;
    
    for (int i=0; i<fNumW; i++) {
        fW[i]=i*w_max/(fNumW-1);
    }
    
}


void GuzinaBaiT::SetTime(int num, double *t)
{
    if (fNumT != num) {
        delete [] ft;
        delete [] fDis;
        delete [] fTrac;
        delete [] fSX;
        delete [] fSY;
        delete [] fSZ;
        
        fNumT=num;
        ft=new double[fNumT];
        fDis=new double[fNumT*9];
        fTrac=new double[fNumT*9];
        
        fSX=new double[fNumT*9];
        fSY=new double[fNumT*9];
        fSZ=new double[fNumT*9];
    }    
    
    MathOperationT::MemCopy(fNumT, t, ft);
}


void GuzinaBaiT::Compute()
{
    long dyn=1;
    long kdgl=1;
    complex<double> ugreen[9], tdgreen[18], tsgreen[18];
    
    complex<double> *dis_temp, *trac_temp; //temporary pointers
    complex<double> *sx_temp, *sy_temp, *sz_temp;
    
    complex<double> SX_static[9], SY_static[9], SZ_static[9];
    complex<double> traction[3];

    //compute Green's function in frequency domain
    for (int i=0; i<fNumW; i++) {
        
        double ww=fW[i];
        if (i==0) {
            ww=0.5;
        }
        
        haldgreen4_(&dyn, &kdgl, &ww, fSrc, fDest, ugreen, tdgreen, tsgreen);
        
        dis_temp=fDis_FD+9*i;
        trac_temp=fTrac_FD+9*i;
        sx_temp=fSX_FD+9*i;
        sy_temp=fSY_FD+9*i;
        sz_temp=fSZ_FD+9*i;
        
        for (int cc=0; cc<9; cc++) {
            dis_temp[cc]=ugreen[cc];
        }
        
        sx_temp[0]=tdgreen[0]; sx_temp[1]=tdgreen[15]; sx_temp[2]=tdgreen[12];
        sx_temp[3]=tdgreen[15]; sx_temp[4]=tdgreen[3]; sx_temp[5]=tdgreen[9];
        sx_temp[6]=tdgreen[12]; sx_temp[7]=tdgreen[9]; sx_temp[8]=tdgreen[6];
        
        sx_temp[0]+=tsgreen[0]; sx_temp[1]+=tsgreen[15]; sx_temp[2]+=tsgreen[12];
        sx_temp[3]+=tsgreen[15]; sx_temp[4]+=tsgreen[3]; sx_temp[5]+=tsgreen[9];
        sx_temp[6]+=tsgreen[12]; sx_temp[7]+=tsgreen[9]; sx_temp[8]+=tsgreen[6];
        
        sy_temp[0]=tdgreen[1]; sy_temp[1]=tdgreen[16]; sy_temp[2]=tdgreen[13];
        sy_temp[3]=tdgreen[16]; sy_temp[4]=tdgreen[4]; sy_temp[5]=tdgreen[10];
        sy_temp[6]=tdgreen[13]; sy_temp[7]=tdgreen[10]; sy_temp[8]=tdgreen[7];
        
        sy_temp[0]+=tsgreen[1]; sy_temp[1]+=tsgreen[16]; sy_temp[2]+=tsgreen[13];
        sy_temp[3]+=tsgreen[16]; sy_temp[4]+=tsgreen[4]; sy_temp[5]+=tsgreen[10];
        sy_temp[6]+=tsgreen[13]; sy_temp[7]+=tsgreen[10]; sy_temp[8]+=tsgreen[7];
        
        sz_temp[0]=tdgreen[2]; sz_temp[1]=tdgreen[17]; sz_temp[2]=tdgreen[14];
        sz_temp[3]=tdgreen[17]; sz_temp[4]=tdgreen[5]; sz_temp[5]=tdgreen[11];
        sz_temp[6]=tdgreen[14]; sz_temp[7]=tdgreen[11]; sz_temp[8]=tdgreen[8];
        
        sz_temp[0]+=tsgreen[2]; sz_temp[1]+=tsgreen[17]; sz_temp[2]+=tsgreen[14];
        sz_temp[3]+=tsgreen[17]; sz_temp[4]+=tsgreen[5]; sz_temp[5]+=tsgreen[11];
        sz_temp[6]+=tsgreen[14]; sz_temp[7]+=tsgreen[11]; sz_temp[8]+=tsgreen[8];
        
        MathOperationT::multAB(3, 3, sx_temp, 3, 1, fNormal, traction);
        for (int i=0; i<3; i++) {
            trac_temp[i*3+0]=traction[i];
        }
        
        MathOperationT::multAB(3, 3, sy_temp, 3, 1, fNormal, traction);
        for (int i=0; i<3; i++) {
            trac_temp[i*3+1]=traction[i];
        }
        
        MathOperationT::multAB(3, 3, sz_temp, 3, 1, fNormal, traction);
        for (int i=0; i<3; i++) {
            trac_temp[i*3+2]=traction[i];
        }
        
        
        //obtain static part
        SX_static[0]=tsgreen[0]; SX_static[1]=tsgreen[15]; SX_static[2]=tsgreen[12];
        SX_static[3]=tsgreen[15]; SX_static[4]=tsgreen[3]; SX_static[5]=tsgreen[9];
        SX_static[6]=tsgreen[12]; SX_static[7]=tsgreen[9]; SX_static[8]=tsgreen[6];
        
        SY_static[0]=tsgreen[1]; SY_static[1]=tsgreen[16]; SY_static[2]=tsgreen[13];
        SY_static[3]=tsgreen[16]; SY_static[4]=tsgreen[4]; SY_static[5]=tsgreen[10];
        SY_static[6]=tsgreen[13]; SY_static[7]=tsgreen[10]; SY_static[8]=tsgreen[7];
        
        SZ_static[0]=tsgreen[2]; SZ_static[1]=tsgreen[17]; SZ_static[2]=tsgreen[14];
        SZ_static[3]=tsgreen[17]; SZ_static[4]=tsgreen[5]; SZ_static[5]=tsgreen[11];
        SZ_static[6]=tsgreen[14]; SZ_static[7]=tsgreen[11]; SZ_static[8]=tsgreen[8];
        
        MathOperationT::multAB(3, 3, SX_static, 3, 1, fNormal, traction);
        for (int i=0; i<3; i++) {
            fTrac_static[i*3+0]=real(traction[i]);
        }
        
        MathOperationT::multAB(3, 3, SY_static, 3, 1, fNormal, traction);
        for (int i=0; i<3; i++) {
            fTrac_static[i*3+1]=real(traction[i]);
        }
        
        MathOperationT::multAB(3, 3, SZ_static, 3, 1, fNormal, traction);
        for (int i=0; i<3; i++) {
            fTrac_static[i*3+2]=real(traction[i]);
        }
    }
    
    //Fourier Synthesis
    for (int it=0; it<fNumT; it++) {
        
        double* U_td=fDis+9*it; //pointer to time-domain displacement
        double* T_td=fTrac+9*it; //pointer to time-domain traction
        double* sx_td=fSX+9*it;
        double* sy_td=fSY+9*it;
        double* sz_td=fSZ+9*it;
        
        
        //initialization of U_td and T_td
        for (int ic=0; ic<9; ic++) {
            U_td[ic]=0;
            T_td[ic]=0;
            sx_td[ic]=0;
            sy_td[ic]=0;
            sz_td[ic]=0;
        }
        
        for (int iw=1; iw<fNumW; iw++) {
            const complex<double>* U_fd_0=fDis_FD+(iw-1)*9;
            const complex<double>* T_fd_0=fTrac_FD+(iw-1)*9;
            const complex<double>* sx_fd_0=fSX_FD+(iw-1)*9;
            const complex<double>* sy_fd_0=fSY_FD+(iw-1)*9;
            const complex<double>* sz_fd_0=fSZ_FD+(iw-1)*9;
            
            const complex<double>* U_fd_1=fDis_FD+iw*9;
            const complex<double>* T_fd_1=fTrac_FD+iw*9;
            const complex<double>* sx_fd_1=fSX_FD+iw*9;
            const complex<double>* sy_fd_1=fSY_FD+iw*9;
            const complex<double>* sz_fd_1=fSZ_FD+iw*9;
            
            complex<double> int_sub=Integral_Sub(fW[iw-1], fW[iw], ft[it]);
            
            //integration
            for (int ic=0; ic<9; ic++) {
                U_td[ic]+=real((U_fd_0[ic]+U_fd_1[ic])*int_sub)/sqrt(2*M_PI);
                T_td[ic]+=real((T_fd_0[ic]+T_fd_1[ic])*int_sub)/sqrt(2*M_PI);
                sx_td[ic]+=real((sx_fd_0[ic]+sx_fd_1[ic])*int_sub)/sqrt(2*M_PI);
                sy_td[ic]+=real((sy_fd_0[ic]+sy_fd_1[ic])*int_sub)/sqrt(2*M_PI);
                sz_td[ic]+=real((sz_fd_0[ic]+sz_fd_1[ic])*int_sub)/sqrt(2*M_PI);
            }
            
        }//end of loop over iw
    }//end of loop over it

}


complex<double> GuzinaBaiT::FTBSpline(double w)
{
    if (w==0){
        return fDT/(4.0*sqrt(2.0*M_PI));
    }
    else{
        double B3=pow(fDT, 3);
        double w4=pow(w,4);
        complex<double> I(0.0, 1.0);
        
        return sqrt(2.0/M_PI)*32.0*exp(-I*w*fDT)*pow(exp(I*w*fDT/4.0)-1.0, 4)/(B3*w4);
    }
}


complex<double> GuzinaBaiT::Integral_Sub(double w1, double w2, double t)
{
    int num_w=20;
    double dw=(w2-w1)/(num_w+0.0);
    
    complex<double> I(0.0, 1.0);
    complex<double> result=0.0;
    
    for (int i=0; i<num_w; i++) {
        double w_temp=w1+dw*(i+0.5);
        result+=FTBSpline(w_temp)*exp(I*w_temp*t)*dw;
    }
    
    return result;
}


void GuzinaBaiT::GetDisplacement(int& num, double** displacement)
{
    num=fNumT;
    *displacement=fDis;
}


void GuzinaBaiT::GetTraction(int& num, double** trac)
{
    num=fNumT;
    *trac=fTrac;
}


double* GuzinaBaiT::GetStaticTraction()
{
    return fTrac_static;
}

void GuzinaBaiT::GetStress(int& num, double ** SX, double ** SY, double ** SZ)
{
    num=fNumT;
    
    *SX=fSX;
    *SY=fSY;
    *SZ=fSZ;
}







