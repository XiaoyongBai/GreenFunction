//
//  PakBaiT.cpp
//  GreenFunction
//
//  Created by Xiaoyong Bai on 5/13/16.
//
//

#include "PakBaiT.h"
#include "MathOperationT.h"

#include "cmath"
#include "complex"
#include <algorithm>


using namespace GreenFunction;
using namespace std;

PakBaiT::PakBaiT()
{
    //Allocate memberss
    fNumT=1;
    ft=new double[fNumT];
    fDis=new double[fNumT*9];
    fTrac=new double[fNumT*9];
    
    fNumW=1000;
    fW=new double[fNumW];
    fPak=new PakT[fNumW];
    
    fFTB=new complex<double> [fNumW];
}


PakBaiT::~PakBaiT()
{
    //Deallocate
    delete [] fDis;
    delete [] fTrac;
    
    delete [] fPak;
}


void PakBaiT::SetMaterial(double E, double nu, double rho)
{
    for (int i=0; i<fNumW; i++) {
        fPak[i].SetMaterial(E, nu, rho);
    }
}


void PakBaiT::SetGeometry(double* src, double* dest, double* normal)
{
    double nn=MathOperationT::VecNormal(3, normal);
    MathOperationT::VecScale(3, normal, 1/nn);
    
    for (int i=0; i<fNumW; i++) {
        fPak[i].SetGeometry(src, dest, normal);
    }
}


void PakBaiT::SetBWidth(double T)
{
    fDT=T;
    
    double w_max=20.0/T;
    
    for (int i=0; i<fNumW; i++) {
        fW[i]=i*w_max/(fNumW-1);
    }
    
}


void PakBaiT::SetTime(int num, double *t)
{
    if (fNumT != num) {
        delete [] ft;
        delete [] fDis;
        delete [] fTrac;
        
        fNumT=num;
        ft=new double[fNumT];
        fDis=new double[fNumT*9];
        fTrac=new double[fNumT*9];
    }    
    
    MathOperationT::MemCopy(fNumT, t, ft);
}


void PakBaiT::Compute()
{
    complex<double> I(0.0, 1.0);
    
    //compute Green's function in frequency domain
    for (int i=0; i<fNumW; i++) {
        fPak[i].Compute(fW[i]);
        
        fFTB[i]=FTBSpline(fW[i]);
    }
    
    double dw=fW[fNumW-1]/fNumW;
    
    //Fourier Synthesis
    for (int it=0; it<fNumT; it++) {
        
        double* U_td=fDis+9*it; //pointer to time-domain displacement
        double* T_td=fTrac+9*it; //pointer to time-domain traction
        
        //initialization of U_td and T_td
        for (int ic=0; ic<9; ic++) {
            U_td[ic]=0;
            T_td[ic]=0;
        }
        
        for (int iw=0; iw<fNumW; iw++) {
            
            const complex<double>* U_fd=fPak[iw].GetDisplacement();
            const complex<double>* T_fd=fPak[iw].GetTraction();
            
            //integration
            for (int ic=0; ic<9; ic++) {
                
                complex<double> EXP=exp(I*fW[iw]*ft[it]);
                
                U_td[ic]+=2*dw*real(U_fd[ic]*fFTB[iw]*EXP)/sqrt(2*M_PI);
                T_td[ic]+=2*dw*real(T_fd[ic]*fFTB[iw]*EXP)/sqrt(2*M_PI);
            }
            
            
        }//end of loop over iw
    }//end of loop over it

}


complex<double> PakBaiT::FTBSpline(double w)
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


void PakBaiT::GetDisplacement(int& num, double** displacement)
{
    num=fNumT;
    *displacement=fDis;
}


void PakBaiT::GetTraction(int& num, double** trac)
{
    num=fNumT;
    *trac=fTrac;
}











