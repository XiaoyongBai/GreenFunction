//
//  KauselT.cpp
//  GreenFunction
//
//  Created by Xiaoyong Bai on 06/15/17.
//
//

#include "KauselT.h"
#include "MathOperationT.h"
#include "cmath"


using namespace GreenFunction;

KauselT::KauselT()
{
    fNumT_BSpline=0;
    fNumT_Heavi=0;
    
    ft_Heavi=NULL;
    ft_BSpline=NULL;
    
    fDis_Heavi=NULL;
    fTrac=NULL;
    
    fDis_BSpline=NULL;
    
    fDT=0.0;
}

KauselT::~KauselT()
{
    if (ft_Heavi) delete [] ft_Heavi;
    if (ft_BSpline) delete [] ft_BSpline;
    
    if (fDis_Heavi) delete [] fDis_Heavi;
    if (fTrac) delete [] fTrac;
    
    if (fDis_BSpline) delete [] fDis_BSpline;
}

void KauselT::SetBWidth(double T)
{
    fDT=T;
}

void KauselT::SetTime(int num, double *t)
{
    if (fNumT_BSpline != num) {
        delete [] ft_Heavi;
        delete [] ft_BSpline;
        
        delete [] fDis_Heavi;
        delete [] fTrac;
        delete [] fDis_BSpline;
        
        fNumT_BSpline=num;
        ft_BSpline=new double[fNumT_BSpline];
        fTrac=new double[fNumT_BSpline*9];
        fDis_BSpline=new double[fNumT_BSpline*9];
        
        fNumT_Heavi=100*fNumT_BSpline;
        ft_Heavi=new double[fNumT_Heavi];
        fDis_Heavi=new double[fNumT_Heavi*9];
    }
    
    MathOperationT::MemCopy(fNumT_BSpline, t, ft_BSpline);
    
    //create time slots for Heaviside loading
    double time_upper=ft_BSpline[0];
    for (int i=0; i<fNumT_BSpline; i++) {
        if (ft_BSpline[i]>time_upper) {
            time_upper=ft_BSpline[i];
        }
    }
    
    double t_inc=time_upper/(fNumT_Heavi-1);
    
    for (int i=0; i<fNumT_Heavi; i++) {
        ft_Heavi[i]=i*t_inc;
    }
    

    for (int i=0; i<fNumT_BSpline*9; i++) {
        fTrac[i]=0;
        fDis_BSpline[i]=0;
    }
    
    for (int i=0; i<fNumT_Heavi*9; i++) {
        fDis_Heavi[i]=0;
    }
    
    //MathOperationT::PrintVector(fNumT_BSpline, ft_BSpline, "time to be computed");
    //MathOperationT::PrintVector(fNumT_Heavi, ft_Heavi, "Time for Heaviside loading");
}


void KauselT::SetGeometry(double *src, double *dest)
{
    if (src[2] != 0 || dest[2] != 0)
    {
        cout<< "Kausel's solution is only valid for surface source and receivers!!!"<<endl;
    }
    
    fX=dest[0]-src[0];
    fY=dest[1]-src[1];
    
    fr=sqrt(fX*fX+fY*fY);
    
    fsTheta=fY/fr;
    fcTheta=fX/fr;
    
}



void KauselT::SetMaterial(double E, double nu, double density)
{
    fE=E;
    fNu=nu;
    fRho=density;
    
    fMu=fE/(2*(1+nu));
    
    double lamda=fMu*2*nu/(1-2*nu);
    Cp=sqrt((lamda+2*fMu)/fRho);
    Cs=sqrt(fMu/fRho);
    
    L_Mu= lamda/fMu;
    cscd2=(1-2*nu)/(2*(1-nu)); //(C2/C1)^2
    fa=sqrt(cscd2);
    
    RayleighWave();
}

void KauselT::RayleighWave()
{
    //compute coefficients of Rayleigh equations
    double a=16*(1-cscd2);
    double b=-8*(3-2*cscd2);
    double c=8.0;
    double d=-1.0;
    
    
    //compute intermidient parameters
    double e=(c-b*b/(3*a))/a;
    double f=(d+2*b*b*b/(27*a*a)-b*c/(3*a))/a;
    
    double D=f*f+4*pow(e, 3.0)/27.0;
    
    std::complex<double> w1=(-f+sqrt(std::complex<double>(D)))/2.0;
    
    std::complex<double> I(0.0, 1.0);
    
    std::complex<double> z1=pow(w1, 1/3.0);
    std::complex<double> z2=z1*exp(I*2.0*M_PI/3.0);
    std::complex<double> z3=z1*exp(I*4.0*M_PI/3.0);
    
    std::complex<double>  s=-e/3.0;
    std::complex<double> y1=z1+s/z1;
    std::complex<double> y2=z2+s/z2;
    std::complex<double> y3=z3+s/z3;
    
    std::complex<double> K1_s_complex=y1-b/(3.0*a);
    K2_s=y2-b/(3.0*a);
    K3_s=y3-b/(3.0*a);

    K1_s=real(K1_s_complex);
    
    if (abs(imag(K2_s))<1e-10) {
        K2_s=real(K2_s);
    }
    
    if (abs(imag(K3_s))<1e-10) {
        K3_s=real(K3_s);
    }
    
    gamma=sqrt(K1_s);
}


void KauselT::Compute()
{
    //Auxilliary parameters
    std::complex<double> D1=(K1_s-K2_s)*(K1_s-K3_s);
    std::complex<double> D2=(K2_s-K1_s)*(K2_s-K3_s);
    std::complex<double> D3=(K3_s-K1_s)*(K3_s-K2_s);
    
    
    std::complex<double> A1=pow(K1_s-0.5,2)*sqrt(std::complex<double>(cscd2-K1_s))/D1;
    std::complex<double> A2=pow(K2_s-0.5,2)*sqrt(cscd2-K2_s)/D2;
    std::complex<double> A3=pow(K3_s-0.5,2)*sqrt(cscd2-K3_s)/D3;
    
    std::complex<double> B1=(1.0-2.0*K1_s)*(1.0-K1_s)/D1;
    std::complex<double> B2=(1.0-2.0*K2_s)*(1.0-K2_s)/D2;
    std::complex<double> B3=(1.0-2.0*K3_s)*(1.0-K3_s)/D3;
    
    std::complex<double> C1=(1.0-K1_s)*sqrt(std::complex<double>(cscd2-K1_s))/D1;
    std::complex<double> C2=(1.0-K2_s)*sqrt(cscd2-K2_s)/D2;
    std::complex<double> C3=(1.0-K3_s)*sqrt(cscd2-K3_s)/D3;
    
    std::complex<double> m1=(1.0-cscd2)/(cscd2-K1_s);
    std::complex<double> m2=(1.0-cscd2)/(cscd2-K2_s);
    std::complex<double> m3=(1.0-cscd2)/(cscd2-K3_s);
    std::complex<double> D=pow(2.0*K1_s-1.0,3)/(8.0*(1.0-cscd2)*pow(K1_s,3)-4.0*K1_s+1.0);
    
    
    double t0=fr/Cs; // reference time
    
    
    std::complex<double> U_H[3];
    std::complex<double> U_X[3];
    std::complex<double> U_Y[3];
    std::complex<double> U_Z[3];

    
    for (int i=0; i<fNumT_Heavi; i++) {
        double tau=ft_Heavi[i]/t0;
        double tau2=tau*tau;
        
        std::complex<double> n=sqrt((tau2-cscd2)/(1.0-cscd2));
        std::complex<double> n2=n*n;
        
        for (int j=0; j<3; j++) {
            U_H[j]=0;
            U_X[j]=0;
            U_Y[j]=0;
            U_Z[j]=0;
        }
        
        if (tau <= fa) {
        
            // do nothing
        
        }else if (tau>fa && tau<1){
            
            U_H[0]=C1/sqrt(std::complex<double>(tau2-K1_s))+C2/sqrt(tau2-K2_s)+C3/sqrt(tau2-K3_s);
            U_H[0]=U_H[0]*(1-fNu)*tau2/2.0;
            
            U_H[1]=1.0-C1*sqrt(std::complex<double>(tau2-K1_s))-C2*sqrt(tau2-K2_s)-C3*sqrt(tau2-K3_s);
            U_H[1]=U_H[1]/2.0;
            
            
            /*std::complex<double> b1=Ellip_K(n);
            std::complex<double> b2=Ellip_PI(n2*m1, n);
            std::complex<double> b3=Ellip_PI(n2*m2, n);
            std::complex<double> b4=Ellip_PI(n2*m3, n);*/
            
            U_Z[0]=2.0*Ellip_K(n)-B1*Ellip_PI(n2*m1, n)-B2*Ellip_PI(n2*m2, n)-B3*Ellip_PI(n2*m3, n);
            U_Z[0]=U_Z[0]/(M_PI*pow(1.0-cscd2, 3.0/2.0));
            
            U_Z[2]=1.0-A1/sqrt(std::complex<double>(tau2-K1_s))-A2/sqrt(tau2-K2_s)-A3/sqrt(tau2-K3_s);
            U_Z[2]=U_Z[2]/2.0;
            
        }else if (tau>=1 && tau<gamma){
            
            std::complex<double> temp=sqrt(std::complex<double>(tau2-K1_s));
            
            U_H[0]=1.0+(1.0-fNu)*tau2*C1/temp;
            
            U_H[1]=1.0-C1*temp;
            
            U_Z[0]=2.0*Ellip_K(1.0/n)-B1*Ellip_PI(m1, 1.0/n)-B2*Ellip_PI(m2, 1.0/n)-B3*Ellip_PI(m3, 1.0/n);
            U_Z[0]=U_Z[0]/(n*M_PI*pow(1.0-cscd2, 3.0/2.0));
            
            U_Z[2]=1.0-A1/temp;
            
        }else{
            
            U_H[0]=1.0;
            
            U_H[1]=1.0;
            
            U_Z[0]=2.0*Ellip_K(1.0/n)-B1*Ellip_PI(m1, 1.0/n)-B2*Ellip_PI(m2, 1.0/n)-B3*Ellip_PI(m3, 1.0/n);
            U_Z[0]=U_Z[0]/( n*M_PI*pow(1.0-cscd2, 3.0/2.0) )+2.0*D/sqrt(std::complex<double>(tau2-K1_s));
            
            U_Z[2]=1.0;
            
        }
        
        U_Z[0]=-U_Z[0]*tau/(8.0*M_PI*fMu*fr);
        U_Z[2]=U_Z[2]*(1.0-fNu)/(2.0*M_PI*fMu*fr);
        
        U_X[0]=U_H[0]*fcTheta/(2.0*M_PI*fMu*fr);
        U_X[1]=U_H[1]*(1.0-fNu)*(-fsTheta)/(2.0*M_PI*fMu*fr);
        U_X[2]=U_Z[0]*(-fcTheta);
        
        U_Y[0]=U_H[0]*fsTheta/(2.0*M_PI*fMu*fr);
        U_Y[1]=U_H[1]*(1.0-fNu)*(fcTheta)/(2.0*M_PI*fMu*fr);
        U_Y[2]=U_Z[0]*(-fsTheta);
        
        double* dis_ptr=fDis_Heavi+9*i;
        
        dis_ptr[0*3+0]=real(U_X[0]*fcTheta-U_X[1]*fsTheta);
        dis_ptr[1*3+0]=real(U_X[0]*fsTheta+U_X[1]*fcTheta);
        dis_ptr[2*3+0]=real(U_X[2]);
        
        dis_ptr[0*3+1]=real(U_Y[0]*fcTheta-U_Y[1]*fsTheta);
        dis_ptr[1*3+1]=real(U_Y[0]*fsTheta+U_Y[1]*fcTheta);
        dis_ptr[2*3+1]=real(U_Y[2]);
        
        dis_ptr[0*3+2]=real(U_Z[0]*fcTheta);
        dis_ptr[1*3+2]=real(U_Z[0]*fsTheta);
        dis_ptr[2*3+2]=real(U_Z[2]);
        
        
        //MathOperationT::PrintVector(9, dis_ptr, "displacement");
        
    }//end of for loop
 
    
    //compute the response to BSpline loading using convolution
    if (fDT==0.0) {
        cout << "!!!!Width of B-Spline loading is 0!!!" <<endl;
    }
    
    for (int i=0; i<fNumT_BSpline*9; i++) {
        fDis_BSpline[i]=0;
    }
    
    for (int i=0; i<fNumT_BSpline; i++) {
        
        double* dis_BSpline=fDis_BSpline+i*9;
        
        double B1, DB1;
        double t_temp=ft_BSpline[i];
        
        int sum_num=0;
        
        for (int j=0; j<fNumT_Heavi; j++) {
            
            if (ft_Heavi[j]>=t_temp-fDT && ft_Heavi[j]<=t_temp) {
                sum_num += 1;
                CubicB1(B1, DB1, t_temp-ft_Heavi[j], fDT);
                
                double* dis_Heavi=fDis_Heavi+j*9;
                
                for (int d=0; d<9; d++) {
                    dis_BSpline[d] += dis_Heavi[d]*DB1;
                }
                
            }
            
        }
        
        if (sum_num==0) {
            cout <<"!!!Number of points in the time convolution is 0.!!!" <<endl;
        }
        
        for (int d=0; d<9; d++) {
            
            if (t_temp<fDT) {
                dis_BSpline[d] = dis_BSpline[d] * t_temp /sum_num;
            }else{
                dis_BSpline[d] = dis_BSpline[d] * fDT/sum_num;
            }
            
        }
        
    }
}


std::complex<double> KauselT::Ellip_K( std::complex<double> n)
{
    
    int size=50;
    double d_phi=(M_PI/2.0)/size;
    
    std::complex<double> K=0.0;
    
    for (int i=0; i<size; i++) {
        double phi=(i+0.5)*d_phi;
        K=K+1.0/sqrt( 1.0-n*n*pow(sin(phi), 2) )*d_phi;
    }
 
    return K;
}


std::complex<double> KauselT::Ellip_PI( std::complex<double> m,std::complex<double> n)
{
    int size=50;
    double d_phi=(M_PI/2.0)/size;
    
    std::complex<double> P=0.0;
    for (int i=0; i<size; i++) {
        double phi=(i+0.5)*d_phi;
        P=P+1.0/( (1.0+m*pow(sin(phi),2) )*sqrt(1.0-n*n*pow(sin(phi),2)))*d_phi;
    }
    
    return P;
}

double KauselT::Heaviside(double t)
{
    if (t>=0)
        return 1.0;
    else
        return 0.0;
}


void KauselT::CubicB1(double& B1, double& DB1, double t, double T)
{
    t=t/T;
    double tp2=pow(t,2);
    double tp3=pow(t,3);
    
    double p1, p2, p3, p4;
    
    p1=32.0/3.0;
    B1=(Heaviside(t)-Heaviside(t-0.25))*p1*tp3;
    
    p1=-32;
    p2=32;
    p3=-8;
    p4=2.0/3.0;
    B1+=(Heaviside(t-0.25)-Heaviside(t-0.5))*(p1*tp3+p2*tp2+p3*t+p4);
    
    p1=32;
    p2=-64;
    p3=40;
    p4=-22.0/3.0;
    B1+=(Heaviside(t-0.5)-Heaviside(t-0.75))*(p1*tp3+p2*tp2+p3*t+p4);
    
    p1=-32.0/3.0;
    p2=32;
    p3=-32;
    p4=32.0/3.0;
    B1+=(Heaviside(t-0.75)-Heaviside(t-1))*(p1*tp3+p2*tp2+p3*t+p4);
    
    DB1 =(Heaviside(t)-Heaviside(t-0.25))*32*tp2;
    DB1+=(Heaviside(t-0.25)-Heaviside(t-0.5))*(-96*tp2+64*t-8);
    DB1+=(Heaviside(t-0.5)-Heaviside(t-0.75))*(96*tp2-128*t+40);
    DB1+=(Heaviside(t-0.75)-Heaviside(t-1))*(-32*tp2+64*t-32);
    DB1=DB1/T;
    
}


void KauselT::GetDisplacement(int& num, double** displacement)
{
    num=fNumT_BSpline;
    *displacement=fDis_BSpline;
}

void KauselT::GetDisplacement_Heavi(int& num, double** time, double** displacement)
{
    num=fNumT_Heavi;
    *time=ft_Heavi;
    *displacement=fDis_Heavi;
}


void KauselT::GetTraction(int& num, double** trac)
{
    num=fNumT_BSpline;
    *trac=fTrac;
}










