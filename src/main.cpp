#include <iostream> 
#include "fstream"
#include "iomanip"
#include "BesselT.h"
#include "MindlinT.h"
#include "PakT.h"
#include "PakBaiT.h"
#include "GuzinaBiMatT.h"
#include "GuzinaBaiT.h"
#include "KauselT.h"
#include "DehoopBaiT.h"
#include <complex>

using namespace std;
using namespace GreenFunction;

void TestBessel();

void TestMindlin();

void TestPak();

void TestPakBai();

void TestGuzinaBiMat();

void TestGuzinaBai();

void TestKausel();

void TestDehoopBai();

int main()
{
    //TestBessel();
    
    //TestMindlin();
    
    TestPak();
    
    //TestPakBai();
    
    //TestGuzinaBiMat();
    
    //TestGuzinaBai();
    
    //TestKausel();
    
    //TestDehoopBai();
    
    return 0;
}


void TestBessel()
{
    int size=1001;
    double start=0.0, end=100.0;
    
    const complex<double> I(0.0,1.0);
    
    complex<double>* x=new complex<double>[size];
    complex<double>* J=new complex<double>[size];
    
    for (int a=0; a<size; a++) {
        x[a]=start+a*end/(size-1.0)*(1.0+0.1*I);

        J[a]=BesselT::BesselJ(2,x[a]);
    }
    
    ofstream J_file("J0.txt");
    
    J_file<<"%" <<setw(9)<<"x" << setw(20) << "real J(0,x)" << setw(20) << "iamg J(0,x)" << endl;
    
    cout.precision(5);
    for (int i=0; i<size; i++) {
        J_file<< setprecision(10) << setw(10) <<real(x[i])<<setw(20)<<real(J[i]) << setw(20) << imag(J[i])<<endl;
        
    }
    
    J_file.close();
    
    cout << "Bessel function is tested \n";
}


void TestMindlin()
{
    double src[]={0,0,0.5};
    double normal[]={1,0,0};
    double E=5e7;
    double nu=0.25;
    
    MindlinT Mind;
    Mind.SetMaterial(E, nu);
    
    int number=1000;
    double x[number];
    
    ofstream U_file("Mindlin_U_CPP.txt");
    U_file<<"%" <<setw(9)<<"x" << setw(15) << "U11" << setw(15) << "U12" << setw(15) << "U13"
                               << setw(15) << "U21" << setw(15) << "U22" << setw(15) << "U23"
                               << setw(15) << "U31" << setw(15) << "U32" << setw(15) << "U33" << endl;

    ofstream SX_file("Mindlin_SX_CPP.txt");
    SX_file<<"%" <<setw(9)<<"x" << setw(15) << "SX11" << setw(15) << "SX12" << setw(15) << "SX13"
                                << setw(15) << "SX21" << setw(15) << "SX22" << setw(15) << "SX23"
                                << setw(15) << "SX31" << setw(15) << "SX32" << setw(15) << "SX33" << endl;
    
    ofstream SY_file("Mindlin_SY_CPP.txt");
    SY_file<<"%" <<setw(9)<<"x" << setw(15) << "SY11" << setw(15) << "SY12" << setw(15) << "SY13"
                                << setw(15) << "SY21" << setw(15) << "SY22" << setw(15) << "SY23"
                                << setw(15) << "SY31" << setw(15) << "SY32" << setw(15) << "SY33" << endl;

    ofstream SZ_file("Mindlin_SZ_CPP.txt");
    SY_file<<"%" <<setw(9)<<"x" << setw(15) << "SZ11" << setw(15) << "SZ12" << setw(15) << "SZ13"
                                << setw(15) << "SZ21" << setw(15) << "SZ22" << setw(15) << "SZ23"
                                << setw(15) << "SZ31" << setw(15) << "SZ32" << setw(15) << "SZ33" << endl;
    
    for (int i=0; i<1000; i++) {
        x[i]=(i+0)*0.01;
        double dest[]={0.5*x[i], -4*x[i], 2.5*x[i]};
        
        Mind.SetGeometry(src, dest, normal);
        Mind.Compute();
        
        const double* U=Mind.GetDisplacement();
        double* SX;
        double* SY;
        double* SZ;
        
        Mind.GetStress(&SX, &SY, &SZ);
        
        U_file<<setprecision(5)<<setw(10)<<x[i]<<setw(15) << U[0] << setw(15) << U[1] << setw(15) << U[2]
                                                << setw(15) << U[3] << setw(15) << U[4] << setw(15) << U[5]
                                                << setw(15) << U[6] << setw(15) << U[7] << setw(15) << U[8] << endl;
        
        SX_file<<setprecision(5)<<setw(10)<<x[i]<<setw(15) << SX[0] << setw(15) << SX[1] << setw(15) << SX[2]
                                                <<setw(15) << SX[3] << setw(15) << SX[4] << setw(15) << SX[5]
                                                <<setw(15) << SX[6] << setw(15) << SX[7] << setw(15) << SX[8] << endl;
        
        SY_file<<setprecision(5)<<setw(10)<<x[i]<<setw(15) << SY[0] << setw(15) << SY[1] << setw(15) << SY[2]
                                                <<setw(15) << SY[3] << setw(15) << SY[4] << setw(15) << SY[5]
                                                <<setw(15) << SY[6] << setw(15) << SY[7] << setw(15) << SY[8] << endl;
        
        SZ_file<<setprecision(5)<<setw(10)<<x[i]<<setw(15) << SZ[0] << setw(15) << SZ[1] << setw(15) << SZ[2]
                                                <<setw(15) << SZ[3] << setw(15) << SZ[4] << setw(15) << SZ[5]
                                                <<setw(15) << SZ[6] << setw(15) << SZ[7] << setw(15) << SZ[8] << endl;
        
    }
    
    U_file.close();
    SX_file.close();
    SY_file.close();
    SZ_file.close();
    
    cout<<"test of MindlinT is over"<<endl;
    
}



void TestPak()
{
    double src[]={0,0,1};
    double dest[]={10, 0,  4};
    double normal[]={1,2,3};
    double E=5e7;
    double nu=0.25;
    double rho=1730;
    
    PakT Pak;
    Pak.SetMaterial(E, nu, rho);
    Pak.SetGeometry(src, dest, normal);

    int number=4000;
    double w[number];
    
    ofstream U_real_file("Pak_U_real_CPP.txt");
    U_real_file<<"%" <<setw(9)<<"w" << setw(15) << "U11" << setw(15) << "U12" << setw(15) << "U13"
                               << setw(15) << "U21" << setw(15) << "U22" << setw(15) << "U23"
                               << setw(15) << "U31" << setw(15) << "U32" << setw(15) << "U33" << endl;
    
    
    ofstream U_imag_file("Pak_U_imag_CPP.txt");
    U_imag_file<<"%" <<setw(9)<<"w" << setw(15) << "U11" << setw(15) << "U12" << setw(15) << "U13"
    << setw(15) << "U21" << setw(15) << "U22" << setw(15) << "U23"
    << setw(15) << "U31" << setw(15) << "U32" << setw(15) << "U33" << endl;
    
    
    ofstream SX_file("Pak_SX_CPP.txt");
    SX_file<<"%" <<setw(9)<<"w" << setw(15) << "SX11" << setw(15) << "SX12" << setw(15) << "SX13"
                                << setw(15) << "SX21" << setw(15) << "SX22" << setw(15) << "SX23"
                                << setw(15) << "SX31" << setw(15) << "SX32" << setw(15) << "SX33" << endl;
    
    ofstream SY_file("Pak_SY_CPP.txt");
    SY_file<<"%" <<setw(9)<<"w" << setw(15) << "SY11" << setw(15) << "SY12" << setw(15) << "SY13"
                                << setw(15) << "SY21" << setw(15) << "SY22" << setw(15) << "SY23"
                                << setw(15) << "SY31" << setw(15) << "SY32" << setw(15) << "SY33" << endl;
    
    ofstream SZ_file("Pak_SZ_CPP.txt");
    SZ_file<<"%" <<setw(9)<<"w" << setw(15) << "SZ11" << setw(15) << "SZ12" << setw(15) << "SZ13"
                                << setw(15) << "SZ21" << setw(15) << "SZ22" << setw(15) << "SZ23"
                                << setw(15) << "SZ31" << setw(15) << "SZ32" << setw(15) << "SZ33" << endl;
    
    for (int i=0; i<number; i++) {
        w[i]=i*1;
        
        Pak.Compute(w[i]);
        
        const complex<double>* U=Pak.GetDisplacement();
        complex<double>* SX;
        complex<double>* SY;
        complex<double>* SZ;
        
        Pak.GetStress(&SX, &SY, &SZ);
        
        U_real_file<<setprecision(5)<<setw(10)<< w[i]<<setw(15) << real(U[0]) << setw(15) << real(U[1]) << setw(15) << real(U[2])
                               <<setw(15) << real(U[3]) << setw(15) << real(U[4]) << setw(15) << real(U[5])
                               <<setw(15) << real(U[6]) << setw(15) << real(U[7]) << setw(15) << real(U[8]) << endl;
        
        U_imag_file<<setprecision(5)<<setw(10)<< w[i]<<setw(15) << imag(U[0]) << setw(15) << imag(U[1]) << setw(15) << imag(U[2])
        <<setw(15) << imag(U[3]) << setw(15) << imag(U[4]) << setw(15) << imag(U[5])
        <<setw(15) << imag(U[6]) << setw(15) << imag(U[7]) << setw(15) << imag(U[8]) << endl;
        
        
        SX_file<<setprecision(5)<<setw(10)<< w[i]<<setw(15) << real(SX[0]) << setw(15) << real(SX[1]) << setw(15) << real(SX[2])
                                <<setw(15) << real(SX[3]) << setw(15) << real(SX[4]) << setw(15) << real(SX[5])
                                <<setw(15) << real(SX[6]) << setw(15) << real(SX[7]) << setw(15) << real(SX[8]) << endl;
        
        SY_file<<setprecision(5)<<setw(10)<<w[i]<<setw(15) << real(SY[0]) << setw(15) << real(SY[1]) << setw(15) << real(SY[2])
                                <<setw(15) << real(SY[3]) << setw(15) << real(SY[4]) << setw(15) << real(SY[5])
                                <<setw(15) << real(SY[6]) << setw(15) << real(SY[7]) << setw(15) << real(SY[8]) << endl;
        
        SZ_file<<setprecision(5)<<setw(10)<<w[i]<<setw(15) << real(SZ[0]) << setw(15) << real(SZ[1]) << setw(15) << real(SZ[2])
                                <<setw(15) << real(SZ[3]) << setw(15) << real(SZ[4]) << setw(15) << real(SZ[5])
                                <<setw(15) << real(SZ[6]) << setw(15) << real(SZ[7]) << setw(15) << real(SZ[8]) << endl;
        
        cout << "computing omega " << w[i] << endl;
        
    }
    
    U_real_file.close();
    U_imag_file.close();
    
    SX_file.close();
    SY_file.close();
    SZ_file.close();
    
    cout<<"test Pak is over"<<endl;
    
}



void TestPakBai()
{
    double src[]={0,0,0};
    double dest[]={20,0,0};
    double normal[]={1,2,3};
    double E=5e7;
    double nu=0.25;
    double rho=1730;
    double DT=0.02;
    
    PakBaiT PakBai;
    
    PakBai.SetBWidth(DT);
    PakBai.SetMaterial(E, nu, rho);
    PakBai.SetGeometry(src, dest, normal);

    int number=1000;
    double time[number];
    for (int i=0; i<number; i++) {
        time[i]=i*3e-4;
    }
    
    ofstream U_file("PakBai_U_CPP.txt");
    U_file<<"%" <<setw(9)<<"time" << setw(15) << "U11" << setw(15) << "U12" << setw(15) << "U13"
                                  << setw(15) << "U21" << setw(15) << "U22" << setw(15) << "U23"
                                  << setw(15) << "U31" << setw(15) << "U32" << setw(15) << "U33" << endl;
    
    ofstream T_file("PakBai_T_CPP.txt");
    T_file<<"%" <<setw(9)<<"time" << setw(15) << "T11" << setw(15) << "T12" << setw(15) << "T13"
                                  << setw(15) << "T21" << setw(15) << "T22" << setw(15) << "T23"
                                  << setw(15) << "T31" << setw(15) << "T32" << setw(15) << "T33" << endl;
    
    PakBai.SetTime(number, time);
    PakBai.Compute();
    
    double* dis, *trac;
    int num;
    
    PakBai.GetDisplacement(num, &dis);
    PakBai.GetTraction(num, &trac);
    
    for (int it=0; it<num; it++) {
        double* U_t=dis+9*it;
        double* T_t=trac+9*it;
        
        U_file<<setprecision(5)<<setw(10)<< time[it] <<setw(15) << U_t[0] << setw(15) << U_t[1] << setw(15) << U_t[2]
                                                     <<setw(15) << U_t[3] << setw(15) << U_t[4] << setw(15) << U_t[5]
                                                     <<setw(15) << U_t[6] << setw(15) << U_t[7] << setw(15) << U_t[8]<< endl;
        
        T_file<<setprecision(5)<<setw(10)<< time[it] <<setw(15) << T_t[0] << setw(15) << T_t[1] << setw(15) << T_t[2]
                                                     <<setw(15) << T_t[3] << setw(15) << T_t[4] << setw(15) << T_t[5]
                                                     <<setw(15) << T_t[6] << setw(15) << T_t[7] << setw(15) << T_t[8]<< endl;
    }
    
    cout<<"test of PakBai is over"<<endl;
    
}



void TestGuzinaBiMat()
{
    double src[]={0.0,0.0,0.5};
    double normal[]={0.0,0.0,1.0};
    double E1=0;
    double nu1=0.25;
    double E2=5e7;
    double nu2=0.25;
    
    GuzinaBiMatT BiMat;
    BiMat.SetMaterial(E1, nu1, E2, nu2);
    
    int number=1000;
    double x[number];
    
    ofstream U_file("BiMat_U_CPP.txt");
    U_file<<"%" <<setw(9)<<"x" << setw(15) << "U11" << setw(15) << "U12" << setw(15) << "U13"
    << setw(15) << "U21" << setw(15) << "U22" << setw(15) << "U23"
    << setw(15) << "U31" << setw(15) << "U32" << setw(15) << "U33" << endl;
    
    ofstream SX_file("BiMat_SX_CPP.txt");
    SX_file<<"%" <<setw(9)<<"x" << setw(15) << "SX11" << setw(15) << "SX12" << setw(15) << "SX13"
    << setw(15) << "SX21" << setw(15) << "SX22" << setw(15) << "SX23"
    << setw(15) << "SX31" << setw(15) << "SX32" << setw(15) << "SX33" << endl;
    
    ofstream SY_file("BiMat_SY_CPP.txt");
    SY_file<<"%" <<setw(9)<<"x" << setw(15) << "SY11" << setw(15) << "SY12" << setw(15) << "SY13"
    << setw(15) << "SY21" << setw(15) << "SY22" << setw(15) << "SY23"
    << setw(15) << "SY31" << setw(15) << "SY32" << setw(15) << "SY33" << endl;
    
    ofstream SZ_file("BiMat_SZ_CPP.txt");
    SY_file<<"%" <<setw(9)<<"x" << setw(15) << "SZ11" << setw(15) << "SZ12" << setw(15) << "SZ13"
    << setw(15) << "SZ21" << setw(15) << "SZ22" << setw(15) << "SZ23"
    << setw(15) << "SZ31" << setw(15) << "SZ32" << setw(15) << "SZ33" << endl;
    
    for (int i=0; i<1000; i++) {
        x[i]=i*0.01;
        double dest[]={0.5*x[i], -4*x[i],2.5*x[i]};
        
        BiMat.SetGeometry(src, dest, normal);
        BiMat.Compute();
        
        const double* U=BiMat.GetDisplacement();
        double* SX;
        double* SY;
        double* SZ;
        
        BiMat.GetStress(&SX, &SY, &SZ);
        
        U_file<<setprecision(5)<<setw(10)<<x[i]<<setw(15) << U[0] << setw(15) << U[1] << setw(15) << U[2]
        << setw(15) << U[3] << setw(15) << U[4] << setw(15) << U[5]
        << setw(15) << U[6] << setw(15) << U[7] << setw(15) << U[8] << endl;
        
        SX_file<<setprecision(5)<<setw(10)<<x[i]<<setw(15) << SX[0] << setw(15) << SX[1] << setw(15) << SX[2]
        <<setw(15) << SX[3] << setw(15) << SX[4] << setw(15) << SX[5]
        <<setw(15) << SX[6] << setw(15) << SX[7] << setw(15) << SX[8] << endl;
        
        SY_file<<setprecision(5)<<setw(10)<<x[i]<<setw(15) << SY[0] << setw(15) << SY[1] << setw(15) << SY[2]
        <<setw(15) << SY[3] << setw(15) << SY[4] << setw(15) << SY[5]
        <<setw(15) << SY[6] << setw(15) << SY[7] << setw(15) << SY[8] << endl;
        
        SZ_file<<setprecision(5)<<setw(10)<<x[i]<<setw(15) << SZ[0] << setw(15) << SZ[1] << setw(15) << SZ[2]
        <<setw(15) << SZ[3] << setw(15) << SZ[4] << setw(15) << SZ[5]
        <<setw(15) << SZ[6] << setw(15) << SZ[7] << setw(15) << SZ[8] << endl;
        
    }
    
    U_file.close();
    SX_file.close();
    SY_file.close();
    SZ_file.close();
    
    cout<<"test of BiMatT is over"<<endl;
    
}



void TestGuzinaBai()
{
    double src[]={0,0,2};
    double dest[]={0,0,4};
    double normal[]={1,2,3};
    double DT=0.02;
    
    GuzinaBaiT GuzinaBai;
    
    GuzinaBai.SetBWidth(DT);
    GuzinaBai.SetGeometry(src, dest, normal);
    
    int number=4000;
    double time[number];
    for (int i=0; i<number; i++) {
        time[i]=i*2e-4;
    }
    
    ofstream U_file("GuzinaBai_U_CPP.txt");
    U_file<<"%" <<setw(9)<<"time" << setw(15) << "U11" << setw(15) << "U12" << setw(15) << "U13"
    << setw(15) << "U21" << setw(15) << "U22" << setw(15) << "U23"
    << setw(15) << "U31" << setw(15) << "U32" << setw(15) << "U33" << endl;
    
    ofstream T_file("GuzinaBai_T_CPP.txt");
    T_file<<"%" <<setw(9)<<"time" << setw(15) << "T11" << setw(15) << "T12" << setw(15) << "T13"
    << setw(15) << "T21" << setw(15) << "T22" << setw(15) << "T23"
    << setw(15) << "T31" << setw(15) << "T32" << setw(15) << "T33" << endl;
    
    ofstream SX_file("GuzinaBai_SX_CPP.txt");
    SX_file<<"%" <<setw(9)<<"time" << setw(15) << "SX11" << setw(15) << "SX12" << setw(15) << "SX13"
    << setw(15) << "SX22" << setw(15) << "SX23" << setw(15) << "SX33" << endl;
 
    ofstream SY_file("GuzinaBai_SY_CPP.txt");
    SY_file<<"%" <<setw(9)<<"time" << setw(15) << "SY11" << setw(15) << "SY12" << setw(15) << "SY13"
    << setw(15) << "SY22" << setw(15) << "SY23" << setw(15) << "SY33" << endl;
    
    ofstream SZ_file("GuzinaBai_SZ_CPP.txt");
    SZ_file<<"%" <<setw(9)<<"time" << setw(15) << "SZ11" << setw(15) << "SZ12" << setw(15) << "SZ13"
    << setw(15) << "SZ22" << setw(15) << "SZ23" << setw(15) << "SZ33" << endl;
    
    
    GuzinaBai.SetTime(number, time);
    GuzinaBai.Compute();
    
    double* dis, *trac;
    double *sx, *sy, *sz;
    int num;
    
    GuzinaBai.GetDisplacement(num, &dis);
    GuzinaBai.GetTraction(num, &trac);
    GuzinaBai.GetStress(num, &sx, &sy, &sz);
    
    for (int it=0; it<num; it++) {
        double* U_t=dis+9*it;
        double* T_t=trac+9*it;
        double* sx_t=sx+9*it;
        double* sy_t=sy+9*it;
        double* sz_t=sz+9*it;
        
        U_file<<setprecision(5)<<setw(10)<< time[it] <<setw(15) << U_t[0] << setw(15) << U_t[1] << setw(15) << U_t[2]
        <<setw(15) << U_t[3] << setw(15) << U_t[4] << setw(15) << U_t[5]
        <<setw(15) << U_t[6] << setw(15) << U_t[7] << setw(15) << U_t[8]<< endl;
        
        T_file<<setprecision(5)<<setw(10)<< time[it] <<setw(15) << T_t[0] << setw(15) << T_t[1] << setw(15) << T_t[2]
        <<setw(15) << T_t[3] << setw(15) << T_t[4] << setw(15) << T_t[5]
        <<setw(15) << T_t[6] << setw(15) << T_t[7] << setw(15) << T_t[8]<< endl;
        
        SX_file<<setprecision(5)<<setw(10)<<time[it]<<setw(15)<<sx_t[0]<<setw(15)<<sx_t[1]<<setw(15)<<sx_t[2]
        <<setw(15)<<sx_t[4]<<setw(15)<<sx_t[5]<<setw(15)<<sx_t[8]<<endl;
        
        SY_file<<setprecision(5)<<setw(10)<<time[it]<<setw(15)<<sy_t[0]<<setw(15)<<sy_t[1]<<setw(15)<<sy_t[2]
        <<setw(15)<<sy_t[4]<<setw(15)<<sy_t[5]<<setw(15)<<sy_t[8]<<endl;
        
        SZ_file<<setprecision(5)<<setw(10)<<time[it]<<setw(15)<<sz_t[0]<<setw(15)<<sz_t[1]<<setw(15)<<sz_t[2]
        <<setw(15)<<sz_t[4]<<setw(15)<<sz_t[5]<<setw(15)<<sz_t[8]<<endl;
    }
    
    cout<<"test of GuzinaBai is over"<<endl;
    
    U_file.close();
    T_file.close();
    SX_file.close();
    SY_file.close();
    SZ_file.close();
    
}



void TestKausel()
{
    double src[]={0,0,0};
    double dest[]={20,0,0};
    
    KauselT Kausel;
    
    Kausel.SetMaterial(5e7, 0.25, 1730);
    Kausel.SetGeometry(src, dest);
    
    int number=1000;
    double time[number];
    for (int i=0; i<number; i++) {
        time[i]=i*3e-4;
    }
    
    Kausel.SetTime(number, time);
    Kausel.SetBWidth(0.02);
    
    Kausel.Compute();
    
    double* dis_BSpline;
    int num_BSpline;
    
    Kausel.GetDisplacement(num_BSpline, &dis_BSpline);
    
    ofstream U_file("Kausel_BSpline_U_CPP.txt");
    U_file<<"%" <<setw(9)<<"time" << setw(15) << "U11" << setw(15) << "U12" << setw(15) << "U13"
    << setw(15) << "U21" << setw(15) << "U22" << setw(15) << "U23"
    << setw(15) << "U31" << setw(15) << "U32" << setw(15) << "U33" << endl;
    
    for (int it=0; it<num_BSpline; it++) {
        double* U_t=dis_BSpline+9*it;
        
        U_file<<std::scientific;
        
        U_file<<setprecision(5)<<setw(10)<< time[it] <<setw(15) << U_t[0] << setw(15) << U_t[1] << setw(15) << U_t[2]
        <<setw(15) << U_t[3] << setw(15) << U_t[4] << setw(15) << U_t[5]
        <<setw(15) << U_t[6] << setw(15) << U_t[7] << setw(15) << U_t[8]<< endl;
        
    }
    U_file.close();
    
    double* dis_Heavi;
    double* time_Heavi;
    int num_Heavi;
    
    Kausel.GetDisplacement_Heavi(num_Heavi, &time_Heavi, &dis_Heavi);
    
    ofstream Heavi_file("Kausel_Heavi_U_CPP.txt");
    Heavi_file<<"%" <<setw(9)<<"time" << setw(15) << "U11" << setw(15) << "U12" << setw(15) << "U13"
    << setw(15) << "U21" << setw(15) << "U22" << setw(15) << "U23"
    << setw(15) << "U31" << setw(15) << "U32" << setw(15) << "U33" << endl;
    
    for (int it=0; it<num_Heavi; it++) {
        double* U_t=dis_Heavi+9*it;
        
        Heavi_file<<std::scientific;
        
        Heavi_file<<setprecision(5)<<setw(10)<< time_Heavi[it] <<setw(15) << U_t[0] << setw(15) << U_t[1] << setw(15)
        << U_t[2]   <<setw(15) << U_t[3] << setw(15) << U_t[4] << setw(15) << U_t[5]
        <<setw(15) << U_t[6] << setw(15) << U_t[7] << setw(15) << U_t[8]<< endl;
        
    }
    Heavi_file.close();
    
    cout<<"test of Kausel is over"<<endl;
    
    
}


void TestDehoopBai()
{
    double src[]={0,0,2};
    double dest[]={0,0,4.0};
    double normal[3]={1,2,3};
    
    DehoopBaiT DehoopBai;
    
    /*int nl=1;
    double E[1]={3e8};
    double nu[1]={0.25};
    double rho[1]={2000};
    double z_interface[1]={0};*/
    
    int nl=2;
    double E[2]={5e7, 10e7};
    double nu[2]={0.25, 0.25};
    double rho[2]={1730, 1730};
    double z_interface[2]={0, 10};
    
    /*int nl=5;
    double E[5]={176.0913, 352.1826, 3*176.0913, 4*176.0913, 5*176.0913};
    double nu[5]={0.25, 0.25, 0.25, 0.25, 0.25};
    double rho[5]={3.3, 3.3, 3.3, 3.3, 3.3};
    double z_interface[5]={0, 10, 20, 30 ,40};*/
    
    DehoopBai.SetLayers(nl, z_interface, E, nu, rho);
    DehoopBai.SetGeometry(src, dest, normal);
    
    
    int number=2000;
    double time[number];
    for (int i=0; i<number; i++) {
        time[i]=i*2e-4;
    }
    
    DehoopBai.SetTime(number, time);
    DehoopBai.SetBWidth(0.02);
    
    DehoopBai.Compute();
    
    double* dis_BSpline;
    int num_BSpline;
    
    DehoopBai.GetDisplacement(num_BSpline, &dis_BSpline);
    
    ofstream U_file("DehoopBai_BSpline_U_CPP.txt");
    U_file<<"%" <<setw(9)<<"time" << setw(15) << "U11" << setw(15) << "U12" << setw(15) << "U13"
    << setw(15) << "U21" << setw(15) << "U22" << setw(15) << "U23"
    << setw(15) << "U31" << setw(15) << "U32" << setw(15) << "U33" << endl;
    
    for (int it=0; it<num_BSpline; it++) {
        double* U_t=dis_BSpline+9*it;
        
        U_file<<std::scientific;
        
        U_file<<setprecision(5)<<setw(10)<< time[it] <<setw(15) << U_t[0] << setw(15) << U_t[1] << setw(15) << U_t[2]
        <<setw(15) << U_t[3] << setw(15) << U_t[4] << setw(15) << U_t[5]
        <<setw(15) << U_t[6] << setw(15) << U_t[7] << setw(15) << U_t[8]<< endl;
        
    }
    U_file.close();
    
    
    double* Traction_BSpline;
    
    DehoopBai.GetTraction(num_BSpline, &Traction_BSpline);
    
    ofstream T_file("DehoopBai_BSpline_T_CPP.txt");
    T_file<<"%" <<setw(9)<<"time" << setw(15) << "T11" << setw(15) << "T12" << setw(15) << "T13"
    << setw(15) << "T21" << setw(15) << "T22" << setw(15) << "T23"
    << setw(15) << "T31" << setw(15) << "T32" << setw(15) << "T33" << endl;
    
    for (int it=0; it<num_BSpline; it++) {
        double* T_t=Traction_BSpline+9*it;
        
        T_file<<std::scientific;
        
        T_file<<setprecision(5)<<setw(10)<< time[it] <<setw(15) << T_t[0] << setw(15) << T_t[1] << setw(15) << T_t[2]
        <<setw(15) << T_t[3] << setw(15) << T_t[4] << setw(15) << T_t[5]
        <<setw(15) << T_t[6] << setw(15) << T_t[7] << setw(15) << T_t[8]<< endl;
        
    }
    T_file.close();
    
    
    double *SX, *SY, *SZ;
    double *time_BSpline;
    
    DehoopBai.GetStress_BSpline(num_BSpline, &time_BSpline, &SX, &SY, &SZ);
    
    ofstream SX_file("DehoopBai_BSpline_SX_CPP.txt");
    SX_file<<"%" <<setw(9)<<"time" << setw(15) << "U11" << setw(15) << "U12" << setw(15) << "U13"
    << setw(15) << "U21" << setw(15) << "U22" << setw(15) << "U23"
    << setw(15) << "U31" << setw(15) << "U32" << setw(15) << "U33" << endl;
    
    for (int it=0; it<num_BSpline; it++) {
        double* SX_t=SX+9*it;
        
        SX_file<<std::scientific;
        
        SX_file<<setprecision(5)<<setw(10)<< time_BSpline[it] <<setw(15) << SX_t[0] << setw(15) << SX_t[1] << setw(15)
        << SX_t[2]   <<setw(15) << SX_t[4] << setw(15) << SX_t[5] <<setw(15)  << SX_t[8]<< endl;
        
    }
    SX_file.close();
    
    
    ofstream SY_file("DehoopBai_BSpline_SY_CPP.txt");
    SY_file<<"%" <<setw(9)<<"time" << setw(15) << "U11" << setw(15) << "U12" << setw(15) << "U13"
    << setw(15) << "U21" << setw(15) << "U22" << setw(15) << "U23"
    << setw(15) << "U31" << setw(15) << "U32" << setw(15) << "U33" << endl;
    
    for (int it=0; it<num_BSpline; it++) {
        double* SX_t=SY+9*it;
        
        SY_file<<std::scientific;
        
        SY_file<<setprecision(5)<<setw(10)<< time_BSpline[it] <<setw(15) << SX_t[0] << setw(15) << SX_t[1] << setw(15)
        << SX_t[2]   <<setw(15) << SX_t[4] << setw(15) << SX_t[5] <<setw(15)  << SX_t[8]<< endl;
        
    }
    SY_file.close();
    
    
    
    ofstream SZ_file("DehoopBai_BSpline_SZ_CPP.txt");
    SZ_file<<"%" <<setw(9)<<"time" << setw(15) << "U11" << setw(15) << "U12" << setw(15) << "U13"
    << setw(15) << "U21" << setw(15) << "U22" << setw(15) << "U23"
    << setw(15) << "U31" << setw(15) << "U32" << setw(15) << "U33" << endl;
    
    for (int it=0; it<num_BSpline; it++) {
        double* SX_t=SZ+9*it;
        
        SZ_file<<std::scientific;
        
        SZ_file<<setprecision(5)<<setw(10)<< time_BSpline[it] <<setw(15) << SX_t[0] << setw(15) << SX_t[1] << setw(15)
        << SX_t[2]   <<setw(15) << SX_t[4] << setw(15) << SX_t[5] <<setw(15)  << SX_t[8]<< endl;
        
    }
    SZ_file.close();
    
    
    
    int num_Heavi;
    double *time_Heavi, *U_Heavi;
    
    DehoopBai.GetDisplacement_Heavi(num_Heavi, &time_Heavi, &U_Heavi);
    
    ofstream U_Heavi_file("DehoopBai_Heavi_U_CPP.txt");
    SX_file<<"%" <<setw(9)<<"time" << setw(15) << "U11" << setw(15) << "U12" << setw(15) << "U13"
    << setw(15) << "U21" << setw(15) << "U22" << setw(15) << "U23"
    << setw(15) << "U31" << setw(15) << "U32" << setw(15) << "U33" << endl;
    
    for (int it=0; it<num_BSpline; it++) {
        double* u_temp=U_Heavi+9*it;
        
        U_Heavi_file<<std::scientific;
        
        U_Heavi_file<<setprecision(5)<<setw(10)<< time_Heavi[it] <<setw(15) << u_temp[0] << setw(15) << u_temp[1] << setw(15)
        << u_temp[2]   <<setw(15) << u_temp[3] << setw(15) << u_temp[4] << setw(15) << u_temp[5]
        <<setw(15) << u_temp[6] << setw(15) << u_temp[7] << setw(15) << u_temp[8]<< endl;
        
    }
    U_Heavi_file.close();

    
    cout<<"test of DehoopBaiT is over"<<endl;
    
    
}





