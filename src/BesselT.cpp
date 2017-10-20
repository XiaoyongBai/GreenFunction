//
//  BesselT.cpp
//  GreenFunction
//
//  Created by Xiaoyong Bai on 5/13/16.
//
//

#include <iostream>

#include "BesselT.h"
#include <cmath>


using namespace GreenFunction;
using namespace std;

//algorithms for numerical recipe in c

double BesselT::BesselJ0(double x)
{
    double J0;
    
    if (x < 8.0)
    {
        double y=x*x;
        double ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
        double ans2=57568490411.0+y*(1029532985.0+y*(9494680.718 +y*(59272.64853+y*(267.8532712+y*1.0))));
        J0=ans1/ans2;
    }
    else
    {
        double z=8.0/x;
        double y=z*z;
        double xx=x-0.785398164;
        double ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4+y*(-0.2073370639e-5+y*0.2093887211e-6)));
        double ans2 = -0.1562499995e-1+y*(0.1430488765e-3+y*(-0.6911147651e-5+y*(0.7621095161e-6 -y*0.934945152e-7)));
        J0=sqrt(0.636619772/x)*(cos(xx)*ans1-z*sin(xx)*ans2);
    }
    
    return J0;
}


complex<double> BesselT::BesselJ0(complex<double> x)
{
    complex<double> xx, y, ans1, ans2, z, J0;
    
    double ax=abs(x);
    
    if (ax < 8.0)
    {
        y=x*x;
        ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
        ans2=57568490411.0+y*(1029532985.0+y*(9494680.718 +y*(59272.64853+y*(267.8532712+y*1.0))));
        J0=ans1/ans2;
    }
    else
    {
        z=8.0/x;
        y=z*z;
        xx=x-0.785398164;
        ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4+y*(-0.2073370639e-5+y*0.2093887211e-6)));
        ans2 = -0.1562499995e-1+y*(0.1430488765e-3+y*(-0.6911147651e-5+y*(0.7621095161e-6 -y*0.934945152e-7)));
        J0=sqrt(0.636619772/x)*(cos(xx)*ans1-z*sin(xx)*ans2);
    }
    
    return J0;
}



double BesselT::BesselJ1(double x)
{
    double xx, y, ans1, ans2, z, J1;
    
    if (x < 8.0)
    {
        y=x*x;
        ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1 +y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
        ans2=144725228442.0+y*(2300535178.0+y*(18583304.74 +y*(99447.43394+y*(376.9991397+y*1.0))));
        J1=ans1/ans2;
    }
    else
    {
        z=8.0/x;
        y=z*z;
        xx=x-2.356194491;
        ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4+y*(0.2457520174e-5+y*(-0.240337019e-6))));
        ans2=0.04687499995+y*(-0.2002690873e-3+y*(0.8449199096e-5+y*(-0.88228987e-6+y*0.105787412e-6)));
        J1=sqrt(0.636619772/x)*(cos(xx)*ans1-z*sin(xx)*ans2);
    }
    
    return J1;
}

complex<double> BesselT::BesselJ1(complex<double> x)
{
    complex<double> xx, y, ans1, ans2, z, J1;
    
    if (abs(x) < 8.0)
    {
        y=x*x;
        ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1 +y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
        ans2=144725228442.0+y*(2300535178.0+y*(18583304.74 +y*(99447.43394+y*(376.9991397+y*1.0))));
        J1=ans1/ans2;
    }
    else
    {
        z=8.0/x;
        y=z*z;
        xx=x-2.356194491;
        ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4+y*(0.2457520174e-5+y*(-0.240337019e-6))));
        ans2=0.04687499995+y*(-0.2002690873e-3+y*(0.8449199096e-5+y*(-0.88228987e-6+y*0.105787412e-6)));
        J1=sqrt(0.636619772/x)*(cos(xx)*ans1-z*sin(xx)*ans2);
    }
    
    return J1;
}



double BesselT::BesselJ(int n, double x)
{
    const double IACC = 40;
    const double BIGNO = 1e10,  BIGNI = 1e-10;
    
    double TOX,BJM,BJ,BJP,SUM,TMP;
    int J, JSUM, M;
    
    if (n == 0) return BesselJ0(x);
    if (n == 1) return BesselJ1(x);
    if (x == 0.0) return 0.0;
    
    TOX = 2.0/x;
    if (x > 1.0*n) {
        BJM = BesselJ0(x);
        BJ  = BesselJ1(x);
        for (J=1; J<n; J++) {
            BJP = J*TOX*BJ-BJM;
            BJM = BJ;
            BJ  = BJP;
        }
        return BJ;
    }
    else {
        M = (int) (2*((n+floor(sqrt(1.0*(IACC*n))))/2));
        TMP = 0.0;
        JSUM = 0;
        SUM = 0.0;
        BJP = 0.0;
        BJ  = 1.0;
        for (J=M; J>0; J--) {
            BJM = J*TOX*BJ-BJP;
            BJP = BJ;
            BJ  = BJM;
            if (fabs(BJ) > BIGNO) {
                BJ  = BJ*BIGNI;
                BJP = BJP*BIGNI;
                TMP = TMP*BIGNI;
                SUM = SUM*BIGNI;
            }
            if (JSUM != 0)  SUM += BJ;
            JSUM = 1-JSUM;
            if (J == n)  TMP = BJP;
        }
        SUM = 2.0*SUM-BJ;
        return (TMP/SUM);
    }
}



complex<double> BesselT::BesselJ(int n, complex<double> x)
{
    const double IACC = 40;
    const double BIGNO = 1e10,  BIGNI = 1e-10;
    
    complex<double> TOX,BJM,BJ,BJP,SUM,TMP;
    int J, JSUM, M;
    
    if (n == 0) return BesselJ0(x);
    if (n == 1) return BesselJ1(x);
    if (abs(x) == 0.0) return 0.0;
    
    TOX = 2.0/x;
    if (abs(x) > 1.0*n) {
        BJM = BesselJ0(x);
        BJ  = BesselJ1(x);
        for (J=1; J<n; J++) {
            BJP = (J+0.0)*TOX*BJ-BJM;
            BJM = BJ;
            BJ  = BJP;
        }
        return BJ;
    }
    else {
        M = (int) (2*((n+floor(sqrt(1.0*(IACC*n))))/2));
        TMP = 0.0;
        JSUM = 0;
        SUM = 0.0;
        BJP = 0.0;
        BJ  = 1.0;
        for (J=M; J>0; J--) {
            BJM = (J+0.0)*TOX*BJ-BJP;
            BJP = BJ;
            BJ  = BJM;
            if (abs(BJ) > BIGNO) {
                BJ  = BJ*BIGNI;
                BJP = BJP*BIGNI;
                TMP = TMP*BIGNI;
                SUM = SUM*BIGNI;
            }
            if (JSUM != 0)  SUM += BJ;
            JSUM = 1.0-JSUM;
            if (J == n)  TMP = BJP;
        }
        SUM = 2.0*SUM-BJ;
        return (TMP/SUM);
    }
}