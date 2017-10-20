//
//  Bessel.h
//  GreenFunction
//
//  Created by Xiaoyong Bai on 5/14/16.
//
//


#ifndef BesselT_hpp
#define BesselT_hpp

#include <complex>


namespace GreenFunction {

class BesselT
{
public:
    //J0
    static double BesselJ0(double x);
    
    static std::complex<double> BesselJ0(std::complex<double> x);
    
    //J1
    static double BesselJ1(double x);
    
    static std::complex<double> BesselJ1(std::complex<double> x);
    
    //J2 and higher order
    static double BesselJ(int n, double x);
    
    static std::complex<double> BesselJ(int n, std::complex<double> x);

    
}; //end of class


} //end of namespace

#endif