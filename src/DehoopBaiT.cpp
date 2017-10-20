//
//  DehoopBaiT.cpp
//  GreenFunction
//
//  Created by Xiaoyong Bai on 06/15/17.
//

#include "DehoopBaiT.h"
#include "MathOperationT.h"
#include "cmath"


using namespace GreenFunction;

DehoopBaiT::DehoopBaiT()
{
    fNumLayer=0;
    fZ_Interface=NULL;
    
    fE=NULL;
    fNu=NULL;
    fMu=NULL;
    fRho=NULL;
    fCp=NULL;
    fCs=NULL;
    
    fNumT_BSpline=0;
    fNumT_Heavi=0;
    
    ft_Heavi=NULL;
    ft_BSpline=NULL;
    
    fDis_Heavi=NULL;
    fSX_Heavi_p1=NULL; fSY_Heavi_p1=NULL; fSZ_Heavi_p1=NULL; //Stresses for Heaviside loading
    fSX_Heavi_p2=NULL; fSY_Heavi_p2=NULL; fSZ_Heavi_p2=NULL;
    
    fSX_BSpline=NULL; fSY_BSpline=NULL; fSZ_BSpline=NULL;
    
    fDis_BSpline=NULL;
    fTrac_BSpline=NULL;
    
    fDT=0.0;
    
    fZ_layer=-1;
    fS_layer=-1;
    
    //initialize wave groups
    fMinLevel=0;
    fMaxLevel=0;
    
    fLevel=NULL;
    fN=NULL;
    fID=NULL;
    fType=NULL;
    fDirection=NULL;
    fLayer=NULL;
    fC=NULL;
    fd=NULL;
    
    fArrival=NULL;
    
    //initialized normalized variables
    fa=1.0;
    
    fMuBar=NULL;
    fLamdaBar=NULL;
    fRhoBar=NULL;
    fCpBar=NULL;
    fCsBar=NULL;
    
    fdBar=NULL;
    fCBar=NULL;
    fArrivalBar=NULL;
}

DehoopBaiT::~DehoopBaiT()
{
    if (fE) delete [] fE;
    if (fNu) delete [] fNu;
    if (fRho) delete [] fRho;
    if (fMu) delete [] fMu;
    if (fCp) delete [] fCp;
    if (fCs) delete [] fCs;
    
    if (ft_Heavi) delete [] ft_Heavi;
    if (ft_BSpline) delete [] ft_BSpline;
    
    if (fDis_Heavi) delete [] fDis_Heavi;
    if (fDis_BSpline) delete [] fDis_BSpline;
    
    if (fSX_Heavi_p1) delete [] fSX_Heavi_p1;
    if (fSY_Heavi_p1) delete [] fSY_Heavi_p1;
    if (fSZ_Heavi_p1) delete [] fSZ_Heavi_p1;
    if (fSX_Heavi_p2) delete [] fSX_Heavi_p2;
    if (fSY_Heavi_p2) delete [] fSY_Heavi_p2;
    if (fSZ_Heavi_p2) delete [] fSZ_Heavi_p2;
    
    if (fSX_BSpline) delete [] fSX_BSpline;
    if (fSY_BSpline) delete [] fSY_BSpline;
    if (fSZ_BSpline) delete [] fSZ_BSpline;
    
    if (fTrac_BSpline) delete [] fTrac_BSpline;
    
    if (fZ_Interface) delete [] fZ_Interface;
    
    if (fMuBar) delete [] fMuBar;
    if (fLamdaBar) delete [] fLamdaBar;
    if (fRhoBar) delete [] fRhoBar;
    if (fCpBar) delete [] fCpBar;
    if (fCsBar) delete [] fCsBar;
    
    DeleteWaves();
}


void DehoopBaiT::SetBWidth(double T)
{
    fDT=T;
}

void DehoopBaiT::SetTime(int num, double *t)
{
    if (fNumT_BSpline != num) {
        delete [] ft_Heavi;
        delete [] ft_BSpline;
        
        delete [] fDis_BSpline;
        delete [] fDis_Heavi;

        if (fSX_Heavi_p1) delete [] fSX_Heavi_p1;
        if (fSY_Heavi_p1) delete [] fSY_Heavi_p1;
        if (fSZ_Heavi_p1) delete [] fSZ_Heavi_p1;
        if (fSX_Heavi_p2) delete [] fSX_Heavi_p2;
        if (fSY_Heavi_p2) delete [] fSY_Heavi_p2;
        if (fSZ_Heavi_p2) delete [] fSZ_Heavi_p2;
        
        if (fSX_BSpline) delete [] fSX_BSpline;
        if (fSY_BSpline) delete [] fSY_BSpline;
        if (fSZ_BSpline) delete [] fSZ_BSpline;
        
        delete [] fTrac_BSpline;
        
        fNumT_BSpline=num;
        ft_BSpline=new double[fNumT_BSpline];
        
        fDis_BSpline=new double[fNumT_BSpline*9];
    
        fSX_BSpline=new double[fNumT_BSpline*9];
        fSY_BSpline=new double[fNumT_BSpline*9];
        fSZ_BSpline=new double[fNumT_BSpline*9];
        fTrac_BSpline=new double[fNumT_BSpline*9];
        
        fNumT_Heavi=2*fNumT_BSpline;
        ft_Heavi=new double[fNumT_Heavi];
        
        fDis_Heavi=new double[fNumT_Heavi*9];
        
        fSX_Heavi_p1=new double[fNumT_Heavi*9];
        fSY_Heavi_p1=new double[fNumT_Heavi*9];
        fSZ_Heavi_p1=new double[fNumT_Heavi*9];
        fSX_Heavi_p2=new double[fNumT_Heavi*9];
        fSY_Heavi_p2=new double[fNumT_Heavi*9];
        fSZ_Heavi_p2=new double[fNumT_Heavi*9];
        
    }
    
    MathOperationT::MemCopy(fNumT_BSpline, t, ft_BSpline);
    
    //create time slots for Heaviside loading
    ft_max=ft_BSpline[0];
    for (int i=0; i<fNumT_BSpline; i++) {
        if (ft_BSpline[i]>ft_max) {
            ft_max=ft_BSpline[i];
        }
    }
    
    double t_inc=ft_max/(fNumT_Heavi-1);
    
    for (int i=0; i<fNumT_Heavi; i++) {
        ft_Heavi[i]=i*t_inc;
    }
    
    
    for (int i=0; i<fNumT_BSpline*9; i++) {
        fTrac_BSpline[i]=0;
        fDis_BSpline[i]=0;
    }
    
    for (int i=0; i<fNumT_Heavi*9; i++) {
        fDis_Heavi[i]=0;
    }
    
    //MathOperationT::PrintVector(fNumT_BSpline, ft_BSpline, "time to be computed");
    //MathOperationT::PrintVector(fNumT_Heavi, ft_Heavi, "Time for Heaviside loading");
}

void DehoopBaiT::SetLayers(int nl, double *z_interface, double *E, double *nu, double *density)
{
    if (fNumLayer != nl) {
        fNumLayer=nl;
        
        
        if (fE) delete [] fE;
        if (fNu) delete [] fNu;
        if (fRho) delete [] fRho;
        if (fMu) delete [] fMu;
        if (fCp) delete [] fCp;
        if (fCs) delete [] fCs;
        if (fZ_Interface) delete [] fZ_Interface;
        
        fE=new double[fNumLayer];
        fNu=new double[fNumLayer];
        fRho=new double[fNumLayer];
        fMu=new double[fNumLayer];
        fCp=new double[fNumLayer];
        fCs=new double[fNumLayer];
        
        fZ_Interface=new double[fNumLayer];
        
        if (fMuBar) delete [] fMuBar;
        if (fLamdaBar) delete [] fLamdaBar;
        if (fRhoBar) delete [] fRhoBar;
        if (fCpBar) delete [] fCpBar;
        if (fCsBar) delete [] fCsBar;
        fMuBar=new double[fNumLayer];
        fLamdaBar=new double[fNumLayer];
        fRhoBar=new double[fNumLayer];
        fCpBar=new double[fNumLayer];
        fCsBar=new double[fNumLayer];
    }
    
    MathOperationT::MemCopy(fNumLayer, E, fE);
    MathOperationT::MemCopy(fNumLayer, nu, fNu);
    MathOperationT::MemCopy(fNumLayer, density, fRho);
    MathOperationT::MemCopy(fNumLayer, z_interface, fZ_Interface);
    
    for (int i=0; i<fNumLayer; i++) {
        fMu[i]=fE[i]/(2.0*(1.0+fNu[i]));
        
        double lamda= fMu[i]*2.0*fNu[i]/(1.0-2.0*fNu[i]);
        
        fCp[i]=sqrt((lamda+2.0*fMu[i])/fRho[i]);
        fCs[i]=sqrt(fMu[i]/fRho[i]);
    }
    
    //Normalization
    fMu_ref=fMu[0];
    fRho_ref=fRho[0];
    fCs_ref=fCs[0];
    ft_ref=fa/fCs_ref;
    
    for (int i=0; i<fNumLayer; i++) {
        fMuBar[i]=fMu[i]/fMu_ref;
        fLamdaBar[i]= (fMu[i]*2.0*fNu[i]/(1.0-2.0*fNu[i]))/fMu_ref;
        fRhoBar[i]=fRho[i]/fRho_ref;
        fCpBar[i]=fCp[i]/fCs_ref;
        fCsBar[i]=fCs[i]/fCs_ref;
    }
    
    //MathOperationT::PrintVector(fNumLayer, fCp, "cp");
    //MathOperationT::PrintVector(fNumLayer, fCs, "cs");
    //MathOperationT::PrintVector(fNumLayer, fMuBar, "MuBar");
    //MathOperationT::PrintVector(fNumLayer, fLamdaBar, "LambdaBar");
    //MathOperationT::PrintVector(fNumLayer, fRhoBar, "RhoBar");
    //MathOperationT::PrintVector(fNumLayer, fCpBar, "cpBar");
    //MathOperationT::PrintVector(fNumLayer, fCsBar, "csBar");
    
    //Initialize the layer of source and receiver point.
    fZ_layer=-1;
    fS_layer=-1;
    
}


void DehoopBaiT::SetGeometry(double *src, double *dest, double* normal)
{
    if (fNumLayer==0) {
        throw "SetGeometry must be called after SetLayers";
    }
    
    if (src[2] < 0 || dest[2] < 0)
    {
        cout<< "depth of points must be >=0 !!!"<<endl;
    }
    
    fX=dest[0]-src[0];
    fY=dest[1]-src[1];
    fZ=dest[2];
    fS=src[2];
    
    fr=sqrt(fX*fX+fY*fY);
    frBar=fr/fa;
    
    if (fr==0) {
        fsTheta=0;
        fcTheta=1.0;
    }else{
        fsTheta=fY/fr;
        fcTheta=fX/fr;
    }

    MathOperationT::MemCopy(3, normal, fNormal);
    
    double nn=MathOperationT::VecNormal(3, fNormal);
    MathOperationT::VecScale(3, fNormal, 1/nn);
    
    
    int z_layer=WhichLayer(fZ);
    int s_layer=WhichLayer(fS);
    
    if (z_layer != fZ_layer || s_layer != fS_layer) {
        fZ_layer=z_layer;
        fS_layer=s_layer;
        
        GenerateWaves();
    }
    
    //generate d and compute arrival times
    for (int i=0; i<fMaxLevel; i++) {
        for (int j=0; j<fN[i]; j++) {
            
            double* d_ptr=fd[i]+j*fLevel[i];
            double* c_ptr=fC[i]+j*fLevel[i];
            int* layer_ptr=fLayer[i]+j*fLevel[i];
            int* id_ptr=fID[i]+j*fLevel[i];
            
            double* dBar_ptr=fdBar[i]+j*fLevel[i];
            
            if (fLevel[i]==1) {
                d_ptr[0]=abs(fZ-fS);
                dBar_ptr[0]=d_ptr[0]/fa;
            }else{
                
                //d for the first segment
                int id=id_ptr[0];
                int layer=layer_ptr[0];
                
                if (id==1 || id==2) {
                    d_ptr[0]=fS-fZ_Interface[layer];
                }else if (id==3 || id==4){
                    d_ptr[0]=fZ_Interface[layer+1]-fS;
                }
                dBar_ptr[0]=d_ptr[0]/fa;
                
                //d for the last segment
                id=id_ptr[fLevel[i]-1];
                layer=layer_ptr[fLevel[i]-1];
                
                if (id==1 || id==2) {
                    d_ptr[fLevel[i]-1]=fZ_Interface[layer+1]-fZ;
                }else if (id==3 || id==4){
                    d_ptr[fLevel[i]-1]=fZ-fZ_Interface[layer];
                }
                dBar_ptr[fLevel[i]-1]=d_ptr[fLevel[i]-1]/fa;
                
                for (int k=1; k<fLevel[i]-1; k++) {
                    layer=layer_ptr[k];
                    d_ptr[k]=fZ_Interface[layer+1]-fZ_Interface[layer];
                    dBar_ptr[k]=d_ptr[k]/fa;

                }
                
            }
            
            
            //Compute arrival time
            fArrival[i][j]=WaveArrival(fLevel[i], fr, d_ptr, c_ptr, id_ptr, layer_ptr);
            fArrivalBar[i][j]=fArrival[i][j]/ft_ref;
        }
        
        //MathOperationT::PrintMatrix(fN[i], fLevel[i], fd[i], "d of wave 1");
        //MathOperationT::PrintMatrix(fN[i], fLevel[i], fdBar[i], "dBar of wave 1");
        //MathOperationT::PrintMatrix(fN[i], fLevel[i], fC[i], "WaveSpeed");
        //MathOperationT::PrintMatrix(fN[i], fLevel[i], fCBar[i], "normalized WaveSpeed");
        //MathOperationT::PrintVector(fN[i], fArrival[i], "Arrival time");
        //MathOperationT::PrintVector(fN[i], fArrivalBar[i], "normalized Arrival time");
    }
    
}


void DehoopBaiT::GenerateWaves()
{
    if (fMaxLevel>0) {
        DeleteWaves();
    }
    
    fMinLevel=abs(fZ_layer-fS_layer)+1;
    fMaxLevel=fMinLevel+2;
    
    fLevel=new int[fMaxLevel];
    fN=new int[fMaxLevel];
    
    fID=new int*[fMaxLevel];
    fType=new int*[fMaxLevel];
    fDirection=new int*[fMaxLevel];
    fLayer=new int*[fMaxLevel];
    fC=new double*[fMaxLevel];
    fd=new double*[fMaxLevel];
    
    fArrival=new double*[fMaxLevel];
    
    fdBar=new double* [fMaxLevel];
    fCBar=new double* [fMaxLevel];
    
    fArrivalBar=new double*[fMaxLevel];
    
    //generate IDs of the waves
    for (int i=0; i<fMaxLevel; i++) {
        fID[i]=NULL;
        fType[i]=NULL;
        fDirection[i]=NULL;
        fLayer[i]=NULL;
        fC[i]=NULL;
        fd[i]=NULL;
        fArrival[i]=NULL;
        
        fCBar[i]=NULL;
        fdBar[i]=NULL;
        fArrivalBar[i]=NULL;
    }
    
    fID[0]=new int[4];
    fLayer[0]=new int[4];
    
    fLevel[0]=1;
    fN[0]=4;
    fID[0][0]=1; fID[0][1]=2; fID[0][2]=3; fID[0][3]=4;

    fLayer[0][0]=fS_layer; fLayer[0][1]=fS_layer;
    fLayer[0][2]=fS_layer; fLayer[0][3]=fS_layer;
    
    //MathOperationT::PrintMatrix(N[0], Level[0], ID[0], "wave 1");
    //MathOperationT::PrintMatrix(N[0], Level[0], Layer[0], "layer of wave 1");
    
    for (int i=1; i<fMaxLevel; i++) {
        WaveExpand(fLevel[i-1], fN[i-1], fID[i-1], fLayer[i-1], fLevel[i], fN[i], fID[i], fLayer[i]);
        
        /*cout <<"level is " << i <<endl;
         MathOperationT::PrintMatrix(N[i], Level[i], ID[i], "wave 1");
         MathOperationT::PrintMatrix(N[i], Level[i], Layer[i], "layer of wave 1");*/
    }
    
    //MathOperationT::PrintVector(max_level, N, "N");
    
    //extract waves that arrive at the receiver point
    
    int* N_add=new int[fMaxLevel]; //Number of waves arriving at the receiver point;
    
    for (int i=0; i<fMaxLevel; i++) {
        N_add[i]=0;
        for (int j=0; j<fN[i]; j++) {
            
            bool add_flag=0;
            
            int* layer_ptr=fLayer[i]+j*fLevel[i];
            int* id_ptr=fID[i]+j*fLevel[i];
            
            if (fLevel[i]==1) {//Speical treatement to direct waves
                if (layer_ptr[0]==fZ_layer) {
                    if (fZ>=fS) {
                        if (id_ptr[0]==3 || id_ptr[0]==4) {
                            add_flag=1;
                        }
                    }else if(fZ<fS){
                        if (id_ptr[0]==1 || id_ptr[0]==2) {
                            add_flag=1;
                        }
                    }
                }
            }else if (fLevel[i]>1){
                if (layer_ptr[fLevel[i]-1]==fZ_layer) {
                    add_flag=1;
                }
            }
            
            if (add_flag==1) {
                
                if (j>N_add[i]) {
                    MathOperationT::MemCopy(fLevel[i], fID[i]+fLevel[i]*j, fID[i]+fLevel[i]*N_add[i]);
                    MathOperationT::MemCopy(fLevel[i], fLayer[i]+fLevel[i]*j, fLayer[i]+fLevel[i]*N_add[i]);
                }
                
                N_add[i]++;
            }
            
        }
        
        fN[i]=N_add[i];
        
        //cout <<"level is " << i+1 << " N_add is " << N_add[i] << endl;
        //MathOperationT::PrintMatrix(fN[i], fLevel[i], fID[i], "wave 1");
        //MathOperationT::PrintMatrix(fN[i], fLevel[i], fLayer[i], "layer of wave 1");
    }
    
    delete [] N_add;
    
    
    //determine wave type, wave direction, wave speeds
    for (int i=0; i<fMaxLevel; i++) {
        fType[i]=new int[fN[i]*fLevel[i]];
        fDirection[i]=new int[fN[i]*fLevel[i]];
        
        fC[i]=new double[fN[i]*fLevel[i]];
        fd[i]=new double[fN[i]*fLevel[i]];
        fArrival[i]=new double[fN[i]];
        
        fCBar[i]=new double[fN[i]*fLevel[i]];
        fdBar[i]=new double[fN[i]*fLevel[i]];
        fArrivalBar[i]=new double[fN[i]];
        
        for (int j=0; j<fN[i]; j++) {
            
            int* type_ptr=fType[i]+j*fLevel[i];
            int* Direction_ptr=fDirection[i]+j*fLevel[i];
            
            double* c_pointer=fC[i]+j*fLevel[i];
            int* id_pointer=fID[i]+j*fLevel[i];
            int* layer_pointer=fLayer[i]+j*fLevel[i];
            
            double* cBar_ptr=fCBar[i]+j*fLevel[i];
            
            for (int k=0; k<fLevel[i]; k++) {
                
                int id=id_pointer[k];
                int layer=layer_pointer[k];
                
                if (id==1) {
                    //upgoing P-wave
                    type_ptr[k]=1;
                    Direction_ptr[k]=1;
                    c_pointer[k]=fCp[layer];
                }else if (id==2){
                    //upgoing S-wave
                    type_ptr[k]=2;
                    Direction_ptr[k]=1;
                    c_pointer[k]=fCs[layer];
                }else if (id==3){
                    //downgoing p-wave
                    type_ptr[k]=1;
                    Direction_ptr[k]=2;
                    c_pointer[k]=fCp[layer];
                }else if (id==4){
                    //downgoing s-wave
                    type_ptr[k]=2;
                    Direction_ptr[k]=2;
                    c_pointer[k]=fCs[layer];
                }
                
                cBar_ptr[k]=c_pointer[k]/fCs_ref;
            }
        }
        
        //MathOperationT::PrintMatrix(fN[i], fLevel[i], fC[i], "speed of wave 1");

    }
    
    
    /*code to debug WaveExpand*/
    /*
     int Level1=1;
     int N1=4;
     int ID1[4]={1,2,3,4};
     int Layer1[4]={0,0,0,0};
     
     int Level2;
     int N2;
     int* ID2=NULL;
     int* Layer2=NULL;
     
     
     WaveExpand(Level1, N1, ID1, Layer1, Level2, N2, ID2, Layer2);
     
     MathOperationT::PrintMatrix(N1, Level1, ID1, "wave 1");
     MathOperationT::PrintMatrix(N1, Level1, Layer1, "layer of wave 1");
     MathOperationT::PrintMatrix(N2, Level2, ID2, "wave 2");
     MathOperationT::PrintMatrix(N2, Level2, Layer2, "layer of wave 2");
     */
    
}


void DehoopBaiT::DeleteWaves()
{
    if (fLevel) delete [] fLevel;
    if (fN) delete [] fN;
    
    for (int i=0; i<fMaxLevel; i++) {
        
        if (fID[i]) delete [] fID[i];
        if (fType[i]) delete[] fType[i];
        if (fDirection[i]) delete [] fDirection[i];
        if (fLayer[i]) delete [] fLayer[i];
        if (fC[i]) delete [] fC[i];
        if (fd[i]) delete [] fd[i];
        if (fArrival[i]) delete [] fArrival[i];
        
        if (fdBar[i]) delete [] fdBar[i];
        if (fCBar[i]) delete [] fCBar[i];
        if (fArrivalBar[i]) delete [] fArrivalBar[i];
        
        /*delete [] fID[i];
        delete[] fType[i];
        delete [] fDirection[i];
        delete [] fLayer[i];
        delete [] fC[i];
        delete [] fd[i];
        delete [] fArrival[i];
        
        delete [] fdBar[i];
        delete [] fCBar[i];
        delete [] fArrivalBar[i];*/
        
    }
    
    if (fID) delete [] fID;
    if (fType) delete [] fType;
    if (fDirection) delete [] fDirection;
    if (fLayer) delete [] fLayer;
    if (fC) delete [] fC;
    if (fd) delete [] fd;
    if (fArrival) delete [] fArrival;

    if (fdBar) delete [] fdBar;
    if (fCBar) delete [] fCBar;
    if (fArrivalBar) delete [] fArrivalBar;
}


void DehoopBaiT::Compute()
{
    /***********compute the dynamic part*********************/
    
    //Loop over time
    for (int i=0; i<fNumT_Heavi; i++) {
        
        double tBar=ft_Heavi[i]/ft_ref;
        //tBar=20;

        std::complex<double> I5[5]={0,0,0,0,0};
        std::complex<double> J10[10]={0,0,0,0,0,0,0,0,0,0};
        
        //Loop over wave groups
        for (int Li=0; Li<fMaxLevel; Li++) {
            
            int ns=fLevel[Li]; //number of segments
            
            //Loop over wave
            for (int wi=0; wi<fN[Li]; wi++) {
                
                if (tBar>fArrivalBar[Li][wi]) {
                    int* ID=fID[Li]+wi*ns;
                    int* Type=fType[Li]+wi*ns;
                    int* Direction=fDirection[Li]+wi*ns;
                    int* Layer=fLayer[Li]+wi*ns;
                    double* dBar=fdBar[Li]+wi*ns;//normalized vertical distance
                    double* CBar=fCBar[Li]+wi*ns; //normalized wave speed
                    
                    //MathOperationT::PrintVector(ns, ID, "ID");
                    //MathOperationT::PrintVector(ns, Type, "Type");
                    //MathOperationT::PrintVector(ns, Direction, "Direction");
                    //MathOperationT::PrintVector(ns, Layer, "Layer");
                    
                    int sample_size=50;
                    
                    std::complex<double> contour_point(1.0, -0.5);
                    contour_point=contour_point*M_PI/2.0;
                    
                    std::complex<double> dPhi=contour_point/(sample_size+0.0);
                    
                    for (int ii=0; ii<sample_size; ii++) {
                        std::complex<double> phi=(ii+0.5)*dPhi;
                        std::complex<double> cPhi=cos(phi);
                        
                        std::complex<double> xi, dxi_dt;
                     
                        Find_xi_Newton(tBar, frBar, ns, dBar, CBar, cPhi, xi, dxi_dt);
                        
                        std::complex<double> M[5], N[10];
                        WaveIntegrand(xi, ns, ID, Type, Direction, Layer, phi, M, N);
                        
                        //MathOperationT::PrintMatrix(5, 1, M, "M");
                        //MathOperationT::PrintMatrix(10, 1, N, "N");
                        
                        
                        for (int mi=0; mi<5; mi++) {
                            I5[mi]+=2.0*M[mi]*dxi_dt*dPhi;
                        }
                        for (int ni=0; ni<10; ni++) {
                            J10[ni]+=2.0*N[ni]*dxi_dt*dPhi;
                        }
                    }
                    
                    
                }//end if time is larger than the arrival time
                

                
                
            }//end loop over wave
        }//end loop over loop wave group
        
        //MathOperationT::PrintMatrix(5, 1, I5, "I5");
        //MathOperationT::PrintMatrix(10, 1, J10, "J10");
        
        double I5_real[5];
        for (int ii=0; ii<5; ii++) {
            I5_real[ii]=real(I5[ii])/(fMu_ref*fa);
        }
        
        double U_X[3], U_Y[3], U_Z[3];
        
        U_X[0]=fcTheta*I5_real[0];
        U_X[1]=-fsTheta*I5_real[1];
        U_X[2]=fcTheta*I5_real[2];
        
        U_Y[0]=fsTheta*I5_real[0];
        U_Y[1]=fcTheta*I5_real[1];
        U_Y[2]=fsTheta*I5_real[2];
        
        U_Z[0]=I5_real[3];
        U_Z[2]=I5_real[4];
        
        double* U_H_ptr=fDis_Heavi+9*i;
        U_H_ptr[0]=U_X[0]*fcTheta-U_X[1]*fsTheta;
        U_H_ptr[3]=U_X[0]*fsTheta+U_X[1]*fcTheta;
        U_H_ptr[6]=U_X[2];
        
        U_H_ptr[1]=U_Y[0]*fcTheta-U_Y[1]*fsTheta;
        U_H_ptr[4]=U_Y[0]*fsTheta+U_Y[1]*fcTheta;
        U_H_ptr[7]=U_Y[2];
        
        U_H_ptr[2]=U_Z[0]*fcTheta;
        U_H_ptr[5]=U_Z[0]*fsTheta;
        U_H_ptr[8]=U_Z[2];
        
        //MathOperationT::PrintMatrix(3, U_H_ptr, "U_H");
        
        //compute part 1 of the stress
        double J10_real[10];
        
        for (int ji=0; ji<10; ji++) {
            J10_real[ji]=real(J10[ji])/(fCs_ref*fa);
        }
        
        
        double SX_cylin[9], SY_cylin[9], SZ_cylin[9];

        //loading in X-direction
        SX_cylin[0]=fcTheta*J10_real[3];
        SX_cylin[4]=fcTheta*J10_real[4];
        SX_cylin[8]=fcTheta*J10_real[2];
        SX_cylin[1]=fsTheta*J10_real[5];
        SX_cylin[2]=fcTheta*J10_real[0];
        SX_cylin[5]=-fsTheta*J10_real[1];
        SX_cylin[3]=SX_cylin[1];
        SX_cylin[6]=SX_cylin[2];
        SX_cylin[7]=SX_cylin[5];
        
        //loading in Y-direction
        SY_cylin[0]=fsTheta*J10_real[3];
        SY_cylin[4]=fsTheta*J10_real[4];
        SY_cylin[8]=fsTheta*J10_real[2];
        SY_cylin[1]=-fcTheta*J10_real[5];
        SY_cylin[2]=fsTheta*J10_real[0];
        SY_cylin[5]=fcTheta*J10_real[1];
        SY_cylin[3]=SY_cylin[1];
        SY_cylin[6]=SY_cylin[2];
        SY_cylin[7]=SY_cylin[5];
        
        //loading in Z-direction
        SZ_cylin[0]=J10_real[8];
        SZ_cylin[4]=J10_real[9];
        SZ_cylin[8]=J10_real[7];
        SZ_cylin[2]=J10_real[6];
        SZ_cylin[6]=SZ_cylin[2];
        SZ_cylin[1]=SZ_cylin[3]=SZ_cylin[5]=SZ_cylin[7]=0.0;
        
        
        //MathOperationT::PrintMatrix(3, SX_cylin, "SX_Cylin");
        //MathOperationT::PrintMatrix(3, SY_cylin, "SY_Cylin");
        //MathOperationT::PrintMatrix(3, SZ_cylin, "SZ_Cylin");
        
        double R[9]={fcTheta, -fsTheta, 0.0, fsTheta,  fcTheta, 0.0, 0.0, 0.0, 1.0};
        
        double* SX_1=fSX_Heavi_p1+9*i;
        double* SY_1=fSY_Heavi_p1+9*i;
        double* SZ_1=fSZ_Heavi_p1+9*i;
        
        MathOperationT::multABAT(3, R, SX_cylin, SX_1);
        MathOperationT::multABAT(3, R, SY_cylin, SY_1);
        MathOperationT::multABAT(3, R, SZ_cylin, SZ_1);

        //MathOperationT::PrintMatrix(3, SX_1, "SX_1");
        //MathOperationT::PrintMatrix(3, SY_1, "SY_1");
        //MathOperationT::PrintMatrix(3, SZ_1, "SZ_1");
        
        //compute part 2 of the stress
        
        double mu=fMu[fZ_layer];
        
        double* SX_2=fSX_Heavi_p2+9*i;
        double* SY_2=fSY_Heavi_p2+9*i;
        double* SZ_2=fSZ_Heavi_p2+9*i;
        
        if (fr>0){
            //loading in X-direction
            SX_cylin[0]=-(2.0*mu/fr)*fcTheta*(I5_real[0]-I5_real[1]);
            SX_cylin[4]=-SX_cylin[0];
            SX_cylin[1]=-(2.0*mu/fr)*fsTheta*(I5_real[0]-I5_real[1]);
            SX_cylin[3]=SX_cylin[1];
            SX_cylin[2]=SX_cylin[6]=SX_cylin[5]=SX_cylin[7]=SX_cylin[8]=0.0;
            
            
            //loading in Y-direction
            SY_cylin[0]=-(2.0*mu/fr)*fsTheta*(I5_real[0]-I5_real[1]);
            SY_cylin[4]=-SY_cylin[0];
            SY_cylin[1]=(2.0*mu/fr)*fcTheta*(I5_real[0]-I5_real[1]);
            SY_cylin[3]=SY_cylin[1];
            SY_cylin[2]=SY_cylin[6]=SY_cylin[5]=SY_cylin[7]=SY_cylin[8]=0.0;

            //loading in Z-direction
            SZ_cylin[0]=-(2.0*mu/fr)*I5_real[3];
            SZ_cylin[4]=-SZ_cylin[0];
            SZ_cylin[1]=SZ_cylin[2]=SZ_cylin[3]=SZ_cylin[5]=SZ_cylin[6]=SZ_cylin[7]=SZ_cylin[8]=0.0;
            
            MathOperationT::multABAT(3, R, SX_cylin, SX_2);
            MathOperationT::multABAT(3, R, SY_cylin, SY_2);
            MathOperationT::multABAT(3, R, SZ_cylin, SZ_2);
            
            //MathOperationT::PrintMatrix(3, SX_2, "SX_2");
            //MathOperationT::PrintMatrix(3, SY_2, "SY_2");
            //MathOperationT::PrintMatrix(3, SZ_2, "SZ_2");
            
        }else if (fr==0){
            MathOperationT::VecSet(9, SX_2, 0.0);
            MathOperationT::VecSet(9, SY_2, 0.0);
            MathOperationT::VecSet(9, SZ_2, 0.0);
        }

        
    }//end loop over time for Heaviside loading
    
    
    Heavi2BSpline();
    
    /***********compute the static part*********************/
    ComputeStaticPart();
}


void DehoopBaiT::ComputeStaticPart()
{
    if (fNumLayer==1) {
        //case of half-space
        double source[3]={0,0,fS};
        double dest[3]={fX,fY,fZ};
        
        fGuzinaBiMat.SetMaterial(0, fNu[0], fE[fS_layer], fNu[fS_layer]);
        fGuzinaBiMat.SetGeometry(source, dest, fNormal);
    }else if (fNumLayer>1){
        //case of Multi-layered half-space
        
        if (fS_layer==0) {
            //Source is in the top layer of the half-space
            if (fS<0.5*fZ_Interface[fS_layer]) {
                double source[3]={0,0,fS};
                double dest[3]={fX,fY,fZ};
                
                fGuzinaBiMat.SetMaterial(0, fNu[0], fE[fS_layer], fNu[fS_layer]);
                fGuzinaBiMat.SetGeometry(source, dest, fNormal);
            }else if (fS>=0.5*fZ_Interface[fS_layer]){
                double source[3]={0,0,fS-fZ_Interface[fS_layer+1]};
                double dest[3]={fX,fY,fZ-fZ_Interface[fS_layer+1]};
                
                fGuzinaBiMat.SetMaterial(fE[fS_layer+1], fNu[fS_layer+1], fE[fS_layer], fNu[fS_layer]);
                fGuzinaBiMat.SetGeometry(source, dest, fNormal);
            }else{
                throw "!!! Source is not in the top layer !!!";
            }
            
        }else if (fS_layer==fNumLayer-1){
            //Source is in the bottome layer
            double source[3]={0,0,fS-fZ_Interface[fS_layer]};
            double dest[3]={fX,fY,fZ-fZ_Interface[fS_layer]};
            
            fGuzinaBiMat.SetMaterial(fE[fS_layer-1], fNu[fS_layer-1], fE[fS_layer], fNu[fS_layer]);
            fGuzinaBiMat.SetGeometry(source, dest, fNormal);
        }else{
            //Source is not in either top nor bottom layers
            if (fS<0.5*(fZ_Interface[fS_layer]+fZ_Interface[fS_layer+1])) {
                //source is in the upper half of the layer
                double source[3]={0,0,fS-fZ_Interface[fS_layer]};
                double dest[3]={fX,fY,fZ-fZ_Interface[fS_layer]};
                
                fGuzinaBiMat.SetMaterial(fE[fS_layer-1], fNu[fS_layer-1], fE[fS_layer], fNu[fS_layer]);
                fGuzinaBiMat.SetGeometry(source, dest, fNormal);
            }else if (fS>=0.5*(fZ_Interface[fS_layer]+fZ_Interface[fS_layer+1])){
                //source is in the lower half of the layer
                double source[3]={0,0,fS-fZ_Interface[fS_layer+1]};
                double dest[3]={fX,fY,fZ-fZ_Interface[fS_layer+1]};
                
                fGuzinaBiMat.SetMaterial(fE[fS_layer+1], fNu[fS_layer+1], fE[fS_layer], fNu[fS_layer]);
                fGuzinaBiMat.SetGeometry(source, dest, fNormal);
            }else{
                throw "!!! source is not in the middle layers !!!";
            }
        }
    }else{
        throw "!!! fNumLayer < 1 !!!!!!";
    }
    
    //computation
    fGuzinaBiMat.Compute();
    
}


void DehoopBaiT::Heavi2BSpline()
{
    if (fDT==0.0) {
        throw "!!!!Width of B-Spline loading is 0!!!";
    }
    
    for (int i=0; i<fNumT_BSpline*9; i++) {
        fDis_BSpline[i]=0;
        
        fSX_BSpline[i]=0;
        fSY_BSpline[i]=0;
        fSZ_BSpline[i]=0;
    }
    
    for (int i=0; i<fNumT_BSpline; i++) {
        
        double* dis_BSpline=fDis_BSpline+i*9;
        
        double* SX_BSpline=fSX_BSpline+i*9;
        double* SY_BSpline=fSY_BSpline+i*9;
        double* SZ_BSpline=fSZ_BSpline+i*9;
        
        double* Trac_BSpline=fTrac_BSpline+i*9;
        
        double B1, DB1, DDB1;
        double t_temp=ft_BSpline[i];
        
        int sum_num=0;
        
        for (int j=0; j<fNumT_Heavi; j++) {
            
            if (ft_Heavi[j]>=t_temp-fDT && ft_Heavi[j]<=t_temp) {
                sum_num += 1;
                CubicB1(B1, DB1, DDB1, t_temp-ft_Heavi[j], fDT);
                
                double* dis_Heavi=fDis_Heavi+j*9;
                
                double* SX_Heavi_p1=fSX_Heavi_p1+j*9;
                double* SY_Heavi_p1=fSY_Heavi_p1+j*9;
                double* SZ_Heavi_p1=fSZ_Heavi_p1+j*9;
                
                double* SX_Heavi_p2=fSX_Heavi_p2+j*9;
                double* SY_Heavi_p2=fSY_Heavi_p2+j*9;
                double* SZ_Heavi_p2=fSZ_Heavi_p2+j*9;
                
                for (int d=0; d<9; d++) {
                    dis_BSpline[d] += dis_Heavi[d]*DB1;
                    
                    SX_BSpline[d] += SX_Heavi_p1[d]*DDB1+SX_Heavi_p2[d]*DB1;
                    SY_BSpline[d] += SY_Heavi_p1[d]*DDB1+SY_Heavi_p2[d]*DB1;
                    SZ_BSpline[d] += SZ_Heavi_p1[d]*DDB1+SZ_Heavi_p2[d]*DB1;
                }
                
            }
            
        }
        
        if (sum_num==0) {
            cout <<"!!!Number of points in the time convolution is 0.!!!" <<endl;
        }
        
        for (int d=0; d<9; d++) {
            
            if (t_temp<fDT) {
                dis_BSpline[d] = dis_BSpline[d] * t_temp/sum_num;
                
                SX_BSpline[d] = SX_BSpline[d] * t_temp/sum_num;
                SY_BSpline[d] = SY_BSpline[d] * t_temp/sum_num;
                SZ_BSpline[d] = SZ_BSpline[d] * t_temp/sum_num;
                
            }else{
                dis_BSpline[d] = dis_BSpline[d] * fDT/sum_num;
                
                SX_BSpline[d] = SX_BSpline[d] * fDT/sum_num;
                SY_BSpline[d] = SY_BSpline[d] * fDT/sum_num;
                SZ_BSpline[d] = SZ_BSpline[d] * fDT/sum_num;
            }
            
        }
        
        //compute Traction based on stress
        double trac_temp[3];
        
        //if (abs(SX_BSpline[0])>0) {
            //MathOperationT::PrintMatrix(3, SX_BSpline, "SX_BSpline");
        //}
        
        MathOperationT::multAB(3, 3, SX_BSpline, 3, 1, fNormal, trac_temp);
        Trac_BSpline[0]=trac_temp[0];
        Trac_BSpline[3]=trac_temp[1];
        Trac_BSpline[6]=trac_temp[2];
        
        MathOperationT::multAB(3, 3, SY_BSpline, 3, 1, fNormal, trac_temp);
        Trac_BSpline[1]=trac_temp[0];
        Trac_BSpline[4]=trac_temp[1];
        Trac_BSpline[7]=trac_temp[2];
        
        MathOperationT::multAB(3, 3, SZ_BSpline, 3, 1, fNormal, trac_temp);
        Trac_BSpline[2]=trac_temp[0];
        Trac_BSpline[5]=trac_temp[1];
        Trac_BSpline[8]=trac_temp[2];
        
    }

}



double DehoopBaiT::Heaviside(double t)
{
    if (t>=0)
        return 1.0;
    else
        return 0.0;
}


void DehoopBaiT::CubicB1(double& B1, double& DB1, double& DDB1, double t, double T)
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
    
    DDB1=    (Heaviside(t)-Heaviside(t-0.25))*64.0*t;
    DDB1+=(Heaviside(t-0.25)-Heaviside(t-0.5))*(-192.0*t+64.0);
    DDB1+=(Heaviside(t-0.5)-Heaviside(t-0.75))*(192.0*t-128.0);
    DDB1+=(Heaviside(t-0.75)-Heaviside(t-1.0))*(-64.0*t+64.0);
    DDB1=DDB1/(T*T);
    
}


void DehoopBaiT::GetDisplacement(int& num, double** displacement)
{
    num=fNumT_BSpline;
    *displacement=fDis_BSpline;
}

void DehoopBaiT::GetDisplacement_Heavi(int& num, double** time, double** displacement)
{
    num=fNumT_Heavi;
    *time=ft_Heavi;
    *displacement=fDis_Heavi;
}


void DehoopBaiT::GetTraction(int& num, double** trac)
{
    num=fNumT_BSpline;
    *trac=fTrac_BSpline;
}

const double* DehoopBaiT::GetStaticTraction()
{
    return fGuzinaBiMat.GetTraction();
}

void DehoopBaiT::GetStress_BSpline(int &num, double** time, double **SX, double **SY, double **SZ)
{
    num=fNumT_BSpline;
    
    *time=ft_BSpline;
    
    *SX=fSX_BSpline;
    *SY=fSY_BSpline;
    *SZ=fSZ_BSpline;
    
}

void DehoopBaiT::GetStress_Heavi(int& num, double** time, double** SX_1, double** SX_2, double** SY_1, double** SY_2,
                     double** SZ_1, double** SZ_2)
{
    num=fNumT_Heavi;
    
    *time=ft_Heavi;
    
    *SX_1=fSX_Heavi_p1; *SX_2=fSX_Heavi_p2;
    *SY_1=fSY_Heavi_p1; *SY_2=fSY_Heavi_p2;
    *SZ_1=fSZ_Heavi_p1; *SZ_2=fSZ_Heavi_p2;
}



//Level1, Level2=level of the wave. Each wave has L1 or L2 segements.
//N1, N2=number of waves in L1, i.e., waves contained in Id1.
//ID1, ID2=ID of wave groups. They are arranged into 1D arrays.
// 1=P, 2=S, 3=p, 4=s;
//Layer1, Layer2=layers in which each segment of the wave propagates.

void DehoopBaiT::WaveExpand(int Level1, int N1, int *ID1, int* Layer1, int& Level2, int& N2, int *& ID2, int*& Layer2)
{
    Level2=Level1+1;
    
    int N2_temp=N1*4;
    N2=0;
    
    if (ID2) delete[] ID2;
    ID2=new int[N2_temp*Level2];
    
    if (Layer2) delete [] Layer2;
    Layer2=new int[N2_temp*Level2];
    
    for (int i=0; i<N1; i++) {
        
        int* ID1_ptr=ID1+i*Level1;
        int* Layer1_ptr=Layer1+i*Level1;
        
        int layer_old=Layer1_ptr[Level1-1]; //layer of the last segment
        
        //No expansion for downwardly-wave in the half-space (i.e. layer fNumLayer-1)
        if (layer_old==fNumLayer-1) {
            if (ID1_ptr[Level1-1]==3 || ID1_ptr[Level1-1]==4) {
                continue;
            }
        }
        
        for (int j=1; j<=4; j++) {
            
            //determine the layer of the new segment
            
            int layer_new=-1;

            if (ID1_ptr[Level1-1]==1 || ID1_ptr[Level1-1]==2) {
                if (j==1 || j==2) {
                    //upwardly transmission
                    layer_new=layer_old-1;
                }else{
                    //reflection
                    layer_new=layer_old;
                }
            }else{
                if (j==1 || j==2) {
                    //reflection
                    layer_new=layer_old;
                }else{
                    //downwardly transmission
                    layer_new=layer_old+1;
                }
            }
            
            //wave is not admittable if it propagates outside the layered half-space
            if (layer_new>=fNumLayer || layer_new<0) {
                continue;
            }
            
            //Expand the ID
            int* ID2_ptr=ID2+N2*Level2;
            MathOperationT::MemCopy(Level1, ID1_ptr, ID2_ptr);
            ID2_ptr[Level1]=j;
            
            //Expand the layer
            int* Layer2_ptr=Layer2+N2*Level2;
            MathOperationT::MemCopy(Level1, Layer1_ptr, Layer2_ptr);
            Layer2_ptr[Level1]=layer_new;
            
            N2 ++;
            
        }
    }
    
}



int DehoopBaiT::WhichLayer(double zcoord)
{
    if (zcoord< fZ_Interface[0])
        throw "WhichLayer: point is not in the half-space";
    
    int layer_id=fNumLayer-1;
    
    for (int i=0; i<=fNumLayer-2; i++) {
        if (zcoord>=fZ_Interface[i] && zcoord<fZ_Interface[i+1]) {
            layer_id=i;
        }
    }

    return layer_id;
}



double DehoopBaiT::WaveArrival(int ns, double r, double *d, double *c, int* id, int* layer)
{
    //special case of direct wave
    if (ns==1) {
        return sqrt(r*r+d[0]*d[0])/c[0];
    }
    
    //special case when source and receiver alligned vertically
    if (r==0) {
        double arrival=0;
        
        for (int i=0; i<ns; i++) {
            arrival += d[i]/c[i];
        }
        
        return arrival;
    }
    
    //************general cases************************
    
    //find head wave speed
    //--step 1: find all the layers involved
    //--setp 2: choose the maximum cd in all involved layers
    int L_min=layer[0];
    int L_max=layer[0];
    
    for (int i=1; i<ns; i++) {
        if (layer[i]<L_max) L_min=layer[i];
        if (layer[i]>L_max) L_max=layer[i];
    }
    
    for (int i=0; i<ns-1; i++) {
        if ( layer[i]==L_min && ((id[i]==1 || id[i]==2)) ) {
            if (L_min != 0) {
                L_min--;
            }
        }
        
        if ( layer[i]==L_max && ((id[i]==3 || id[i]==4)) ) {
            L_max++;
        }
    }
    
    double head_c=fCp[L_min]; //potential head wave speed
    
    for (int i=1; i<=L_max; i++) {
        if (head_c<fCp[i]) {
            head_c=fCp[i];
        }
    }
    
    double c_max=c[0];
    
    for (int i=1; i<ns; i++) {
        if (c_max<c[i]) c_max=c[i];
    }
    
    //case of Head wave
    if (head_c > c_max) {
        double r_critical=0;
        
        for (int i=0; i<ns; i++) {
            double theta=asin(c[i]/head_c);
            if (theta<M_PI/2) {
                r_critical += d[i]*tan(theta);
            }
        }

        if (r>=r_critical) {
            double head_t=r/head_c;
            
            for (int i=0; i<ns; i++) {
                head_t += d[i]*sqrt(1.0/pow(c[i], 2.0)-1.0/pow(head_c,2));
            }
            
            //MathOperationT::PrintVector(ns, c, "wave speed");
            
            return head_t;
        }
        
    }

    //case of normal wave
    //----special case of source, or receiver, or both is on the layer interface
    int head_flag=0;
    head_c=0;
    double r_critical=0;
    
    if (d[0]==0 && d[ns-1]>0) {
        if (c[0]==c_max) {
            
            head_flag=1;
            head_c=c[0];
            
            for (int i=1; i<ns; i++) {
                if (c[i]>=c[0]) {
                    head_flag=0;
                    break;
                }
            }
        }
    }else if (d[0]>0 && d[ns-1]==0){
        if (c[ns-1]==c_max) {
            head_flag=1;
            head_c=c[ns-1];
            
            for (int i=0; i<ns-1; i++) {
                if (c[i]>=c[ns-1]) {
                    head_flag=0;
                    break;
                }
            }
        }
   }else if (d[0]==0 && d[ns-1]==0){
       
       double c_temp=c[0];
       if (c_temp<c[ns-1]) c_temp=c[ns-1];
       
       if (c_temp==c_max) {
           
           if (ns==2) {
               head_flag=1;
               head_c=c_temp;
           }else if (ns>2){
               head_flag=1;
               head_c=c_temp;
               
               for (int i=1; i<ns-1; i++) {
                   if (c[i]>=c_temp) {
                       head_flag=0;
                       break;
                   }
               }
           }else{
               head_flag=0;
           }
       }
        
    }
        
    
    if ( head_flag==1) {
        for (int i=0; i<ns; i++) {
            double theta=asin(c[i]/head_c);
            if (theta<M_PI/2) {
                r_critical += d[i]*tan(theta);
            }
            
        }
    }
    
    // compute arrival time
    double arrival=0;
    
    if (head_flag==1 && r>=r_critical){
        arrival=r/head_c;
        for (int i=0; i<ns; i++) {
            arrival+=d[i]*sqrt(1.0/pow(c[i],2)-1.0/pow(head_c,2));
        }
        
    }else{
        //--use Newton-Raphson method to solve the incidece angle, and then the
        //  travle time
        //Theta is the incidence angle;
        //x=sin(theta);
        //The key is to find x
        double x_max=c[0]/c_max;
        
        double x=0.5*x_max; //initial guess of x
        double f=r;
        
        for (int i=0; i<ns; i++) {
            double sin_theta=x*c[i]/c[0]; //Snell's law
            double cos_theta=sqrt(1.0-sin_theta*sin_theta);
            
            f -= d[i]*sin_theta/cos_theta;
        }
        
        int iter_n=1;
        
        while (abs(f/r)>1e-8) {
            double df=0;
            
            for (int i=0; i<ns; i++) {
                double sin_theta=x*c[i]/c[0]; //Snell's law
                double cos_theta=sqrt(1.0-sin_theta*sin_theta);
                
                 df=df-d[i]*(c[i]/c[0])/pow(cos_theta,3);
            }
            
            double dx=-f/df;
            
            if (x+dx>=x_max){
                x=(x+x_max)/2;
            }else{
                x += dx;
            }
            
            f=r;
            for (int i=0; i<ns; i++) {
                double sin_theta=x*c[i]/c[0]; //Snell's law
                double cos_theta=sqrt(1.0-sin_theta*sin_theta);
                
                f -= d[i]*sin_theta/cos_theta;
            }
            
            iter_n++;
            if (iter_n >=50){
                throw "!!!!Compute arrival time::slow convergence!!!!!";
            }

        }//end of while
        
        arrival=0;
        
        for (int i=0; i<ns; i++) {
            double sin_theta=x*c[i]/c[0]; //Snell's law
            double cos_theta=sqrt(1.0-sin_theta*sin_theta);
            
            arrival += d[i]/(cos_theta*c[i]);
        }
        

    }//end if wave is normal
        
    return arrival;
    
}


void DehoopBaiT::Find_xi_Newton(double t, double r, int ns, double *d, double *c, std::complex<double> cPhi,
                                std::complex<double> &xi, std::complex<double> &dxi_dt)
{

    //MathOperationT::PrintVector(ns, d, "d");
    //MathOperationT::PrintVector(ns, c, "c");
    
    std::complex<double>* alpha=new std::complex<double>[ns];
    double* K=new double[ns];
    double K_min=1.0/c[0];
    
    for (int i=0; i<ns; i++) {
        K[i]=1/c[i];
        
        if (K[i] < K_min) {
            K_min=K[i];
        }
        
        K[i]=K[i]*K[i];
    }
    
    std::complex<double> I(0.0, 1.0);
    
    xi=0.5*K_min*(1.0+I); //initial guess of xi
    
    std::complex<double> f_xi=-t-I*xi*r*cPhi;
    
    for (int i=0; i<ns; i++) {
        alpha[i]=sqrt(complex<double>(xi*xi+K[i]));
        f_xi=f_xi+d[i]*alpha[i];
    }

    
    int iter_n=1;
    
    while (abs(f_xi) > 1e-8*t) {
        
        std::complex<double> df=-I*r*cPhi;
        for (int i=0; i<ns; i++) {
            df=df+d[i]*xi/alpha[i];
        }
        
        std::complex<double> dxi=-f_xi/df;
        
        xi=xi+dxi;
        
        f_xi=-t-I*xi*r*cPhi;
        for (int i=0; i<ns; i++) {
            alpha[i]=sqrt(complex<double>(xi*xi+K[i]));
            f_xi=f_xi+d[i]*alpha[i];
        }
        
        iter_n=iter_n+1;
        if (iter_n>20){
            throw "--Find_xi_Newton:: slow convergence!!!";
        }

    }
    
    std::complex<double> dt_dxi=-I*r*cPhi;
    for (int i=0; i<ns; i++) {
        dt_dxi += d[i]*xi/alpha[i];
    }
    
    dxi_dt=1.0/dt_dxi;
    
    
    delete [] alpha;
    delete [] K;
}


void DehoopBaiT::WaveIntegrand(std::complex<double> xi, int ns, int *ID, int* Type, int* Direction, int *Layer,
                               std::complex<double> phi, std::complex<double> *M, std::complex<double> *N)
{
    /*MathOperationT::PrintVector(ns, ID, "ID");
    MathOperationT::PrintVector(ns, Type, "Type");
    MathOperationT::PrintVector(ns, Direction, "Direction");
    MathOperationT::PrintVector(ns, Layer, "Layer");*/
    
    int SH_flag=1;
    for (int i=0; i<ns; i++) {
        if (Type[i]==1) {
            SH_flag=0;
            break;
        }
    }
    
    //compute M(1), M(2) and M(3), which are related to loading in X,Y-direction
    std::complex<double> w_sv, w_sh;
    Source_X(xi, fMuBar[fS_layer], fRhoBar[fS_layer], fCpBar[fS_layer], fCsBar[fS_layer], ID[0], w_sv, w_sh);
    
    for (int i=0; i<ns-1; i++) {
        
        //find the two layers near the interface
        int L1, L2;
        
        if (Direction[i]==1) {
            L1=Layer[i]-1;
        }else{
            L1=Layer[i];
        }
        L2=L1+1;
        
        std::complex<double> SV_factor, SH_factor;
        
        Propagator(xi, L1, L2, Direction[i], Direction[i+1], Type[i], Type[i+1], SV_factor, SH_factor);
        
        w_sv=SV_factor*w_sv;
        
        if (SH_flag==1) {
            w_sh=SH_factor*w_sh;
        }
        
    }
    
    std::complex<double> W_SV[4]={0.0, 0.0, 0.0, 0.0};
    std::complex<double> W_SH[2]={0.0, 0.0};

    
    if (ID[ns-1]==1) {
        //P-wave
        W_SV[2]=w_sv;
    }else if (ID[ns-1]==2) {
        //S-wave
        W_SV[3]=w_sv;
        
        if (SH_flag==1) {
            W_SH[1]=w_sh;
        }
        
    }else if (ID[ns-1]==3) {
        //p-wave
        W_SV[0]=w_sv;
        
    }else if (ID[ns-1]==4) {
        //s-wave
        W_SV[1]=w_sv;
        
        if (SH_flag==1) {
            W_SH[0]=w_sh;
        }
    }
    
    std::complex<double> cc_sv[24];
    std::complex<double> cc_sh[6];
    
    CC(xi, cc_sv, cc_sh);
    
    std::complex<double> I(0.0, 1.0);
    std::complex<double> V_sv[6];
    std::complex<double> V_sh[3];
    
    
    MathOperationT::multAB(6, 4, cc_sv, 4, 1, W_SV, V_sv);
    
    M[0]=-V_sv[0]*(1.0+cos(2.0*phi));
    M[1]=-V_sv[0]*(1.0-cos(2.0*phi));
    M[2]=-2.0*I*V_sv[1]*cos(phi);
    
    N[0]=-V_sv[2]*(1.0+cos(2.0*phi));
    N[1]=-V_sv[2]*(1.0-cos(2.0*phi));
    N[2]=-2.0*I*V_sv[3]*cos(phi);
    N[3]=-2.0*I*V_sv[5]*cos(phi);
    N[4]=-2.0*I*V_sv[4]*cos(phi);
    
    if (SH_flag==1){
        MathOperationT::multAB(3, 2, cc_sh, 2, 1, W_SH, V_sh);
        
        M[0]=M[0]+V_sh[0]*(1.0-cos(2.0*phi));
        M[1]=M[1]+V_sh[0]*(1.0+cos(2.0*phi));
        
        N[0]=N[0]+V_sh[1]*(1.0-cos(2.0*phi));
        N[1]=N[1]+V_sh[1]*(1.0+cos(2.0*phi));
        N[5]=2.0*V_sh[2]*cos(phi);
    }

    /*MathOperationT::PrintMatrix(6, 4, cc_sv, "cc_sv");
    MathOperationT::PrintMatrix(4, 1, W_SV, "W_sv");
    MathOperationT::PrintMatrix(6, 1, V_sv, "V_sv");
    MathOperationT::PrintMatrix(5, 1, M, "M");
    MathOperationT::PrintMatrix(10, 1, N, "N");*/
    
    
    //compute M(4) and M(5), which are related to loading in Z-direction
    Source_Z(xi, fMuBar[fS_layer], fRhoBar[fS_layer], fCpBar[fS_layer], fCsBar[fS_layer], ID[0], w_sv);

    for (int i=0; i<ns-1; i++) {
        
        //find the two layers near the interface
        int L1, L2;
        
        if (Direction[i]==1) {
            L1=Layer[i]-1;
        }else{
            L1=Layer[i];
        }
        L2=L1+1;
        
        std::complex<double> SV_factor, SH_factor;
        
        Propagator(xi, L1, L2, Direction[i], Direction[i+1], Type[i], Type[i+1], SV_factor, SH_factor);
        
        w_sv=SV_factor*w_sv;
    }
    
    if (ID[ns-1]==1) {
        //P-wave
        W_SV[2]=w_sv;
    }else if (ID[ns-1]==2) {
        //S-wave
        W_SV[3]=w_sv;
        
    }else if (ID[ns-1]==3) {
        //p-wave
        W_SV[0]=w_sv;
        
    }else if (ID[ns-1]==4) {
        //s-wave
        W_SV[1]=w_sv;
    }
    
    MathOperationT::multAB(6, 4, cc_sv, 4, 1, W_SV, V_sv);
    
    M[3]=-I*V_sv[0]*cos(phi);
    M[4]=V_sv[1];
    
    N[6]=-I*V_sv[2]*cos(phi);
    N[7]=V_sv[3];
    N[8]=V_sv[5];
    N[9]=V_sv[4];
    
    for (int i=0; i<5; i++) {
        M[i]=M[i]*xi/M_PI;
    }
    for (int i=0; i<10; i++) {
        N[i]=N[i]*xi/M_PI;
    }

    /*MathOperationT::PrintMatrix(4, 1, W_SV, "W_sv");
    MathOperationT::PrintMatrix(6, 1, V_sv, "V_sv");
    MathOperationT::PrintMatrix(5, 1, M, "M");
    MathOperationT::PrintMatrix(10, 1, N, "N");*/
}

void DehoopBaiT::Source_X(std::complex<double> xi, double mu, double rho, double cd, double cs, int ID,
                          std::complex<double>& Source_sv, std::complex<double>& Source_sh)
{
    std::complex<double> alpha=sqrt(complex<double>(xi*xi+1.0/pow(cd,2.0)));
    std::complex<double> beta =sqrt(complex<double>(xi*xi+1.0/pow(cs,2.0)));
    
    std::complex<double> I(0.0, 1.0);
    
    if (ID==1) {
        //upgoing P-wave
        Source_sv=(1.0/(2.0*rho))*( (-xi/(2.0*alpha)) * (1.0/(2.0*M_PI)) );;
        Source_sh=0.0;
        
    }else if (ID==2){
        //upgoing S-wave
        Source_sv=(1.0/(2.0*rho))*( (1.0/2.0) * (1.0/(2.0*M_PI)) );
        Source_sh=(-I/2.0) * ( 1.0/(2.0*mu*xi*beta) * (1.0/(2.0*M_PI)) ) ;
        
    }else if (ID==3){
        //downgoing p-wave
        Source_sv=(-1.0/(2.0*rho))*( (xi/(2.0*alpha)) * (1.0/(2.0*M_PI)) );
        Source_sh=0.0;
        
    }else if (ID==4){
        //downgoing s-wave
        Source_sv=(-1.0/(2.0*rho))*( (1.0/2.0) * (1.0/(2.0*M_PI)) );
        Source_sh=(-I/2.0) * ( 1.0/(2.0*mu*xi*beta) * (1.0/(2.0*M_PI)) ) ;
        
    }
    
}


void DehoopBaiT::Source_Z(std::complex<double> xi, double mu, double rho, double cd, double cs, int ID, std::complex<double>& Source_sv)
{
    std::complex<double> beta =sqrt(complex<double>(xi*xi+1.0/pow(cs,2.0)));

    if (ID==1) {
        //upgoing P-wave
        Source_sv=(1.0/(2.0*rho))*(1.0/(2.0*M_PI));
        
    }else if (ID==2){
        //upgoing S-wave
        Source_sv=(1.0/(2.0*rho))*( (-xi/beta) * (1.0/(2.0*M_PI)) );
        
    }else if (ID==3){
        //downgoing p-wave
        Source_sv=(-1.0/(2.0*rho))* (1.0/(2.0*M_PI));
        
    }else if (ID==4){
        //downgoing s-wave
        Source_sv=(-1.0/(2.0*rho))*( (xi/beta) * (1.0/(2.0*M_PI)) );
    }

}



void DehoopBaiT::Propagator(std::complex<double>xi, int L1, int L2, int D1, int D2, int Type1, int Type2,
                std::complex<double>& SV_factor, std::complex<double>& SH_factor)
{
    double mu1, mu2;
    double rho1, rho2;
    double cd1, cd2;
    double cs1, cs2;
    std::complex<double> alpha1, alpha2;
    std::complex<double> beta1, beta2;
    
    if (L1==-1) {
        mu1=0;
        rho1=0;
        cd1=1;
        cs1=1;
    }else{
        mu1=fMuBar[L1];
        rho1=fRhoBar[L1];
        cd1=fCpBar[L1];
        cs1=fCsBar[L1];
    }
    
    mu2=fMuBar[L2];
    rho2=fRhoBar[L2];
    cd2=fCpBar[L2];
    cs2=fCsBar[L2];
    
    std::complex<double>  xi_s=xi*xi;
    alpha1=sqrt( std::complex<double>( xi_s+1.0/pow(cd1, 2.0) ) );
    beta1=sqrt( std::complex<double>( xi_s+1.0/pow(cs1, 2.0) ) );
    alpha2=sqrt( std::complex<double>( xi_s+1.0/pow(cd2, 2.0) ) );
    beta2=sqrt( std::complex<double>( xi_s+1.0/pow(cs2, 2.0) ) );
    
    std::complex<double> S;
    S  = 4.0*xi_s*pow(mu2-mu1,2.0)*(xi_s-alpha1*beta1)*(xi_s-alpha2*beta2);
    S +=-4.0*xi_s*(mu2-mu1)*( rho1*(xi_s-alpha2*beta2)-rho2*(xi_s-alpha1*beta1) );
    S += xi_s*pow(rho2-rho1,2.0)-(rho1*beta2+rho2*beta1)*(rho1*alpha2+rho2*alpha1);
    S *= -1;
    
    
    //computation of entries of transmission/reflection matrix
    SH_factor=0.0;
    
    if (D1==1 && D2==1) {
        //up-up transmission
        
        if (Type1==1 && Type2==1) {
            //P to P; TuN(1,1) is used;
            SV_factor=2.0*xi_s*(mu2-mu1)*(beta2-beta1)-(rho1*beta2+rho2*beta1);
            SV_factor=SV_factor*(-2.0*rho2*alpha2)/S;
            
        }else if (Type1==1 && Type2==2){
            //P to S; TuN(2,1) is used;
            SV_factor=2.0*xi*(mu2-mu1)*(xi_s-alpha1*beta2)+xi*(rho2-rho1);
            SV_factor=SV_factor*(-2.0*rho2*alpha2)/S;
            
        }else if (Type1==2 && Type2==1){
            //S to P; TuN(1,2) is used;
            SV_factor=2.0*xi*(mu2-mu1)*(xi_s-alpha2*beta1)+xi*(rho2-rho1);
            SV_factor=SV_factor*(-2.0*rho2*beta2)/S;
            
        }else if (Type1==2 && Type2==2){
            //S to S; TuN(2,2) is used;
            SV_factor=2.0*xi_s*(mu2-mu1)*(alpha2-alpha1)-(rho1*alpha2+rho2*alpha1);
            SV_factor=SV_factor*(-2.0*rho2*beta2)/S;
            
            SH_factor=2.0*mu2*beta2/(mu1*beta1+mu2*beta2);
        }
        
    }else if (D1==2 && D2==2) {
        //down-down transmission
        
        if (Type1==1 && Type2==1) {
            //P to P; TdN(1,1) is used;
            SV_factor=2.0*xi_s*(mu2-mu1)*(beta2-beta1)-(rho1*beta2+rho2*beta1);
            SV_factor=SV_factor*(-2*rho1*alpha1)/S;
            
        }else if (Type1==1 && Type2==2){
            //P to S; TdN(2,1) is used;
            SV_factor=2.0*xi*(mu2-mu1)*(xi_s-alpha2*beta1)+xi*(rho2-rho1);
            SV_factor=SV_factor*(-2.0*rho1*alpha1)/S;

        }else if (Type1==2 && Type2==1){
            //S to P; TdN(1,2) is used;
            SV_factor=2.0*xi*(mu2-mu1)*(xi_s-alpha1*beta2)+xi*(rho2-rho1);
            SV_factor=SV_factor*(-2.0*rho1*beta1)/S;

        }else if (Type1==2 && Type2==2){
            //S to S; TdN(2,2) is used;
            SV_factor=2.0*xi_s*(mu2-mu1)*(alpha2-alpha1)-(rho1*alpha2+rho2*alpha1);
            SV_factor=SV_factor*(-2.0*rho1*beta1)/S;
            
            SH_factor=2.0*mu1*beta1/(mu1*beta1+mu2*beta2);
        }
        
    }else if (D1==1 && D2==2) {
        //up-down reflection
        if (mu1==0) {
            //Reflection of the free surface
            std::complex<double> Rm=pow(beta2*beta2+xi_s,2.0)-4.0*xi_s*alpha2*beta2;
            std::complex<double> Rp=pow(beta2*beta2+xi_s,2.0)+4.0*xi_s*alpha2*beta2;
            
            if (Type1==1 && Type2==1) {
                //P to P; RuN(1,1) is used;
                SV_factor=-Rp/Rm;
                
            }else if (Type1==1 && Type2==2){
                //P to S; RuN(2,1) is used;
                SV_factor=-4.0*xi*alpha2*(beta2*beta2+xi_s)/Rm;
                
            }else if (Type1==2 && Type2==1){
                //S to P; RuN(1,2) is used;
                SV_factor=-4.0*xi*beta2*(beta2*beta2+xi_s)/Rm;
                
            }else if (Type1==2 && Type2==2){
                //S to S; RuN(2,2) is used;
                SV_factor=-Rp/Rm;
                
                SH_factor=1.0;
            }
            
            
        }else{
            //Reflection of other interface
            
            if (Type1==1 && Type2==1) {
                //P to P; RuN(1,1) is used;
                SV_factor=          4.0*xi_s*pow(mu2-mu1,2.0)*(xi_s-alpha1*beta1)*(xi_s+alpha2*beta2);
                SV_factor=SV_factor-4.0*xi_s*(mu2-mu1)*(rho1*(xi_s+alpha2*beta2)-rho2*(xi_s-alpha1*beta1));
                SV_factor=SV_factor+xi_s*pow(rho2-rho1,2.0)+(rho1*beta2+rho2*beta1)*(rho1*alpha2-rho2*alpha1);
                SV_factor=SV_factor/S;
                
            }else if (Type1==1 && Type2==2){
                //P to S; RuN(2,1) is used;
                SV_factor=          4.0*xi_s*pow(mu2-mu1,2)*(xi_s-alpha1*beta1);
                SV_factor=SV_factor-2.0*(mu2-mu1)*(2.0*rho1*xi_s-rho2*(xi_s-alpha1*beta1)-rho1*(rho2-rho1));
                SV_factor=SV_factor*2.0*xi*alpha2/S;
                
            }else if (Type1==2 && Type2==1){
                //S to P; RuN(1,2) is used;
                SV_factor=          4.0*xi_s*pow(mu2-mu1,2.0)*(xi_s-alpha1*beta1);
                SV_factor=SV_factor-2.0*(mu2-mu1)*(2.0*rho1*xi_s-rho2*(xi_s-alpha1*beta1)-rho1*(rho2-rho1));
                SV_factor=SV_factor*2.0*xi*beta2/S;
                
            }else if (Type1==2 && Type2==2){
                //S to S; RuN(2,2) is used;
                SV_factor=          4.0*xi_s*pow(mu2-mu1,2.0)*(xi_s-alpha1*beta1)*(xi_s+alpha2*beta2);
                SV_factor=SV_factor-4.0*xi_s*(mu2-mu1)*(rho1*(xi_s+alpha2*beta2)-rho2*(xi_s-alpha1*beta1));
                SV_factor=SV_factor+xi_s*pow(rho2-rho1,2.0)+(rho1*beta2-rho2*beta1)*(rho1*alpha2+rho2*alpha1);
                SV_factor=SV_factor/S;
                
                SH_factor=(mu2*beta2-mu1*beta1)/(mu1*beta1+mu2*beta2);
            }
        }
        
        
    }else if (D1==2 && D2==1) {
        //down-up reflection
        if (Type1==1 && Type2==1) {
            //P to P; RdN(1,1) is used;

            SV_factor=          4.0*xi_s*pow(mu2-mu1,2.0)*(xi_s+alpha1*beta1)*(xi_s-alpha2*beta2);
            SV_factor=SV_factor-4.0*xi_s*(mu2-mu1)*(rho1*(xi_s-alpha2*beta2)-rho2*(xi_s+alpha1*beta1));
            SV_factor=SV_factor+xi_s*pow(rho2-rho1,2.0)-(rho1*beta2+rho2*beta1)*(rho1*alpha2-rho2*alpha1);
            SV_factor=SV_factor/S;

        }else if (Type1==1 && Type2==2){
            //P to S; RdN(2,1) is used;
            
            SV_factor=         -4.0*xi_s*pow(mu2-mu1,2.0)*(xi_s-alpha2*beta2);
            SV_factor=SV_factor-2.0*(mu2-mu1)*(2.0*rho2*xi_s-rho1*(xi_s-alpha2*beta2)-rho2*(rho2-rho1));
            SV_factor=SV_factor*2.0*xi*alpha1/S;
            
        }else if (Type1==2 && Type2==1){
            //S to P; RdN(1,2) is used;
            
            SV_factor=         -4.0*xi_s*pow(mu2-mu1,2.0)*(xi_s-alpha2*beta2);
            SV_factor=SV_factor-2.0*(mu2-mu1)*(2.0*rho2*xi_s-rho1*(xi_s-alpha2*beta2)-rho2*(rho2-rho1));
            SV_factor=SV_factor*2.0*xi*beta1/S;

        }else if (Type1==2 && Type2==2){
            //S to S; RdN(2,2) is used;

            SV_factor=          4.0*xi_s*pow(mu2-mu1,2.0)*(xi_s+alpha1*beta1)*(xi_s-alpha2*beta2);
            SV_factor=SV_factor-4.0*xi_s*(mu2-mu1)*(rho1*(xi_s-alpha2*beta2)-rho2*(xi_s+alpha1*beta1));
            SV_factor=SV_factor+xi_s*pow(rho2-rho1,2.0)-(rho1*beta2-rho2*beta1)*(rho1*alpha2+rho2*alpha1);
            SV_factor=SV_factor/S;
            
            SH_factor=(mu1*beta1-mu2*beta2)/(mu1*beta1+mu2*beta2);
        }
    }

    
}

void DehoopBaiT::CC(std::complex<double> xi, std::complex<double>* cc_sv, std::complex<double>* cc_sh)
{
    double lambda=fLamdaBar[fZ_layer];
    double mu=fMuBar[fZ_layer];
    double cd=fCpBar[fZ_layer];
    double cs=fCsBar[fZ_layer];
    
    std::complex<double> xi_s=xi*xi;
    
    std::complex<double> alpha=sqrt(std::complex<double>(xi_s+1/pow(cd,2.0)));
    std::complex<double> beta =sqrt(std::complex<double>(xi_s+1/pow(cs,2.0)));
    std::complex<double> alpha_s=alpha*alpha;
    std::complex<double> beta_s=beta*beta;
    
    
    //////////////////////////////////////////
    /////////////////////////////////////////
    cc_sv[0]=-xi;
    cc_sv[1]=beta;
    cc_sv[2]=-xi;
    cc_sv[3]=-beta;
    
    cc_sv[4]=-alpha;
    cc_sv[5]=xi;
    cc_sv[6]=alpha;
    cc_sv[7]=xi;
    
    cc_sv[8]=2*mu*xi*alpha;
    cc_sv[9]=-mu*(beta_s+xi_s);
    cc_sv[10]=-cc_sv[8];
    cc_sv[11]=cc_sv[9];
    
    cc_sv[12]=mu*(beta_s+xi_s);
    cc_sv[13]=-2.0*mu*xi*beta;
    cc_sv[14]=cc_sv[12];
    cc_sv[15]=-cc_sv[13];
    
    cc_sv[16]=lambda*(alpha_s-xi_s);
    cc_sv[17]=0;
    cc_sv[18]=cc_sv[16];
    cc_sv[19]=0;
    
    cc_sv[20]=lambda*alpha_s-(lambda+2.0*mu)*xi_s;
    cc_sv[21]=2.0*mu*xi*beta;
    cc_sv[22]=cc_sv[20];
    cc_sv[23]=-cc_sv[21];
    
    //////////////////////////////////////////
    //////////////////////////////////////////
    std::complex<double> I(0.0, 1.0);
    cc_sh[0]=I*xi;
    cc_sh[1]=cc_sh[0];
    cc_sh[2]=-I*mu*xi*beta;
    cc_sh[3]=-cc_sh[2];
    cc_sh[4]=mu*xi_s;
    cc_sh[5]=cc_sh[4];
}


