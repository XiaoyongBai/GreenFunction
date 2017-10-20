!----------------------------------------------------------------------------------------------------------------------------------!
MODULE GFLA
!----------------------------------------------------------------------------------------------------------------------------------! 

USE KINDS 

integer(kind=k9) :: nl, ml, nrec, ns, ms, ls, lo, ka, kb,homog      
real   (kind=dp) :: rho0, cc(8,8), ks0, mu0            

integer(kind=k9), parameter :: nint=9  
real   (kind=dp), parameter :: o=1.0e-32_dp, eps=1.0e-9_dp, hbmin=8.0e-3_dp, magfac=1.0e9_dp     
complex(kind=dp), parameter :: pi=(0.07957747154594768_dp,0.0_dp), ic=(0.0_dp,1.0_dp), zero=(0.0_dp,0.0_dp) 
complex(kind=dp), parameter :: xpy(3,-1:1)=magfac*reshape((/ pi, ic*pi,zero, zero,zero,zero, pi,-ic*pi,zero/),(/3,3/)) 
complex(kind=dp), parameter :: xmy(3,-1:1)=magfac*reshape((/-pi,-ic*pi,zero, zero,zero,zero, pi,-ic*pi,zero/),(/3,3/))   
complex(kind=dp), parameter :: zzz(3,-1:1)=magfac*reshape((/zero,zero,zero, zero,zero,2.0_dp*pi, zero,zero,zero/),(/3,3/))

integer(kind=k9), allocatable :: flag(:)
real   (kind=dp), allocatable :: rho(:), hh(:), hb(:), zh(:), zzcoord(:)
complex(kind=dp), allocatable :: mu(:), lam(:) , pv(:)
real   (kind=dp), allocatable :: v(:) 
complex(kind=dp), allocatable :: zeta(:), zetas(:), kp(:), ks(:), alfa(:,:), beta(:,:), cmps(:,:,:,:), cmsh(:,:,:,:) 
complex(kind=dp), allocatable :: difa(:,:), difb(:,:), diab(:,:), dfaa(:,:), dfab(:,:), dfbb(:,:), difc(:,:), difd(:,:) 

complex(kind=dp), allocatable :: tups0(:,:,:,:), tdps0(:,:,:,:), rups0(:,:,:,:), rdps0(:,:,:,:)  
complex(kind=dp), allocatable :: tush0(:,:)    , tdsh0(:,:)    , rush0(:,:)    , rdsh0(:,:)  
complex(kind=dp), allocatable :: tups(:,:,:,:) , tdps(:,:,:,:) , rups(:,:,:,:) , rdps(:,:,:,:)  
complex(kind=dp), allocatable :: tush(:,:)     , tdsh(:,:)     , rush(:,:)     , rdsh(:,:)  
complex(kind=dp), allocatable :: thups(:,:,:,:), thdps(:,:,:,:), rhups(:,:,:,:), rhdps(:,:,:,:)  
complex(kind=dp), allocatable :: thush(:,:)    , thdsh(:,:)    , rhush(:,:)    , rhdsh(:,:)           

complex(kind=dp), allocatable :: wups(:,:,:,:), wdps(:,:,:,:), wush(:,:)  , wdsh(:,:)       
complex(kind=dp), allocatable :: sups(:,:,:,:), sdps(:,:,:,:), sush(:,:,:), sdsh(:,:,:)           
complex(kind=dp), allocatable :: mups(:,:,:,:), mdps(:,:,:,:), mush(:,:,:), mdsh(:,:,:), usmp(:,:,:,:), tsmp(:,:,:,:)              
complex(kind=dp), allocatable :: usmp_global(:,:,:,:)
!----------------------------------------------------------------------------------------------------------------------------------!
END MODULE GFLA 
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
MODULE GFLB
!----------------------------------------------------------------------------------------------------------------------------------!

USE KINDS 

integer(kind=k9) :: sing  
real   (kind=dp) :: sst, rst, tst, zst, stat0(3), statc(3), stats(3)  
complex(kind=dp) :: expp, expm, mst1, mst2, vst1, vst2                  

!----------------------------------------------------------------------------------------------------------------------------------!
END MODULE GFLB 
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
MODULE GFLC
!----------------------------------------------------------------------------------------------------------------------------------! 

USE KINDS  

real   (kind=dp)            :: dcont, hcont, lpath           
real   (kind=dp), parameter :: lcont=1.6_dp, lpmin=400_dp, lpmax=2.0e3_dp, hcont_min=0.01_dp, dcont_min=0.01_dp, hcont_max=0.30_dp

!----------------------------------------------------------------------------------------------------------------------------------!
END MODULE GFLC 
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
MODULE GFLJ
!----------------------------------------------------------------------------------------------------------------------------------! 

USE KINDS 

integer(kind=k9), parameter :: bn(0:7,0:7)=reshape((/1,1,1,1,1,1,1,1,0,1,2,3,4,5,6,7,0,0,1,3,6,10,15,21,0,0,0,1,4,10,20,35,0,0,0, &  
                    0,1,5,15,35,0,0,0,0,0,1,6,21,0,0,0,0,0,0,1,7,0,0,0,0,0,0,0,1/),(/8,8/)), n(0:2)=(/0,1,2/), jmax=200, kmax=10000  

!----------------------------------------------------------------------------------------------------------------------------------!
END MODULE GFLJ
!----------------------------------------------------------------------------------------------------------------------------------! 

!----------------------------------------------------------------------------------------------------------------------------------!
MODULE GFLR
!----------------------------------------------------------------------------------------------------------------------------------! 

USE KINDS  

integer(kind=k9) :: rflag  !-indicates whether series expansion of bessel fcs applies throughout the integration path (1=yes, 2=no) 
real   (kind=dp) :: rfac   !-normalization factor used for calculation of t_rr, t_tt, t_rt (rflag=1 => rfac=1, rflag=2 = rfac=1/r)

!----------------------------------------------------------------------------------------------------------------------------------!
END MODULE GFLR
!----------------------------------------------------------------------------------------------------------------------------------! 

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE HALDGREEN4 (dyngreen, kdgl, omega, xyzs, xyzo, ugreen, tdgreen, tsgreen)  
!----------------------------------------------------------------------------------------------------------------------------------!

!===purpose===!

!-Dynamic point-load Green's functions for the viscoelastic layered half-space 
!-Version 4.80 (limit r->0, complex v, SINPA debugged, new XINC and JINT)  May 23, 2001   
!-BY FARZAD IT IS MODIFIED FOR MULTI DOMAIN CALLS ALSO IN OBTAINING SING PARAMETER IN SUBROUTINE 
!GFL_SINPA : the small value of -.00000001*(zh(ls)-zh(ls-1)) is deleted(compare with haldgreen3 and halsgreen3)
!dyngreen=0 if only singular part is needed and dyngreen =1 if both singular and total green's function is needed. 
!  uses same normalization ks0 for singualr part and dynamic Green's function.
!===external variables===!

!--omega--= circular frequency of excitation [rad/sec]  
!---xyzs--= cartesian coordinates of the source point [m]
!---xyzo--= cartesian coordinates of the observation point [m]
!--ugreen-= displacement green's functions [m/MN]
!-tdgreen-= stress green's functions ==> regular  part [MPa/MN]      
!-tsgreen-= stress green's functions ==> singular part [MPa/MN]      
 
!===additional notes===! 

!-layer properties and integration parameters are red by the subroutine gfl_init from file "layers.dat" 
!-green's functions are evaluated in their original dimensional format (i.e. they are not normalized) 
!-first index of the green's functions refers to the source direction 

!~Copyright Bojan Guzina, 2001      

USE GFLA 

IMPLICIT NONE 

real   (kind=dp), intent(in ) :: omega, xyzs(3), xyzo(3)  
complex(kind=dp), intent(out) :: ugreen(3,3), tdgreen(3,6), tsgreen(3,6)
integer(kind=k9)              :: init=0, lsold=-1, kdgl, dyngreen  
real   (kind=dp)              :: xs, ys, zs, xo, yo, zo, s, r, t, z, omegaold=-1.0_dp, sold=-1.0_dp   
complex(kind=dp)              :: uhat(3,3), tdhat(3,6), tshat(3,6)
real   (kind=dp)              :: radius

integer :: iii;

!~internal spatial coordinates 

xs = xyzs(1)
ys = xyzs(2) 
zs = xyzs(3)  

xo = xyzo(1)
yo = xyzo(2) 
zo = xyzo(3)   

radius=sqrt((xs-xo)**2+(ys-yo)**2)

!~initialization procedure at the start of the program

if (init/=kdgl) then
  if (init/=0) call GFL_DEALLOC  !-if 2nd domain or greater, deallocate variables before reallocating them
  call GFL_INIT   !-memory allocation, initialization of featuring arrays 

  init = kdgl 

  omegaold = -1.0_dp
  sold = -1.0_dp
  lsold = -1
  
end if  

!~evaluation of the transmission/reflection matrices for each frequency   
 
if (omegaold/=omega) then

  call GFL_FDPAR (omega, radius)   !-evaluation of the frequency-dependent parameters

  call GFL_CONTO           !-evaluation of the integration variable zeta and radicals alfa, beta

  call GFL_AUXAR           !-evaluation of the auxiliary arrays 

  call GFL_CFMAT           !-evaluation of the coeffficient matrices

  call GFL_RFMAT           !-evaluation of the reference matrices

  call GFL_TRMAT           !-evaluation of the transmission/reflection matrices


end if  

!~conversion from cartesian to dimensionless cylindrical coordinates

call GFL_COORD (xs, ys, zs, xo, yo, zo, ks0, s, r, t, z) 

!~calculate number of the layer (ls) containing the source 

if (sold/=s) call GFL_LSNUM (zs)   

!~evaluation of propagation matrices for each ls   

if (omegaold/=omega.or.lsold/=ls) then 

  call GFL_PRMAT    !-evaluation of propagation matrices 

  call GFL_FCMAT    !-evaluation of factorized matrices 

  lsold = ls  

end if  

!~evaluation of the loading coefficients for each source 

if (omegaold/=omega.or.sold/=s) then

  call GFL_SORC1 (s)  !-evaluation of the loading coefficients    

  call GFL_SORC2      !-evaluation of the loading coefficients    

  omegaold = omega;  sold = s 

end if   

!~calculate number of the layer (lo) containing the observation point  

call GFL_LONUM (zo)  

!~calculate parameters needed for the singularity treatment 

call GFL_SINPA (s, r, t, z, zo, zs)  

!~evaluation of the auxiliary exponential/trigonomteric functions of angular coordinate 

call GFL_THETA (t) 

!~adaptive integration procedure (numerical inverse hankel transform)  

call GFL_XHAT (s, r, z, uhat, tdhat, tshat) 

!~transformation of green's functions from cylindrical to cartesian coordinates 

call GFL_TRANS (t, uhat, tdhat, ugreen, tdgreen)

call GFL_TRANS (t, uhat, tshat, ugreen, tsgreen)



!open(unit=101, file='./usmp.txt', action='write')
!    loop2: Do iii=1, ms
!        write(unit=101, fmt='(E15.5,E15.5,E15.5,E15.5,E15.5)'), real(zeta(iii)), real(usmp_global(iii,1,1,0)), &
!                             aimag(usmp_global(iii,1,1,0)),real(usmp_global(iii,1,1,2)), aimag(usmp_global(iii,1,1,2))
!    end Do loop2
!close(unit=101)

!----------------------------------------------------------------------------------------------------------------------------------!
CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_INIT 
!----------------------------------------------------------------------------------------------------------------------------------!

!~initialization routine (input of the layer properties and integration parameters, dynamic memory allocation) 

USE GFLA 
USE GFLC 

integer(kind=k9) :: j      

!open (unit=kdgl, file='layers'//char(48+kdgl)//'.dat', action='read')   !~open layers1.dat, layers2.dat, etc...
open (unit=kdgl, file='layerscoord'//char(48+kdgl)//'.dat', action='read')   !~open layerscoord1.dat, layerscoord2.dat,etc..
!~data input: number of layers, maximium recursion depth, separation and height of the contour integration path  

read(kdgl,*) nl, nrec, dcont, hcont   

!~total number of layers including half-space (ml>1) 

ml = nl+1 

!~maximum number of sampling points along the integration path  

if (nrec<1) nrec = 1

ns = 8
loop1: do j = 1,nrec 
  ns = 2*ns - 1
end do loop1 
ms = nint*(ns-1) + 1 

!~dynamic memory allocation

call GFL_ALLOC

!~data input: material and geometric layer properties    
!~mu(j) = complex shear modulus [MPa]; v(j) = poisson's ratio; rho(j) = mass density [Gg/m^3]; h(j) = layer thickness [m] 
!read from layerscoord.dat files. Note the data is in terms of mu,poisson's ratio, density and the z coordinate of 
!bottom of the layer. The top layer z coordinate starts from zero
zzcoord(0) =0.0_dp
loop2: do j = 1,ml 
   !  read(kdgl,*) mu(j), v(j), rho(j), hh(j)
   read(kdgl,*) mu(j), v(j), rho(j), zzcoord(j)
   hh(j)=zzcoord(j)-zzcoord(j-1)   
end do loop2

!~frequency-independent material constants used for normalization 

mu0  = real(mu(sum(minloc(real(mu)/rho))))    
rho0 = rho(sum(minloc(real(mu)/rho))) 

!~normalization of material parameters  

mu  = mu/mu0
rho = rho/rho0  

!~auxiliary material parameters  

lam = 2.0_dp*mu*v/(1.0_dp-2.0_dp*v) 
ks  = sqrt(rho/mu)
kp  = sqrt(rho/(lam+2.0_dp*mu)) 
pv  = 2.0_dp-4.0_dp*v

!~initialization of remaining arrays 

cc = reshape((/0.0000000000000000e-0_dp, 0.0000000000000000e-0_dp, 0.0000000000000000e-0_dp, 1.0000000000000000e-0_dp, &
               0.0000000000000000e-0_dp, 0.0000000000000000e-0_dp, 0.0000000000000000e-0_dp, 0.0000000000000000e-0_dp, &
              -9.5238095238095238e-3_dp, 1.0000000000000000e-1_dp,-0.6000000000000000e-0_dp,-2.5000000000000000e-1_dp, &
               1.0000000000000000e-0_dp,-3.0000000000000000e-1_dp, 6.6666666666666667e-2_dp,-7.1428571428571429e-3_dp, &
               5.5555555555555556e-3_dp,-7.5000000000000000e-2_dp, 7.5000000000000000e-1_dp,-1.3611111111111111e-0_dp, &
               7.5000000000000000e-1_dp,-7.5000000000000000e-2_dp, 5.5555555555555556e-3_dp, 0.0000000000000000e-0_dp, &
               1.1111111111111111e-2_dp,-9.8611111111111111e-2_dp, 6.6666666666666667e-2_dp, 3.4027777777777778e-1_dp, &
              -6.1111111111111111e-1_dp, 3.7083333333333333e-1_dp,-8.8888888888888889e-2_dp, 9.7222222222222222e-3_dp, &
              -6.9444444444444444e-3_dp, 8.3333333333333333e-2_dp,-2.7083333333333333e-1_dp, 3.8888888888888889e-1_dp, &
              -2.7083333333333333e-1_dp, 8.3333333333333333e-2_dp,-6.9444444444444444e-3_dp, 0.0000000000000000e-0_dp, &
              -1.3888888888888889e-3_dp,-2.7777777777777778e-3_dp, 3.7500000000000000e-2_dp,-9.7222222222222222e-2_dp, &
               1.1805555555555556e-1_dp,-7.5000000000000000e-2_dp, 2.3611111111111111e-2_dp,-2.7777777777777778e-3_dp, &
               1.3888888888888889e-3_dp,-8.3333333333333333e-3_dp, 2.0833333333333333e-2_dp,-2.7777777777777778e-2_dp, &
               2.0833333333333333e-2_dp,-8.3333333333333333e-3_dp, 1.3888888888888889e-3_dp, 0.0000000000000000e-0_dp, &
              -1.9841269841269841e-4_dp, 1.3888888888888889e-3_dp,-4.1666666666666667e-3_dp, 6.9444444444444444e-3_dp, &
              -6.9444444444444444e-3_dp, 4.1666666666666667e-3_dp,-1.3888888888888889e-3_dp, 1.9841269841269841e-4_dp/), (/8,8/)) 

!~check if the half-space is homogeneous (yes- homog=1; no- homog=0; eps=1.0e-9)  

homog = 1

loop3: do j = 2,ml 
  if (abs(mu(j)-mu(j-1))>eps*abs(mu(j)).or.abs(v(j)-v(j-1))>eps*abs(v(j)).or.abs(rho(j)-rho(j-1))>eps*rho(j)) then
    homog=0;  exit 
  end if  
end do loop3 

close(kdgl)   

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_INIT 
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_ALLOC   
!----------------------------------------------------------------------------------------------------------------------------------!

!~dynamic memory allocation

USE GFLA 

!~arrays containing material properties, thicknesses and p/s wave numbers of the layers 

allocate (mu(ml), v(ml), pv(ml), lam(ml), rho(ml), hh(ml), hb(ml), zh(0:ml), kp(ml), ks(ml), zzcoord(0:ml) ) 

!~arrays containing coordinates/flags of the sampling points, alfas, betas and coefficient matrices   

allocate (zeta(ms), zetas(ms), flag(ns), alfa(ms,ml), beta(ms,ml), cmps(ms,ml,6,4), cmsh(ms,ml,3,2))

!~arrays containing reference matrices 

allocate (tups0(ms,nl,2,2), tdps0(ms,nl,2,2), rups0(ms,0:nl,2,2), rdps0(ms,nl,2,2))  
allocate (tush0(ms,nl)    , tdsh0(ms,nl)    , rush0(ms,0:nl)    , rdsh0(ms,nl)) 

!~arrays containing transmission/reflection matrices 

allocate (tups(ms,nl,2,2), tdps(ms,nl,2,2), rups(ms,0:nl,2,2), rdps(ms,nl,2,2))  
allocate (tush(ms,nl)    , tdsh(ms,nl)    , rush(ms,0:nl)    , rdsh(ms,nl)) 

!~arrays containing generalized propagation matrices 

allocate (thups(ms,nl,2,2), thdps(ms,nl,2,2), rhups(ms,0:nl,2,2), rhdps(ms,ml,2,2))  
allocate (thush(ms,nl)    , thdsh(ms,nl)    , rhush(ms,0:nl)    , rhdsh(ms,ml)) 

!~arrays containing factorized propagation matrices 

allocate (wups(ms,ml,2,2), wdps(ms,ml,2,2), wush(ms,ml), wdsh(ms,ml))  

!~arrays containing the normalized loading coefficients due to the point source   

allocate (sups(ms,3,-1:1,2), sdps(ms,3,-1:1,2), sush(ms,3,-1:1), sdsh(ms,3,-1:1))  

!~arrays containing the auxiliary loading coefficients due to the point source   

allocate (mups(ms,3,-1:1,2), mdps(ms,3,-1:1,2), mush(ms,3,-1:1), mdsh(ms,3,-1:1))  

!~arrays containing the integrand evaluated at the sampling points    

allocate (usmp(ns,3,3,0:2), tsmp(ns,3,8,0:2))  

allocate (usmp_global(ms,3,3,0:2))
!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_ALLOC   
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_DEALLOC   
!----------------------------------------------------------------------------------------------------------------------------------!

!~dynamic memory deallocation before handling domains 2 through ndom

USE GFLA 

!~arrays containing material properties, thicknesses and p/s wave numbers of the layers 

deallocate (mu, v, pv, lam, rho, hh, hb, zh, kp, ks, zzcoord) 

!~arrays containing coordinates/flags of the sampling points, alfas, betas and coefficient matrices   

deallocate (zeta, zetas, flag, alfa, beta, cmps, cmsh)

!~arrays containing reference matrices 

deallocate (tups0, tdps0, rups0, rdps0)  
deallocate (tush0, tdsh0, rush0, rdsh0) 

!~arrays containing transmission/reflection matrices 

deallocate (tups, tdps, rups, rdps)  
deallocate (tush, tdsh, rush, rdsh) 

!~arrays containing generalized propagation matrices 

deallocate (thups, thdps, rhups, rhdps)  
deallocate (thush, thdsh, rhush, rhdsh) 

!~arrays containing factorized propagation matrices 

deallocate (wups, wdps, wush, wdsh)  

!~arrays containing the normalized loading coefficients due to the point source   

deallocate (sups, sdps, sush, sdsh)  

!~arrays containing the auxiliary loading coefficients due to the point source   

deallocate (mups, mdps, mush, mdsh)  

!~arrays containing the integrand evaluated at the sampling points    

deallocate (usmp, tsmp)  

deallocate (usmp_global)
!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_DEALLOC   
!----------------------------------------------------------------------------------------------------------------------------------!


!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_FDPAR (omega, radius)
!----------------------------------------------------------------------------------------------------------------------------------!

!~evaluation of the frequency-dependent parameters

USE GFLA
USE GFLC

real   (kind=dp), intent(in) :: omega, radius
integer(kind=k9)             :: j        

!~reference wave number used for normalization

ks0 = omega*sqrt(rho0/mu0)  

!~normalized layer thicknesses 

hb = ks0 * hh  

!~limit on the minimum dimensionless layer thickness 

if (minval(hb)<hbmin) then 
  print *, 'Minimum dimensionless layer thickness too small' 
  print *, 'Increase the layer thicknesses or the frequency';  stop
end if  

!~dimensionless coordinates of the layer interfaces 

zh(0) = 0.0_dp
zh(1:ml) = ks0*zzcoord(1:ml)

!loop1: do j = 1,ml
!   zh(j) = zh(j-1) + hb(j)  !--this still introduces error in summing layer thicknesses to obtain layer depth
!end do loop1 


!~adjust the contour taking int account the omega and radius
if (radius*ks0>20.0_dp*(hcont_max/hcont)) then
    hcont=19.9_dp*hcont_max/(radius*ks0)
    dcont=hcont
end if

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_FDPAR     
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_CONTO   
!----------------------------------------------------------------------------------------------------------------------------------!

!~evaluation of the integration variable zeta and radicals alfa, beta along the contour path 

USE GFLA 
USE GFLC 

integer(kind=k9) :: j, k        
complex(kind=dp) :: contour(nint+1), zmax, dz       

!~adaptive length of the integration path (zmax) as measured along the real axis    

if (homog==0) then  
  zmax  = max(lpmin, lpmax*hbmin/minval(hb))  !-heterogeneous half-space (lpmin < zmax < lpmax)
else  
  zmax  = lpmin  !-homogeneous half-space (zmax = lpmin = 100)
end if 

lpath = real(zmax) 

!~check the height the contour path  (input parameter hcont)  

if (hcont>hcont_max) then 
  hcont = hcont_max;  print 100, 'Height of the contour: hcont>hcont_max ==> set hcont=', hcont_max  
else if (hcont<hcont_min) then 
  hcont = hcont_min;  print 100, 'Height of the contour: hcont<hcont_min ==> set hcont=', hcont_min  
end if   

!~check the initial slope of the contour path  (input parameter dcont)  
       
if (dcont>hcont) then 
  dcont = hcont;  print 100, 'Initial contour slope: theta_init<Pi/4 ==> set dcont=', hcont    
elseif (dcont<dcont_min) then 
  dcont = dcont_min;  print 100, 'Initial contour slope: theta_init~Pi/2 ==> set dcont=', dcont_min    
end if   

!~nodes of the contour path 

contour = (/zero, 0.50_dp*cmplx(dcont,hcont,dp), cmplx(dcont,hcont,dp),                                & 
            cmplx(dcont+0.33_dp*(lcont-dcont),hcont,dp), cmplx(dcont+0.66_dp*(lcont-dcont),hcont,dp),  & 
            cmplx(lcont,hcont,dp), cmplx(2.00_dp*lcont,0.00_dp,dp), 0.10_dp*zmax, 0.40_dp*zmax, zmax/) 

!~coordinates of the sampling points (integration variable zeta)  

loop1: do j = 1,nint 
  dz = (contour(j+1)-contour(j))/(ns-1)
  loop2: do k = 1,ns
    zeta((j-1)*(ns-1)+k) = contour(j) + (k-1)*dz 
  end do loop2
end do loop1 

!~avoid evaluation of the integrand exactly at zeta=0 

zeta(1) = (0.0_dp,1.0e-6_dp)  

!~square of the integration variable 

zetas = zeta*zeta

!~evaluation of radicals alfa, beta at the sampling points    

loop3: do j = 1,ml  
  alfa(:,j) = sqrt(zetas-kp(j)*kp(j))
  beta(:,j) = sqrt(zetas-ks(j)*ks(j)) 
end do loop3 

!~transiton point between the standard and modified representation of p-sv waves 

ka = min(sum(minloc(abs(zeta-2.0_dp*minval(hb)))),ms-1)
kb = ka + 1       

100 format (a, f4.2)

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_CONTO   
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_AUXAR      
!----------------------------------------------------------------------------------------------------------------------------------!

!~regular/asymptotic evaluation of the auxiliary arrays: difa, difb, diab, dfaa, dfbb, dfab, difc, difd  
!~asymptotic expressions are used for abs(zeta)>320*abs(ks(j))        

USE GFLA 

integer(kind=k9)              :: j, k, kc, kd, iflag=0
real   (kind=dp)              :: zt        
complex(kind=dp)              :: a, b      
complex(kind=dp), allocatable :: x(:)     

!~dynamic memory allocation 

if (iflag==0) then  
  iflag = 1
else 
  deallocate (difa, difb, diab, dfaa, dfab, dfbb, difc, difd) 
end if  

allocate (difa(kb:ms,ml), difb(kb:ms,ml), diab(kb:ms,ml), dfaa(kb:ms,ml))
allocate (dfab(kb:ms,ml), dfbb(kb:ms,ml), difc(kb:ms,ml), difd(kb:ms,ml)) 

loop1: do j = 1,ml  !-loop over the (layers + half-space) 

!~~~transition point between standard and asymptotic expressions  

  zt = 320.0_dp*abs(ks(j)) 
  kc = ms 

  loop2: do k = kb,ms 
    if (abs(zeta(k))>zt) then 
      kc = k 
      exit loop2
    end if  
  end do loop2 

  kd = kc+1 

!~~~standard evaluation of the auxiliary arrays    

  a = kp(j)*kp(j) 
  b = ks(j)*ks(j) 

  difa(kb:kc,j) = zeta(kb:kc) - alfa(kb:kc,j) 
  difb(kb:kc,j) = zeta(kb:kc) - beta(kb:kc,j) 
  diab(kb:kc,j) = alfa(kb:kc,j) - beta(kb:kc,j) 

  dfaa(kb:kc,j) = a  
  dfbb(kb:kc,j) = b   

  dfab(kb:kc,j) = 1.0_dp / (zetas(kb:kc) - alfa(kb:kc,j)*beta(kb:kc,j)) 

  difc(kb:kc,j) = zeta(kb:kc) * (2.0_dp+pv(j) - (1.0_dp+pv(j))*dfbb(kb:kc,j)*dfab(kb:kc,j))
  difd(kb:kc,j) = zeta(kb:kc) * (1.0_dp - (1.0_dp+pv(j))/pv(j)*dfaa(kb:kc,j)*dfab(kb:kc,j)) 

!~~~asymptotic evaluation of the auxiliary arrays  

  if (kc<ms) then 

    allocate (x(kd:ms))
    x = 1.0_dp/(zetas(kd:ms))

    difa(kd:ms,j) = 0.5_dp/zeta(kd:ms)*(a+x*(0.25_dp*a*a*(1.0_dp+x*0.5_dp*a))) 
    difb(kd:ms,j) = 0.5_dp/zeta(kd:ms)*(b+x*(0.25_dp*b*b*(1.0_dp+x*0.5_dp*b))) 
    diab(kd:ms,j) = 0.5_dp/zeta(kd:ms)*(b-a+x*(0.25_dp*(b*b-a*a)+x*0.125_dp*(b**3-a**3))) 

    dfaa(kd:ms,j) = a  
    dfbb(kd:ms,j) = b  
    dfab(kd:ms,j) = 2.0_dp/(a+b+x*(0.25_dp*(a-b)**2+x*(0.125_dp*(a-b)*(a*a-b*b)-x*0.03125_dp*a*b*(2.0_dp*(a*a+b*b)+a*b))))
 
    difc(kd:ms,j) = zeta(kd:ms) * (2.0_dp+pv(j) - (1.0_dp+pv(j))*dfbb(kd:ms,j)*dfab(kd:ms,j)) 
    difd(kd:ms,j) = zeta(kd:ms) * (1.0_dp - (1.0_dp+pv(j))/pv(j)*dfaa(kd:ms,j)*dfab(kd:ms,j)) 

    deallocate (x) 

  end if   

end do loop1 

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_AUXAR      
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_CFMAT      
!----------------------------------------------------------------------------------------------------------------------------------!

!~evaluation of the coefficient matrices: cmps, cmsh 

USE GFLA 
   
integer(kind=k9) :: j  
complex(kind=dp) :: pvm2, pvm4      

loop1: do j = 1,ml  !-loop over the (layers + half-space) 

!~~~jth layer, p-sv waves, standard representation for moderate values of zeta

  cmps(1:ka,j,1,1) = -zeta(1:ka)
  cmps(1:ka,j,1,2) =  beta(1:ka,j) 
  cmps(1:ka,j,1,3) =  cmps(1:ka,j,1,1)
  cmps(1:ka,j,1,4) = -cmps(1:ka,j,1,2) 

  cmps(1:ka,j,2,1) = -alfa(1:ka,j)
  cmps(1:ka,j,2,2) =  zeta(1:ka)
  cmps(1:ka,j,2,3) = -cmps(1:ka,j,2,1)
  cmps(1:ka,j,2,4) =  cmps(1:ka,j,2,2) 

  cmps(1:ka,j,3,1) =  2.0_dp*mu(j)*zeta(1:ka)*alfa(1:ka,j)
  cmps(1:ka,j,3,2) = -mu(j)*(beta(1:ka,j)*beta(1:ka,j)+zetas(1:ka)) 
  cmps(1:ka,j,3,3) = -cmps(1:ka,j,3,1)
  cmps(1:ka,j,3,4) =  cmps(1:ka,j,3,2) 

  cmps(1:ka,j,4,1) =  mu(j)*(beta(1:ka,j)*beta(1:ka,j)+zetas(1:ka))
  cmps(1:ka,j,4,2) = -2.0_dp*mu(j)*zeta(1:ka)*beta(1:ka,j) 
  cmps(1:ka,j,4,3) =  cmps(1:ka,j,4,1)
  cmps(1:ka,j,4,4) = -cmps(1:ka,j,4,2) 

  cmps(1:ka,j,5,1) =  lam(j)*(alfa(1:ka,j)*alfa(1:ka,j)-zetas(1:ka))  
  cmps(1:ka,j,5,2) =  0.0_dp  
  cmps(1:ka,j,5,3) =  cmps(1:ka,j,5,1)
  cmps(1:ka,j,5,4) =  cmps(1:ka,j,5,2) 

  cmps(1:ka,j,6,1) =  lam(j)*alfa(1:ka,j)*alfa(1:ka,j) - (lam(j)+2.0_dp*mu(j))*zetas(1:ka) 
  cmps(1:ka,j,6,2) =  2.0_dp*mu(j)*zeta(1:ka)*beta(1:ka,j)  
  cmps(1:ka,j,6,3) =  cmps(1:ka,j,6,1)
  cmps(1:ka,j,6,4) = -cmps(1:ka,j,6,2) 

!~~~jth layer, p-sv waves, modified representation for large values of zeta   

  pvm2 = pv(j)-2.0_dp 
  pvm4 = pv(j)-4.0_dp

  cmps(kb:ms,j,1,1) =  1.0_dp
  cmps(kb:ms,j,1,2) = -1.0_dp 
  cmps(kb:ms,j,1,3) =  1.0_dp
  cmps(kb:ms,j,1,4) =  1.0_dp 

  cmps(kb:ms,j,2,1) = -pv(j)
  cmps(kb:ms,j,2,2) = -1.0_dp
  cmps(kb:ms,j,2,3) =  pv(j)
  cmps(kb:ms,j,2,4) = -1.0_dp 

  cmps(kb:ms,j,3,1) =  mu(j)*(pvm2*zeta(kb:ms) + difc(:,j) + difa(:,j)*dfbb(:,j)*dfab(:,j))
  cmps(kb:ms,j,3,2) =  mu(j)*(2.0_dp*zeta(kb:ms) - difa(:,j)*dfbb(:,j)*dfab(:,j)) 
  cmps(kb:ms,j,3,3) = -cmps(kb:ms,j,3,1)
  cmps(kb:ms,j,3,4) =  cmps(kb:ms,j,3,2)  

  cmps(kb:ms,j,4,1) =  mu(j)*(pv(j)*zeta(kb:ms) - difc(:,j) - pv(j)*(difb(:,j)*dfbb(:,j)*dfab(:,j)))
  cmps(kb:ms,j,4,2) =  mu(j)*(2.0_dp*zeta(kb:ms) - difb(:,j)*dfbb(:,j)*dfab(:,j))
  cmps(kb:ms,j,4,3) =  cmps(kb:ms,j,4,1) 
  cmps(kb:ms,j,4,4) = -cmps(kb:ms,j,4,2)  

  cmps(kb:ms,j,5,1) = -mu(j)*(pvm2*zeta(kb:ms) - pvm2*(difd(:,j)+difb(:,j)*dfaa(:,j)*dfab(:,j)))    
  cmps(kb:ms,j,5,2) =  mu(j)*(pvm2/pv(j)*difb(:,j)*dfaa(:,j)*dfab(:,j))  
  cmps(kb:ms,j,5,3) =  cmps(kb:ms,j,5,1) 
  cmps(kb:ms,j,5,4) = -cmps(kb:ms,j,5,2)  

  cmps(kb:ms,j,6,1) = -mu(j)*(pvm4*zeta(kb:ms) - pvm2*(difd(:,j)+difb(:,j)*dfaa(:,j)*dfab(:,j)))   
  cmps(kb:ms,j,6,2) = -mu(j)*(2.0_dp*zeta(kb:ms) - pvm2/pv(j)*difb(:,j)*dfaa(:,j)*dfab(:,j))   
  cmps(kb:ms,j,6,3) =  cmps(kb:ms,j,6,1)  
  cmps(kb:ms,j,6,4) = -cmps(kb:ms,j,6,2) 

!~~~jth layer, sh waves, universal representation 

  cmsh(:,j,1,1) =  ic*zeta 
  cmsh(:,j,1,2) =  cmsh(:,j,1,1)

  cmsh(:,j,2,1) = -ic*mu(j)*zeta*beta(:,j)
  cmsh(:,j,2,2) = -cmsh(:,j,2,1)

  cmsh(:,j,3,1) =  mu(j)*zetas
  cmsh(:,j,3,2) =  cmsh(:,j,3,1)

end do loop1 

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_CFMAT      
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_RFMAT      
!----------------------------------------------------------------------------------------------------------------------------------!

!~evaluation of the reference transmission/reflection matrices/coefficients: tups0, tdps0, rups0, rdps0, tush0, tdsh0, rush0, rdsh0

USE GFLA 

integer(kind=k9) :: j       
real   (kind=dp) :: dr     
complex(kind=dp) :: dm, shdt(ms), sdet(1:ka), zab0(1:ka), zab1(1:ka), dela(1:ka), delb(1:ka)        
complex(kind=dp) :: rhap(1:ka), rhbp(1:ka), rham(1:ka), rhbm(1:ka), sdat(kb:ms), v0, v1, opv0, omv0, opv1, omv1  

!~free surface, p-sv waves, standard representation for moderate values of zeta

zab1 = zetas(1:ka) - alfa(1:ka,1)*beta(1:ka,1) 

sdet = -1.0_dp/(4.0_dp*zetas(1:ka)*zab1-ks(1)*ks(1)*(4.0_dp*zetas(1:ka)-ks(1)*ks(1)))  
   
!-reference reflection matrix  Rou_0  

rups0(1:ka,0,1,1) = sdet * ((beta(1:ka,1)*beta(1:ka,1)+zetas(1:ka))**2+4.0_dp*zetas(1:ka)*alfa(1:ka,1)*beta(1:ka,1)) 
rups0(1:ka,0,1,2) = sdet * 4.0_dp*zeta(1:ka)*beta(1:ka,1)*(beta(1:ka,1)*beta(1:ka,1) + zetas(1:ka)) 
rups0(1:ka,0,2,1) = sdet * 4.0_dp*zeta(1:ka)*alfa(1:ka,1)*(beta(1:ka,1)*beta(1:ka,1) + zetas(1:ka))
rups0(1:ka,0,2,2) = sdet * ((beta(1:ka,1)*beta(1:ka,1)+zetas(1:ka))**2+4.0_dp*zetas(1:ka)*alfa(1:ka,1)*beta(1:ka,1)) 

!~free surface, p-sv waves, modified representation for large values of zeta

sdat = 1.0_dp/(cmps(kb:ms,1,4,2)*cmps(kb:ms,1,3,1)-cmps(kb:ms,1,4,1)*cmps(kb:ms,1,3,2))

!-reference reflection matrix  Rou_0  

rups0(kb:ms,0,1,1) = sdat * (cmps(kb:ms,1,4,1)*cmps(kb:ms,1,3,2) + cmps(kb:ms,1,4,2)*cmps(kb:ms,1,3,1))
rups0(kb:ms,0,1,2) = sdat * (-2.0_dp*cmps(kb:ms,1,4,2)*cmps(kb:ms,1,3,2)) 
rups0(kb:ms,0,2,1) = sdat * (-2.0_dp*cmps(kb:ms,1,4,1)*cmps(kb:ms,1,3,1))  
rups0(kb:ms,0,2,2) = sdat * (cmps(kb:ms,1,4,1)*cmps(kb:ms,1,3,2) + cmps(kb:ms,1,4,2)*cmps(kb:ms,1,3,1))   

!~free surface, sh waves, universal representation 

rush0(:,0) = 1.0_dp   

loop1: do j = 1,nl  !-loop over the layer interfaces 

!~~~jth interface, p-sv waves, standard representation for moderate values of zeta

  zab0 = zetas(1:ka) - alfa(1:ka,j)*beta(1:ka,j) 
  zab1 = zetas(1:ka) - alfa(1:ka,j+1)*beta(1:ka,j+1) 

  rhap = rho(j)*alfa(1:ka,j+1) + rho(j+1)*alfa(1:ka,j) 
  rhbp = rho(j)*beta(1:ka,j+1) + rho(j+1)*beta(1:ka,j) 
  rham = rho(j)*alfa(1:ka,j+1) - rho(j+1)*alfa(1:ka,j) 
  rhbm = rho(j)*beta(1:ka,j+1) - rho(j+1)*beta(1:ka,j) 

  dm   =  mu(j+1)-mu(j)
  dr   =  0.5_dp*(rho(j+1)-rho(j)) 

  sdet = -1.0_dp/(4.0_dp*zetas(1:ka)*(dm*(dm*zab0*zab1+rho(j)*zab1-rho(j+1)*zab0)+dr**2)-rhap*rhbp)  

!---reference transmission matrix Tju_0 

  zab0 = zetas(1:ka) - alfa(1:ka,j)*beta(1:ka,j+1) 
  zab1 = zetas(1:ka) - alfa(1:ka,j+1)*beta(1:ka,j) 

  dela = alfa(1:ka,j+1) - alfa(1:ka,j) 
  delb = beta(1:ka,j+1) - beta(1:ka,j) 
   
  tups0(1:ka,j,1,1) = sdet * 2.0_dp*rho(j+1)*alfa(1:ka,j+1)*(2.0_dp*zetas(1:ka)*dm*delb+rhbp) 
  tups0(1:ka,j,1,2) = sdet * 4.0_dp*rho(j+1)*beta(1:ka,j+1)*(zeta(1:ka)*(dm*zab1-dr)) 
  tups0(1:ka,j,2,1) = sdet * 4.0_dp*rho(j+1)*alfa(1:ka,j+1)*(zeta(1:ka)*(dm*zab0-dr))
  tups0(1:ka,j,2,2) = sdet * 2.0_dp*rho(j+1)*beta(1:ka,j+1)*(2.0_dp*zetas(1:ka)*dm*dela+rhap) 

!---reference transmission matrix Tjd_0  

  tdps0(1:ka,j,1,1) = sdet * 2.0_dp*rho(j)*alfa(1:ka,j)*(2.0_dp*zetas(1:ka)*dm*delb+rhbp)  
  tdps0(1:ka,j,1,2) = sdet * 4.0_dp*rho(j)*beta(1:ka,j)*(zeta(1:ka)*(dm*zab0-dr)) 
  tdps0(1:ka,j,2,1) = sdet * 4.0_dp*rho(j)*alfa(1:ka,j)*(zeta(1:ka)*(dm*zab1-dr))
  tdps0(1:ka,j,2,2) = sdet * 2.0_dp*rho(j)*beta(1:ka,j)*(2.0_dp*zetas(1:ka)*dm*dela+rhap) 

!---reference reflection matrix Rju_0 

  zab0 = zetas(1:ka) - alfa(1:ka,j)*beta(1:ka,j) 
  zab1 = zetas(1:ka) + alfa(1:ka,j+1)*beta(1:ka,j+1) 

  rups0(1:ka,j,1,1) = sdet * (4.0_dp*zetas(1:ka)*(dm*(dm*zab0*zab1+rho(j)*zab1-rho(j+1)*zab0)+dr**2)+rham*rhbp)
  rups0(1:ka,j,1,2) = sdet * 4.0_dp*zeta(1:ka)*beta(1:ka,j+1)*(2.0_dp*zetas(1:ka)*dm*(dm*zab0+rho(j))-dm*rho(j+1)*zab0-rho(j)*dr)
  rups0(1:ka,j,2,1) = sdet * 4.0_dp*zeta(1:ka)*alfa(1:ka,j+1)*(2.0_dp*zetas(1:ka)*dm*(dm*zab0+rho(j))-dm*rho(j+1)*zab0-rho(j)*dr)
  rups0(1:ka,j,2,2) = sdet * (4.0_dp*zetas(1:ka)*(dm*(dm*zab0*zab1+rho(j)*zab1-rho(j+1)*zab0)+dr**2)+rhap*rhbm)
  
!---reference reflection matrix Rjd_0  

  zab0 = zetas(1:ka) + alfa(1:ka,j)*beta(1:ka,j) 
  zab1 = zetas(1:ka) - alfa(1:ka,j+1)*beta(1:ka,j+1) 

  rdps0(1:ka,j,1,1) = sdet * (4.0_dp*zetas(1:ka)*(dm*(dm*zab0*zab1+rho(j)*zab1-rho(j+1)*zab0)+dr**2)-rham*rhbp)
  rdps0(1:ka,j,1,2) = sdet * 4.0_dp*zeta(1:ka)*beta(1:ka,j)*(2.0_dp*zetas(1:ka)*dm*(rho(j+1)-dm*zab1)-dm*rho(j)*zab1-rho(j+1)*dr)
  rdps0(1:ka,j,2,1) = sdet * 4.0_dp*zeta(1:ka)*alfa(1:ka,j)*(2.0_dp*zetas(1:ka)*dm*(rho(j+1)-dm*zab1)-dm*rho(j)*zab1-rho(j+1)*dr)
  rdps0(1:ka,j,2,2) = sdet * (4.0_dp*zetas(1:ka)*(dm*(dm*zab0*zab1+rho(j)*zab1-rho(j+1)*zab0)+dr**2)-rhap*rhbm)
 
!~~~jth interface, p-sv waves, modified representation for large values of zeta     

  v0   = pv(j)
  v1   = pv(j+1)
  opv0 = 1.0_dp + v0
  omv0 = 1.0_dp - v0
  opv1 = 1.0_dp + v1
  omv1 = 1.0_dp - v1  

  sdat = 1.0_dp/(cmps(kb:ms,j+1,4,1)*(opv0*cmps(kb:ms,j+1,3,2)-omv0*cmps(kb:ms,j,3,2)-2*cmps(kb:ms,j,3,1)) +                     &
                 cmps(kb:ms,j+1,4,2)*((v1+v0)*cmps(kb:ms,j,3,2)-opv0*cmps(kb:ms,j+1,3,1)-omv1*cmps(kb:ms,j,3,1)) +               & 
                 cmps(kb:ms,j,4,1)*(opv1*cmps(kb:ms,j,3,2)-omv1*cmps(kb:ms,j+1,3,2)-2*cmps(kb:ms,j+1,3,1)) +                     &
                 cmps(kb:ms,j,4,2)*((v1+v0)*cmps(kb:ms,j+1,3,2)-opv1*cmps(kb:ms,j,3,1)-omv0*cmps(kb:ms,j+1,3,1)))  
 
!---reference transmission matrix Tju_0 

  tups0(kb:ms,j,1,1) = sdat * (cmps(kb:ms,j+1,4,1)*(v1*(cmps(kb:ms,j+1,3,2)+cmps(kb:ms,j,3,2))-2.0_dp*cmps(kb:ms,j+1,3,1)) +     &
                               cmps(kb:ms,j+1,4,2)*(v1*cmps(kb:ms,j,3,2)-cmps(kb:ms,j+1,3,1)) -                                  & 
                               cmps(kb:ms,j,4,2)*(cmps(kb:ms,j+1,3,1)-v1*cmps(kb:ms,j+1,3,2)))*2.0_dp   
  tups0(kb:ms,j,1,2) = sdat * (cmps(kb:ms,j+1,4,1)*(cmps(kb:ms,j+1,3,2)-cmps(kb:ms,j,3,2)) +                                     &
                               cmps(kb:ms,j+1,4,2)*(cmps(kb:ms,j+1,3,1)+omv1*cmps(kb:ms,j+1,3,2)-cmps(kb:ms,j,3,2)) +            &
                               cmps(kb:ms,j,4,2)*(v1*cmps(kb:ms,j+1,3,2)-cmps(kb:ms,j+1,3,1)))*2.0_dp 
  tups0(kb:ms,j,2,1) = sdat * (cmps(kb:ms,j+1,4,1)*(v1*(cmps(kb:ms,j,3,1)-cmps(kb:ms,j+1,3,2))+omv0*cmps(kb:ms,j+1,3,1)) +       &
                               cmps(kb:ms,j+1,4,2)*(v1*cmps(kb:ms,j,3,1)-v0*cmps(kb:ms,j+1,3,1)) +                               & 
                               cmps(kb:ms,j,4,1)*(v1*cmps(kb:ms,j+1,3,2)-cmps(kb:ms,j+1,3,1)))*2.0_dp   
  tups0(kb:ms,j,2,2) = sdat * (cmps(kb:ms,j+1,4,1)*(v0*cmps(kb:ms,j+1,3,2)-cmps(kb:ms,j,3,1)) +                                  &
                               cmps(kb:ms,j+1,4,2)*((v0+v1)*cmps(kb:ms,j+1,3,2)-cmps(kb:ms,j+1,3,1)-cmps(kb:ms,j,3,1)) -         & 
                               cmps(kb:ms,j,4,1)*(cmps(kb:ms,j+1,3,1)-v1*cmps(kb:ms,j+1,3,2)))*2.0_dp    

!---reference transmission matrix Tjd_0 

  tdps0(kb:ms,j,1,1) = sdat * (cmps(kb:ms,j+1,4,2)*(v0*cmps(kb:ms,j,3,2)-cmps(kb:ms,j,3,1)) +                                    &
                               cmps(kb:ms,j,4,1)*(v0*(cmps(kb:ms,j,3,2)+cmps(kb:ms,j+1,3,2))-2.0_dp*cmps(kb:ms,j,3,1)) -         &
                               cmps(kb:ms,j,4,2)*(cmps(kb:ms,j,3,1)-v0*cmps(kb:ms,j+1,3,2)))*2.0_dp        
  tdps0(kb:ms,j,1,2) = sdat * 2.0_dp*(cmps(kb:ms,j+1,4,2)*(cmps(kb:ms,j,3,1)-v0*cmps(kb:ms,j,3,2)) +                             &
                               cmps(kb:ms,j,4,1)*(cmps(kb:ms,j+1,3,2)-cmps(kb:ms,j,3,2)) +                                       &
                               cmps(kb:ms,j,4,2)*(cmps(kb:ms,j+1,3,2)-omv0*cmps(kb:ms,j,3,2)-cmps(kb:ms,j,3,1)))      
  tdps0(kb:ms,j,2,1) = sdat * 2.0_dp*(cmps(kb:ms,j+1,4,1)*(cmps(kb:ms,j,3,1)-v0*cmps(kb:ms,j,3,2)) +                             &
                               cmps(kb:ms,j,4,2)*(v1*cmps(kb:ms,j,3,1)-v0*cmps(kb:ms,j+1,3,1)) +                                 &
                               cmps(kb:ms,j,4,1)*(v0*(cmps(kb:ms,j,3,2)-cmps(kb:ms,j+1,3,1))-omv1*cmps(kb:ms,j,3,1)))        
  tdps0(kb:ms,j,2,2) = sdat * (cmps(kb:ms,j+1,4,1)*(cmps(kb:ms,j,3,2)*v0-cmps(kb:ms,j,3,1)) +                                    &
                               cmps(kb:ms,j,4,2)*(cmps(kb:ms,j,3,2)*(v1+v0)-cmps(kb:ms,j+1,3,1)-cmps(kb:ms,j,3,1)) +             &
                               cmps(kb:ms,j,4,1)*(v1*cmps(kb:ms,j,3,2)-cmps(kb:ms,j+1,3,1)))*2.0_dp
     
!---reference reflection matrix Rju_0  

  rups0(kb:ms,j,1,1) = sdat * (cmps(kb:ms,j+1,4,1)*(omv0*cmps(kb:ms,j,3,2)-opv0*cmps(kb:ms,j+1,3,2)+2.0_dp*cmps(kb:ms,j,3,1)) +  &
                               cmps(kb:ms,j+1,4,2)*(opv1*cmps(kb:ms,j,3,1)-opv0*cmps(kb:ms,j+1,3,1)+(v1-v0)*cmps(kb:ms,j,3,2)) + &
                               cmps(kb:ms,j,4,1)*(opv1*cmps(kb:ms,j+1,3,2)-omv1*cmps(kb:ms,j,3,2)-2.0_dp*cmps(kb:ms,j+1,3,1)) +  &      
                               cmps(kb:ms,j,4,2)*(omv1*cmps(kb:ms,j,3,1)-omv0*cmps(kb:ms,j+1,3,1)+(v1-v0)*cmps(kb:ms,j+1,3,2))) 
  rups0(kb:ms,j,1,2) = sdat * (cmps(kb:ms,j+1,4,2)*(opv0*cmps(kb:ms,j+1,3,2)-(cmps(kb:ms,j,3,1)+cmps(kb:ms,j,3,2))) +            &
                               cmps(kb:ms,j,4,1)*(cmps(kb:ms,j+1,3,2)-cmps(kb:ms,j,3,2)) +                                       &
                               cmps(kb:ms,j,4,2)*(cmps(kb:ms,j,3,1)-v0*cmps(kb:ms,j+1,3,2)))*2.0_dp     
  rups0(kb:ms,j,2,1) = sdat * (cmps(kb:ms,j+1,4,1)*(opv0*cmps(kb:ms,j+1,3,1)-v1*(cmps(kb:ms,j,3,1)+cmps(kb:ms,j,3,2))) -         &
                               cmps(kb:ms,j,4,1)*(cmps(kb:ms,j+1,3,1)-v1*cmps(kb:ms,j,3,2)) +                                    &
                               cmps(kb:ms,j,4,2)*(v0*cmps(kb:ms,j+1,3,1)-v1*cmps(kb:ms,j,3,1)))*2.0_dp      
  rups0(kb:ms,j,2,2) = sdat * (cmps(kb:ms,j+1,4,1)*(opv0*cmps(kb:ms,j,3,2)-opv0*cmps(kb:ms,j+1,3,2)) -                           &
                               cmps(kb:ms,j+1,4,2)*(omv1*cmps(kb:ms,j,3,1)+opv0*cmps(kb:ms,j+1,3,1)-(v1+v0)*cmps(kb:ms,j,3,2)) + & 
                               cmps(kb:ms,j,4,1)*(omv1*cmps(kb:ms,j+1,3,2)-omv1*cmps(kb:ms,j,3,2)) +                             &
                               cmps(kb:ms,j,4,2)*(omv1*cmps(kb:ms,j,3,1)+opv0*cmps(kb:ms,j+1,3,1)-(v1+v0)*cmps(kb:ms,j+1,3,2)))  
     
!---reference reflection matrix Rjd_0  

  rdps0(kb:ms,j,1,1) = sdat * (cmps(kb:ms,j+1,4,1)*(opv0*cmps(kb:ms,j,3,2)-omv0*cmps(kb:ms,j+1,3,2)-2.0_dp*cmps(kb:ms,j,3,1)) +  &   
                               cmps(kb:ms,j+1,4,2)*(omv0*cmps(kb:ms,j+1,3,1)-omv1*cmps(kb:ms,j,3,1)+(v0-v1)*cmps(kb:ms,j,3,2)) - &
                               cmps(kb:ms,j,4,1)*(opv1*cmps(kb:ms,j,3,2)-omv1*cmps(kb:ms,j+1,3,2)-2.0_dp*cmps(kb:ms,j+1,3,1)) +  &       
                               cmps(kb:ms,j,4,2)*(opv0*cmps(kb:ms,j+1,3,1)-opv1*cmps(kb:ms,j,3,1)+(v0-v1)*cmps(kb:ms,j+1,3,2)))  
  rdps0(kb:ms,j,1,2) = sdat *  2.0_dp*(cmps(kb:ms,j+1,4,1)*(cmps(kb:ms,j+1,3,2)-cmps(kb:ms,j,3,2)) +                             &
                               cmps(kb:ms,j+1,4,2)*(v1*cmps(kb:ms,j,3,2)-cmps(kb:ms,j+1,3,1)) +                                  &
                               cmps(kb:ms,j,4,2)*(cmps(kb:ms,j+1,3,1)+cmps(kb:ms,j+1,3,2)-opv1*cmps(kb:ms,j,3,2)))                 
  rdps0(kb:ms,j,2,1) = sdat *  2.0_dp*(cmps(kb:ms,j+1,4,1)*(cmps(kb:ms,j,3,1)-v0*cmps(kb:ms,j+1,3,2)) +                          &
                               cmps(kb:ms,j+1,4,2)*(v0*cmps(kb:ms,j+1,3,1)-v1*cmps(kb:ms,j,3,1)) +                               &
                               cmps(kb:ms,j,4,1)*(v0*(cmps(kb:ms,j+1,3,1)+cmps(kb:ms,j+1,3,2))-opv1*cmps(kb:ms,j,3,1)))          
  rdps0(kb:ms,j,2,2) = sdat * (cmps(kb:ms,j+1,4,1)*(omv0*cmps(kb:ms,j,3,2)-omv0*cmps(kb:ms,j+1,3,2)) +                           &
                               cmps(kb:ms,j+1,4,2)*(omv0*cmps(kb:ms,j+1,3,1)+opv1*cmps(kb:ms,j,3,1)-(v0+v1)*cmps(kb:ms,j,3,2)) - &
                               cmps(kb:ms,j,4,1)*(opv1*cmps(kb:ms,j,3,2)-opv1*cmps(kb:ms,j+1,3,2)) -                             &
                               cmps(kb:ms,j,4,2)*(omv0*cmps(kb:ms,j+1,3,1)+opv1*cmps(kb:ms,j,3,1)-(v0+v1)*cmps(kb:ms,j+1,3,2)))

!~~~jth interface, sh waves, universal representation 

  shdt = 1.0_dp/(mu(j)*beta(:,j)+mu(j+1)*beta(:,j+1))          

  tush0(:,j) = shdt * 2.0_dp*mu(j+1)*beta(:,j+1)   
  tdsh0(:,j) = shdt * 2.0_dp*mu(j)*beta(:,j)   
  rush0(:,j) = shdt * (mu(j+1)*beta(:,j+1)-mu(j)*beta(:,j))   
  rdsh0(:,j) = shdt * (mu(j)*beta(:,j)-mu(j+1)*beta(:,j+1))   

end do loop1 

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_RFMAT      
!----------------------------------------------------------------------------------------------------------------------------------! 

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_TRMAT      
!----------------------------------------------------------------------------------------------------------------------------------!

!~evaluation of the transmission/reflection matrices/coefficients: thps, rhps, thsh, rhsh

USE GFLA 
 
integer(kind=k9) :: k, j    
complex(kind=dp) :: expu(2,2), expd(2,2), exsu, exsd  

loop1: do k = 1,ms  !-loop over the sampling points 

!~~~exponential matrices/coefficients for the free surface

  call GFL_EXMAT (k, 1, 1, zh(0), zh(0), expu, expd, exsu, exsd) 

!~~~p-sv waves, reflection matrix for the free surface  
 
  rups(k,0,:,:) = matmul(rups0(k,0,:,:),expu)

!~~~sh waves, reflection coefficient for the j_th interface  

  rush(k,0) = rush0(k,0)*exsu  

  loop2: do j = 1,nl  !-loop over the layer interfaces  

!~~~~~exponential matrices/coefficients for the j_th interface

    call GFL_EXMAT (k, j+1, j, zh(j), zh(j), expu, expd, exsu, exsd) 

!~~~~~p-sv waves, transmission/reflection matrices for the j_th interface  

    tups(k,j,:,:) = matmul(tups0(k,j,:,:),expu)
    tdps(k,j,:,:) = matmul(tdps0(k,j,:,:),expd)
    rups(k,j,:,:) = matmul(rups0(k,j,:,:),expu)
    rdps(k,j,:,:) = matmul(rdps0(k,j,:,:),expd)

!~~~~~sh waves, transmission/reflection coefficients for the j_th interface  

    tush(k,j) = tush0(k,j)*exsu   
    tdsh(k,j) = tdsh0(k,j)*exsd  
    rush(k,j) = rush0(k,j)*exsu 
    rdsh(k,j) = rdsh0(k,j)*exsd 

  end do loop2 

end do loop1 

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_TRMAT    
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_EXMAT (k, ju, jd, xu, xd, expu, expd, exsu, exsd)      
!----------------------------------------------------------------------------------------------------------------------------------!

!~evaluation of the exponential matrices/coefficients: expu, expd, exsu, exsd 

USE GFLA 

integer(kind=k9), intent(in ) :: k, ju, jd  
real   (kind=dp), intent(in ) :: xu, xd      
complex(kind=dp), intent(out) :: expu(2,2), expd(2,2), exsu, exsd   
complex(kind=dp)              :: epsi 

!~sh waves, exponenetial coefficients

exsu = exp(beta(k,ju)*(xu-zh(ju))) 
exsd = exp(beta(k,jd)*(zh(jd-1)-xd)) 

!~p-sv waves, exponential matrices

if (k<kb) then  !-standard representation for moderate values of zeta

  expu = reshape((/exp(alfa(k,ju)*(xu-zh(ju)))  , zero, zero, exsu/), (/2,2/)) 
  expd = reshape((/exp(alfa(k,jd)*(zh(jd-1)-xd)), zero, zero, exsd/), (/2,2/)) 

else  !-modified representation for large values of zeta

  epsi = dfab(k,ju)/(1.0_dp+pv(ju)) * (exp(diab(k,ju)*(xu-zh(ju))) - 1.0_dp) 

  expu = exsu * reshape((/1.0_dp+epsi*difa(k,ju)*(zeta(k)+pv(ju)*beta(k,ju))          ,   & 
                          epsi*(zeta(k)+pv(ju)*beta(k,ju))*(alfa(k,ju)+pv(ju)*zeta(k)),   & 
                          epsi*difa(k,ju)*difb(k,ju), 1.0_dp+epsi*difb(k,ju)*(alfa(k,ju)+pv(ju)*zeta(k))/), (/2,2/)) 

  epsi = dfab(k,jd)/(1.0_dp+pv(jd)) * (exp(diab(k,jd)*(zh(jd-1)-xd)) - 1.0_dp) 

  expd = exsd * reshape((/1.0_dp+epsi*difa(k,jd)*(zeta(k)+pv(jd)*beta(k,jd))          ,   & 
                         -epsi*(zeta(k)+pv(jd)*beta(k,jd))*(alfa(k,jd)+pv(jd)*zeta(k)),   & 
                         -epsi*difa(k,jd)*difb(k,jd), 1.0_dp+epsi*difb(k,jd)*(alfa(k,jd)+pv(jd)*zeta(k))/), (/2,2/)) 

end if  

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_EXMAT      
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_COORD (xs, ys, zs, xo, yo, zo, ks0, s, r, t, z) 
!----------------------------------------------------------------------------------------------------------------------------------!

!~conversion from cartesian to dimensionless cylindrical coordinates   

real   (kind=dp), intent(in ) :: xs, ys, zs, xo, yo, zo, ks0   
real   (kind=dp), intent(out) :: s, r, t, z     
real   (kind=dp)              :: dx, dy, eps=1.0e-120_dp  

!~auxiliary variables 

dx = xo-xs
dy = yo-ys  

!~dimensionless cylindrical coordinates

s = ks0*zs
r = ks0*sqrt(dx*dx+dy*dy) 
t = atan2(dy,dx+sign(eps,dx)) 
z = ks0*zo  

if (r==0.0_dp.and.(z-s)==0.0_dp) then   
  print  * , 'Source and receiver points coincide' 
  print 100, xs, ys, zs  
  print 200, xo, yo, zo;  stop   
end if     
  
100 format (' x_s = (', es10.3e2, ',', es10.3e2, ',', es10.3e2, ')') 
200 format (' x_r = (', es10.3e2, ',', es10.3e2, ',', es10.3e2, ')') 

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_COORD   
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_LSNUM (zscoord) 
!----------------------------------------------------------------------------------------------------------------------------------! 

!~number of the layer containing the source 

USE GFLA

real   (kind=dp), intent(in) :: zscoord
integer(kind=k9)             :: j   

!~number of the layer containing the source     

loop1: do j = 1,ml
  if (zzcoord(j-1)<=zscoord .and. zscoord<zzcoord(j)) ls = j 
end do loop1 

if (zscoord>zzcoord(ml)) then 
  print *, 'Source depth zscoord > zzcoord(nl+1)' 
  print *, 'Thickness of the half-space must be increased';  stop  
end if  

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_LSNUM  
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_PRMAT      
!----------------------------------------------------------------------------------------------------------------------------------! 

!~evaluation of the generalized propagation matrices/coefficients: thups, thdps, thush, thdsh  

USE GFLA 

integer(kind=k9)            :: j, k   
complex(kind=dp)            :: inmat(2,2) 
complex(kind=dp), parameter :: identity(2,2)=reshape((/(1.0_dp,0.0_dp),zero,zero,(1.0_dp,0.0_dp)/),(/2,2/))

!~upward propagation matrices (interfaces above the source)  

!-propagation matrices for the 0_th interface (free-surface) 

if (ls>1) then  !-standard propagation matrices  

  rhups(:,0,:,:) = rups(:,0,:,:)   !-p-sv waves
  rhush(:,0)     = rush(:,0)       !-sh waves  

else  !-normalized propagation matrices

  rhups(:,0,:,:) = rups0(:,0,:,:)  !-p-sv waves
  rhush(:,0)     = rush0(:,0)      !-sh waves  

end if  

!-propagation matrices for the j_th interface, 1<=j<=(ls-1)    

loop1: do j = 1,ls-1   !-loop over the layer interfaces  
  loop2: do k = 1,ms  !-loop over the sampling points 

    call GFL_INMAT (identity-matmul(rdps(k,j,:,:),rhups(k,j-1,:,:)), inmat) 

    if (j<ls-1) then  !-standard propagation matrices  

      thups(k,j,:,:) = matmul(inmat,tups(k,j,:,:))                                                       !-p-sv waves
      rhups(k,j,:,:) = matmul(matmul(tdps(k,j,:,:),rhups(k,j-1,:,:)),thups(k,j,:,:)) + rups(k,j,:,:)     !-p-sv waves
      thush(k,j)     = tush(k,j)/(1.0_dp-rdsh(k,j)*rhush(k,j-1))                                         !-sh waves 
      rhush(k,j)     = tdsh(k,j)*rhush(k,j-1)*thush(k,j) + rush(k,j)                                     !-sh waves 

    else  !-normalized propagation matrices 

      thups(k,j,:,:) = matmul(inmat,tups0(k,j,:,:))                                                      !-p-sv waves
      rhups(k,j,:,:) = matmul(matmul(tdps(k,j,:,:),rhups(k,j-1,:,:)),thups(k,j,:,:)) + rups0(k,j,:,:)    !-p-sv waves       
      thush(k,j)     = tush0(k,j)/(1.0_dp-rdsh(k,j)*rhush(k,j-1))                                        !-sh waves
      rhush(k,j)     = tdsh(k,j)*rhush(k,j-1)*thush(k,j) + rush0(k,j)                                    !-sh waves

    end if  

  end do loop2
end do loop1  

!~downward propagation matrices (interfaces below the source)  

!-propagation matrices for the ml_th interface (half-space) 

rhdps(:,ml,:,:) = zero   !-p-sv waves  
rhdsh(:,ml)     = zero   !-sh waves 

!-propagation matrices for the j_th interface, ls<=j<=nl   

loop3: do j = nl,ls,-1  !-loop over the layer interfaces  
  loop4: do k = 1,ms   !-loop over the sampling points 

    call GFL_INMAT (identity-matmul(rups(k,j,:,:),rhdps(k,j+1,:,:)), inmat) 

    if (j>ls) then !-standard propagation matrices 

      thdps(k,j,:,:) = matmul(inmat,tdps(k,j,:,:))                                                       !-p-sv waves
      rhdps(k,j,:,:) = matmul(matmul(tups(k,j,:,:),rhdps(k,j+1,:,:)),thdps(k,j,:,:)) + rdps(k,j,:,:)     !-p-sv waves
      thdsh(k,j)     = tdsh(k,j)/(1.0_dp-rush(k,j)*rhdsh(k,j+1))                                         !-sh waves 
      rhdsh(k,j)     = tush(k,j)*rhdsh(k,j+1)*thdsh(k,j) + rdsh(k,j)                                     !-sh waves 

    else  !-normalized propagation matrices 

      thdps(k,j,:,:) = matmul(inmat,tdps0(k,j,:,:))                                                      !-p-sv waves  
      rhdps(k,j,:,:) = matmul(matmul(tups(k,j,:,:),rhdps(k,j+1,:,:)),thdps(k,j,:,:)) + rdps0(k,j,:,:)    !-p-sv waves 
      thdsh(k,j)     = tdsh0(k,j)/(1.0_dp-rush(k,j)*rhdsh(k,j+1))                                        !-sh waves 
      rhdsh(k,j)     = tush(k,j)*rhdsh(k,j+1)*thdsh(k,j) + rdsh0(k,j)                                    !-sh waves
  
    end if  

  end do loop4
end do loop3 
     
!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_PRMAT        
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_INMAT (matrix, inverse)     
!----------------------------------------------------------------------------------------------------------------------------------!

!~inversion of the 2x2 matrix 

complex(kind=dp), intent(in)  :: matrix(2,2) 
complex(kind=dp), intent(out) :: inverse(2,2)
complex(kind=dp)              :: det 

det = matrix(1,1)*matrix(2,2) - matrix(1,2)*matrix(2,1) 

if (abs(det)<1.0e-12_dp) then
  print *, 'Singular matrix detected in gfl_inmat';  stop 
end if  

inverse(1,1) =  matrix(2,2)/det 
inverse(1,2) = -matrix(1,2)/det
inverse(2,1) = -matrix(2,1)/det
inverse(2,2) =  matrix(1,1)/det

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_INMAT     
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_FCMAT      
!----------------------------------------------------------------------------------------------------------------------------------!

!~evaluation of the factorized matrices/coefficients: wups, wdps, wush, wdsh 

USE GFLA 

integer(kind=k9) :: j, k   

!~factorized matrices for the layers above the source 

!-factorized matrices for the (ls-1)_th layer 

if (ls>1) then  

  loop1: do k = 1,ms     

    wups(k,ls-1,:,:) = thups(k,ls-1,:,:)                           !-p-sv waves 
    wdps(k,ls-1,:,:) = matmul(rhups(k,ls-2,:,:),wups(k,ls-1,:,:))  !-p-sv waves

    wush(k,ls-1) = thush(k,ls-1)                                   !-sh waves 
    wdsh(k,ls-1) = rhush(k,ls-2)*wush(k,ls-1)                      !-sh waves   

  end do loop1 

end if  

!-factorized matrices for the j_th layer (j=ls-2,ls-3,...) 

loop2: do j = ls-2,1,-1   !-loop over the layer interfaces  
  loop3: do k = 1,ms     !-loop over the sampling points 

    wups(k,j,:,:) = matmul(thups(k,j,:,:),wups(k,j+1,:,:))         !-p-sv waves 
    wdps(k,j,:,:) = matmul(rhups(k,j-1,:,:),wups(k,j,:,:))         !-p-sv waves

    wush(k,j) = thush(k,j)*wush(k,j+1)                             !-sh waves 
    wdsh(k,j) = rhush(k,j-1)*wush(k,j)                             !-sh waves 

  end do loop3
end do loop2 
 
!~factorized matrices for the layers below the source 

!-factorized matrices for the (ls+1)_th layer 

if (ls<ml) then  

  loop4: do k = 1,ms     

    wdps(k,ls+1,:,:) = thdps(k,ls,:,:)                             !-p-sv waves
    wups(k,ls+1,:,:) = matmul(rhdps(k,ls+1,:,:),wdps(k,ls+1,:,:))  !-p-sv waves 

    wdsh(k,ls+1)     = thdsh(k,ls)                                 !-sh waves
    wush(k,ls+1)     = rhdsh(k,ls+1)*wdsh(k,ls+1)                  !-sh waves 

  end do loop4 

end if  

!-factorized matrices for the j_th layer (j=ls+2,ls+3,...) 

loop5: do j = ls+2,ml  !-loop over the layer interfaces  
  loop6: do k = 1,ms   !-loop over the sampling points 

    wdps(k,j,:,:) = matmul(thdps(k,j-1,:,:),wdps(k,j-1,:,:))       !-p-sv waves 
    wups(k,j,:,:) = matmul(rhdps(k,j,:,:),wdps(k,j,:,:))           !-p-sv waves

    wdsh(k,j) = thdsh(k,j-1)*wdsh(k,j-1)                           !-sh waves 
    wush(k,j) = rhdsh(k,j)*wdsh(k,j)                               !-sh waves 
   
  end do loop6
end do loop5  
   
!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_FCMAT        
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_SORC1 (s) 
!----------------------------------------------------------------------------------------------------------------------------------!

!~calculation of the normalized source matrices/coefficients (routine #1)  

USE GFLA

real   (kind=dp), intent(in) :: s
integer(kind=k9)             :: j, k  
complex(kind=dp)             :: coef(kb:ms), epsi(kb:ms)    

loop1: do j = 1,3     !-loop over the source directions 
  loop2: do k = -1,1  !-loop over the fourier componenets  

    if (((j==1.or.j==2).and.k==0).or.(j==3.and.(k==-1.or.k==1))) then 

!~~~~~~~p-sv waves, zero source vectors 

      sups(:,j,k,:) = 0.0_dp  
      sdps(:,j,k,:) = 0.0_dp  

!~~~~~~~sh waves, zero source coefficients  

      sush(:,j,k) = 0.0_dp   
      sdsh(:,j,k) = 0.0_dp    

    else  

!~~~~~~~p-sv waves, source matrices, standard representation for moderate values of zeta

      sups(1:ka,j,k,1) = -0.5_dp/rho(ls)*exp(alfa(1:ka,ls)*(zh(ls-1)-s)) * (zzz(j,k)-zeta(1:ka)/alfa(1:ka,ls)*xmy(j,k))   
      sups(1:ka,j,k,2) =  0.5_dp/rho(ls)*exp(beta(1:ka,ls)*(zh(ls-1)-s)) * (zzz(j,k)*zeta(1:ka)/beta(1:ka,ls)-xmy(j,k)) 
      sdps(1:ka,j,k,1) =  0.5_dp/rho(ls)*exp(alfa(1:ka,ls)*(s-zh(ls)))   * (zzz(j,k)+zeta(1:ka)/alfa(1:ka,ls)*xmy(j,k))   
      sdps(1:ka,j,k,2) =  0.5_dp/rho(ls)*exp(beta(1:ka,ls)*(s-zh(ls)))   * (zzz(j,k)*zeta(1:ka)/beta(1:ka,ls)+xmy(j,k))   

!~~~~~~~p-sv waves, source matrices, modified representation for large values of zeta

      coef = 0.5_dp/(rho(ls)*(1.0_dp+pv(ls))) * exp(beta(kb:ms,ls)*(zh(ls-1)-s))
      epsi = (zzz(j,k)-zeta(kb:ms)/alfa(kb:ms,ls)*xmy(j,k)) * (exp(diab(:,ls)*(zh(ls-1)-s))-1.0_dp)   

      sups(kb:ms,j,k,1) =  coef * ((zzz(j,k)/beta(kb:ms,ls)-xmy(j,k)/alfa(kb:ms,ls))/dfab(:,ls) + epsi*difa(:,ls))   
      sups(kb:ms,j,k,2) = -coef * ((zzz(j,k)/beta(kb:ms,ls)+pv(ls)*xmy(j,k)/alfa(kb:ms,ls))/dfab(:,ls) - &  
                                   epsi*(alfa(kb:ms,ls)+pv(ls)*zeta(kb:ms)))  

      coef = 0.5_dp/(rho(ls)*(1.0_dp+pv(ls)))*exp(beta(kb:ms,ls)*(s-zh(ls)))
      epsi = (zzz(j,k)+zeta(kb:ms)/alfa(kb:ms,ls)*xmy(j,k)) * (exp(diab(:,ls)*(s-zh(ls)))-1.0_dp)  

      sdps(kb:ms,j,k,1) = -coef * ((zzz(j,k)/beta(kb:ms,ls)+xmy(j,k)/alfa(kb:ms,ls))/dfab(:,ls) + epsi*difa(:,ls))   
      sdps(kb:ms,j,k,2) = -coef * ((zzz(j,k)/beta(kb:ms,ls)-pv(ls)*xmy(j,k)/alfa(kb:ms,ls))/dfab(:,ls) - &  
                                   epsi*(alfa(kb:ms,ls)+pv(ls)*zeta(kb:ms)))  
 
!~~~~~~~sh waves, source coefficients, universal representation 

      sush(:,j,k) = -0.5_dp*ic*exp(beta(:,ls)*(zh(ls-1)-s)) * xpy(j,k)/(mu(ls)*zeta*beta(:,ls))          
      sdsh(:,j,k) = -0.5_dp*ic*exp(beta(:,ls)*(s-zh(ls)))   * xpy(j,k)/(mu(ls)*zeta*beta(:,ls))           

    end if  

  end do loop2
end do loop1 

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_SORC1 
!----------------------------------------------------------------------------------------------------------------------------------!
  
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_SORC2  
!----------------------------------------------------------------------------------------------------------------------------------!

!~calculation of the normalized source matrices/coefficients (routine #2)  

USE GFLA

integer(kind=k9)            :: i, j, k      
complex(kind=dp)            :: inpsu(2,2), inpsd(2,2), inshu, expu(2,2), expd(2,2), exsu, exsd          
complex(kind=dp), parameter :: identity(2,2)=reshape((/(1.0_dp,0.0_dp),zero,zero,(1.0_dp,0.0_dp)/),(/2,2/))

!~p-sv waves 

loop1: do k = 1,ms  !-loop over the sampling points 
 
!~~~exponential matrices/coefficients 

  call GFL_EXMAT (k, ls, ls, zh(ls-1), zh(ls), expu, expd, exsu, exsd)     

!~~~auxiliary matrices/coefficients 

  call GFL_INMAT (identity-matmul(matmul(rhdps(k,ls,:,:),expd),matmul(rhups(k,ls-1,:,:),expu)), inpsu)   !-p-sv waves
  call GFL_INMAT (identity-matmul(matmul(rhups(k,ls-1,:,:),expu),matmul(rhdps(k,ls,:,:),expd)), inpsd)   !-p-sv waves 

  inshu = 1.0_dp/(1.0_dp-rhdsh(k,ls)*exsd*rhush(k,ls-1)*exsu)   !-sh-waves  

  loop2: do i = 1,3     !-loop over the source directions 
    loop3: do j = -1,1  !-loop over the fourier componenets  

      if (((i==1.or.i==2).and.j==0).or.(i==3.and.(j==-1.or.j==1))) then 

!~~~~~~~~~p-sv waves, zero source vectors 

        mups(k,i,j,:) = 0.0_dp
        mdps(k,i,j,:) = 0.0_dp
        sups(k,i,j,:) = 0.0_dp
        sdps(k,i,j,:) = 0.0_dp 

!~~~~~~~~~sh waves, zero source coefficients  

        mush(k,i,j) = 0.0_dp
        mdsh(k,i,j) = 0.0_dp  
        sush(k,i,j) = 0.0_dp
        sdsh(k,i,j) = 0.0_dp  

      else 

!~~~~~~~~~p-sv waves, source vectors 

        mups(k,i,j,:) = matmul(inpsu, matmul(rhdps(k,ls,:,:),sdps(k,i,j,:))   + &
                        matmul(matmul(matmul(rhdps(k,ls,:,:),expd),rhups(k,ls-1,:,:)),sups(k,i,j,:)))
        mdps(k,i,j,:) = matmul(inpsd, matmul(rhups(k,ls-1,:,:),sups(k,i,j,:)) + &
                        matmul(matmul(matmul(rhups(k,ls-1,:,:),expu),rhdps(k,ls,:,:)),sdps(k,i,j,:)))

        sups(k,i,j,:) = matmul(expu,mups(k,i,j,:)) + sups(k,i,j,:)
        sdps(k,i,j,:) = matmul(expd,mdps(k,i,j,:)) + sdps(k,i,j,:)

!~~~~~~~~~sh waves, source coefficients  

        mush(k,i,j) = inshu * (rhdsh(k,ls)*sdsh(k,i,j)   + rhdsh(k,ls)*exsd*rhush(k,ls-1)*sush(k,i,j))
        mdsh(k,i,j) = inshu * (rhush(k,ls-1)*sush(k,i,j) + rhush(k,ls-1)*exsu*rhdsh(k,ls)*sdsh(k,i,j)) 

        sush(k,i,j) = exsu*mush(k,i,j) + sush(k,i,j)
        sdsh(k,i,j) = exsd*mdsh(k,i,j) + sdsh(k,i,j)  

     end if    

    end do loop3
  end do loop2

end do loop1   

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_SORC2 
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_LONUM (zocoord) 
!----------------------------------------------------------------------------------------------------------------------------------!

!~number of the layer containing the observation point     

USE GFLA

real   (kind=dp), intent(in) :: zocoord
integer(kind=k9)             :: j   

loop1: do j = 1,ml 
  if (zzcoord(j-1)<=zocoord .and. zocoord<zzcoord(j)) lo = j 
end do loop1

if (zocoord>zzcoord(ml)) then 
  print *, 'Warning: zocoord > zzcoord(nl+1)' 
  print *, 'Thickness of the half-space must be increased';  stop  
end if  

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_LONUM  
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_SINPA (s, r, t, z, zocoord, zscoord) 
!----------------------------------------------------------------------------------------------------------------------------------!

!~calculation of parameters needed for the singularity treatment      

USE GFLA
USE GFLB

real   (kind=dp), intent(in) :: s, r, t, z , zocoord, zscoord

!~material properties of the 'lower' half space in bimaterial formulation 

mst2 = mu(ls)
vst2 = v(ls) 

!~local coordinate system (rst,tst,zst) and the material properties of the 'upper' half space in bimaterial formulation 

if ((zscoord-zzcoord(ls-1))<=(zzcoord(ls)-zscoord).or.ls==ml) then 
  if(ls>1) then   
    mst1 = mu(ls-1)
    vst1 = v(ls-1)  
  else 
    mst1 = 0.0_dp
    vst1 = 0.0_dp 
  end if 
  sst  =  s-zh(ls-1) 
  rst  =  r 
  tst  =  t 
  zst  =  z-zh(ls-1) 
  sing = 1 
else 
  mst1 =  mu(ls+1)
  vst1 =  v(ls+1) 
  sst  =  zh(ls)-s 
  rst  =  r 
  tst  = -t 
  zst  =  zh(ls)-z 
  sing = -1
end if   

!print 100, 'sst=',sst, ' rst=',rst, ' tst=',tst, ' zst=',zst, ' sing=',sing
!100 format (1x, 4(a, es12.3e3, 1x), a, i2)
!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_SINPA  
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_THETA (t)   
!----------------------------------------------------------------------------------------------------------------------------------!

!~exponential and trigonometric functions of the angular coordinate theta  

USE GFLA 
USE GFLB 

real   (kind=dp), intent(in) :: t

expp = exp( ic*t)
expm = exp(-ic*t)

stat0 = (/0.0_dp, 0.0_dp, 1.0_dp/) 
statc = (/cos(tst), cos(tst-1.570796326794897_dp), 0.0_dp/) 
stats = (/sin(tst), sin(tst-1.570796326794897_dp), 0.0_dp/) 

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_THETA 
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_XHAT (s, r, z, uhat, tdhat, tshat) 
!----------------------------------------------------------------------------------------------------------------------------------!

!~adaptive integration procedure (numerical inverse hankel transform)   

USE GFLA
USE GFLB
USE GFLC 
USE GFLR 

real   (kind=dp), intent(in ) :: s, r, z 
complex(kind=dp), intent(out) :: uhat(3,3), tdhat(3,6), tshat(3,6) 
integer(kind=k9)              :: k, inter, irec, jref, jstp, loc(2)  
real   (kind=dp)              :: difu(3,3), dift(3,6), refu(3,3), reft(3,6) 
complex(kind=dp)              :: uint(3,3), tint(3,6), uold(3,3), told(3,6), uinc(3,3), tinc(3,6), ushat(3,3), that(3,6)    
real   (kind=dp), parameter   :: epsr=2.0e-3_dp, epsa=2.0e-7_dp, tran=10.0_dp            

!~check if the radial distance is too large for accurate evaluation 

if (r>20.0_dp*(hcont_max/hcont)) then
  print *, 'Warning: frequency/radial_distance too large for accurate integration' 
  print *, 'Maximium allowed (r_bar) =', real(20.0_dp*(hcont_max/hcont)) 
  print *, 'Actual           (r_bar) =', real(r) 
  print *, 'Decrease the contour height';  stop
end if  

!~check if series representation of bessel functions applies throughout the integration path      
     
if (r*lpath<tran) then
  rflag=1;  rfac=1.0_dp    !-yes, radius is "small"; no external normalization needed to calculate tdhat(:,k) & tshat(:,k), k=1,2,6
else
  rflag=2;  rfac=1.0_dp/r  !-no, radius is "large"; external normalization needed to calculate tdhat(:,k) & tshat(:,k), k=1,2,6
end if  
    
!~initialization of the green's functions  

call GFL_BIMAT (sing, mst1, vst1, mst2, vst2, sst, rst, tst, zst, mu0, ks0, ushat, tshat)   

if (dyngreen == 0) then
    tdhat=(0.0_dp,0.0_dp)
    uhat=ushat
    return
endif

uhat = magfac*ushat  !-assign the static bi-material values  
that = magfac*tshat  !-magfac is used to improve accuracy of far-field calculations     

!~loop over the contour intervals (1:/, 2:/, 3:--, 4:--, 5:--, 6:\, 7:__, 8:__, 9:__)

loop1: do inter = 1,nint  

!---initialize flags of the sampling points for each interval

flag = 0
                
!~~~recursive integration of the inter_th contour interval (uint, tint) 

  loop2: do irec = 0,nrec  !-irec is the recursion nuber     

    if (irec>0) then 
      uold = uint
      told = tint
    end if  
    uint = 0.0_dp 
    tint = 0.0_dp 

!-----loop over the subintervals (uinc, tinc)  

    jref = 1 + (inter-1)*(ns-1)  !-reference number for the subinterval location  
    jstp = 2**(nrec-irec)        !-subinterval length for the irec_th recursiion     

    loop3: do k = 1,2**irec  
      call GFL_XINC (inter, jref, jstp, s, r, z, tran, uinc, tinc)  !-subinterval integration 
      uint = uint + uinc                                  
      tint = tint + tinc                                  
      jref = jref + 7*jstp     
    end do loop3 

!-----general convergence criterion for recursive integration  

    if (irec>0) then 

      refu = max(epsr*max(abs(uold),abs(uhat)), epsa)  !-treshold value for displacements  
      reft = max(epsr*max(abs(told),abs(that)), epsa)  !-treshold value for stresses  

      if (all(abs(uint-uold)<refu).and.all(abs(tint-told)<reft)) exit loop2  

    end if   

!-----warning statement  

    if (irec==nrec) then 

      difu = 0.0_dp 
      dift = 0.0_dp 

      where (abs(uint-uold)>refu) difu = abs(uint-uold)
      where (abs(tint-told)>reft) dift = abs(tint-told)

      if (maxval(difu)>maxval(dift)) then 
        loc = maxloc(difu) 
        if (abs(uint(loc(1),loc(2)))>epsa) print 100, inter, 'Uint =', abs(uint(loc(1),loc(2))), maxval(difu), loc         
      else  
        loc = maxloc(dift) 
        if (z>0.0_dp.and.abs(tint(loc(1),loc(2)))>epsa) print 100, inter, 'Tint =', abs(tint(loc(1),loc(2))), maxval(dift), loc
      end if  

    end if  

  end do loop2   !-end recursive integration 

!~~~add contributions from the inter_th contour interval

  uhat = uhat + uint    
  that = that + tint     

usmp_global((inter-1)*(ns-1)+1:inter*(ns-1)+1,:,:,:)=usmp(:,:,:,:)

!~~~general convergence criterion for addaptive integration   
  
  if (inter>6.and.all(abs(uint)<max(epsr*abs(uhat),epsa)).and.all(abs(tint)<max(epsr*abs(that),epsa))) exit loop1

end do loop1  !-end loop over the contour intervals   





!~undo the magnification used to aid accuracy of far-field calculations

uhat = uhat/magfac
that = that/magfac 

!~regular (dynamic) part of the stress green's functions 

tdhat = that - tshat 

100 format ('Slow convergence, Int =', i2, 2x, a, es8.1e2, 2x, 'Diff =', es8.1e2, 1x, 2i2) 

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_XHAT   
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_XINC (inter, jref, jstp, s, r, z, tran, uinc, tinc)  
!----------------------------------------------------------------------------------------------------------------------------------!

!~integration routine for a single subinterval   

USE GFLA 
USE GFLR  

integer(kind=k9), intent(in ) :: inter, jref, jstp      
real   (kind=dp), intent(in ) :: s, r, z, tran      
complex(kind=dp), intent(out) :: uinc(3,3), tinc(3,6)   
integer(kind=k9)              :: j, k       
real   (kind=dp)              :: cofa(0:2), cofb(0:2)             
complex(kind=dp)              :: zlow, zupp, dzet, zet0, ztra, jint(0:2,8), cint(8), sint(8), wght(0:2,8)    
complex(kind=dp)              :: fu(3,3,0:2,8), ft(3,8,0:2,8), taux(3,2), pp(0:2,8,8), qq(0:2,8,8)  

!-initialization 

uinc = 0.0_dp;  tinc = 0.0_dp;  taux = 0.0_dp;  pp = 0.0_dp;  qq = 0.0_dp;  wght = 0.0_dp;  cofa = (/1.0_dp,r,r/);  cofb = 1.0_dp

!~interval parameters 

zlow = zeta(jref)          !-lower limit of the interval
zupp = zeta(jref+7*jstp)   !-upper limit of the interval 
dzet = (zupp-zlow)/7.0_dp  !-spacing between the sampling points   
zet0 = zlow + 3.0_dp*dzet  !-reference point for numerical integration  

!~sampling of the integrand to be splined 

loop1: do j = 1,8 
  call GFL_RKER (s, z, jref+(j-1)*jstp, jref+(j-1)*jstp-(inter-1)*(ns-1), fu(:,:,:,j), ft(:,:,:,j))   
end do loop1     

!~semi-analytical integration of Int (zet^k Jn(r*zet) dzet) from zlow to zupp, k=0,7, j=0,2 

if (r*abs(zupp)<=tran) then  !-integration based on series expansion of bessel functions (rflag=1/2) 

  call GFL_JINT (r, zlow, zupp, zet0, dzet, jint)   !-jint(1,:) & jint(2,:) are divided (normalized) by r

  loop2: do k = 0,2
    wght(k,:) = matmul(cc,jint(k,:))
  end do loop2   

  if (rflag==2)  cofb = cofa  !-undo internal normalization of jint when calculating taux for large r 

else !-integration based on asymptotic expansion of bessel functions (rflag=2) 
                     
  call GFL_CINT (r, zlow, zupp, zet0, dzet, cint, sint)  

  loop3: do j = 1,8   
    call GFL_PPQQ (r*zeta(jref+(j-1)*jstp), pp(:,j,j), qq(:,j,j))  
  end do loop3 

  loop4: do k = 0,2
    wght(k,:) = matmul(matmul(pp(k,:,:),cc),cint) + matmul(matmul(qq(k,:,:),cc),sint)
  end do loop4  

  cofa = 1.0_dp  !-normalization parameter for uinc & tinc 
 
end if  

!~assembly loop: interpolating polynomials times the bessel integrals   
  
loop5: do j = 1,3  !-loop over the source directions   
  loop6: do k = 0,2  !-loop over the orders of bessel functions (0,1,2) 

    uinc(j,:) = uinc(j,:) + cofa(k)*matmul(fu(j,1:3,k,:),wght(k,:))  !-cofa/=1 cancels internal normalization of jint(k,:), k=1,2       
    tinc(j,:) = tinc(j,:) + cofa(k)*matmul(ft(j,1:6,k,:),wght(k,:))  !-cofa/=1 cancels internal normalization of jint(k,:), k=1,2   
    taux(j,:) = taux(j,:) + cofb(k)*matmul(ft(j,7:8,k,:),wght(k,:))  !-cofb/=1 cancels internal normalization of jint(k,:), k=1,2  
      
  end do loop6
end do loop5    

!~stress green's function increments: t_rr, t_tt, t_rt    

tinc(:,1) = tinc(:,1) - rfac*taux(:,1)   
tinc(:,2) = tinc(:,2) + rfac*taux(:,1)  
tinc(:,6) = tinc(:,6) - rfac*taux(:,2)      

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_XINC  
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_JINT (r, zl, zu, z0, dz, jint)
!----------------------------------------------------------------------------------------------------------------------------------!

!~semi-analytical evaluation of the integrals:  jint(n,k) = Int (((zet-zet0)/dz)^k J_n(r*zet) dzet), zet=zl,zu, k=0,7, j=0,2   
!~note: integrals jint(1,k) and jint(2,k) are normalized (divided) by r 

USE GFLA  
USE GFLJ 

real   (kind=dp), intent(in ) :: r
complex(kind=dp), intent(in ) :: zl, zu, z0, dz
complex(kind=dp), intent(out) :: jint(0:2,8)  
integer(kind=k9)              :: i, j   
complex(kind=dp)              :: rzls, rzus, trmil(0:2), trmiu(0:2), coeff(0:2), trmjl(0:2), trmju(0:2), yint(0:2,0:7), yold(0:2)

!-initialization of local and global variables 

yint = 0.0_dp
jint = 0.0_dp  

rzls = (0.5_dp*r*zl)**2 
rzus = (0.5_dp*r*zu)**2

trmil = (/zl, 0.50_dp*zl**2, 0.25_dp*r*zl**3/)    
trmiu = (/zu, 0.50_dp*zu**2, 0.25_dp*r*zu**3/)  

!-computation loop: yint(n,i) = integral(z^i)*besselj(n,r*z))

loop1: do i = 0,7  

  coeff = (/1.0_dp,1.0_dp,0.5_dp/)
  trmjl = trmil
  trmju = trmiu  

  loop2: do j = 1,jmax 

    yold = yint(:,i) 
    
    yint(:,i) = yint(:,i) + coeff*(trmju-trmjl)/(n+2*j+i-1) 

    if (all(abs(yint(:,i)-yold)/(abs(yold)+o)<eps)) exit loop2 

    coeff = -coeff/(j*(n+j)) 
    trmjl =  trmjl*rzls 
    trmju =  trmju*rzus   

  end do loop2 

  if (j>=jmax) stop 'No convergence in gfl_jint -> program terminated'

  trmil = trmil*zl
  trmiu = trmiu*zu  

  loop3: do j = 0,i 

    jint(:,i+1) = jint(:,i+1) + bn(i,j)*yint(:,i-j)*(-z0)**j

  end do loop3 

  jint(:,i+1) = jint(:,i+1)/dz**i

end do loop1 

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_JINT   
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_JINT2 (r, zl, zu, z0, dz, jint)
!----------------------------------------------------------------------------------------------------------------------------------!

!===purpose===! 

!-semi-analytical evaluation of the integrals:  jint(k,j) = Int (((zet-zet0)/dz)^j J_k(r*zet) dzet), zet=zl,zu, j=0,7, k=0,1,2

USE KINDS
USE GFLA 
USE GFLJ 

real   (kind=dp), intent(in ) :: r
complex(kind=dp), intent(in ) :: zl, zu, z0, dz 
complex(kind=dp), intent(out) :: jint(0:2,8)
integer(kind=k9)              :: i, j, l, m  
real   (kind=dp)              :: rbs
complex(kind=dp)              :: a(0:2,0:1), bj(0:2), rzls, rzus     
complex(kind=dp)              :: trmil(0:2), trmiu(0:2), coeff(0:2), trmjl(0:2), trmju(0:2), yint(0:2,0:7), yold(0:2) 
integer(kind=k9), allocatable :: binom(:,:), binom1(:,:)
complex(kind=dp), allocatable :: coefa(:), coefa1(:)

!-initialization 

jint = 0.0_dp  

if (abs(dz)<1.0_dp) then

  rbs = (0.5_dp*r)**2 

  loop1: do i = 0,7

    j = 0

    allocate (coefa1(0:2*j+2), binom1(2,0:2*j+2))

    binom1 = 0;  binom1(:,0) = 1;  binom1(1,1) = 1;  binom1(2,1) = 2;  binom1(2,2) = 1

    coefa1(0) =      dz*(4.0_dp**(i+1)-(-3.0_dp)**(i+1)) / (i+1)
    coefa1(1) = (dz**2)*(4.0_dp**(i+2)-(-3.0_dp)**(i+2)) / (i+2)
    coefa1(2) = (dz**3)*(4.0_dp**(i+3)-(-3.0_dp)**(i+3)) / (i+3) 

    a(:,:) = 0.0_dp
    a(0,0) = coefa1(0)
    a(1,0) = z0*coefa1(0) + coefa1(1)
    a(2,0) = (z0**2)*coefa1(0) + 2._dp*z0*coefa1(1) + coefa1(2)

    jint(0,i+1) = a(0,0)
    jint(1,i+1) = 0.5_dp*r*a(1,0)
    jint(2,i+1) = 0.5_dp*rbs*a(2,0)

    bj = jint(:,i+1)

    j = 1
 
    loop2: do 

      allocate (binom(2,0:2*j+2), coefa(0:2*j+2))

      binom = 0;  coefa(0:2*(j-1)+2) = coefa1;  binom(:,0:2*(j-1)+2) = binom1

      deallocate (coefa1, binom1)

      yold = jint(:,i+1)
  
      if (j/=1)  a(:,0) = a(:,1)

      a(:,1) = 0.0_dp

      loop3: do m = 1,2

        binom(1,:) = binom(2,:) 

        loop4: do l = 1,2*j+2
          binom(2,l) = binom(1,l-1) + binom(1,l)
        end do loop4 
     
    end do loop3

      coefa(2*j+1) = (dz**(2*j+2))*(4.0_dp**(i+2*j+2) - (-3.0_dp)**(i+2*j+2))/(i+2*j+2)
      coefa(2*j+2) = (dz**(2*j+3))*(4.0_dp**(i+2*j+3) - (-3.0_dp)**(i+2*j+3))/(i+2*j+3)

      a(0,1) = a(2,0)

      loop5: do l = 0,2*j+2

        if (l<=2*j+1)  a(1,1) = a(1,1) + binom(1,l)*(z0**l)*coefa(2*j+1-l)
   
        a(2,1) = a(2,1) + binom(2,l)*(z0**l)*coefa(2*j+2-l)

      end do loop5

      bj = -rbs*a(:,1)*bj / (j*(j+n)*a(:,0))
    
      jint(:,i+1) = jint(:,i+1) + bj 

      if (all(abs((jint(:,i+1)-yold)/jint(:,i+1))<eps))  exit

      allocate (binom1(2,0:2*j+2), coefa1(0:2*j+2))

      binom1 = binom;  coefa1 = coefa

      j = j + 1

      deallocate (binom, coefa)

    end do loop2

    deallocate (binom, coefa)

  end do loop1 

else  !-case when dz>1 

  yint = 0.0_dp

  rzls = (0.5_dp*r*zl)**2 
  rzus = (0.5_dp*r*zu)**2

  trmil = (/zl, 0.50_dp*r*zl**2, 0.25_dp*r**2*zl**3/)    
  trmiu = (/zu, 0.50_dp*r*zu**2, 0.25_dp*r**2*zu**3/)   

!---computation loop: yint(n,i) = integral(z^i)*besselj(n,r*z))

  loop6: do i = 0,7  

    coeff = (/1.0_dp,1.0_dp,0.5_dp/)
    trmjl = trmil
    trmju = trmiu  

    loop7: do j = 1,jmax

      yold = yint(:,i) 
    
      yint(:,i) = yint(:,i) + coeff*(trmju-trmjl)/(n+2*j+i-1) 

      if (all(abs(yint(:,i)-yold)/(abs(yold)+o)<eps)) exit loop7  

      coeff = -coeff/(j*(n+j)) 
      trmjl =  trmjl*rzls 
      trmju =  trmju*rzus   

    end do loop7 

    if (j>=jmax) stop 'No convergence in gfl_jint ==> program terminated'

    trmil = trmil*zl
    trmiu = trmiu*zu  

    loop8: do j = 0,i 
      jint(:,i+1) = jint(:,i+1) + bn(i,j)*yint(:,i-j)*(-z0)**j
    end do loop8 

    jint(:,i+1) = jint(:,i+1)/dz**i

  end do loop6 

end if 

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_JINT2   
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_CINT (r, zl, zu, z0, dz, cint, sint)
!----------------------------------------------------------------------------------------------------------------------------------!

!~semi-analytical integration of trigonometric functions times powers of zeta  (Int (zet^k sin/cos(r*zet) dzet), k=0,7, j=0,2)

USE GFLA

real   (kind=dp), intent(in ) :: r
complex(kind=dp), intent(in ) :: zl, zu, dz, z0
complex(kind=dp), intent(out) :: cint(8), sint(8)   
integer(kind=k9)              :: k  
complex(kind=dp)              :: lz, dx, dxi, xp, xm, xmp(7), czm, czp, szp, sc0, sc1, sc3, sc4, sc5, sc6, sc7, a(7), cf   
real   (kind=dp), parameter   :: h=0.5_dp  
 
lz  = zu-zl 
dx  = r*lz/7.0_dp  
dxi = 1.0_dp/(dx*dx)
xp  = 0.5_dp*r*(zu+zl) 
xm  = 0.5_dp*r*(zu-zl)
cf  = dx/(r*dz) 
czm = cos(xm) 
czp = cos(xp) 
szp = sin(xp)

loop1 : do k = 1,7
 a(k) = ((xp-r*z0)/dx)**k
end do loop1 

!-auxiliary expressions: sc0, sc1, sc3, sc4, sc5, sc6, sc7

if (abs(xm)>1.0_dp) then  !-closed format of the auxiliary integrals for (r*zm)>1

  sc0 = sin(xm)/xm

  sc1 = sc0 - czm 

  sc3 = (3.675e1_dp+3.0_dp*a(2))*sc0 - & 
        (1.225e1_dp+3.0_dp*a(2))*czm - 0.6e1_dp*sc1*dxi

  sc4 = (1.47e2_dp*a(1)+4.0_dp*a(3))*sc0 - & 
        (4.90e1_dp*a(1)+4.0_dp*a(3))*czm - 2.4e1_dp*a(1)*sc1*dxi 

  sc5 = (7.503125e2_dp+3.675e2_dp*a(2)+5.0_dp*a(4))*sc0 -  &
        (1.500625e2_dp+1.225e2_dp*a(2)+5.0_dp*a(4))*czm - 2.0e1_dp*sc3*dxi

  sc6 = (4.501875e3_dp*a(1)+7.350e2_dp*a(3)+6.0_dp*a(5))*sc0 - &
        (9.003750e2_dp*a(1)+2.450e2_dp*a(3)+6.0_dp*a(5))*czm - 3.0e1_dp*sc4*dxi

  sc7 = (1.2867859375e4_dp+1.57565625e4_dp*a(2)+1.28625e3_dp*a(4)+7.0_dp*a(6))*sc0 - &
        (1.8382656250e3_dp+3.15131250e3_dp*a(2)+4.28750e2_dp*a(4)+7.0_dp*a(6))*czm - 4.2e1_dp*sc5*dxi

else  !-series expansion of the auxiliary integrals for (r*zm)<1

  xmp(1) = xm*xm
  loop2 : do k = 2,7
     xmp(k) = xmp(k-1)*xm*xm 
  end do loop2 

  sc0 =   1.0000000000000000e+00_dp        - 1.6666666666666667e-01_dp*xmp(1) + & 
          8.3333333333333333e-03_dp*xmp(2) - 1.9841269841269841e-04_dp*xmp(3) + &    
          2.7557319223985891e-06_dp*xmp(4) - 2.5052108385441719e-08_dp*xmp(5) + & 
          1.6059043836821615e-10_dp*xmp(6) - 7.6471637318198165e-13_dp*xmp(7)

  sc1 =   3.3333333333333333e-01_dp*xmp(1) - 3.3333333333333333e-02_dp*xmp(2) + & 
          1.1904761904761905e-03_dp*xmp(3) - 2.2045855379188712e-05_dp*xmp(4) + & 
          2.5052108385441719e-07_dp*xmp(5) - 1.9270852604185938e-09_dp*xmp(6) + &
          1.0706029224547743e-11_dp*xmp(7)

  sc3 = ( 2.450000000000000e+00_dp + 1.000000000000000e+00_dp*a(2))*xmp(1) + &
        (-2.916666666666666e-01_dp - 1.000000000000000e-01_dp*a(2))*xmp(2) + &
        ( 1.134259259259259e-02_dp + 3.571428571428571e-03_dp*a(2))*xmp(3) + &
        (-2.209595959595960e-04_dp - 6.613756613756615e-05_dp*a(2))*xmp(4) + &
        ( 2.596747388414056e-06_dp + 7.515632515632516e-07_dp*a(2))*xmp(5) + &
        (-2.045922184811074e-08_dp - 5.781255781255782e-09_dp*a(2))*xmp(6) + &
        ( 1.124133068577513e-10_dp + 3.211808767364323e-11_dp*a(2))*xmp(7)

  sc4 = ( 9.800000000000000e+00_dp*a(1) + 1.333333333333333e+00_dp*a(3))*xmp(1) + &
        (-1.166666666666666e+00_dp*a(1) - 1.333333333333333e-01_dp*a(3))*xmp(2) + &
        ( 4.537037037037037e-02_dp*a(1) + 4.761904761904762e-03_dp*a(3))*xmp(3) + &
        (-8.838383838383840e-04_dp*a(1) - 8.818342151675490e-05_dp*a(3))*xmp(4) + &
        ( 1.038698955365622e-05_dp*a(1) + 1.002084335417669e-06_dp*a(3))*xmp(5) + &
        (-8.183688739244290e-08_dp*a(1) - 7.708341041674374e-09_dp*a(3))*xmp(6) + &
        ( 4.496532274310052e-10_dp*a(1) + 4.282411689819098e-11_dp*a(3))*xmp(7)

  sc5 = ( 2.143750000000000e+01_dp + 2.450000000000001e+01_dp*a(2) + 1.666666666666666e+00_dp*a(4))*xmp(1) + & 
        (-2.778935185185184e+00_dp - 2.916666666666667e+00_dp*a(2) - 1.666666666666667e-01_dp*a(4))*xmp(2) + &
        ( 1.136837121212121e-01_dp + 1.134259259259259e-01_dp*a(2) + 5.952380952380953e-03_dp*a(4))*xmp(3) + &
        (-2.290331196581196e-03_dp - 2.209595959595960e-03_dp*a(2) - 1.102292768959435e-04_dp*a(4))*xmp(4) + &
        ( 2.756880144032921e-05_dp + 2.596747388414056e-05_dp*a(2) + 1.252605419272086e-06_dp*a(4))*xmp(5) + &
        (-2.203300814411925e-07_dp - 2.045922184811074e-07_dp*a(2) - 9.635426302092970e-09_dp*a(4))*xmp(6) + &
        ( 1.265682912690674e-09_dp + 1.157195805888616e-09_dp*a(2) + 5.353014612273873e-11_dp*a(4))*xmp(7)

  sc6 = ( 1.286249999999999e+02_dp*a(1) + 4.900000000000002e+01_dp*a(3) + 2.000000000000000e+00_dp*a(5))*xmp(1) + &
        (-1.667361111111111e+01_dp*a(1) - 5.833333333333334e+00_dp*a(3) - 2.000000000000000e-01_dp*a(5))*xmp(2) + &
        ( 6.821022727272728e-01_dp*a(1) + 2.268518518518519e-01_dp*a(3) + 7.142857142857143e-03_dp*a(5))*xmp(3) + &
        (-1.374198717948718e-02_dp*a(1) - 4.419191919191920e-03_dp*a(3) - 1.322751322751323e-04_dp*a(5))*xmp(4) + &
        ( 1.654128086419753e-04_dp*a(1) + 5.193494776828113e-05_dp*a(3) + 1.503126503126503e-06_dp*a(5))*xmp(5) + &
        (-1.321980488647155e-06_dp*a(1) - 4.091844369622148e-07_dp*a(3) - 1.156251156251156e-08_dp*a(5))*xmp(6) + &
        ( 7.594097476144045e-09_dp*a(1) + 2.314391611777233e-09_dp*a(3) + 6.423617534728647e-11_dp*a(5))*xmp(7)

  sc7 = ( 2.042517361111108e+02_dp      + 4.501875000000003e+02_dp*a(2)           & 
         +8.575000000000000e+01_dp*a(4) + 2.333333333333333e+00_dp*a(6))*xmp(1) + & 
        (-2.785250946969697e+01_dp      - 5.835763888888889e+01_dp*a(2)           & 
         -1.020833333333333e+01_dp*a(4) - 2.333333333333333e-01_dp*a(6))*xmp(2) + & 
        ( 1.178375400641025e+00_dp      + 2.387357954545455e+00_dp*a(2)           & 
         +3.969907407407407e-01_dp*a(4) + 8.333333333333333e-03_dp*a(6))*xmp(3) + & 
        (-2.431568287037037e-02_dp      - 4.809695512820514e-02_dp*a(2)           & 
         -7.733585858585857e-03_dp*a(4) - 1.543209876543210e-04_dp*a(6))*xmp(4) + & 
         (2.975695456164207e-04_dp      + 5.789448302469136e-04_dp*a(2)           & 
         +9.088615859449190e-05_dp*a(4) + 1.753647586980920e-06_dp*a(6))*xmp(5) + & 
        (-2.422441153915189e-06_dp      - 4.643942488611607e-06_dp*a(2)           & 
         -7.160727646838757e-07_dp*a(4) - 1.348959682293016e-08_dp*a(6))*xmp(6) + & 
        ( 1.405389140626560e-08_dp      + 2.663529767422312e-08_dp*a(2)           & 
         +4.050185320610157e-09_dp*a(4) + 7.494220457183420e-11_dp*a(6))*xmp(7)

end if  !-end calculation of the auxiliary integrals

!-integrals: cint(k)=int(((z-z0)/dz)**k cos(r*z)), sint(k)=int(((z-z0)/dz)**k sin(r*z))

cint(1) = lz*czp*sc0
sint(1) = lz*szp*sc0  

cint(2) = lz*cf*(a(1)*czp*sc0 - szp*sc1/dx) 
sint(2) = lz*cf*(a(1)*szp*sc0 + czp*sc1/dx)

cint(3) = lz*cf**2*(h*czp*((2.45e1_dp+2.0_dp*a(2))*sc0-0.4e1_dp*sc1*dxi)-2.0_dp*a(1)*szp*sc1/dx)
sint(3) = lz*cf**2*(h*szp*((2.45e1_dp+2.0_dp*a(2))*sc0-0.4e1_dp*sc1*dxi)+2.0_dp*a(1)*czp*sc1/dx) 

cint(4) = lz*cf**3*(h*czp*((7.35e1_dp*a(1)+2.0_dp*a(3))*sc0-1.2e1_dp*a(1)*sc1*dxi)-szp*sc3/dx)
sint(4) = lz*cf**3*(h*szp*((7.35e1_dp*a(1)+2.0_dp*a(3))*sc0-1.2e1_dp*a(1)*sc1*dxi)+czp*sc3/dx) 

cint(5) = lz*cf**4*(h*czp*((3.00125e2_dp+1.47e2_dp*a(2)+2.0_dp*a(4))*sc0-0.8e1_dp*sc3*dxi)-szp*sc4/dx)
sint(5) = lz*cf**4*(h*szp*((3.00125e2_dp+1.47e2_dp*a(2)+2.0_dp*a(4))*sc0-0.8e1_dp*sc3*dxi)+czp*sc4/dx)

cint(6) = lz*cf**5*(h*czp*((1.500625e3_dp*a(1)+2.450e2_dp*a(3)+2.0_dp*a(5))*sc0-1.0e1_dp*sc4*dxi)-szp*sc5/dx)
sint(6) = lz*cf**5*(h*szp*((1.500625e3_dp*a(1)+2.450e2_dp*a(3)+2.0_dp*a(5))*sc0-1.0e1_dp*sc4*dxi)+czp*sc5/dx)

cint(7) = lz*cf**6*(h*czp*((3.67653125e3_dp+4.501875e3_dp*a(2)+3.675e2_dp*a(4)+2.0_dp*a(6))*sc0-1.2e1_dp*sc5*dxi)-szp*sc6/dx)
sint(7) = lz*cf**6*(h*szp*((3.67653125e3_dp+4.501875e3_dp*a(2)+3.675e2_dp*a(4)+2.0_dp*a(6))*sc0-1.2e1_dp*sc5*dxi)+czp*sc6/dx)  

cint(8) = lz*cf**7*(h*czp*((2.573571875e4_dp*a(1)+1.0504375e4_dp*a(3)+5.145e2_dp*a(5)+2.0_dp*a(7))*sc0-1.4e1_dp*sc6*dxi)-szp*sc7/dx)
sint(8) = lz*cf**7*(h*szp*((2.573571875e4_dp*a(1)+1.0504375e4_dp*a(3)+5.145e2_dp*a(5)+2.0_dp*a(7))*sc0-1.4e1_dp*sc6*dxi)+czp*sc7/dx)
 
!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_CINT   
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_PPQQ (rzet, pp, qq)
!----------------------------------------------------------------------------------------------------------------------------------!

!~polynomials (p,q) featuring the asymptotic expansions of bessel functions    

USE GFLA 

complex(kind=dp), intent(in ) :: rzet  
complex(kind=dp), intent(out) :: pp(0:2), qq(0:2) 
integer(kind=k9)              :: k         
complex(kind=dp)              :: rtzs, coef,  psum(0:2), qsum(0:2), trmp(0:2), trmq(0:2)   
integer(kind=k9), parameter   :: n(0:2)=(/0,4,16/)  
real   (kind=dp), parameter   :: phia(0:2)=(/0.785398163397448_dp,2.356194490192345_dp,3.926990816987241_dp/) 
 
rtzs = 64.0_dp*rzet*rzet 
coef = sqrt(0.6366197723675814_dp/rzet)
trmp = 1.0_dp 
trmq = (n-1)/(8.0_dp*rzet)  
psum = trmp
qsum = trmq
loop1: do k = 1,ceiling(abs(rzet)) 
  trmp = -trmp*(n-(4*k-3)*(4*k-3))*(n-(4*k-1)*(4*k-1))/(rtzs*k*(4*k-2))
  trmq = -trmq*(n-(4*k-1)*(4*k-1))*(n-(4*k+1)*(4*k+1))/(rtzs*k*(4*k+2)) 
  psum =  psum + trmp
  qsum =  qsum + trmq   
end do loop1

pp = coef*(psum*cos(phia)+qsum*sin(phia))
qq = coef*(psum*sin(phia)-qsum*cos(phia)) 

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_PPQQ   
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_RKER (s, z, kz, lz, fu, ft)    
!----------------------------------------------------------------------------------------------------------------------------------!

!~kernel (regular part) of the hankel inversion integrals 

USE GFLA 
USE GFLB  

integer(kind=k9), intent(in ) :: kz, lz  
real   (kind=dp), intent(in ) :: s, z  
complex(kind=dp), intent(out) :: fu(3,3,0:2), ft(3,8,0:2)   
integer(kind=k9)              :: i, j    
complex(kind=dp)              :: psterm(3,-1:1,4), shterm(3,-1:1,2), disigps(3,-1:1,6), disigsh(3,-1:1,3)      
complex(kind=dp)              :: expu(2,2), expd(2,2), exsu, exsd, zeto, alfo, beto, coef, epsi          

zeto = zeta(kz) 
alfo = alfa(kz,lo)         
beto = beta(kz,lo)         

!~check if the integrand has already been evaluated   

if (flag(lz)==1) then
  fu = usmp(lz,:,:,:)
  ft = tsmp(lz,:,:,:)
  return  
else 
  flag(lz) = 1 
end if    

!~evaluation of the exponential matrices   

call GFL_EXMAT (kz, lo, lo, z, z, expu, expd, exsu, exsd) 

loop2: do i = 1,3     !-loop over the source directions
  loop3: do j = -1,1  !-loop over the fourier components

!~~~~~integration constants times the exponential functions   

    if (((i==1.or.i==2).and.j==0).or.(i==3.and.(j==-1.or.j==1))) then 

        psterm(i,j,:) = 0.0_dp   !-zero vector, p-sv waves 
        shterm(i,j,:) = 0.0_dp   !-zero vector, sh waves  

    else  

      if (lo<ls) then  !-receiver in the layer above the source 

!---------p-sv waves, total contribution 

        psterm(i,j,1:2) = matmul(expd,matmul(wdps(kz,lo,:,:),sups(kz,i,j,:)))      
        psterm(i,j,3:4) = matmul(expu,matmul(wups(kz,lo,:,:),sups(kz,i,j,:)))       
 
!---------sh waves, total contribution

        shterm(i,j,1) = exsd*wdsh(kz,lo)*sush(kz,i,j)                            
        shterm(i,j,2) = exsu*wush(kz,lo)*sush(kz,i,j)                            

      elseif (lo>ls) then  !-receiver in the layer below the source

!---------p-sv waves, total contribution

        psterm(i,j,1:2) = matmul(expd,matmul(wdps(kz,lo,:,:),sdps(kz,i,j,:)))      
        psterm(i,j,3:4) = matmul(expu,matmul(wups(kz,lo,:,:),sdps(kz,i,j,:)))      

!---------sh waves, total contribution

        shterm(i,j,1) = exsd*wdsh(kz,lo)*sdsh(kz,i,j)                            
        shterm(i,j,2) = exsu*wush(kz,lo)*sdsh(kz,i,j)                            

      else  !-receiver in the layer containing the source 

!---------p-sv waves, part1

        psterm(i,j,1:2) = matmul(expd,mdps(kz,i,j,:)) 
        psterm(i,j,3:4) = matmul(expu,mups(kz,i,j,:)) 

!---------sh waves, part1 

        shterm(i,j,1) = exsd*mdsh(kz,i,j) 
        shterm(i,j,2) = exsu*mush(kz,i,j) 

        if (z>=s) then  !-receiver below the source

          if(kz<kb) then 

!-------------p-sv waves, part2, standard representation for moderate values of zeta  

            psterm(i,j,1) = psterm(i,j,1) + 0.5_dp/rho(ls)*exp(alfo*(s-z))*(zzz(i,j)+zeto/alfo*xmy(i,j))       
            psterm(i,j,2) = psterm(i,j,2) + 0.5_dp/rho(ls)*exp(beto*(s-z))*(zzz(i,j)*zeto/beto+xmy(i,j))     

          else 

!-------------p-sv waves, part2, modified representation for large values of zeta  

            coef = 0.5_dp/(rho(ls)*(1.0_dp+pv(ls))) * exp(beto*(s-z))
            epsi = (zzz(i,j)+zeto/alfo*xmy(i,j)) * (exp(diab(kz,ls)*(s-z))-1.0_dp)  

            psterm(i,j,1) = psterm(i,j,1) - coef*((zzz(i,j)/beto+xmy(i,j)/alfo)/dfab(kz,ls)+epsi*difa(kz,ls))   
            psterm(i,j,2) = psterm(i,j,2) - coef*((zzz(i,j)/beto-pv(ls)*xmy(i,j)/alfo)/dfab(kz,ls)-epsi*(alfo+pv(ls)*zeto))  

          end if  

!-------------sh waves, part2 

          shterm(i,j,1) = shterm(i,j,1) - 0.5_dp*ic*exp(beto*(s-z))*xpy(i,j)/(mu(ls)*zeto*beto)           

        else  !-receiver above the source

          if(kz<kb) then 

!-------------p-sv waves, part2, standard representation for moderate values of zeta  

            psterm(i,j,3) = psterm(i,j,3) - 0.5_dp/rho(ls)*exp(alfo*(z-s))*(zzz(i,j)-zeto/alfo*xmy(i,j))     
            psterm(i,j,4) = psterm(i,j,4) + 0.5_dp/rho(ls)*exp(beto*(z-s))*(zzz(i,j)*zeto/beto-xmy(i,j))     
         
          else 

!-------------p-sv waves, part2, modified representation for large values of zeta  

            coef = 0.5_dp/(rho(ls)*(1.0_dp+pv(ls))) * exp(beto*(z-s))
            epsi = (zzz(i,j)-zeto/alfo*xmy(i,j)) * (exp(diab(kz,ls)*(z-s))-1.0_dp)

            psterm(i,j,3) = psterm(i,j,3) + coef*((zzz(i,j)/beto-xmy(i,j)/alfo)/dfab(kz,ls)+epsi*difa(kz,ls))   
            psterm(i,j,4) = psterm(i,j,4) - coef*((zzz(i,j)/beto+pv(ls)*xmy(i,j)/alfo)/dfab(kz,ls)-epsi*(alfo+pv(ls)*zeto))  

          end if  

!-------------sh waves, part2 

          shterm(i,j,2) = shterm(i,j,2) - 0.5_dp*ic*exp(beto*(z-s))*xpy(i,j)/(mu(ls)*zeto*beto)           

        end if      

      end if  

    end if   

!~~~~~transformed displacements and stresses  

    disigps(i,j,:) = matmul(cmps(kz,lo,:,:),psterm(i,j,:)) 
    disigsh(i,j,:) = matmul(cmsh(kz,lo,:,:),shterm(i,j,:)) 

  end do loop3
end do loop2 

!~evaluation of the integrand (total kernel) 

!-displacements u_r 

epsi = 0.5_dp*ks0/mu0*zeta(kz)

usmp(lz,:,1,0) =  epsi * ((disigsh(:,-1,1)+disigps(:,-1,1))*expm + (disigsh(:,1,1)-disigps(:,1,1))*expp)       
usmp(lz,:,1,1) =  epsi * (2.0_dp*disigps(:,0,1)) 
usmp(lz,:,1,2) =  epsi * ((disigsh(:,-1,1)-disigps(:,-1,1))*expm + (disigsh(:,1,1)+disigps(:,1,1))*expp)       

!-displacements u_t 

epsi = -0.5_dp*ic*ks0/mu0*zeta(kz)

usmp(lz,:,2,0) =  epsi * ((disigsh(:,-1,1)+disigps(:,-1,1))*expm - (disigsh(:,1,1)-disigps(:,1,1))*expp)       
usmp(lz,:,2,1) =  epsi * (2.0_dp*disigsh(:,0,1)) 
usmp(lz,:,2,2) = -epsi * ((disigsh(:,-1,1)-disigps(:,-1,1))*expm - (disigsh(:,1,1)+disigps(:,1,1))*expp)       

!-displacements u_z 

epsi = ks0/mu0*zeta(kz)

usmp(lz,:,3,0) =  epsi * (disigps(:,0,2)) 
usmp(lz,:,3,1) =  epsi * (disigps(:,1,2)*expp-disigps(:,-1,2)*expm) 
usmp(lz,:,3,2) =  0.0_dp  

!-stresses t_rr 

epsi = ks0*ks0*zeta(kz)

tsmp(lz,:,1,0) =  epsi * (disigps(:,0,6)) 
tsmp(lz,:,1,1) =  epsi * (disigps(:,1,6)*expp-disigps(:,-1,6)*expm) 
tsmp(lz,:,1,2) =  0.0_dp  

!-stresses t_tt 

tsmp(lz,:,2,0) =  epsi * (disigps(:,0,5)) 
tsmp(lz,:,2,1) =  epsi * (disigps(:,1,5)*expp-disigps(:,-1,5)*expm) 
tsmp(lz,:,2,2) =  0.0_dp  

!-stresses t_zz 

tsmp(lz,:,3,0) =  epsi * (disigps(:,0,4)) 
tsmp(lz,:,3,1) =  epsi * (disigps(:,1,4)*expp-disigps(:,-1,4)*expm) 
tsmp(lz,:,3,2) =  0.0_dp  

!-stresses t_zt 

epsi = -0.5_dp*ic*ks0*ks0*zeta(kz)

tsmp(lz,:,4,0) =  epsi * ((disigsh(:,-1,2)+disigps(:,-1,3))*expm - (disigsh(:,1,2)-disigps(:,1,3))*expp)       
tsmp(lz,:,4,1) =  epsi * (2.0_dp*disigsh(:,0,2)) 
tsmp(lz,:,4,2) = -epsi * ((disigsh(:,-1,2)-disigps(:,-1,3))*expm - (disigsh(:,1,2)+disigps(:,1,3))*expp)       

!-stresses t_zr 

epsi = 0.5_dp*ks0*ks0*zeta(kz)

tsmp(lz,:,5,0) =  epsi * ((disigsh(:,-1,2)+disigps(:,-1,3))*expm + (disigsh(:,1,2)-disigps(:,1,3))*expp)       
tsmp(lz,:,5,1) =  epsi * (2.0_dp*disigps(:,0,3)) 
tsmp(lz,:,5,2) =  epsi * ((disigsh(:,-1,2)-disigps(:,-1,3))*expm + (disigsh(:,1,2)+disigps(:,1,3))*expp)       

!-stresses t_rt 

epsi = ks0*ks0*zeta(kz)

tsmp(lz,:,6,0) =  epsi * (disigsh(:,0,3)) 
tsmp(lz,:,6,1) =  epsi * (disigsh(:,1,3)*expp-disigsh(:,-1,3)*expm) 
tsmp(lz,:,6,2) =  0.0_dp    

!-auxiliary term {2mu0*ks0*mu(lo) * (u_r + i*(u_t(1)*expp - u_t(-1)*expm))} 

epsi = 2.0_dp*mu0*ks0*mu(lo) 
coef = 0.5_dp*ks0/mu0*zeta(kz)

tsmp(lz,:,7,0) = 0.0_dp      
tsmp(lz,:,7,1) = epsi * (usmp(lz,:,1,1))   
tsmp(lz,:,7,2) = epsi * (usmp(lz,:,1,2) + coef*((disigsh(:,1,1)+disigps(:,1,1))*expp + (disigsh(:,-1,1)-disigps(:,-1,1))*expm)) 

!-auxiliary term {2mu0*ks0*mu(lo) * (u_t - i*(u_r(1)*expp - u_r(-1)*expm))} 

coef = 0.5_dp*ic*ks0/mu0*zeta(kz)

tsmp(lz,:,8,0) = 0.0_dp      
tsmp(lz,:,8,1) = epsi * (usmp(lz,:,2,1))   
tsmp(lz,:,8,2) = epsi * (usmp(lz,:,2,2) - coef*((disigsh(:,1,1)+disigps(:,1,1))*expp - (disigsh(:,-1,1)-disigps(:,-1,1))*expm))   

!~singularity teratment (arrays fu and ft are used temporary to store the singular parts) 
       
!-singular part of the kernel  

call GFL_SKER (mu0, ks0, zeta(kz), fu, ft)     

!-regular parts of the kernel  

usmp(lz,:,:,:) = usmp(lz,:,:,:) - magfac*fu   
tsmp(lz,:,:,:) = tsmp(lz,:,:,:) - magfac*ft       
     
!~assign the result to output variables 

fu = usmp(lz,:,:,:)
ft = tsmp(lz,:,:,:) 

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_RKER   
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_SKER (mu0, ks0, zet, using, tsing)    
!----------------------------------------------------------------------------------------------------------------------------------!
 
!~kernel (singular part) of the hankel inversion integrals 

USE GFLB  

real   (kind=dp), intent(in ) :: mu0, ks0     
complex(kind=dp), intent(in ) :: zet    
complex(kind=dp), intent(out) :: using(3,3,0:2), tsing(3,8,0:2)     
complex(kind=dp)              :: ome1, ome2, gam1, gam2, gam3, dom1, dom2, dga1, dga2, dga3, muj, laj, factor, coef     

using = 0.0_dp
tsing = 0.0_dp 

!~evaluation of the static bimaterial influence functions

call GFL_OMEGA (zet, ome1, ome2, gam1, gam2, gam3, dom1, dom2, dga1, dga2, dga3)  

!~choice of the appropriate material parameters 

if (zst<0.0_dp.or.(zst==0.0_dp.and.sing<0)) then  
  muj = mst1
  laj = (2.0_dp*vst1*mst1)/(1.0_dp-2.0_dp*vst1)
else 
  muj = mst2
  laj = (2.0_dp*vst2*mst2)/(1.0_dp-2.0_dp*vst2)
end if  

!~kernel for the displacement green's functions  

factor = 0.07957747154594768_dp*ks0/(mu0*mst2)*zet   

using(:,1,0) =  factor*statc*(gam2+gam1)
using(:,1,1) = -factor*stat0*(2.0_dp*gam3) 
using(:,1,2) =  factor*statc*(gam2-gam1)

using(:,2,0) = -factor*stats*(gam2+gam1)
using(:,2,2) =  factor*stats*(gam2-gam1)

using(:,3,0) =  factor*stat0*(2.0_dp*ome2)
using(:,3,1) =  factor*statc*(2.0_dp*ome1)  

!~kernel for the stress green's functions  

factor = 0.1591549430918953_dp*ks0*ks0/mst2*zet      

tsing(:,1,0) =  factor*stat0*(laj*dom2-(laj+2.0_dp*muj)*zet*gam3)
tsing(:,1,1) =  factor*statc*(laj*dom1-(laj+2.0_dp*muj)*zet*gam1) 

tsing(:,2,0) =  factor*stat0*laj*(dom2-zet*gam3)
tsing(:,2,1) =  factor*statc*laj*(dom1-zet*gam1) 

tsing(:,3,0) =  factor*stat0*((laj+2.0_dp*muj)*dom2-laj*zet*gam3)
tsing(:,3,1) =  factor*statc*((laj+2.0_dp*muj)*dom1-laj*zet*gam1) 

tsing(:,4,0) = -factor*stats*muj*0.5_dp*(dga2+dga1+zet*ome1)
tsing(:,4,2) =  factor*stats*muj*0.5_dp*(dga2-dga1-zet*ome1)

tsing(:,5,0) =  factor*statc*muj*0.5_dp*(dga2+dga1+zet*ome1)
tsing(:,5,1) = -factor*stat0*muj*(dga3+zet*ome2) 
tsing(:,5,2) =  factor*statc*muj*0.5_dp*(dga2-dga1-zet*ome1)   

tsing(:,6,1) =  factor*stats*muj*(zet*gam2)  

!~auxiliary kernel for the stress green's functions  

coef   = 2.0_dp*ks0*mu0*muj 
factor = 0.1591549430918953_dp*ks0/(mu0*mst2)*zet   
   
tsing(:,7,1) =  coef*using(:,1,1) 
tsing(:,7,2) =  coef*factor*statc*(gam2-gam1)       
     
tsing(:,8,2) =  coef*factor*stats*(gam2-gam1)      

if (sing<0) then   

  using(2:3,:,:) = -using(2:3,:,:)  !-reversing the x,z source directions 
  tsing(2:3,:,:) = -tsing(2:3,:,:)  !-reversing the x,z source directions 
    
  using(:,2,:) = -using(:,2,:)  !-conversion to the global coord. system  
  using(:,3,:) = -using(:,3,:)  !-conversion to the global coord. system
  tsing(:,5,:) = -tsing(:,5,:)  !-conversion to the global coord. system 
  tsing(:,6,:) = -tsing(:,6,:)  !-conversion to the global coord. system 
  tsing(:,8,:) = -tsing(:,8,:)  !-conversion to the global coord. system 
  
end if     

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_SKER   
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_OMEGA (zet, ome1, ome2, gam1, gam2, gam3, dom1, dom2, dga1, dga2, dga3)    
!----------------------------------------------------------------------------------------------------------------------------------!

!~static bimaterial influence functions: omega_k, gamma_j and their derivatives (k=1,2, j=1,2,3)  

USE GFLB  

complex(kind=dp), intent(in ) :: zet  
complex(kind=dp), intent(out) :: ome1, ome2, gam1, gam2, gam3, dom1, dom2, dga1, dga2, dga3  
real   (kind=dp)              :: szs   
complex(kind=dp)              :: exd1, exd2, mm1, mm2, vv1, vv2, ww1, ww2, cfd1, cfd2, trm0, trm1, trm2 

!~auxiliary parameters 

vv1  = 3.0_dp-4.0_dp*vst1
vv2  = 3.0_dp-4.0_dp*vst2
ww1  = 1.0_dp-2.0_dp*vst1
ww2  = 1.0_dp-2.0_dp*vst2
mm1  = (mst1+vv1*mst2)
mm2  = (mst2+vv2*mst1)   

!~influence functions and their derivatives

if (zst<0.0_dp.or.(zst==0.0_dp.and.sing<0)) then  !-upper layer (region1, zst<0)

  exd1 =  exp(-zet*abs(zst-sst)) 
  cfd1 =  mst2*exd1/(2.0_dp*mm1*mm2)  
  
  trm1 =  zet*(zst*mm2-sst*mm1) 
  trm0 =  mst1*ww1*vv2-mst2*ww2*vv1

  ome1 =  cfd1*(trm1-trm0)/zet 
  gam3 = -cfd1*(trm1+trm0)/zet 
  dom1 =  cfd1*(trm1-trm0+mm2) 
  dga3 = -cfd1*(trm1+trm0+mm2) 

  trm0 =  mst1*(2.0_dp-2.0_dp*vst1)*vv2+mst2*(2.0_dp-2.0_dp*vst2)*vv1

  ome2 = -cfd1*(trm1-trm0)/zet 
  gam1 =  cfd1*(trm1+trm0)/zet 
  dom2 = -cfd1*(trm1-trm0+mm2) 
  dga1 =  cfd1*(trm1+trm0+mm2) 

  gam2 =  mst2*exd1/(zet*(mst1+mst2))
  dga2 =  mst2*exd1/(mst1+mst2)

else  !-lower layer (region2/3, zst>0) 

  exd1 =  exp(-zet*abs(zst-sst)) 
  exd2 =  exp(-zet*(zst+sst)) 
  cfd1 =  exd1/(8.0_dp*(1.0_dp-vst2)*mm1*mm2) 
  cfd2 =  exd2/(8.0_dp*(1.0_dp-vst2)*mm1*mm2) 

  trm2 =  2.0_dp*zet*zet*(mst1-mst2)*mm1*zst*sst 
  trm1 =  zet*(mst1-mst2)*vv2*mm1*(zst-sst) 
  trm0 =  4.0_dp*mst2*(1.0_dp-vst2)*(mst1*ww1*vv2-mst2*ww2*vv1)  

  ome1 =  (cfd1*zet*(zst-sst)*mm1*mm2 + cfd2*(trm2-trm1-trm0))/zet  
  gam3 =  (cfd1*zet*(sst-zst)*mm1*mm2 + cfd2*(trm2+trm1-trm0))/zet  

  trm1 =  zet*(mst1-mst2)*mm1*(zst*vv2-sst*(1.0_dp-4.0_dp*vst2)) 
  trm0 =  mst1**2*vv2+mst2**2*vv1*(1.0_dp-8.0_dp*vst2+8.0_dp*vst2**2)-mst1*mst2*2.0*ww1*vv2*ww2   

  dom1 = -cfd1*(zet*abs(zst-sst)-1.0_dp)*mm1*mm2 - cfd2*(trm2-trm1+trm0)  

  trm1 =  zet*(mst1-mst2)*mm1*(zst*vv2-sst*(5.0_dp-4.0_dp*vst2)) 
  trm0 =  mst1**2*vv2-mst2**2*vv1*(7.0_dp-16.0_dp*vst2+8.0_dp*vst2**2)+mst1*mst2*2.0*ww1*vv2*(2.0_dp+ww2)  

  dga3 =  cfd1*(zet*abs(zst-sst)-1.0_dp)*mm1*mm2 - cfd2*(trm2+trm1-trm0)  

  trm1 =  zet*(mst1-mst2)*vv2*mm1*(zst+sst) 
  trm0 =  mst1**2*vv2**2-mst2**2*vv1*(5.0_dp-12.0_dp*vst2+8.0_dp*vst2**2)+mst1*mst2*2.0*ww1*vv2*ww2   

  ome2 = (cfd1*(vv2+zet*abs(zst-sst))*mm1*mm2 - cfd2*(trm2+trm1+trm0))/zet  
  gam1 = (cfd1*(vv2-zet*abs(zst-sst))*mm1*mm2 - cfd2*(trm2-trm1+trm0))/zet  

  trm1 =  zet*(mst1-mst2)*mm1*(zst*vv2+sst*(1.0_dp-4.0_dp*vst2)) 
  trm0 =  2.0_dp*ww2*(mst1**2*vv2-mst2**2*vv1*ww2)-4.0_dp*mst1*mst2*vst2*vv2*ww1  

  if (sing>0) then 
    szs =  sign(1.0_dp,zst-sst) 
  else 
    szs = -sign(1.0_dp,sst-zst) 
  end if  

  dom2 = -cfd1*(2.0_dp-4.0_dp*vst2+zet*abs(zst-sst))*szs*mm1*mm2 + cfd2*(trm2+trm1+trm0) 

  trm1 =  zet*(mst1-mst2)*mm1*(zst*vv2+sst*(5.0_dp-4.0_dp*vst2)) 
  trm0 =  4.0_dp*(1.0_dp-vst2)*(mst1**2*vv2-mst2**2*vv1*(2.0_dp-2.0_dp*vst2)+mst1*mst2*vv2*ww1)

  dga1 = -cfd1*(4.0_dp-4.0_dp*vst2-zet*abs(zst-sst))*szs*mm1*mm2 + cfd2*(trm2-trm1+trm0) 

  gam2 = (exd1*0.5_dp - exd2*0.5_dp*(mst1-mst2)/(mst1+mst2))/zet 
  dga2 = -exd1*0.5_dp*szs + exd2*0.5_dp*(mst1-mst2)/(mst1+mst2) 

end if   

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_OMEGA     
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_BIMAT (sing, mu1, v1, mu2, v2, s, r, t, z, mu0, ks0, ushat, tshat)
!----------------------------------------------------------------------------------------------------------------------------------!

!~author  : bojan guzina, university of colorado at boulder, april 20, 1994; version 1.0  
!~purpose : static point-load green's functions for the joined half-spaces 

!~input  : sing       => integer defining whether the source point is below or above the relevant interface
!~input  : mu1, v1    => dimensionless shear_modulus/poisson's_ratio of the upper medium (z<0)
!~input  : mu2, v2    => dimensionless shear_modulus/poisson's_ratio of the lower medium (z>0)
!~input  : (0, 0, s)  => dimensionless cylindrical coordinates of the source point  
!~input  : (r, t, z)  => dimensionless cylindrical coordinates of the observation point  
!~input  : mu0, ks0   => reference shear modulus [kPa] and wave number [0.0316228/m] used for normalization   
!~output : ushat(i,j) => displ. green's functions due to the (i=1=>x,i=2=>y,i=3=>z) unit point-load [m/kN]  
!~output : tshat(i,j) => stress green's functions due to the (i=1=>x,i=2=>y,i=3=>z) unit point-load [1/m^2] 

USE GFLR 

integer(kind=k9), intent(in ) :: sing 
real   (kind=dp), intent(in ) :: s, r, t, z, mu0, ks0 
complex(kind=dp), intent(in ) :: mu1, v1, mu2, v2 
complex(kind=dp), intent(out) :: ushat(3,3), tshat(3,6)
complex(kind=dp)              :: ustat(3,3), tstat(3,8)
complex(kind=dp)              :: rcof


!~evaluation of the green's functions (in terms of dimensional variables)      

call GFB_DRIVER (sing, mu0*mu1, v1, mu0*mu2, v2, s/ks0, r/ks0, t, z/ks0, ustat, tstat)   

!~conversion from local to global cylindrical coordinates       

if (sing<0) then 

  ustat(2:3,:) = -ustat(2:3,:)  !-reversing the source directions, u^2_*, u^3_*
  tstat(2:3,:) = -tstat(2:3,:)  !-reversing the source directions, tau^2_*, tau^3_* 

  ustat(:,2) = -ustat(:,2)  !-transofrmation from local to global coord. system, u^*_t
  ustat(:,3) = -ustat(:,3)  !-transofrmation from local to global coord. system, u^*_z
  tstat(:,5) = -tstat(:,5)  !-transofrmation from local to global coord. system, tau^*_{zr}
  tstat(:,6) = -tstat(:,6)  !-transofrmation from local to global coord. system, tau^*_{rt}
  tstat(:,8) = -tstat(:,8)  !-transofrmation from local to global coord. system, tau^*_{rt} (auxiliary term) 

end if    

!~conversion to the format used by dynamic solution  

if (rflag==2) then;  rcof=ks0*rfac;  else; rcof=rfac;  end if  

ushat(:,:) = ustat(:,:) 
tshat(:,1) = tstat(:,1) - rcof*tstat(:,7)   
tshat(:,2) = tstat(:,2) + rcof*tstat(:,7)  
tshat(:,3) = tstat(:,3)
tshat(:,4) = tstat(:,4)
tshat(:,5) = tstat(:,5) 
tshat(:,6) = tstat(:,6) - rcof*tstat(:,8)      

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_BIMAT 
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFB_DRIVER (sing, mu1, v1, mu2, v2, s, r, t, z, ustat, tstat)  
!----------------------------------------------------------------------------------------------------------------------------------!

!~static point-load green's functions for the joined half-spaces 

USE GFLR 

!~input  : sing       => integer defining whether the source point is below or above the relevant interface
!~input  : mu1, v1    => dimensional shear_modulus/poisson's_ratio of the upper medium (z<0)
!~input  : mu2, v2    => dimensional shear_modulus/poisson's_ratio of the lower medium (z>0)
!~input  : (0, 0, s)  => dimensional cylindrical coordinates of the source point  
!~input  : (r, t, z)  => dimensional cylindrical coordinates of the observation point  
!~output : ustat(i,j) => displ. green's functions due to the (i=1=>x,i=2=>y,i=3=>z) point-load
!~output : tstat(i,j) => stress green's functions due to the (i=1=>x,i=2=>y,i=3=>z) point-load 

!~notation : ur =ustat(i,1), ut =ustat(i,2), uz =ustat(i,3), where subscript "t" stands for "t"
!~notation : trr=tstat(i,1), ttt=tstat(i,2), tzz=tstat(i,3), tzt=tstat(i,4), tzr=tstat(i,5), trt=tstat(i,6)

integer(kind=k9), intent(in ) :: sing 
real   (kind=dp), intent(in ) :: s, r, t, z   
complex(kind=dp), intent(in ) :: mu1, v1, mu2, v2 
complex(kind=dp), intent(out) :: ustat(3,3), tstat(3,8)
real   (kind=dp)              :: ccc, sss      
complex(kind=dp)              :: ogi(5,3,0:2), muj, laj, lmj, fac 

!~choice of the appropriate material parameters 

if (z<0.0_dp.or.(z==0.0_dp.and.sing<0)) then  
  muj = mu1
  laj = (2.0_dp*v1*mu1)/(1.0_dp-2.0_dp*v1)
else 
  muj = mu2
  laj = (2.0_dp*v2*mu2)/(1.0_dp-2.0_dp*v2)
end if  

!~functions of the angular coordinate   

ccc = cos(t);  sss = sin(t)   

!~auxiliary parameters

fac = 0.1591549430918953_dp/mu2 
lmj = laj+2.0_dp*muj 

!~normalized hankel inversion integrals ogi(i,j,k)  (ogi(i,j,1) and ogi(i,j,2) are divided by r)   

call OGINT1 (sing, mu1, v1, mu2, v2, s, r, z, ogi)

!~green's functions that require normalization for small r 

if (rflag==2)  ogi(:,:,1:2) = r * ogi(:,:,1:2)  !-undo the normalization for moderate to large r   

tstat(1,7) =  (2.0*muj*fac) * ccc*(ogi(4,1,2)-ogi(3,1,2)) 
tstat(1,8) =  (2.0*muj*fac) * sss*(ogi(4,1,2)-ogi(3,1,2))
 
tstat(2,7) =  (2.0*muj*fac) * sss*(ogi(4,1,2)-ogi(3,1,2)) 
tstat(2,8) = -(2.0*muj*fac) * ccc*(ogi(4,1,2)-ogi(3,1,2))      

tstat(3,7) = -(2.0*muj*fac) * (ogi(5,1,1))
tstat(3,8) =  (    0.0    )  

!~green's functions that do not require any normalization 

if (rflag==1)  ogi(:,:,1:2) = r * ogi(:,:,1:2)  !-undo the normalization for small values of r     

ustat(1,1) =  (0.5*fac) * ccc*(ogi(4,1,0)+ogi(3,1,0)+ogi(4,1,2)-ogi(3,1,2))
ustat(1,2) = -(0.5*fac) * sss*(ogi(4,1,0)+ogi(3,1,0)-ogi(4,1,2)+ogi(3,1,2))
ustat(1,3) =  (    fac) * ccc*(ogi(1,1,1))

tstat(1,1) =  (        fac) * ccc*(laj*ogi(1,3,1)-lmj*ogi(3,2,1))
tstat(1,2) =  (    laj*fac) * ccc*(ogi(1,3,1)-ogi(3,2,1)) 
tstat(1,3) =  (        fac) * ccc*(lmj*ogi(1,3,1)-laj*ogi(3,2,1))
tstat(1,4) = -(0.5*muj*fac) * sss*(ogi(4,3,0)+ogi(3,3,0)+ogi(1,2,0)-ogi(4,3,2)+ogi(3,3,2)+ogi(1,2,2))
tstat(1,5) =  (0.5*muj*fac) * ccc*(ogi(4,3,0)+ogi(3,3,0)+ogi(1,2,0)+ogi(4,3,2)-ogi(3,3,2)-ogi(1,2,2))
tstat(1,6) =  (    muj*fac) * sss*(ogi(4,2,1))   

ustat(2,1) =  (0.5*fac) * sss*(ogi(4,1,0)+ogi(3,1,0)+ogi(4,1,2)-ogi(3,1,2))
ustat(2,2) =  (0.5*fac) * ccc*(ogi(4,1,0)+ogi(3,1,0)-ogi(4,1,2)+ogi(3,1,2))
ustat(2,3) =  (    fac) * sss*(ogi(1,1,1))

tstat(2,1) =  (        fac) * sss*(laj*ogi(1,3,1)-lmj*ogi(3,2,1))
tstat(2,2) =  (    laj*fac) * sss*(ogi(1,3,1)-ogi(3,2,1)) 
tstat(2,3) =  (        fac) * sss*(lmj*ogi(1,3,1)-laj*ogi(3,2,1))
tstat(2,4) =  (0.5*muj*fac) * ccc*(ogi(4,3,0)+ogi(3,3,0)+ogi(1,2,0)-ogi(4,3,2)+ogi(3,3,2)+ogi(1,2,2))
tstat(2,5) =  (0.5*muj*fac) * sss*(ogi(4,3,0)+ogi(3,3,0)+ogi(1,2,0)+ogi(4,3,2)-ogi(3,3,2)-ogi(1,2,2))
tstat(2,6) = -(    muj*fac) * ccc*(ogi(4,2,1))

ustat(3,1) = -(fac) * (ogi(5,1,1))
ustat(3,2) =  (0.0)   
ustat(3,3) =  (fac) * (ogi(2,1,0)) 

tstat(3,1) =  (    fac) * (laj*ogi(2,3,0)-lmj*ogi(5,2,0))
tstat(3,2) =  (laj*fac) * (    ogi(2,3,0)-    ogi(5,2,0))
tstat(3,3) =  (    fac) * (lmj*ogi(2,3,0)-laj*ogi(5,2,0))
tstat(3,4) =  (  0.0  )
tstat(3,5) = -(muj*fac) * (    ogi(5,3,1)+    ogi(2,2,1))
tstat(3,6) =  (  0.0  )

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFB_DRIVER
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE OGINT1 (sing, mu1, v1, mu2, v2, s, r, z, ogi)
!----------------------------------------------------------------------------------------------------------------------------------!

!~evaluation of the generalized integrals: ogi(i,j,k) = Int (omegamma(xi) xi J_k(xi*r) dxi), xi=0,infinity (routine #1) 
!~note: integrals ogi(i,j,1) and ogi(i,j,2) are normalized (divided) by r   

!~notation : i=1 => omegamma(xi)=xxx*omega_1(xi), ... i=5 => omegamma(xi)=xxx*gamma_3(xi)
!~notation : j=1 => xxx=(1.0), j=2 => xxx=(xi), j=3 => xxx=(d/dz) 
!~notation : j_k  is the bessel function of integer order k 

integer(kind=k9), intent(in ) :: sing 
real   (kind=dp), intent(in ) :: s, r, z 
complex(kind=dp), intent(in ) :: mu1, v1, mu2, v2 
complex(kind=dp), intent(out) :: ogi(5,3,0:2) 
integer(kind=k9)              :: k      
real   (kind=dp)              :: j1(0:2,0:3), j2(0:2,0:3), d1, d2, d3, szs 
complex(kind=dp)              :: mx1, mx2, o1v1, o1v2, o2v1, o2v2, t4v1, t4v2, t2v2, o4v2, f4v2, c(0:9) 

!~definition of the auxiliary variables 

d1 = abs(z-s)
d2 =    (z+s) 
d3 =    (z-s)   

szs = sign(1.0_dp,z-s) 
    
o1v1 = 1.0-1.0*v1
o1v2 = 1.0-1.0*v2
o2v1 = 1.0-2.0*v1
o2v2 = 1.0-2.0*v2
t4v1 = 3.0-4.0*v1
t4v2 = 3.0-4.0*v2
t2v2 = 3.0-2.0*v2
o4v2 = 1.0-4.0*v2
f4v2 = 5.0-4.0*v2
mx1  = mu1+t4v1*mu2
mx2  = mu2+t4v2*mu1

!~initialize the array of auxiliary coefficients

c    = 0.0_dp 
c(0) = mu1*o2v1*t4v2-mu2*o2v2*t4v1  

!~selection of proper omegas and gammas depending on the source/receiver location 

if (z<0.0_dp.or.(z==0.0_dp.and.sing<0)) then  !-receiver and source are in different media

  call BASINT1 (r, d1, j1)

  c(1) = mu2/(2.0*mx1*mx2)
  c(2) = mu2/(mu1+mu2)
  c(3) = mu1*2.0*o1v1*t4v2+mu2*2.0*o1v2*t4v1 

else  !-receiver and source are in the same medium 

  call BASINT1 (r, d1, j1)
  call BASINT1 (r, d2, j2)

  c(1) = 1.0/(8.0*o1v2*mx1*mx2)
  c(2) = 1.0/(2.0*(mu1+mu2))
  c(3) = 2.0*mx1*(mu1-mu2)*z*s
  c(4) = mx1*(mu1-mu2)
  c(5) = mu1**2*t4v2**2-mu2**2*t4v1*(5.0-12.0*v2+8.0*v2**2)+2.0*mu1*mu2*o2v1*t4v2*o2v2
  c(6) = mu1**2*t4v2+mu2**2*t4v1*(1.0-8.0*v2+8.0*v2**2)-2.0*mu1*mu2*o2v1*t4v2*o2v2
  c(7) = 2.0*mu1**2*t4v2*o2v2-2.0*mu2**2*t4v1*o2v2**2-4.0*mu1*mu2*v2*t4v2*o2v1
  c(8) = 4.0*o1v2*(mu1**2*t4v2-2.0*mu2**2*t4v1*o1v2+mu1*mu2*t4v2*o2v1)
  c(9) = mu1**2*t4v2-mu2**2*t4v1*(7.0-16.0*v2+8.0*v2**2)+2.0*mu1*mu2*o2v1*t4v2*t2v2

end if 

loop2: do k = 0,2
  call  OGINT2 (sing, k, mu1, mu2, mx1, mx2, o1v2, o2v2, t4v2, o4v2, f4v2, c, s, z, d1, d2, d3, szs, j1, j2, ogi(:,:,k))
end do loop2

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE OGINT1
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE OGINT2 (sing, k, mu1, mu2, mx1, mx2, o1v2, o2v2, t4v2, o4v2, f4v2, c, s, z, d1, d2, d3, szs, j1, j2, ogi)
!----------------------------------------------------------------------------------------------------------------------------------!

!~evaluation of the generalized integrals: ogi(i,j,k) = Int (omegamma(xi) xi J_k(xi*r) dxi), xi=0,infinity (routine #2) 
!~note: integrals ogi(i,j,1) and ogi(i,j,2) are normalized (divided) by r   

integer(kind=k9), intent(in ) :: sing, k 
real   (kind=dp), intent(in ) :: z, s, d1, d2, d3, szs, j1(0:2,0:3), j2(0:2,0:3) 
complex(kind=dp), intent(in ) :: mu1, mu2, mx1, mx2, o1v2, o2v2, t4v2, o4v2, f4v2, c(0:9) 
complex(kind=dp), intent(out) :: ogi(5,3)  
integer(kind=k9)              :: stab 
complex(kind=dp)              :: wx2, w1v2

!~stabilization procedure for ogi(3,3,2) and ogi(4,3,2) which appear in the form of a difference 
!~terms containing j1(2,1) may become very large for r->0 and should be subtracted prior to multiplication by j1(2,1)    

if (k==2) then ; stab = 1 ; else ; stab = 0 ; end if 

wx2  = mx2  - stab*c(2)/c(1)
w1v2 = o1v2 - stab*c(2)/c(1)*(mu1+mu2)/(4.0*mx1*mx2)

if (z<0.0_dp.or.(z==0.0_dp.and.sing<0)) then  !-observation point and the source are in different media

  ogi(1,1) =  c(1) * ((z*mx2-s*mx1)*j1(k,1) - c(0)*j1(k,0))
  ogi(2,1) = -c(1) * ((z*mx2-s*mx1)*j1(k,1) - c(3)*j1(k,0))
  ogi(3,1) =  c(1) * ((z*mx2-s*mx1)*j1(k,1) + c(3)*j1(k,0))
  ogi(4,1) =  c(2) * j1(k,0)
  ogi(5,1) = -c(1) * ((z*mx2-s*mx1)*j1(k,1) + c(0)*j1(k,0))
  ogi(1,2) =  c(1) * ((z*mx2-s*mx1)*j1(k,2) - c(0)*j1(k,1))
  ogi(2,2) = -c(1) * ((z*mx2-s*mx1)*j1(k,2) - c(3)*j1(k,1))
  ogi(3,2) =  c(1) * ((z*mx2-s*mx1)*j1(k,2) + c(3)*j1(k,1))
  ogi(4,2) =  c(2) * j1(k,1)
  ogi(5,2) = -c(1) * ((z*mx2-s*mx1)*j1(k,2) + c(0)*j1(k,1))
  ogi(1,3) =  c(1) * ((z*mx2-s*mx1)*j1(k,2) - (c(0)-mx2)*j1(k,1))
  ogi(2,3) = -c(1) * ((z*mx2-s*mx1)*j1(k,2) - (c(3)-mx2)*j1(k,1))
  ogi(3,3) =  c(1) * ((z*mx2-s*mx1)*j1(k,2) + (c(3)+wx2)*j1(k,1))
  ogi(4,3) =  c(2) * (1-stab)*j1(k,1) 
  ogi(5,3) = -c(1) * ((z*mx2-s*mx1)*j1(k,2) + (c(0)+mx2)*j1(k,1))

else  !-observation point and the source are in the same medium 

  ogi(1,1) =  c(1) * ( mx1*mx2*d3*j1(k,1) + c(3)*j2(k,2) - c(4)*t4v2*d3*j2(k,1) - 4.0*mu2*o1v2*c(0)*j2(k,0))
  ogi(2,1) =  c(1) * ( mx1*mx2*(d1*j1(k,1)+t4v2*j1(k,0)) - c(3)*j2(k,2) - c(4)*t4v2*d2*j2(k,1) - c(5)*j2(k,0))
  ogi(3,1) =  c(1) * (-mx1*mx2*(d1*j1(k,1)-t4v2*j1(k,0)) - c(3)*j2(k,2) + c(4)*t4v2*d2*j2(k,1) - c(5)*j2(k,0))
  ogi(4,1) =  c(2) * ((mu1+mu2)*j1(k,0)-(mu1-mu2)*j2(k,0)) 
  ogi(5,1) =  c(1) * (-mx1*mx2*d3*j1(k,1) + c(3)*j2(k,2) + c(4)*t4v2*d3*j2(k,1) - 4.0*mu2*o1v2*c(0)*j2(k,0)) 
  ogi(1,2) =  c(1) * ( mx1*mx2*d3*j1(k,2) + c(3)*j2(k,3) - c(4)*t4v2*d3*j2(k,2) - 4.0*mu2*o1v2*c(0)*j2(k,1))
  ogi(2,2) =  c(1) * ( mx1*mx2*(d1*j1(k,2)+t4v2*j1(k,1)) - c(3)*j2(k,3) - c(4)*t4v2*d2*j2(k,2) - c(5)*j2(k,1))
  ogi(3,2) =  c(1) * (-mx1*mx2*(d1*j1(k,2)-t4v2*j1(k,1)) - c(3)*j2(k,3) + c(4)*t4v2*d2*j2(k,2) - c(5)*j2(k,1))           
  ogi(4,2) =  c(2) * ((mu1+mu2)*j1(k,1)-(mu1-mu2)*j2(k,1))      
  ogi(5,2) =  c(1) * (-mx1*mx2*d3*j1(k,2) + c(3)*j2(k,3) + c(4)*t4v2*d3*j2(k,2) - 4.0*mu2*o1v2*c(0)*j2(k,1))  
  ogi(1,3) =  c(1) * (-mx1*mx2*(d1*j1(k,2)-j1(k,1)) - c(3)*j2(k,3) + c(4)*(z*t4v2-s*o4v2)*j2(k,2) - c(6)*j2(k,1))
  ogi(2,3) =  c(1) * (-mx1*mx2*szs*(d1*j1(k,2)+2.0*o2v2*j1(k,1)) + c(3)*j2(k,3) + c(4)*(z*t4v2+s*o4v2)*j2(k,2) + c(7)*j2(k,1))
  ogi(3,3) =  c(1) * ( mx1*mx2*szs*(d1*j1(k,2)-4.0*w1v2*j1(k,1)) + c(3)*j2(k,3) - c(4)*(z*t4v2+s*f4v2)*j2(k,2) + c(8)*j2(k,1))
  ogi(4,3) =  c(2) * (-(1-stab)*(mu1+mu2)*szs*j1(k,1)+(mu1-mu2)*j2(k,1)) 
  ogi(5,3) =  c(1) * ( mx1*mx2*(d1*j1(k,2)-j1(k,1)) - c(3)*j2(k,3) - c(4)*(z*t4v2-s*f4v2)*j2(k,2) + c(9)*j2(k,1))

end if 

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE OGINT2
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE BASINT1 (r, d, jx)
!----------------------------------------------------------------------------------------------------------------------------------!

!~evaluation of the basic integrals: jx(n,k) = Int (e^(-xi*d) xi^k J_n(xi*r) dxi), xi=0,infinity 
!~note: integrals jx(1,k) and jx(2,k) are normalized (divided) by r  

!~notation : J_n  is the bessel function of integer order n

real   (kind=dp), intent(in ) :: r, d
real   (kind=dp), intent(out) :: jx(0:2,0:3)
integer(kind=k9)              :: n, k
real   (kind=dp)              :: rs, ds, dist, coef

!~definition of the auxiliary variables 

rs   = r*r
ds   = d*d
dist = sqrt(ds+rs)  
coef = 1.0_dp/(dist+d)

!~calculation of the normalized integrals: jx(0,k), jx(1,k)/r, jx(2,k)/r    

loop1: do n = 0,2 
  loop2 : do k = 0,3 
    call BASINT2 (n, k, r, rs, d, ds, dist, coef, jx(n,k))
  end do loop2
end do loop1  

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE BASINT1
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE BASINT2 (n, k, r, rs, d, ds, dist, coef, jx)
!----------------------------------------------------------------------------------------------------------------------------------!

!~evaluation of the normalized integrals: jx(n=0,k), jx(n=1,k)/r, jx(n=2,k)/r    

integer(kind=k9), intent(in ) :: n, k
real   (kind=dp), intent(in ) :: r, rs, d, ds, dist, coef
real   (kind=dp), intent(out) :: jx 

!~evaluation of jx(0,k), jx(1,k)/r, jx(2,k)/r^2    

select case (k)
case(0)
  jx = coef**n/dist
case(1)
  jx = coef**n/dist**3 * (n*dist+d)
case(2)
  jx = coef**n/dist**5 * (n*dist*(n*dist+3.0_dp*d) + 2.0_dp*ds - rs)
case(3) 
  jx = coef**n/dist**7 * (n*dist*((n*dist)**2+11.0_dp*ds-4.0*rs) + 3.0_dp*d*(2.0_dp*(n*dist)**2 + 2.0_dp*ds - 3.0_dp*rs)) 
end select 

!~conversion jx(2,k)/r = (jx(2,k)/r^2) * r    

if (n==2)  jx = r * jx     

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE BASINT2 
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE GFL_TRANS (t, ul, tl, ug, tg)
!----------------------------------------------------------------------------------------------------------------------------------!

!~transformation of displacement/stress green's functions from cylindrical to cartesian coordinate system

!~notation : t => angular cylindrical coordinate theta 
!~notation : ul, tl => green's functions in cylindrical coordinates
!~notation : ug, tg => green's functions in cartesian coordinates 

real   (kind=dp), intent(in ) :: t 
complex(kind=dp), intent(in ) :: ul(3,3), tl(3,6)
complex(kind=dp), intent(out) :: ug(3,3), tg(3,6)
integer(kind=k9)              :: k
real   (kind=dp)              :: q(3,3) 

!~compose the transformation matrix

q = reshape((/cos(t), -sin(t), 0.0_dp, sin(t), cos(t), 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/), (/3, 3/)) 

!~transformation of the displacement green's functions

loop1: do k = 1,3  !-loop over the source directions  

  ug(k,1) = q(1,1)*ul(k,1) + q(2,1)*ul(k,2) + q(3,1)*ul(k,3) 
  ug(k,2) = q(1,2)*ul(k,1) + q(2,2)*ul(k,2) + q(3,2)*ul(k,3) 
  ug(k,3) = q(1,3)*ul(k,1) + q(2,3)*ul(k,2) + q(3,3)*ul(k,3) 

end do loop1

!~transformation of the stress green's functions

loop2: do k = 1,3  !-loop over the source directions 

  tg(k,1) = q(1,1)*q(1,1)*tl(k,1) + q(2,1)*q(2,1)*tl(k,2) + q(3,1)*q(3,1)*tl(k,3) +  &
            2.0_dp * (q(3,1)*q(2,1)*tl(k,4) + q(3,1)*q(1,1)*tl(k,5) + q(1,1)*q(2,1)*tl(k,6)) 

  tg(k,2) = q(1,2)*q(1,2)*tl(k,1) + q(2,2)*q(2,2)*tl(k,2) + q(3,2)*q(3,2)*tl(k,3) +  &
            2.0_dp * (q(3,2)*q(2,2)*tl(k,4) + q(3,2)*q(1,2)*tl(k,5) + q(1,2)*q(2,2)*tl(k,6)) 

  tg(k,3) = q(1,3)*q(1,3)*tl(k,1) + q(2,3)*q(2,3)*tl(k,2) + q(3,3)*q(3,3)*tl(k,3) +  &
            2.0_dp * (q(3,3)*q(2,3)*tl(k,4) + q(3,3)*q(1,3)*tl(k,5) + q(1,3)*q(2,3)*tl(k,6)) 

  tg(k,4) = q(1,3)*q(1,2)*tl(k,1) + q(2,3)*q(2,2)*tl(k,2) + q(3,3)*q(3,2)*tl(k,3) + (q(3,3)*q(2,2)+q(2,3)*q(3,2))*tl(k,4) + &
            (q(3,3)*q(1,2)+q(1,3)*q(3,2))*tl(k,5) +  (q(1,3)*q(2,2)+q(2,3)*q(1,2))*tl(k,6)

  tg(k,5) = q(1,3)*q(1,1)*tl(k,1) + q(2,3)*q(2,1)*tl(k,2) + q(3,3)*q(3,1)*tl(k,3) + (q(3,3)*q(2,1)+q(2,3)*q(3,1))*tl(k,4) + & 
            (q(3,3)*q(1,1)+q(1,3)*q(3,1))*tl(k,5) + (q(1,3)*q(2,1)+q(2,3)*q(1,1))*tl(k,6)

  tg(k,6) = q(1,1)*q(1,2)*tl(k,1) + q(2,1)*q(2,2)*tl(k,2) + q(3,1)*q(3,2)*tl(k,3) + (q(3,1)*q(2,2)+q(2,1)*q(3,2))*tl(k,4) + &
            (q(3,1)*q(1,2)+q(1,1)*q(3,2))*tl(k,5) + (q(1,1)*q(2,2)+q(2,1)*q(1,2))*tl(k,6) 

end do loop2 

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE GFL_TRANS  
!----------------------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE HALDGREEN4   
!----------------------------------------------------------------------------------------------------------------------------------!

