MODULE Ho_basis

  PRIVATE

  PUBLIC :: RadHO

CONTAINS
  ! log(Gamma(xx)) From numerical recipies
  FUNCTION gammln(xx)
    DOUBLE PRECISION gammln,xx
    ! Returns the value ln[GAM(xx)] for xx > 0.
    INTEGER j
    DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
    SAVE cof,stp
    DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,&
          24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,&
          -.5395239384953d-5,2.5066282746310005d0/
    x=xx
    y=x
    tmp=x+5.5d0
    tmp=(x+0.5d0)*log(tmp)-tmp
    ser=1.000000000190015d0
    do j=1,6
       y=y+1.d0
       ser=ser+cof(j)/y
    enddo
    gammln=tmp+log(stp*ser/x)
    return
  END FUNCTION gammln

  
!  L = L_(n)^(l+1/2) ! Using recursion relations

! Originally returned sign(L) and log(abs(L)) /Daniel W
  Subroutine LaguerreL2(n, l, RVEC, FVEC, FVEC_S, dimr) 
    IMPLICIT NONE
    INTEGER :: n,l,Ik,II, dimr
    DOUBLE PRECISION :: RVEC(dimr), FVEC(dimr),FVEC_S(dimr),bin,alpha
    DOUBLE PRECISION :: L0(dimr), L1(dimr)

    IF(n.lt.0) THEN
       FVEC   = 0.D0
       FVEC_S = 0.D0
       RETURN
    END IF

    alpha = l + 0.5d0
    
    L1 = 0.d0; FVEC = 1.d0;
    DO ii = 1 , n 
       L0 = L1
       L1 = FVEC
       FVEC = ((2.d0 * ii- 1.d0 + alpha- RVEC)* L1- (ii- 1.d0 + alpha)* L0)/ dble(ii)
    END DO
  
    FVEC_S = sign(1.D0,FVEC)
    !Commented away statement below I want the value not the log /Daniel W
    !FVEC   = log(abs(FVEC)) 
    !End Comment /Daniel W
    RETURN 
  END Subroutine LaguerreL2

! Calculates a vector of values of g_nl(r), where g_nl(r) is the radial part of a HO wave function with size parameter b as defined
! on page 49 of Jouni Sohonen, From Nucleons to Nucleus, Springer
  Subroutine RadHO(n, l, b, RVEC, FVEC, dimr)
    IMPLICIT NONE
    INTEGER :: n,l,dimr
    DOUBLE PRECISION :: nR,lR,b, RVEC(dimr), FVEC(dimr), FVEC_S(dimr), FVECtmp(dimr), lnfac
    nR = DBLE(n) 
    lR = DBLE(l)
    CALL LaguerreL2(n, l, (RVEC/b)**2, FVEC, FVEC_S, dimr)
     
    ! gamma(n+1) = n!
    IF (n == 0) THEN
      lnfac = 0d0
    ELSE
      lnfac = gammln(nR+1d0)
    END IF

   FVEC = SQRT(2d0/b**3)* EXP(0.5d0* ( lnfac - gammln(nR+lR+1.5d0) ) )* (RVEC/b)**l* EXP(-RVEC**2/2/b**2)* FVEC 
   
  END Subroutine RadHO

! Computes overlap certain of radial HO functions < R_n0^(bN) | R_00^(ba) > with a closed expression
! n is nubmer of nodes of radial ho-wave function with l=0 and oscillator length bN, ba is oscillator length of n=0,l=0 wave function
  FUNCTION olrho(n,bN,ba)
    ! Returns the value of < R_n0^(bN) | R_00^(ba) >
    IMPLICIT NONE
    INTEGER :: n
    DOUBLE PRECISION :: nDbl, olrho, bN, ba  
    IF (n<0 .OR. bN<=0d0 .OR. ba<=0d0) THEN
      olrho = 0d0
      RETURN
    END IF

    IF (n==0) THEN
      olrho = ( 2d0*bN*ba/(bN**2+ba**2) )**(3d0/2d0)
      RETURN
    END IF
    
    nDbl = DBLE(n)
    olrho = EXP( 0.5d0*gammln(2*nDbl + 2) - gammln(nDbl + 1) )*2d0**(3d0/2d0 - nDbl)*(bN*ba)**(3d0/2d0)*(bN**2 - ba**2)**nDbl/(bN**2 + ba**2)**(3d0/2d0+nDbl)
    

  END FUNCTION olrho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION ylm(l,m,theta, phi )

  INTEGER :: l, m
  DOUBLE PRECISION :: theta, phi, lD, mD 
  COMPLEX :: ylm
  DOUBLE PRECISION , PARAMETER ::  pi = 2d0*ACOS(0.0)
  
  IF(ABS(m)>l) THEN
    ylm = 0d0
    RETURN
  END IF
  
  lD = DBLE(l)
  IF(m >= 0) THEN
    mD = DBLE(m)
    ylm = sqrt( (2*l+1d0)/4d0/pi ) * exp( .5d0*( gammln(lD-mD+1d0) - gammln(lD+mD+1d0) ) )*plgndr(l,m,COS(theta))*exp( (0d0,1d0)*m*phi )
    RETURN
  ELSE !Y_l(-m) = (-1)^m Y^*_lm
    mD = DBLE(-1d0*m)
    m = -1*m
    ylm = (-1d0)**m*sqrt( (2*l+1d0)/4d0/pi ) * exp( .5d0*( gammln(lD-mD+1d0) - gammln(lD+mD+1d0) ) )*plgndr(l,m,COS(theta))*exp( (0d0,-1d0)*m*phi )
  END IF  

  END FUNCTION ylm



!FROM NUMERICAL RECIPES IN FORTRAN77 and FORTRAN 90 PressW.H. et al. Cambridge University Publishing 1997
!page 247
! CHANGED REAL TO DOUBLE PRECISION and removed do labels

  FUNCTION plgndr(l,m,x)
  INTEGER l,m
  DOUBLE PRECISION :: plgndr,x
  !Computes the associated Legendre polynomial Plm (x). Here m and l are integers satisfying
  !0 <= m <= l, while x lies in the range -1 <= x <= 1.
  INTEGER i,ll
  DOUBLE PRECISION :: fact,pll,pmm,pmmp1,somx2
  if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.) THEN !pause ?bad arguments in plgndr?
    plgndr = 0d0
    RETURN
  end if
  
  pmm=1.0d0  !Compute P_m^m
  if(m.gt.0) then
      somx2=sqrt((1.0d0-x)*(1.0d0+x))
      fact=1.0d0
      do i=1,m
	pmm=-pmm*fact*somx2
	fact=fact+2.d0
      enddo 
  endif
  if(l.eq.m) then
      plgndr=pmm
  else
    pmmp1=x*(2d0*m+1d0)*pmm !Compute P_m+1^m 
    if(l.eq.m+1) then
      plgndr=pmmp1
    else !Compute P_l^m , l > m + 1.
      do ll=m+2,l
	pll=(x*(2d0*ll-1d0)*pmmp1-(ll+m-1d0)*pmm)/(ll-m)
	pmm=pmmp1
	pmmp1=pll
      enddo 
      plgndr=pll
    endif
  endif
  return
  END FUNCTION plgndr
    

END MODULE
