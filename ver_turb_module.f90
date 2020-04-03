MODULE VTURB_MOD

! A random displacement model is implemented to simulate sub-grid scale
!   turbulent particle motion in the vertical (z) direction. 
!
! Created by:             Elizabeth North
! Modified by:            Zachary Schlag
! Created on:             2003
! Last Modified on:       18 Aug 2008

IMPLICIT NONE
PUBLIC

CONTAINS

!      ********************** Vertical Turbulence *************************
!     ***********************************************************************
!    ** Random Displacement Model (Visser 1997 MEPS 158:275-281) for        **
!    ** simulating displacement due to turbulent diffusion                  **
!    ** (vertical direction)                                                **
!    **    z(t+1)= z + K'(z)dt + R{ 2/r K[z+0.5K'(z)dt]dt}**0.5             **
!    **    where z = particle vertical location at time t                   **
!    **      K' = dK/dz (Kprime) and K = vertical diffusivity (KH from ROMS)**
!    **      dt = time step of RDM (deltat)                                 **
!    **      R  = random process with mean = 0 and standard deviation = r.  **
!    **                                                                     **
!    ** Programmed by EW North February 2003 UMCES HPL enorth@hpl.umces.edu **
!    *************************************************************************

  SUBROUTINE VTurb(Xpar,Ypar,Zpar,ets,ex,ix,ng,TurbV)
    USE PARAM_MOD,  ONLY: s_w,idt,t_b,t_c,t_f,serr,smth,sub,deltat,AKSback
    use grid_mod, only: getWlevel
    USE INT_MOD,    ONLY: linint,polintd,getInterp3d,getInterp2D,getInterpAKs
    USE PDF_MOD,   ONLY: norm

  !  USE TENSION_MOD, ONLY: TSPSI,HVAL,HPVAL
    IMPLICIT NONE
	real ( kind = 8 ) smooth,ppvalu

    INTEGER, INTENT(IN) :: ets
    DOUBLE PRECISION, INTENT(IN) :: Xpar,Ypar,Zpar,ex(3),ix(3)
    DOUBLE PRECISION, INTENT(OUT) :: TurbV

   
	
    DOUBLE PRECISION :: DEV,r,zetab,zetac,zetaf,zeta,zfit(4),KHfit(4),ast,v(100)
    INTEGER :: i,j,k,jlo,loop,itop,ibot
    DOUBLE PRECISION :: tdepth,slopem,ParZc,Kprimec,KprimeZc,newZc,KH3rdc,Z3rdc, &
                        thisyc,ey(3),Pwc_KHb,Pwc_KHc,Pwc_KHf,zb,zc,zf,m,dz
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::fitz,fitKH,     &
      Pwc_KH,newxb,newyb,newxc,newyc,newxf,newyf,z

    !TSPACK Variables 
    INTEGER :: IER,SigErr
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: dy
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: a,coefsm

    !Number of values to proliferate to
    INTEGER :: np,ng
	DOUBLE PRECISION :: ttime(10),dTz
	
	! sub = 4
    np =int(dble(s_w(ng))*sub)
    !ALLOCATE VARIABLES
    ALLOCATE(Pwc_KH(s_w(ng)))     
    ALLOCATE(z(s_w(ng)))
    ALLOCATE(dy(np))
    ALLOCATE(a(np,4))
    ALLOCATE(coefsm(4,np))
    ALLOCATE(fitz(np))
    ALLOCATE(fitKH(np))

    !       *********************************************************
    !       *           Find kh in water column profile             *
    !       *********************************************************
	
	! background = 1.0D-6
	! dy=.0001
	! smth=6.0E-1
	v=0.0
	a=0.0
	coefsm=0.0
	fitz=0.0
	fitKH=0.0
	 tdepth = DBLE(-1.0)* getInterp2D("depth",ng,Xpar,Ypar,t_c)
     zetab =  getInterp2D("zeta",ng,Xpar,Ypar,t_b)
     zetac =  getInterp2D("zeta",ng,Xpar,Ypar,t_c)
     zetaf =  getInterp2D("zeta",ng,Xpar,Ypar,t_f)
    !i. find KH in water column profile at particle locationv(i)
		
    ey(1) = zetab
    ey(2) = zetac
    ey(3) = zetaf
    zeta = polintd(ex,ey,3,ix(2))
	!write(*,"(F10.4,F10.4,F10.4,F10.4)") zetab,zetac,zetaf,zeta

	
    do i=1,s_w(ng)
	    zb=getWlevel(zetab,tdepth,ng,i)
	    zc=getWlevel(zetac,tdepth,ng,i)
	    zf=getWlevel(zetaf,tdepth,ng,i)
		
		ey(1) = zb
		ey(2) = zc
		ey(3) = zf
		z(i) = polintd(ex,ey,3,ix(2))
		! Pwc_KHb(i) = getInterp3d("AKs",ng,Xpar,Ypar,zb(i),t_b,2,zetab,tdepth)
		! Pwc_KHc(i) = getInterp3d("AKs",ng,Xpar,Ypar,zc(i),t_c,2,zetac,tdepth)
		! Pwc_KHf(i) = getInterp3d("AKs",ng,Xpar,Ypar,zf(i),t_f,2,zetaf,tdepth)
		call getInterpAKs(ng,Xpar,Ypar,t_b,i,Pwc_KHb)
		call getInterpAKs(ng,Xpar,Ypar,t_c,i,Pwc_KHc)
		call getInterpAKs(ng,Xpar,Ypar,t_f,i,Pwc_KHf)
		
		ey(1) = Pwc_KHb
		ey(2) = Pwc_KHc
		ey(3) = Pwc_KHf
		Pwc_KH(i) = polintd(ex,ey,3,ix(2))
		
	
		!WRITE(*,*) i,z(i),ZPar
		
    enddo
	
	dz=(z(s_w(ng))-z(1))/dble(np-1)
    fitz(1)=z(1)
	fitz(np)=z(s_w(ng))
    fitKH(1)=Pwc_KH(1)
	fitKH(np)=Pwc_KH(s_w(ng))
	do i=2,np-1
		fitz(i)=fitz(i-1)+dz
		
		call linint(z,Pwc_KH,s_w(ng),fitz(i),fitKH(i),m)
		dy(i)=serr
	end do
	
	
	
	 ! do i = 1, s_w(ng)
	

        ! write ( *, '(2e15.7)' ) z(i),Pwc_KH(i)
     ! end do
	! write(*,*) '---------------'
	! do i = 1,np
	

        ! write ( *, '(2e15.7)' ) fitz(i),fitKH(i)
     ! end do
	
	
	! do i=1,s_w(ng)
	  		! if (ZPar.LT.z(i)) then
				! itop=i+1
				! ibot=i-2
				! exit
		! endif
		
    ! enddo
	! if (ibot.eq.0) then
		! ibot=ibot+1
		! itop=itop+1
	! endif
	! if (itop.eq.s_w(ng)+1) then
		! ibot=ibot-1
		! itop=itop-1
	! endif
		
     ! write(*,*) '---------------'
	 ! do i=1,4
		 ! zfit(i)=z(ibot+i-1)
		 ! KHfit(i)=Pwc_KH(ibot+i-1)
	 ! enddo	
	
	! ! dy=n

	! write(*,*) j
	! ! write(*,*) np
	! write(*,*) v
    
	ast=smooth ( fitz,fitKH, dy,np,smth, v, a )
		
    
	! write(*,*) j
	coefsm=0.0
	do i = 1, np
       coefsm(1:4,i) = a(i,1:4)
    end do
	
	 ! do i = 1, np
	   ! do j = 1, 4
	     
         ! a(i,j) = ppvalu ( fitz, coefsm, np-1, 4, fitz(i), j-1 )
       ! end do
       ! ! write ( *, '(6e15.7)' ) fitz(i),a(i,1:4)
     ! end do
	
	! ! ifitx(1)=z(1)
	! do i=2,p2	
		! ifitx(i)=ifitx(i-1)+0.3
	! end do
	
    ! write(*,*) '---------------'
	! do i=1,p2	
	  ! do j = 1, 4
		! KHfit(j)= ppvalu ( z, coefsm, s_w(ng)-1, 4, ifitx(i), j-1 )
		
	 ! end do
		
		! write ( *, '(5e15.7)' ) ifitx(i),KHfit
	! end do
		
    !  
    ! enddo

    ! !       *********************************************************
    ! !       *          Create Extra Points for Moving Average       *
    ! !       *********************************************************

    
	

    ! !vii. fit a tension spline to water column profile of KH using TSPACK
    ! !SigErr=0
  

    ! !viii. Initialize. Set deltat, number of iterations, initial z-coordinate
     ! deltat=1.0  ! dt= 1 sec

     loop= idt/int(deltat) ! number of iterations of RDM loop  
     ParZc = Zpar  !set initial particle depth  
	 !loop=1

    !       *********************************************************
    !       *           Random Displacement Model Loop              *
    !       *********************************************************

    !ix. Begin iterations
    do i=1,loop   
    !  a. Determine the second term in the equation: K'(z)deltat
    !    1. Find Kprime & solve for second term in RDM equation
    !       (K'(z)deltat = KprimeZ)
      Kprimec=0.0
      KprimeZc=0.0
	  
      if (ParZc.LE.tdepth .OR. ParZc.GE.zeta) then
        Kprimec=0.0                                                
      else
      !  CALL linint(fitx,ifity,p2,ParZc,thisyc,Kprimec)  
	     
         Kprimec = ppvalu ( fitz, coefsm, np-1, 4, ParZc, 1 )
		
      endif
	  
	  
      KprimeZc=Kprimec*deltat   

    !  b. Determine the 3rd term in the RDM equation: 
    !     R{ 2/r K[z+0.5K'(z)dt]dt}**0.5
    !    i. Find K at location of [z+0.5K'(z)dt] = Z3rd
    !      1. calculate Z3rd and make sure within boudaries
      Z3rdc = ParZc + DBLE(0.5)*KprimeZc

    !      2. Find KH at the location Z3rd
      KH3rdc=0.0
      if (Z3rdc.LT.tdepth .OR. Z3rdc.GT.zeta) then
        KH3rdc=AKSback                                
      else
         !! CALL linint(ifitx,ifity,p2,Z3rdc,KH3rdc,slopem)
         KH3rdc = ppvalu ( fitz, coefsm, np-1, 4, Z3rdc, 0 )
        if (KH3rdc.LT.AKSback) KH3rdc=AKSback    
      endif

    !  c. Solve the entire equation
      DEV=norm()        ! the random deviate
      r=1.                ! the standard deviation of the random deviate
	  

	!write(*,"(F10.4,F10.4,F10.4)") KprimeZc,zeta
	dTz=KprimeZc + DEV* (DBLE(2.0)/r * KH3rdc*deltat)**0.5
    newZc = ParZc + dTz
    !x. update particle z-coordinate
    ParZc = newZc
	  
	  
		if (ParZc.LE.tdepth) then
		ParZc = tdepth +  ABS(ParZc-tdepth)
	 endif
	 if (ParZc.GE.zeta) then
		ParZc = zeta - ABS(ParZc-zeta)
	endif
   enddo   !End RDM iterations

!           *********************************************************
!           *       Calculate displacement due to vertical turb     *
!           *********************************************************

	! call CPU_TIME(ttime(7))	
	! write(*,*)'----'
	! write(*,*) ttime(2)-ttime(1)
	! write(*,*) ttime(3)-ttime(2)
	! write(*,*) ttime(4)-ttime(3)
	! write(*,*) ttime(5)-ttime(4)
	! write(*,*) ttime(6)-ttime(5)
	! write(*,*) ttime(7)-ttime(6)
	! write(*,*) ttime(7)-ttime(10)
    !xi. find vertical displacement of particle due to turbulent diffusion
    TurbV = ParZc-Zpar
    !TurbV = 0.0


! **************** End Vertical Turbulence (RDM) *************************
! ************************************************************************

    !DEALLOCATE VARIABLES


    DEALLOCATE(fitz)
    DEALLOCATE(fitKH)
    DEALLOCATE(Pwc_KH)
    DEALLOCATE(z)
    DEALLOCATE(dy)
    DEALLOCATE(a)
    DEALLOCATE(coefsm)
	


  END SUBROUTINE VTurb

END MODULE VTURB_MOD