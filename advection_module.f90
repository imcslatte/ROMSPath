MODULE ADVECTION_MOD


! Particle Advection module separated from main particel update progema.  
!
! Created by:             Elias Hunter
! Modified by:            Elias Hunter
! Created on:             8/6/2019
! Last Modified on:       8/6/2019

IMPLICIT NONE
PUBLIC

CONTAINS


!    *************************************************************************

  SUBROUTINE RKAdvect(Xpar,Ypar,Zpar,ex,ix,pm,pn,ng,ets,AdvectX,AdvectY,AdvectZ)
    USE PARAM_MOD,  ONLY: idt
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: Xpar,Ypar,ex(3),ix(3),pm,pn
	INTEGER, INTENT(IN) :: ets
    DOUBLE PRECISION, INTENT(OUT) :: AdvectX,AdvectY,AdvectZ

    DOUBLE PRECISION :: kn1_u,kn1_v,kn1_w,kn2_u,kn2_v,kn2_w,kn3_u,kn3_v,kn3_w,kn4_u,kn4_v,kn4_w, &
     UAD,VAD,WAD,x1,x2,x3,y1,y2,y3,z1,z2,z3,P_U,P_V,P_W
	
    INTEGER :: i,j,k,jlo,loop,itop,ibot
	DOUBLE PRECISION :: P_depth,P_zeta,Zpar


    
    INTEGER :: ng


	

    !       *********************************************************
    !       *       Runga Kutta            *
    !       *********************************************************
		
	   CALL find_currents(Xpar,Ypar,Zpar,ex,ix,ng,ets,Uad,Vad,Wad,P_depth,P_zeta)
		
     
      ! !Store advection currents at original coordinates
       kn1_u = Uad*pm
       kn1_v = Vad*pn
       kn1_w = Wad

      ! !Estimate new coordinates for next RK position
       x1 = Xpar + Uad*pm*DBLE(idt)/DBLE(2) 
       y1 = Ypar + Vad*pn*DBLE(idt)/DBLE(2) 
       z1 = Zpar + Wad*DBLE(idt)/DBLE(2) 
      ! if(z1 .GT. minpartdepth) z1 = minpartdepth - DBLE(0.000001)
      ! if(z1 .LT. maxpartdepth) z1 = maxpartdepth + DBLE(0.000001)
			! OPEN(1,FILE='testdata2')
			! write(1,"(F10.4,F10.4)") P_depth,P_zeta
			! write(1,"(E14.5,E14.5,E14.5,E14.5)") Xpar,Ypar,Zpar
			! write(1,"(E14.5,E14.5,E14.5,E14.5)") x1,y1,z1
			! CLOSE(1)
      ! !Find advection currents at estimated next RK position
       CALL find_currents(x1,y1,z1,ex,ix,ng,ets,Uad,Vad,Wad,P_depth,P_zeta)
      ! CALL find_currents(x1,y1,z1,Pwc_zb,Pwc_zc,Pwc_zf,Pwc_wzb,Pwc_wzc,        &
                         ! Pwc_wzf,P_zb,P_zc,P_zf,ex,ix,p,2,Uad,Vad,Wad)

      ! !Store advection currents at 2nd RK position
       kn2_u = Uad*pm
       kn2_v = Vad*pn
       kn2_w = Wad

      ! !Estimate new coordinates for next RK position
       x2 = Xpar + Uad*pm*DBLE(idt)/DBLE(2) 
       y2 = Ypar + Vad*pn*DBLE(idt)/DBLE(2) 
       z2 = Zpar + Wad*DBLE(idt)/DBLE(2) 
      ! if(z2 .GT. minpartdepth) z2 = minpartdepth - DBLE(0.000001)
      ! if(z2 .LT. maxpartdepth) z2 = maxpartdepth + DBLE(0.000001)

      ! !Find advection currents at estimated next RK position
	  
		
		
       CALL find_currents(x2,y2,z2,ex,ix,ng,ets,Uad,Vad,Wad,P_depth,P_zeta)
      ! CALL find_currents(x2,y2,z2,Pwc_zb,Pwc_zc,Pwc_zf,Pwc_wzb,Pwc_wzc,        &
                         ! Pwc_wzf,P_zb,P_zc,P_zf,ex,ix,p,2,Uad,Vad,Wad)

      ! !Store advection currents at 3rd RK position
      kn3_u = Uad*pm
      kn3_v = Vad*pn
      kn3_w = Wad
	  	
      ! !Calculate the coordinates at the final position

       x3 = Xpar + Uad*pm*DBLE(idt)/DBLE(2) 
       y3 = Ypar + Vad*pn*DBLE(idt)/DBLE(2) 
       z3 = Zpar + Wad*DBLE(idt)/DBLE(2) 
      ! if(z3 .GT. minpartdepth) z3 = minpartdepth - DBLE(0.000001)
      ! if(z3 .LT. maxpartdepth) z3 = maxpartdepth + DBLE(0.000001)

      ! !Find advection currents at the final position
       CALL find_currents(x3,y3,z3,ex,ix,ng,ets,Uad,Vad,Wad,P_depth,P_zeta)
      ! CALL find_currents(x3,y3,z3,Pwc_zb,Pwc_zc,Pwc_zf,Pwc_wzb,Pwc_wzc,        &
                         ! Pwc_wzf,P_zb,P_zc,P_zf,ex,ix,p,3,Uad,Vad,Wad)

      ! !Store advection currents at final position
      kn4_u = Uad*pm
      kn4_v = Vad*pn
      kn4_w = Wad

      ! !Use the RK formula to get the final Advection values
      P_U = (kn1_u + DBLE(2.0)*kn2_u + DBLE(2.0)*kn3_u + kn4_u)/DBLE(6.0)
      P_V = (kn1_v + DBLE(2.0)*kn2_v + DBLE(2.0)*kn3_v + kn4_v)/DBLE(6.0)
      P_W = (kn1_w + DBLE(2.0)*kn2_w + DBLE(2.0)*kn3_w + kn4_w)/DBLE(6.0)

      AdvectX = idt*(P_U)
      AdvectY = idt*(P_V)
      AdvectZ = idt*P_W
      
	
	


  END SUBROUTINE RKAdvect
  
  

 SUBROUTINE find_currents(Xpar,Ypar,Zpar,ex,ix,ng,ets,Uad,Vad,Wad,tdepth,zeta)
    !This Subroutine calculates advection currents at the particle's 
    !  location in space and time
	USE PARAM_MOD, ONLY: t_b,t_c,t_f,s_rho,s_w,zob
    USE GRID_MOD,  ONLY: GRIDS
    USE INT_MOD,    ONLY: getInterp2D,getInterp3D,polintd
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: Xpar,Ypar,ex(3),ix(3)
	DOUBLE PRECISION, INTENT(INOUT) ::Zpar
    DOUBLE PRECISION, INTENT(OUT) :: Uad,Vad,Wad,tdepth,zeta

	INTEGER, INTENT(IN) :: ets

    INTEGER :: i,ii,iii,ng,version
    DOUBLE PRECISION :: P_Ub,P_Uc,P_Uf,P_Vb,P_Vc,P_Vf,P_Wb,P_Wc,P_Wf,ey(3),    &
      Pwc_ub,Pwc_uc,Pwc_uf,Pwc_vb,Pwc_vc,Pwc_vf,Pwc_wb,Pwc_wc,Pwc_wf

     DOUBLE PRECISION :: pm,pn,zetab,zetac,zetaf,&
						zb,zc,zf
	  
	  tdepth = DBLE(-1.0)* getInterp2D("depth",ng,Xpar,Ypar,t_c)
      zetab =  getInterp2D("zeta",ng,Xpar,Ypar,t_b)
      zetac =  getInterp2D("zeta",ng,Xpar,Ypar,t_c)
      zetaf =  getInterp2D("zeta",ng,Xpar,Ypar,t_f)
	  
	  
	!  write(*,"(F10.4,F10.4,F10.4)") zetab,zetac,zetaf
	  !*******************************
	  !*******************************
	  !*****KLUDGE WARNIN******************
	  !*******************************
	  version=2

	  !*******************************
	  !*******************************
	  !*******************************
	  !*******************************
      ! ! !Check if particle location above or below boundary, If so, place
      ! ! !  just within boundary (1 mm)
      if (Zpar.LT.tdepth) then
        Zpar = tdepth + DBLE(0.001)
        !IF(TrackCollisions) hitBottom(n) = hitBottom(n) + 1
      endif
     

	  
      ey(1) = zetab
      ey(2) = zetac
      ey(3) = zetaf
      zeta = polintd(ex,ey,3,ix(2))

      if (Zpar.GT.zeta) Zpar = zeta - DBLE(0.001) 
	  
	  

  
      ! ! *********************************************************
      ! ! *                                                       *
      ! ! *             Create Matrix of Z-coordinates            *
      ! ! *                                                       *
      ! ! *********************************************************

      ! !Create matrix of z-coordinates at particle and at each node for
      ! !  back, center, forward times
	  !write(*,*) '------------------'

    !           *********************************************************
    !           *                                                       *
    !           *       Calculate U,V,W in Water Column Profile         *
    !           *                                                       *
    !           *********************************************************


    !i. Determine if particle is deep enough that velocities are affected by 
    !  the bottom.
    ! If so, apply log layer between deepest current velocity predicitons 
    ! (deepest rho s-level for u,v and deepest w s-level for w) and bottom.
    ! ! OR, if below z0, set advection velocities to 0.0
    ! if ((Zpar .LT. Pwc_wzb(1)+z0) .OR. &
        ! (Zpar .LT. Pwc_wzc(1)+z0) .OR. &
        ! (Zpar .LT. Pwc_wzf(1)+z0)      ) then

      ! Uad = 0.0
      ! Vad = 0.0
      ! Wad = 0.0

    ! elseif ((Zpar .LT. Pwc_zb(1)) .OR. &
            ! (Zpar .LT. Pwc_zc(1)) .OR. &
            ! (Zpar .LT. Pwc_zf(1))      ) then
	 ! write(*,*) '----------------'
       Pwc_Uc = getInterp3d("u",ng,Xpar,Ypar,Zpar,t_c,1,zeta,tdepth)
       Pwc_Uf = getInterp3d("u",ng,Xpar,Ypar,Zpar,t_f,1,zeta,tdepth)
       Pwc_Ub = getInterp3d("u",ng,Xpar,Ypar,Zpar,t_b,1,zeta,tdepth)
	   
       Pwc_Vb = getInterp3d("v",ng,Xpar,Ypar,Zpar,t_b,1,zeta,tdepth)
       Pwc_Vc = getInterp3d("v",ng,Xpar,Ypar,Zpar,t_c,1,zeta,tdepth)
       Pwc_Vf = getInterp3d("v",ng,Xpar,Ypar,Zpar,t_f,1,zeta,tdepth)
       Pwc_Wb = getInterp3d("w",ng,Xpar,Ypar,Zpar,t_b,2,zeta,tdepth)
       Pwc_Wc = getInterp3d("w",ng,Xpar,Ypar,Zpar,t_c,2,zeta,tdepth)
       Pwc_Wf = getInterp3d("w",ng,Xpar,Ypar,Zpar,t_f,2,zeta,tdepth)
	   
	   ! write(*,*) zeta,tdepth
	  ! ! write(*,*) Xpar,Ypar,Zpar
	   ! write(*,*)Pwc_Ub,Pwc_Uc,Pwc_Uf
		
      !  u(z)  = [ u(z1) / (log(z1/zo) ] * (log (z/zo) 
      !where:
      !  u is current velocity
      !  z1 is height of first sigma level above bottom
      !  z0 is roughness height of model
      !  z is height of desired velocity
      !
      !  Note that Pwc_wzb(1) = P_depth = Depth at particle location

      ! P_Ub=Pwc_Ub*log10((Zpar-Pwc_wzb(1))/z0)/log10((Pwc_zb(1) -Pwc_wzb(1))/z0)
      ! P_Uc=Pwc_Uc*log10((Zpar-Pwc_wzb(1))/z0)/log10((Pwc_zc(1) -Pwc_wzb(1))/z0)
      ! P_Uf=Pwc_Uf*log10((Zpar-Pwc_wzb(1))/z0)/log10((Pwc_zf(1) -Pwc_wzb(1))/z0)
      ! P_Vb=Pwc_Vb*log10((Zpar-Pwc_wzb(1))/z0)/log10((Pwc_zb(1) -Pwc_wzb(1))/z0)
      ! P_Vc=Pwc_Vc*log10((Zpar-Pwc_wzb(1))/z0)/log10((Pwc_zc(1) -Pwc_wzb(1))/z0)
      ! P_Vf=Pwc_Vf*log10((Zpar-Pwc_wzb(1))/z0)/log10((Pwc_zf(1) -Pwc_wzb(1))/z0)
      ! P_Wb=Pwc_Wb*log10((Zpar-Pwc_wzb(1))/z0)/log10((Pwc_wzb(2)-Pwc_wzb(1))/z0)
      ! P_Wc=Pwc_Wc*log10((Zpar-Pwc_wzb(1))/z0)/log10((Pwc_wzc(2)-Pwc_wzb(1))/z0)
      ! P_Wf=Pwc_Wf*log10((Zpar-Pwc_wzb(1))/z0)/log10((Pwc_wzf(2)-Pwc_wzb(1))/z0)
	  
	  
      P_Ub=Pwc_Ub
      P_Uc=Pwc_Uc
      P_Uf=Pwc_Uf
      P_Vb=Pwc_Vb
      P_Vc=Pwc_Vc
      P_Vf=Pwc_Vf
      P_Wb=Pwc_Wb
      P_Wc=Pwc_Wc
      P_Wf=Pwc_Wf


      !     *********************************************************
      !     *        Find Internal b,c,f and Advection Values       *
      !     *********************************************************
      !
      ! ii. fit polynomial to hydrodynamic model output and find internal 
      !     b,c,f values

      !a. U velocity
      ! 1. Prepare external time step values
      if (ets .EQ. 1) then
        ey=0.0
        ey(1) = P_Ub
        ey(2) = P_Ub
        ey(3) = P_Uc
      else
        ey=0.0
        ey(1) = P_Ub
        ey(2) = P_Uc
        ey(3) = P_Uf
      endif
	
      ! 2. Get Advection value
      if(version .EQ. 1) then
        Uad = polintd(ex,ey,3,ix(1))
      elseif (version .EQ. 2) then
        Uad = polintd(ex,ey,3,ix(2))
      else
        Uad = polintd(ex,ey,3,ix(3))
      endif

      !b. V velocity
      ! 1. Prepare external time step values
      if (ets .EQ. 1) then
        ey=0.0
        ey(1) = P_Vb
        ey(2) = P_Vb
        ey(3) = P_Vc
      else
        ey=0.0
        ey(1) = P_Vb
        ey(2) = P_Vc
        ey(3) = P_Vf
      endif
		
      ! 2. Get Advection value
      if(version .EQ. 1) then
        Vad = polintd(ex,ey,3,ix(1))
      elseif (version .EQ. 2) then
        Vad = polintd(ex,ey,3,ix(2))
      else
        Vad = polintd(ex,ey,3,ix(3))
      endif


      !c. W velocity
      ! 1. Prepare external time step values
      if (ets .EQ. 1) then
        ey=0.0
        ey(1) = P_Wb
        ey(2) = P_Wb
        ey(3) = P_Wc
      else
        ey=0.0
        ey(1) = P_Wb
        ey(2) = P_Wc
        ey(3) = P_Wf
      endif

      ! 2. Get Advection value
      if(version .EQ. 1) then
        Wad = polintd(ex,ey,3,ix(1))
      elseif (version .EQ. 2) then
        Wad = polintd(ex,ey,3,ix(2))
      else
        Wad = polintd(ex,ey,3,ix(3))
      endif


	  
	
	   
    RETURN
  END SUBROUTINE find_currents
	

END MODULE ADVECTION_MOD