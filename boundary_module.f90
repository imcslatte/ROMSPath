MODULE BOUNDARY_MOD

! This module contains variables and subroutines associated with the 
! creation of the land/sea boundaries.  The main purpose of this module 
! is to create the land/sea boundaries from a given masked rho grid.  
! The main subroutine in the module, createBounds, determines the number
! of boundary points, allocates an array of that size, and fills it with 
! the boundary points in order.
! ROMSPath Version: 1.0.1

IMPLICIT NONE
PRIVATE
SAVE

!*****************************************************************
!*                          VARIABLES                            *
!*****************************************************************


  !final boundary variables, after reformatting from bnds
  !DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) ::
  !INTEGER, ALLOCATABLE, DIMENSION (:) :: 
  !DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: bnd_x,bnd_y

  !TRUE if the boundary is land, FALSE if it is open ocean
  !LOGICAL, ALLOCATABLE, DIMENSION(:) :: land

  INTEGER :: i,j,ng
  PUBLIC :: bounds,zbounds

CONTAINS

!*****************************************************************
!*                    FUNCTIONS & SUBROUTINES                    *
!*****************************************************************


    !*********************************************************
    !*                     Boundaries                   *
    !*********************************************************


SUBROUTINE bounds(ng,Ipar,Jpar,nmask,ingrid,obound)
!This subroutine is for adding boundary points to bnds
! INPUT:
  use grid_mod,   	only: GRIDS 
  use param_mod,    only: xi_rho,eta_rho	
  use int_mod,    only: inside	
  IMPLICIT NONE
  DOUBLE PRECISION,  INTENT(IN) :: Ipar,Jpar
  DOUBLE PRECISION :: X,Y,ecut,fcut
  DOUBLE PRECISION :: tX(4),tY(4),m(2,2)
  INTEGER :: I,J,it,jt,n
  LOGICAL :: incell
  INTEGER, INTENT(IN) :: ng 
  INTEGER, INTENT(OUT) :: nmask 
  LOGICAL, INTENT(OUT) :: ingrid,obound 
	ecut=0.1 ! POSSIBLY ADD TO PARAM_MOD	
	fcut=1.0-ecut! POSSIBLY ADD TO PARAM_MOD	
	ingrid=.FALSE.
	obound=.FALSE.
	I=floor(Ipar)
	J=floor(Jpar)
	X=Ipar-dble(I)
	Y=Jpar-dble(J)
	
	
	tX=0.0d0
	tY=0.0d0
	
	if ((I.GT.1).and.(J.GT.1).and.(I.LT.(xi_rho(ng))).and.(J.LT.(eta_rho(ng))))  then
		m(1,1) = GRIDS(ng)%mask_rho(I,J)
		m(2,1) = GRIDS(ng)%mask_rho(I+1,J)
		m(2,2) = GRIDS(ng)%mask_rho(I+1,J+1)
		m(1,2) = GRIDS(ng)%mask_rho(I,J+1)
		nmask=m(1,1)+m(2,1)+m(2,2)+m(1,2)
		if (nmask.eq.4) then
			ingrid=.TRUE.
			! if ((I.EQ.1).or.(J.EQ.1).or.(I.EQ.xi_rho(ng)-1).or.(J.EQ.eta_rho(ng)-1))  then
				! obound=.TRUE.
			! endif	
		endif
		
		if (nmask.eq.3) then
			n=1
			do it=1,2
				do jt=1,2
					if (m(it,jt).eq.1)  then
						tX(n)=dble(it)-1.0
						tY(n)=dble(jt)-1.0
						n=n+1
					endif
				enddo
			enddo
			
			incell=inside (3, tX, tY, X, Y)
			if (incell) then
				ingrid=.TRUE.
				! if ((I.EQ.1).or.(J.EQ.1).or.(I.EQ.xi_rho(ng)-1).or.(J.EQ.eta_rho(ng)-1))  then
					! obound=.TRUE.
				! endif
			endif				
		endif
		
		if (nmask.eq.2) then
			if (m(1,1).eq.1)  then
				if ((X.lt.ecut).and.(Y.lt.ecut)) then
					ingrid=.TRUE.
				endif						
			endif
			if (m(2,1).eq.1)  then
				if ((X.gt.fcut).and.(Y.lt.ecut)) then
					ingrid=.TRUE.
				endif						
			endif
			if (m(2,2).eq.1)  then
				if ((X.gt.fcut).and.(Y.gt.fcut)) then
					ingrid=.TRUE.
				endif						
			endif
			if (m(1,2).eq.1)  then
				if ((X.lt.ecut).and.(Y.gt.fcut)) then
					ingrid=.TRUE.
				endif						
			endif	
			
			
		endif
		
		
		if ((I.LE.1).or.(J.LE.1).or.(I.GE.xi_rho(ng)-1).or.(J.GE.eta_rho(ng)-1))  then
			obound=.TRUE.
		endif
		
	else
	
		obound=.TRUE.
	endif
	

	
	
END SUBROUTINE bounds

SUBROUTINE zbounds(ng,Ipar,Jpar,Zpar,ingrid,t)
!This subroutine is for adding boundary points to bnds
! INPUT:
  use grid_mod,   	only: GRIDS 
  use param_mod,    only: xi_rho,eta_rho	
  use INT_MOD,   only:   getInterp2D
  IMPLICIT NONE
  DOUBLE PRECISION,  INTENT(IN) :: Ipar,Jpar
  DOUBLE PRECISION,  INTENT(INOUT) :: Zpar
  INTEGER :: I,J
  INTEGER, INTENT(IN) :: ng ,t
  LOGICAL, INTENT(OUT) :: ingrid
  DOUBLE PRECISION :: tdepth,tzeta
  
	ingrid=.FALSE.
	I=floor(Ipar)
	J=floor(Jpar)
	if (I.GT.0) then
		tdepth = DBLE(-1.0)* getInterp2D("depth",ng,Ipar,Jpar,1)
		tzeta =  getInterp2D("zeta",ng,Ipar,Jpar,t)
		
		
		if ((Zpar.LT.tzeta).and.(Zpar.GT.tdepth)) then
			ingrid=.TRUE.
		endif
		
		if (Zpar.GE.tzeta) then
		    Zpar=tzeta-0.01D0 !set particle depth to 1 less than 
			ingrid=.TRUE.
		endif
		
		
	endif
	
END SUBROUTINE zbounds



 
	


  ! This subroutine calculates the intersection between the particle
  ! trajectory and the boundary line in a grid cell, and then calculates
  ! the reflection, returning the new particle location
  ! subroutine intersect_reflect(Xpos,Ypos,nXpos,nYpos,fintersectX,fintersectY,  &
    ! freflectX,freflectY,intersectf,skipbound,isWater)
    ! IMPLICIT NONE
    ! INTEGER, INTENT(OUT) :: intersectf
    ! INTEGER, INTENT(INOUT) :: skipbound
    ! DOUBLE PRECISION, INTENT(IN) :: Xpos,Ypos,nXpos,nYpos
    ! DOUBLE PRECISION, INTENT(OUT) :: fintersectX,fintersectY,freflectX,freflectY
    ! LOGICAL, OPTIONAL, INTENT(OUT) :: isWater
    ! INTEGER :: i,intersect,skipboundi
    ! DOUBLE PRECISION :: crossk,dPBC,mBCperp,rx1,rx2,ry1,ry2,Bp,distBC,dist1,   &
      ! dist2,intersctx,interscty,rPxyzX,rPxyzY,Mbc,Bbc,Mp,bcx1,bcy1,bcx2,bcy2,  &
      ! bBCperp,xhigh,xlow,yhigh,ylow,d_Pinter,dtest,bxhigh,bxlow,byhigh,bylow

    ! distBC=0.0
    ! Mbc = 0.0
    ! Bbc = 0.0
    ! Mp = 0.0
    ! Bp = 0.0
    ! intersect=0
    ! intersectf=0
    ! skipboundi = skipbound
    ! fintersectX = -999999.
    ! fintersectY = -999999.
    ! freflectX = -999999.
    ! freflectY = -999999.
    ! dtest = 999999.
    ! isWater = .FALSE.

    ! if (Xpos.GE.nXpos) then
      ! xhigh = Xpos
      ! xlow = nXpos
    ! else
      ! xhigh = nXpos
      ! xlow = Xpos
    ! endif  

    ! if (Ypos.GE.nYpos) then
      ! yhigh = Ypos
      ! ylow = nYpos
    ! else
      ! yhigh = nYpos
      ! ylow = Ypos
    ! endif

    ! do i=1,nbounds

      ! if (i == skipbound) cycle

        ! intersect = 0
        ! bcx1=bnd_x(1,i)
        ! bcy1=bnd_y(1,i)
        ! bcx2=bnd_x(2,i)
        ! bcy2=bnd_y(2,i)

        ! !If the boundary segment end points are both east, west, north, or 
        ! !  south of the particle's previous or new location, cycle to next 
        ! !  boundary
        ! if( ((bcx1 > xhigh) .AND. (bcx2 > xhigh)) .OR. &
            ! ((bcx1 < xlow ) .AND. (bcx2 < xlow )) .OR. &
            ! ((bcy1 > yhigh) .AND. (bcy2 > yhigh)) .OR. &
            ! ((bcy1 < ylow ) .AND. (bcy2 < ylow ))      ) cycle
        
        ! if (bcx1.GE.bcx2) then
          ! bxhigh = bcx1
          ! bxlow = bcx2
        ! else
          ! bxhigh = bcx2
          ! bxlow = bcx1
        ! endif  

        ! if (bcy1.GE.bcy2) then
          ! byhigh = bcy1
          ! bylow = bcy2
        ! else
          ! byhigh = bcy2
          ! bylow = bcy1
        ! endif

        ! !First determine if an undefined denominator is possible
        ! if (bcx1.EQ.bcx2 .OR. nXpos.EQ.Xpos ) then
          ! !test if they both vertical, if so cycle because they cannot intersect
          ! if (bcx1.EQ.bcx2 .AND. nXpos.EQ.Xpos ) cycle
          ! !test if perpendicular and parrallel to coordinate axes
          ! if (bcx1.EQ.bcx2 .AND. nYpos.EQ.Ypos ) then
            ! !undefined denominator, perp. & || to axes
            ! intersctx = bcx1
            ! interscty = nYpos
            ! if (intersctx.LE.xhigh  .AND. intersctx.GE.xlow .AND.              &
              ! interscty.LE.yhigh  .AND. interscty.GE.ylow .AND.                &
              ! intersctx.LE.bxhigh .AND. intersctx.GE.bxlow .AND.               &
              ! interscty.LE.byhigh .AND. interscty.GE.bylow  ) then
              ! dPBC=sqrt((intersctx-nXpos)**2+(interscty-nYpos)**2)
              ! rx1=nXpos+(DBLE(2.0)*dPBC)
              ! ry1=nYpos
              ! rx2=nXpos-(DBLE(2.0)*dPBC)
              ! ry2=nYpos
              ! dist1 = sqrt( (intersctx-rx1)**2 + (interscty-ry1)**2 )
              ! dist2 = sqrt( (intersctx-rx2)**2 + (interscty-ry2)**2 )
              ! if(dist1.LT.dist2) then
                ! rPxyzX= rx1
                ! rPxyzY= ry1
              ! elseif(dist1.GT.dist2) then
                ! rPxyzX= rx2
                ! rPxyzY= ry2
              ! endif
              ! intersect=1
            ! endif
          ! elseif (nXpos.EQ.Xpos .AND. bcy1.EQ.bcy2 ) then
            ! !undefined denominator, perp. & || to axes
            ! intersctx = nXpos
            ! interscty = bcy1
            ! if (intersctx.LE.xhigh .AND.  intersctx.GE.xlow .AND.              &
              ! interscty.LE.yhigh  .AND. interscty.GE.ylow .AND.                &
              ! intersctx.LE.bxhigh .AND. intersctx.GE.bxlow .AND.               &
              ! interscty.LE.byhigh .AND. interscty.GE.bylow  ) then
              ! dPBC=sqrt((intersctx-nXpos)**2+(interscty-nYpos)**2)
              ! rx1=nXpos
              ! ry1=nYpos+(DBLE(2.0)*dPBC)
              ! rx2=nXpos
              ! ry2=nYpos-(DBLE(2.0)*dPBC)
              ! dist1 = sqrt( (intersctx-rx1)**2 + (interscty-ry1)**2 )
              ! dist2 = sqrt( (intersctx-rx2)**2 + (interscty-ry2)**2 )
              ! if(dist1.LT.dist2) then
                ! rPxyzX= rx1
                ! rPxyzY= ry1
              ! elseif(dist1.GT.dist2) then
                ! rPxyzX= rx2
                ! rPxyzY= ry2
              ! endif
              ! intersect=1
            ! endif
          ! elseif (bcx1.EQ.bcx2 .AND. nYpos.NE.Ypos ) then
            ! !undefined denominator, not perpendicular
            ! Mp = (nYpos-Ypos)/(nXpos-Xpos)
            ! Bp = Ypos - Mp*Xpos
            ! intersctx = bcx1
            ! interscty = Mp*intersctx + Bp
            ! if (intersctx.LE.xhigh .AND.  intersctx.GE.xlow .AND.              &
                ! interscty.LE.yhigh  .AND. interscty.GE.ylow .AND.              &
                ! intersctx.LE.bxhigh .AND. intersctx.GE.bxlow .AND.             &
                ! interscty.LE.byhigh .AND. interscty.GE.bylow  ) then
              ! dPBC = nXpos-intersctx
              ! rx1=nXpos+(DBLE(2.0)*dPBC)
              ! ry1=nYpos
              ! rx2=nXpos-(DBLE(2.0)*dPBC)
              ! ry2=nYpos
              ! dist1 = sqrt( (intersctx-rx1)**2 + (interscty-ry1)**2 )
              ! dist2 = sqrt( (intersctx-rx2)**2 + (interscty-ry2)**2 )
              ! if(dist1.LT.dist2) then
                ! rPxyzX= rx1
                ! rPxyzY= ry1
              ! elseif(dist1.GT.dist2) then
                ! rPxyzX= rx2
                ! rPxyzY= ry2
              ! endif
              ! intersect=1
            ! endif
          ! elseif (nXpos.EQ.Xpos .AND. bcy1.NE.bcy2  ) then
            ! !undefined denominator, not perpendicular
            ! Mbc = (bcy2-bcy1)/(bcx2-bcx1)
            ! Bbc = bcy2 - Mbc*bcx2
            ! intersctx = nXpos
            ! interscty = Mbc*intersctx + Bbc
            ! if (intersctx.LE.xhigh .AND.  intersctx.GE.xlow .AND.              &
                ! interscty.LE.yhigh  .AND. interscty.GE.ylow .AND.              &
                ! intersctx.LE.bxhigh .AND. intersctx.GE.bxlow .AND.             &
                ! interscty.LE.byhigh .AND. interscty.GE.bylow ) then
              ! !Now use cross product to determine the distance of the particle 
              ! !  from the boundary
              ! distBC = sqrt((bcx1-bcx2)**2+(bcy1-bcy2)**2)
              ! crossk= ((nXpos-bcx1)*(bcy2-bcy1)) - ((bcx2-bcx1)*(nYpos-bcy1))
              ! dPBC = sqrt(crossk**2)/distBC
              ! !find line perpendicular to boundary
              ! mBCperp = DBLE(-1.0)/Mbc
              ! bBCperp = nYpos - mBCperp*nXpos
              ! !find two potential reflection points
              ! rx1 = sqrt( ((DBLE(2.0)*dPBC)**2)/(DBLE(1.0)+mBCperp**2) ) +nXpos
              ! ry1 = mBCperp*rx1 + bBCperp
              ! rx2 = sqrt( ((DBLE(2.0)*dPBC)**2)/(DBLE(1.0)+mBCperp**2) )       &
                  ! * DBLE(-1.0) + nXpos
              ! ry2 = mBCperp*rx2 + bBCperp
              ! !point closest to intersection of boundary and particle trajectory
              ! !  is the right one
              ! dist1 = sqrt( (intersctx-rx1)**2 + (interscty-ry1)**2 )
              ! dist2 = sqrt( (intersctx-rx2)**2 + (interscty-ry2)**2 )
              ! if(dist1.LT.dist2) then
                ! rPxyzX= rx1
                ! rPxyzY= ry1
              ! elseif(dist1.GT.dist2) then
                ! rPxyzX= rx2
                ! rPxyzY= ry2
              ! endif
              ! intersect=1
            ! endif
          ! endif
        ! else

          ! if(intersect == 0)then

            ! Mbc = (bcy2-bcy1)/(bcx2-bcx1)
            ! Bbc = bcy2 - Mbc*bcx2
            ! Mp = (nYpos-Ypos)/(nXpos-Xpos)
            ! Bp = Ypos - Mp*Xpos
            ! intersctx = (Bbc - Bp)/(Mp - Mbc)
            ! interscty = Mp*intersctx + Bp

            ! !when bc parallel with x-axis, byhigh=bylow=intersecty
            ! if (Mbc.EQ.0.0) interscty = byhigh
        
            ! if (intersctx.LE.xhigh .AND.  intersctx.GE.xlow .AND.              &
                ! interscty.LE.yhigh  .AND. interscty.GE.ylow .AND.              &
                ! intersctx.LE.bxhigh .AND. intersctx.GE.bxlow .AND.             &
                ! interscty.LE.byhigh .AND. interscty.GE.bylow  ) then

              ! if (Mbc.EQ.0.0) then  !inverse slope denominator not OK
                ! dPBC = nYpos-bcy1
                ! rx1=nXpos
                ! ry1=nYpos+(DBLE(2.0)*dPBC)
                ! rx2=nXpos
                ! ry2=nYpos-(DBLE(2.0)*dPBC)
                ! dist1 = sqrt( (intersctx-rx1)**2 + (interscty-ry1)**2 )
                ! dist2 = sqrt( (intersctx-rx2)**2 + (interscty-ry2)**2 ) 
                ! if(dist1.LT.dist2) then
                  ! rPxyzX= rx1
                  ! rPxyzY= ry1
                ! elseif(dist1.GT.dist2) then
                  ! rPxyzX= rx2
                  ! rPxyzY= ry2
                ! endif
                  ! intersect=1
                ! endif  

              ! if(intersect == 0)then

                ! !Now use cross product to determine the distance of the
                ! !  particle from the boundary
                ! distBC = sqrt((bcx1-bcx2)**2+(bcy1-bcy2)**2)
                ! crossk= ((nXpos-bcx1)*(bcy2-bcy1)) - ((bcx2-bcx1)*(nYpos-bcy1))
                ! dPBC = sqrt(crossk**2)/distBC
                ! !find line perpendicular to boundary
                ! mBCperp = DBLE(-1.0)/Mbc
                ! bBCperp = nYpos - mBCperp*nXpos
                ! !find two potential reflection points
                ! rx1 = sqrt(((DBLE(2.0)*dPBC)**2)/(DBLE(1.0)+mBCperp**2)) +nXpos
                ! ry1 = mBCperp*rx1 + bBCperp
                ! rx2 = sqrt(((DBLE(2.0)*dPBC)**2)/(DBLE(1.0)+mBCperp**2))       &
                    ! * DBLE(-1.0) + nXpos
                ! ry2 = mBCperp*rx2 + bBCperp
                ! !point closest to intersection of boundary and particle 
                ! !  trajectory is the right one
                ! dist1 = sqrt( (intersctx-rx1)**2 + (interscty-ry1)**2 )
                ! dist2 = sqrt( (intersctx-rx2)**2 + (interscty-ry2)**2 )
                ! if(dist1.LT.dist2) then
                  ! rPxyzX= rx1
                  ! rPxyzY= ry1
                ! elseif(dist1.GT.dist2) then
                  ! rPxyzX= rx2
                  ! rPxyzY= ry2
                ! endif
                ! intersect=1
                
              ! endif
            ! endif  
          ! endif
        ! endif


        ! d_Pinter = sqrt( (Xpos-intersctx)**2 + (Ypos-interscty)**2 )
        ! if( (intersect .EQ. 1) .AND. (d_Pinter .LT. dtest) ) then
          ! fintersectX = intersctx
          ! fintersectY = interscty
          ! freflectX = rPxyzX
          ! freflectY = rPxyzY
          ! intersectf = 1
          ! dtest = d_Pinter
          ! skipboundi = i
          ! isWater = .NOT. land(i)
        ! endif

    ! enddo

    ! skipbound = skipboundi
  ! END SUBROUTINE intersect_reflect




END MODULE
