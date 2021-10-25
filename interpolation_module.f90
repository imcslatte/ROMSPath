	MODULE INT_MOD

!  The Interpolation Module contains two procedures that interpolate data:
!  Subroutine linint uses linear interpolation.
!  Subroutine polintd uses polynomial interpolation.
!
!  Created by:            Elizabeth North
!  Modified by:           Elias Hunter
!  Created on:            2003
!  Last Modified on:      26 Apr 2019
! ROMSPath Version: 1.0.1

IMPLICIT NONE
PUBLIC

CONTAINS


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~                                                               ~~
! ~~                     SUBROUTINE linint                         ~~
! ~~                                                               ~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE linint(xa,ya,n,x,y,m)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    DOUBLE PRECISION, INTENT(IN) :: xa(n),ya(n),x
    DOUBLE PRECISION, INTENT(OUT) :: y,m

    INTEGER :: jlo,jhi,k
    DOUBLE PRECISION :: b

    jlo=1  
    jhi=n

    !determine which two xa location that x lies between
    do
      k=(jhi+jlo)/2

      if(xa(k).gt.x)then
        jhi=k
      else
        jlo=k
      endif

      if (jhi-jlo == 1) exit
    enddo

    !calculate the slope and intersect
    m = ( ya(jlo) - ya(jhi) )/ ( xa(jlo) - xa(jhi) )
    b = ya(jlo) - m*xa(jlo)

    !linearly interpolate y
    y = m*x + b

    return

  END SUBROUTINE linint


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~                                                               ~~
! ~~                     FUNCTION polintd                          ~~
! ~~                                                               ~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  DOUBLE PRECISION FUNCTION polintd(xa,ya,n,x)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    DOUBLE PRECISION, INTENT(IN) :: xa(n),ya(n),x

    INTEGER :: i,ns
    DOUBLE PRECISION :: dif,dift,a,b,c

    ns=1
    dif=abs(x-xa(1))

    !Here we find the index, ns, of the closest table entry
    do i=1,n
      dift=abs(x-xa(i))
      if (dift.lt.dif) then
        ns=i
        dif=dift
      endif
    enddo 

    !calculate the value of c
    c = (xa(2)-x) * ( (ya(3)-ya(2)) / (xa(2)-xa(3)) )
    c = c - (xa(2)-x) * ( (ya(2)-ya(1)) / (xa(1)-xa(2)) )
    c = c / (xa(1) - xa(3))

    !calculate the values of a and b
    if (ns .EQ. 3) then
      a = (ya(3)-ya(2))/(xa(2)-xa(3))
      b = xa(3)-x
    else
      a = (ya(2)-ya(1))/(xa(1)-xa(2))
      b = xa(1)-x
    endif

    !calculate the value of polintd via polynomial interpolation
    polintd = ya(ns) + (xa(ns)-x)*a + b*c

  END FUNCTION polintd
  
  
  
  
  subroutine LL2ij(longrd,latgrd,angle,inlon,inlat,numpar,	&
				   xi_rho,eta_rho,Ipar,Jpar)
	use param_mod, only: Eradius,deg2rad,rad2deg
	integer :: numpar,xi_rho,eta_rho
    double precision,intent(in) :: longrd(1:xi_rho,1:eta_rho),latgrd(1:xi_rho,1:eta_rho),angle(1:xi_rho,1:eta_rho)				
	double precision,intent(in) :: inlon(1:numpar),inlat(1:numpar)	
	double precision,intent(out) :: Ipar(1:numpar),Jpar(1:numpar)	
	integer :: i,Imin,Imax,Jmin,Jmax,i0,j0,j					
	logical :: found, foundi, foundj      
	double precision :: aa2, ang, bb2, diag2, dx, dy, phi
    double precision :: xfac, xpp, yfac, ypp
	
!    OPEN(10,FILE='TESTSTUFF.dat',STATUS='REPLACE')
	do i=1,numpar
		Imin=1
		Imax=xi_rho
		Jmin=1
		Jmax=eta_rho
	
		

		found=inbox(Imin,Imax, Jmin,Jmax,    &
     &                      longrd, latgrd,  &
     &                      inlon(i), inlat(i))
		

		
		IF (found) THEN
			DO while (((Imax-Imin).gt.1).or.((Jmax-Jmin).gt.1))
                IF ((Imax-Imin).gt.1) THEN
                  i0=(Imin+Imax)/2
                  found=inbox(Imin, i0, Jmin,Jmax,    &
     &                      longrd, latgrd,  &
     &                      inlon(i), inlat(i))
                  IF (found) THEN
                    Imax=i0
                  ELSE
                    Imin=i0
                  END IF
                END IF
                IF ((Jmax-Jmin).gt.1) THEN
                  j0=(Jmin+Jmax)/2
                  found=inbox(Imin, Imax, Jmin,j0,    &
     &                      longrd, latgrd,  &
     &                      inlon(i), inlat(i))
                  IF (found) THEN
                    Jmax=j0
                  ELSE
                    Jmin=j0
                  END IF
                END IF
				
            END DO
		
		endif
	
	
!		write(10,"(F10.4,F10.4,F10.4,F10.4,F10.4,F10.4)") inlon(i), inlat(i),longrd(Imin,Jmin),latgrd(Imin,Jmin),longrd(Imax,Jmax),latgrd(Imax,Jmax)
!		write(10,"(F10.4,F10.4,I4,I4,I4,I4)") inlon(i), inlat(i),Imin,Imax,Jmin,Jmax
				found=inbox(Imin,Imax, Jmin,Jmax,    &
     &                      longrd, latgrd,  &
     &                      inlon(i), inlat(i))
	 if (found) then
			yfac=Eradius*deg2rad
			xfac=yfac*COS(inlat(i)*deg2rad)
			xpp=(inlon(i)-longrd(Imin,Jmin))*xfac
			ypp=(inlat(i)-latgrd(Imin,Jmin))*yfac
			
			diag2=((longrd(Imin+1,Jmin)-longrd(Imin,Jmin+1))*xfac)**2+      &
		 &            ((latgrd(Imin+1,Jmin)-latgrd(Imin,Jmin+1))*yfac)**2
				aa2=((longrd(Imin,Jmin)-longrd(Imin+1,Jmin))*xfac)**2+          &
		 &          ((latgrd(Imin,Jmin)-latgrd(Imin+1,Jmin))*yfac)**2
				bb2=((longrd(Imin,Jmin)-longrd(Imin,Jmin+1))*xfac)**2+          &
		 &          ((latgrd(Imin,Jmin)-latgrd(Imin,Jmin+1))*yfac)**2
				phi=ASIN((diag2-aa2-bb2)/(2.0*SQRT(aa2*bb2)))
			
			
			ang=angle(Imin,Jmin)
			dx=xpp*COS(ang)+ypp*SIN(ang)
			dy=ypp*COS(ang)-xpp*SIN(ang)
	
			dx=dx+dy*TAN(phi)
			dy=dy/COS(phi)
	
			dx=MIN(MAX(0.0,dx/SQRT(aa2)),1.0)
			dy=MIN(MAX(0.0,dy/SQRT(bb2)),1.0)
			Ipar(i)=DBLE(Imin)+dx
			Jpar(i)=DBLE(Jmin)+dy
		else 
			
			Ipar(i)=-1.0
			Jpar(i)=-1.0
		endif
	enddo
!	CLOSE(10)
  end subroutine LL2ij
  
  
  LOGICAL FUNCTION inbox(Li,Ui, Lj,Uj,                     &
     &                      longrd, latgrd,  &
     &                      inlon, inlat)
  integer, intent(in) :: Li,Ui, Lj,Uj
  double precision,intent(in) :: longrd(:,:),latgrd(:,:)
  double precision,intent(in) :: inlon, inlat
  integer ::  Nb, i, j, shft, ic
  double precision, dimension(2*(Uj-Lj+Ui-Li)+1) :: Lonb, Latb
  
      Nb=2*(Uj-Lj+Ui-Li)
      shft=1-Li
      DO i=Li,Ui-1
        Lonb(i+shft)=longrd(i,Lj)
        Latb(i+shft)=latgrd(i,Lj)
      END DO
      shft=1-Lj+Ui-Li
      DO j=Lj,Uj-1
        Lonb(j+shft)=longrd(Ui,j)
        Latb(j+shft)=latgrd(Ui,j)
      END DO
      shft=1+Uj-Lj+2*Ui-Li
      DO i=Ui,Li+1,-1
        Lonb(shft-i)=longrd(i,Uj)
        Latb(shft-i)=latgrd(i,Uj)
      END DO
      shft=1+2*Uj-Lj+2*(Ui-Li)
      DO j=Uj,Lj+1,-1
        Lonb(shft-j)=longrd(Li,j)
        Latb(shft-j)=latgrd(Li,j)
      END DO

	  inbox=inside(Nb, Lonb, Latb, inlon, inlat)

  RETURN
  END FUNCTION inbox
  LOGICAL FUNCTION inside (Nb, Xb, Yb, Xo, Yo)
!
!=======================================================================
!                                                                      !
!  Given the vectors Xb and Yb of size Nb, defining the coordinates    !
!  of a closed polygon,  this function find if the point (Xo,Yo) is    !
!  inside the polygon.  If the point  (Xo,Yo)  falls exactly on the    !
!  boundary of the polygon, it still considered inside.                !
!                                                                      !
!  This algorithm does not rely on the setting of  Xb(Nb)=Xb(1) and    !
!  Yb(Nb)=Yb(1).  Instead, it assumes that the last closing segment    !
!  is (Xb(Nb),Yb(Nb)) --> (Xb(1),Yb(1)).                               !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Reid, C., 1969: A long way from Euclid. Oceanography EMR,         !
!      page 174.                                                       !
!                                                                      !
!  Algorithm:                                                          !
!                                                                      !
!  The decision whether the point is  inside or outside the polygon    !
!  is done by counting the number of crossings from the ray (Xo,Yo)    !
!  to (Xo,-infinity), hereafter called meridian, by the boundary of    !
!  the polygon.  In this counting procedure,  a crossing is counted    !
!  as +2 if the crossing happens from "left to right" or -2 if from    !
!  "right to left". If the counting adds up to zero, then the point    !
!  is outside.  Otherwise,  it is either inside or on the boundary.    !
!                                                                      !
!  This routine is a modified version of the Reid (1969) algorithm,    !
!  where all crossings were counted as positive and the decision is    !
!  made  based on  whether the  number of crossings is even or odd.    !
!  This new algorithm may produce different results  in cases where    !
!  Xo accidentally coinsides with one of the (Xb(k),k=1:Nb) points.    !
!  In this case, the crossing is counted here as +1 or -1 depending    !
!  of the sign of (Xb(k+1)-Xb(k)).  Crossings  are  not  counted if    !
!  Xo=Xb(k)=Xb(k+1).  Therefore, if Xo=Xb(k0) and Yo>Yb(k0), and if    !
!  Xb(k0-1) < Xb(k0) < Xb(k0+1),  the crossing is counted twice but    !
!  with weight +1 (for segments with k=k0-1 and k=k0). Similarly if    !
!  Xb(k0-1) > Xb(k0) > Xb(k0+1), the crossing is counted twice with    !
!  weight -1 each time.  If,  on the other hand,  the meridian only    !
!  touches the boundary, that is, for example, Xb(k0-1) < Xb(k0)=Xo    !
!  and Xb(k0+1) < Xb(k0)=Xo, then the crossing is counted as +1 for    !
!  segment k=k0-1 and -1 for segment k=k0, resulting in no crossing.   !
!                                                                      !
!  Note 1: (Explanation of the logical condition)                      !
!                                                                      !
!  Suppose  that there exist two points  (x1,y1)=(Xb(k),Yb(k))  and    !
!  (x2,y2)=(Xb(k+1),Yb(k+1)),  such that,  either (x1 < Xo < x2) or    !
!  (x1 > Xo > x2).  Therefore, meridian x=Xo intersects the segment    !
!  (x1,y1) -> (x2,x2) and the ordinate of the point of intersection    !
!  is:                                                                 !
!                                                                      !
!                 y1*(x2-Xo) + y2*(Xo-x1)                              !
!             y = -----------------------                              !
!                          x2-x1                                       !
!                                                                      !
!  The mathematical statement that point  (Xo,Yo)  either coinsides    !
!  with the point of intersection or lies to the north (Yo>=y) from    !
!  it is, therefore, equivalent to the statement:                      !
!                                                                      !
!         Yo*(x2-x1) >= y1*(x2-Xo) + y2*(Xo-x1),   if   x2-x1 > 0      !
!  or                                                                  !
!         Yo*(x2-x1) <= y1*(x2-Xo) + y2*(Xo-x1),   if   x2-x1 < 0      !
!                                                                      !
!  which, after noting that  Yo*(x2-x1) = Yo*(x2-Xo + Xo-x1) may be    !
!  rewritten as:                                                       !
!                                                                      !
!        (Yo-y1)*(x2-Xo) + (Yo-y2)*(Xo-x1) >= 0,   if   x2-x1 > 0      !
!  or                                                                  !
!        (Yo-y1)*(x2-Xo) + (Yo-y2)*(Xo-x1) <= 0,   if   x2-x1 < 0      !
!                                                                      !
!  and both versions can be merged into  essentially  the condition    !
!  that (Yo-y1)*(x2-Xo)+(Yo-y2)*(Xo-x1) has the same sign as x2-x1.    !
!  That is, the product of these two must be positive or zero.         !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: Nb
      double precision, intent(in) :: Xo, Yo
      DOUBLE PRECISION, intent(inout) :: Xb(:), Yb(:)
!
!  Local variable declarations.
!
      integer, parameter :: Nstep =128
      integer :: crossings, i, inc, k, kk, nc
      integer, dimension(Nstep) :: Sindex
      DOUBLE PRECISION :: dx1, dx2, dxy
!
!-----------------------------------------------------------------------
!  Find intersections.
!-----------------------------------------------------------------------
!
!  Set crossings counter and close the contour of the polygon.
!

	
      crossings=0
      Xb(Nb+1)=Xb(1)
      Yb(Nb+1)=Yb(1)

!  The search is optimized.  First select the indices of segments
!  where Xb(k) is different from Xb(k+1) and Xo falls between them.
!  Then, further investigate these segments in a separate loop.
!  Doing it in two stages takes less time because the first loop is
!  pipelined.
!
      DO kk=0,Nb-1,Nstep
        nc=0
        DO k=kk+1,MIN(kk+Nstep,Nb)
          IF (((Xb(k+1)-Xo)*(Xo-Xb(k)).ge.0.0).and.                  &
     &        (Xb(k).ne.Xb(k+1))) THEN
            nc=nc+1
            Sindex(nc)=k
          END IF
        END DO
        DO i=1,nc
          k=Sindex(i)
          IF (Xb(k).ne.Xb(k+1)) THEN
            dx1=Xo-Xb(k)
            dx2=Xb(k+1)-Xo
            dxy=dx2*(Yo-Yb(k))-dx1*(Yb(k+1)-Yo)
            inc=0
            IF ((Xb(k).eq.Xo).and.(Yb(k).eq.Yo)) THEN
              crossings=1
              goto 10
            ELSE IF (((dx1.eq.0.0).and.(Yo.ge.Yb(k  ))).or.          &
     &              ((dx2.eq.0.0).and.(Yo.ge.Yb(k+1)))) THEN
              inc=1
            ELSE IF ((dx1*dx2.gt.0.0).and.                           &
     &              ((Xb(k+1)-Xb(k))*dxy.ge.0.0)) THEN  ! see note 1
              inc=2
            END IF
            IF (Xb(k+1).gt.Xb(k)) THEN
              crossings=crossings+inc
            ELSE
              crossings=crossings-inc
            END IF
          END IF
        END DO
      END DO
!
!  Determine if point (Xo,Yo) is inside of closed polygon.
!
  10  IF (crossings.eq.0) THEN
        inside=.FALSE.
      ELSE
        inside=.TRUE.
      END IF
      RETURN
      END FUNCTION inside
	  DOUBLE PRECISION FUNCTION getInterp2d(var,ng,Ipar,Jpar,t)
		!This Function returns the interpolated value at the particle's location 
		!  
		use grid_mod, only: GRIDS
		use hydro_mod, only: HYDRODATA
		IMPLICIT NONE
		CHARACTER(LEN=*), INTENT(IN) :: var
		INTEGER, INTENT(IN) :: ng
		DOUBLE PRECISION, INTENT(IN) :: Ipar,Jpar
		INTEGER :: Inode,Jnode,t,i
		DOUBLE PRECISION :: v(4),m(4),cff(4)
		DOUBLE PRECISION :: cff4,X,Y
		DOUBLE PRECISION :: Ival,nwater

	    Inode=floor(Ipar)
		Jnode=floor(Jpar)
		X=Ipar-Inode
		Y=Jpar-Jnode
		
		m(1) = GRIDS(ng)%mask_rho(Inode,Jnode)
		m(2) = GRIDS(ng)%mask_rho(Inode+1,Jnode)
		m(3) = GRIDS(ng)%mask_rho(Inode+1,Jnode+1)
		m(4) = GRIDS(ng)%mask_rho(Inode,Jnode+1)
	!	 f(x,y)=f(0,0)(1-x)(1-y)+f(1,0)x(1-y)+f(0,1)(1-x)y+f(1,1)xy,}
	
	    cff(1)=(1-X)*(1-Y)
	    cff(2)=X*(1-Y)
	    cff(3)=X*Y
	    cff(4)=(1-X)*Y
		!Determine which data to interpolate from
		SELECT CASE(var)
		  CASE("depth")
			v(1) = GRIDS(ng)%H(Inode,Jnode)
			v(2) = GRIDS(ng)%H(Inode+1,Jnode)
			v(3) = GRIDS(ng)%H(Inode+1,Jnode+1)
			v(4) = GRIDS(ng)%H(Inode,Jnode+1)
		  CASE("pm")
			v(1) = GRIDS(ng)%pm(Inode,Jnode)
			v(2) = GRIDS(ng)%pm(Inode+1,Jnode)
			v(3) = GRIDS(ng)%pm(Inode+1,Jnode+1)
			v(4) = GRIDS(ng)%pm(Inode,Jnode+1)
		  CASE("pn")
			v(1) = GRIDS(ng)%pn(Inode,Jnode)
			v(2) = GRIDS(ng)%pn(Inode+1,Jnode)
			v(3) = GRIDS(ng)%pn(Inode+1,Jnode+1)
			v(4) = GRIDS(ng)%pn(Inode,Jnode+1)
		  CASE("lon")
			v(1) = GRIDS(ng)%lon_rho(Inode,Jnode)
			v(2) = GRIDS(ng)%lon_rho(Inode+1,Jnode)
			v(3) = GRIDS(ng)%lon_rho(Inode+1,Jnode+1)
			v(4) = GRIDS(ng)%lon_rho(Inode,Jnode+1)
		  CASE("lat")
			v(1) = GRIDS(ng)%lat_rho(Inode,Jnode)
			v(2) = GRIDS(ng)%lat_rho(Inode+1,Jnode)
			v(3) = GRIDS(ng)%lat_rho(Inode+1,Jnode+1)
			v(4) = GRIDS(ng)%lat_rho(Inode,Jnode+1)
		  CASE("zeta")
				
			v(1) = HYDRODATA(ng)%zeta(Inode,Jnode,t)
			v(2) = HYDRODATA(ng)%zeta(Inode+1,Jnode,t)
			v(3) = HYDRODATA(ng)%zeta(Inode+1,Jnode+1,t)
			v(4) = HYDRODATA(ng)%zeta(Inode,Jnode+1,t)
			
		
			call getcoeff(X,Y,m,cff)
			
			
		  CASE DEFAULT
			write(*,*) 'Problem interpolating ',var
			write(*,*) ' '
			write(*,*) 'The Program Cannot Continue and Will Terminate'
			stop
		END SELECT

		
		! if(tOK == 1) then 
		!Ival=cff1*v(1)+cff2*v(2)+cff3*v(3)+cff4*v(4)
		!		write(*,*) v(1),v(2),v(3),v(4),Ival
		

		
		getInterp2D = cff(1)*v(1)+cff(2)*v(2)+cff(3)*v(3)+cff(4)*v(4)
		


	  END FUNCTION getInterp2D
	  
	DOUBLE PRECISION FUNCTION getInterp3d(var,ng,Ipar,Jpar,Zpar,t,ztype,zeta,depth)
		!This Function returns the interpolated value at the particle's location 
		!  using the interpolation variables stored from function setInterp, and
		!  the hydrodynamic variables that have been read in

		use grid_mod, only: GRIDS,getSlevel,getWlevel
		use hydro_mod, only: HYDRODATA
		use param_mod, only: s_rho,s_w,xi_rho,eta_rho,xi_u,eta_u,xi_v,eta_v
		
		IMPLICIT NONE
		CHARACTER(LEN=*), INTENT(IN) :: var
		INTEGER, INTENT(IN) :: ng,t,ztype
		DOUBLE PRECISION, INTENT(IN) :: Ipar,Jpar,Zpar,zeta,depth
		INTEGER :: Inode,Jnode,nz,i,Itop,Ibot
		DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: zprof
		DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: v
		DOUBLE PRECISION :: m(4),cff(4),zcff(2)
		DOUBLE PRECISION :: X,Y,uvIpar,uvJpar,Z
		DOUBLE PRECISION :: Ival,zlev(2),nwater
	    SELECT CASE(ztype)
		  CASE(int(1))
			  nz=s_rho(ng)
			  itop=nz
			  ibot=nz			  
		      do i=1,nz
				
				zlev(2)=getSlevel(zeta,depth,ng,i)
				if (zlev(2).gt.Zpar) then
					Itop=i
					Ibot=i-1
					if (i.eq.1) then
						zlev(1)=depth
					else
						zlev(1)=getSlevel(zeta,depth,ng,i-1)
					endif
					exit
				endif
				
			  enddo
			  
		    if (i.gt.nz) then
				Itop=nz
				Ibot=nz
				zcff(1)=0.5
				zcff(2)=0.5
			elseif (Itop.eq.1) then
				Itop=1
				Ibot=1
				zcff(1)=0.5
				zcff(2)=0.5
			else
				Z=(Zpar-zlev(1))/(zlev(2)-zlev(1))
				zcff(1)=1.0-Z
				zcff(2)=Z
				
			endif
		  CASE (int(2))
			  nz=s_w(ng)
			  itop=nz
			  ibot=nz		
		      do i=1,nz
				zlev(2)=getWlevel(zeta,depth,ng,i)
				if (zlev(2).gt.Zpar) then
					zlev(1)=getWlevel(zeta,depth,ng,i-1)
					Itop=i
					Ibot=i-1
					exit
				endif
			enddo
			
			Z=(Zpar-zlev(1))/(zlev(2)-zlev(1))
			zcff(1)=1.0-Z
			zcff(2)=Z
			
			
		  CASE DEFAULT
			write(*,*) 'NOT A VALID DEPTH TYPE'
			write(*,*) ' '
			write(*,*) 'The Program Cannot Continue and Will Terminate'
			stop
		END SELECT
		
		ALLOCATE(v(4,nz))
		ALLOCATE(zprof(nz))
		
		Inode=floor(Ipar)
		Jnode=floor(Jpar)
		X=Ipar-Inode
		Y=Jpar-Jnode
		

	!	 f(x,y)=f(0,0)(1-x)(1-y)+f(1,0)x(1-y)+f(0,1)(1-x)y+f(1,1)xy,}
		cff(1)=(1-X)*(1-Y)
	    cff(2)=X*(1-Y)
	    cff(3)=X*Y
	    cff(4)=(1-X)*Y
		!Determine which data to interpolate from
		
		SELECT CASE(var)
		  CASE("u")
			uvIpar=Ipar-0.5d0
			Inode=floor(uvIpar)
			
			
			
			X=uvIpar-Inode
			cff(1)=(1-X)*(1-Y)
			cff(2)=X*(1-Y)
			cff(3)=X*Y
			cff(4)=(1-X)*Y
			
		  
			
			if (Inode.GE.xi_u(ng)) then
				
				v(1,:) = HYDRODATA(ng)%U(xi_u(ng),Jnode,:,t)
				v(2,:) = HYDRODATA(ng)%U(xi_u(ng),Jnode,:,t)
				v(3,:) = HYDRODATA(ng)%U(xi_u(ng),Jnode+1,:,t)	
				v(4,:) = HYDRODATA(ng)%U(xi_u(ng),Jnode+1,:,t)	
#ifdef WETDRY
				m(1) = HYDRODATA(ng)%wetdry_mask_u(xi_u(ng),Jnode,t)
				m(2) = 0
				m(3) = 0
				m(4) = HYDRODATA(ng)%wetdry_mask_u(xi_u(ng),Jnode+1,t)
#else			
				m(1) = GRIDS(ng)%mask_u(xi_u(ng),Jnode)
				m(2) = 0
				m(3) = 0
				m(4) = GRIDS(ng)%mask_u(xi_u(ng),Jnode+1)	
#endif					
			
			elseif  (Inode.LT.1) then
				v(1,:) = HYDRODATA(ng)%U(1,Jnode,:,t)
				v(2,:) = HYDRODATA(ng)%U(1,Jnode,:,t)
				v(3,:) = HYDRODATA(ng)%U(1,Jnode+1,:,t)	
				v(4,:) = HYDRODATA(ng)%U(1,Jnode+1,:,t)	
#ifdef WETDRY		
				m(1) = 0
				m(2) = HYDRODATA(ng)%wetdry_mask_u(1,Jnode,t)
				m(3) = HYDRODATA(ng)%wetdry_mask_u(1,Jnode+1,t)
				m(4) = 0
#else	
				m(1) = 0
				m(2) = GRIDS(ng)%mask_u(1,Jnode)
				m(3) = GRIDS(ng)%mask_u(1,Jnode+1)
				m(4) = 0	
#endif			
			else
				v(1,:) = HYDRODATA(ng)%U(Inode,Jnode,:,t)
				v(2,:) = HYDRODATA(ng)%U(Inode+1,Jnode,:,t)
				v(3,:) = HYDRODATA(ng)%U(Inode+1,Jnode+1,:,t)	
				v(4,:) = HYDRODATA(ng)%U(Inode,Jnode+1,:,t)
#ifdef WETDRY	
		        m(1) = HYDRODATA(ng)%wetdry_mask_u(Inode,Jnode,t)
				m(2) = HYDRODATA(ng)%wetdry_mask_u(Inode+1,Jnode,t)
				m(3) = HYDRODATA(ng)%wetdry_mask_u(Inode+1,Jnode+1,t)
				m(4) = HYDRODATA(ng)%wetdry_mask_u(Inode,Jnode+1,t)
#else			
				m(1) = GRIDS(ng)%mask_u(Inode,Jnode)
				m(2) = GRIDS(ng)%mask_u(Inode+1,Jnode)
				m(3) = GRIDS(ng)%mask_u(Inode+1,Jnode+1)
				m(4) = GRIDS(ng)%mask_u(Inode,Jnode+1)
#endif		
			endif
			
			nwater=m(1)+m(2)+m(3)+m(4)
			
			if ((nwater.LT.4.0).and.(nwater.GT.0.0)) then
			
				do i=1,nz
					if (m(1).eq.0.0d0) v(1,i)=0.0
					if (m(2).eq.0.0d0) v(2,i)=0.0
					if (m(3).eq.0.0d0) v(3,i)=0.0
					if (m(4).eq.0.0d0) v(4,i)=0.0
					
				enddo
			elseif (nwater.EQ.0.0) then
				cff=0
				v=0
			endif
			
				
				
		  CASE("v")
			uvJpar=Jpar-0.5
			Jnode=floor(uvJpar)
			if (Jnode.lt.1) Jnode=1
			if (Jnode.gt.eta_v(ng)) Inode=eta_v(ng)
			Y=uvJpar-Jnode
			cff(1)=(1-X)*(1-Y)
			cff(2)=X*(1-Y)
			cff(3)=X*Y
			cff(4)=(1-X)*Y
			
			if (Jnode.GE.eta_v(ng)) then
				v(1,:) = HYDRODATA(ng)%V(Inode,eta_v(ng),:,t)
				v(2,:) = HYDRODATA(ng)%V(Inode+1,eta_v(ng),:,t)
				v(3,:) = HYDRODATA(ng)%V(Inode+1,eta_v(ng),:,t)	
				v(4,:) = HYDRODATA(ng)%V(Inode,eta_v(ng),:,t)
				
				
#ifdef WETDRY				
				m(1) = HYDRODATA(ng)%wetdry_mask_v(Inode,eta_v(ng),t)
				m(2) = HYDRODATA(ng)%wetdry_mask_v(Inode+1,eta_v(ng),t)
				m(3) = 0
				m(4) = 0
#else			
				m(1) = GRIDS(ng)%mask_v(Inode,eta_v(ng))
				m(2) = GRIDS(ng)%mask_v(Inode+1,eta_v(ng))
				m(3) = 0
				m(4) = 0
#endif			
			elseif (Jnode.LT.1) then
			
				v(1,:) = HYDRODATA(ng)%V(Inode,1,:,t)
				v(2,:) = HYDRODATA(ng)%V(Inode+1,1,:,t)
				v(3,:) = HYDRODATA(ng)%V(Inode+1,1,:,t)	
				v(4,:) = HYDRODATA(ng)%V(Inode,1,:,t)
#ifdef WETDRY				
				m(1) = 0
				m(2) = 0
				m(3) = HYDRODATA(ng)%wetdry_mask_v(Inode+1,1,t)
				m(4) = HYDRODATA(ng)%wetdry_mask_v(Inode,1,t)
#else				
				m(1) = 0
				m(2) = 0
				m(3) = GRIDS(ng)%mask_v(Inode+1,1)
				m(4) = GRIDS(ng)%mask_v(Inode,1)
#endif
			else
				v(1,:) = HYDRODATA(ng)%V(Inode,Jnode,:,t)
				v(2,:) = HYDRODATA(ng)%V(Inode+1,Jnode,:,t)
				v(3,:) = HYDRODATA(ng)%V(Inode+1,Jnode+1,:,t)	
				v(4,:) = HYDRODATA(ng)%V(Inode,Jnode+1,:,t)
#ifdef WETDRY	
				m(1) = HYDRODATA(ng)%wetdry_mask_v(Inode,Jnode,t)
				m(2) = HYDRODATA(ng)%wetdry_mask_v(Inode+1,Jnode,t)
				m(3) = HYDRODATA(ng)%wetdry_mask_v(Inode+1,Jnode+1,t)
				m(4) = HYDRODATA(ng)%wetdry_mask_v(Inode,Jnode+1,t)
#else
				m(1) = GRIDS(ng)%mask_v(Inode,Jnode)
				m(2) = GRIDS(ng)%mask_v(Inode+1,Jnode)
				m(3) = GRIDS(ng)%mask_v(Inode+1,Jnode+1)
				m(4) = GRIDS(ng)%mask_v(Inode,Jnode+1)
#endif
			endif
			nwater=m(1)+m(2)+m(3)+m(4)
			if ((nwater.LT.4.0).and.(nwater.GT.0.0)) then
				!call getcoeff(X,Y,m,cff)
				do i=1,nz
					if (m(1).eq.0.0d0) v(1,i)=0.0
					if (m(2).eq.0.0d0) v(2,i)=0.0
					if (m(3).eq.0.0d0) v(3,i)=0.0
					if (m(4).eq.0.0d0) v(4,i)=0.0
				enddo
			elseif (nwater.EQ.0.0) then
				cff=0
				v=0
			endif
			
		  CASE("salt")
			v(1,:) = HYDRODATA(ng)%salt(Inode,Jnode,:,t)
			v(2,:) = HYDRODATA(ng)%salt(Inode+1,Jnode,:,t)
			v(3,:) = HYDRODATA(ng)%salt(Inode+1,Jnode+1,:,t)	
			v(4,:) = HYDRODATA(ng)%salt(Inode,Jnode+1,:,t)	
#ifdef WETDRY	
			m(1) = HYDRODATA(ng)%wetdry_mask_rho(Inode,Jnode,t)
			m(2) = HYDRODATA(ng)%wetdry_mask_rho(Inode+1,Jnode,t)
			m(3) = HYDRODATA(ng)%wetdry_mask_rho(Inode+1,Jnode+1,t)
			m(4) = HYDRODATA(ng)%wetdry_mask_rho(Inode,Jnode+1,t)	
#else			
			m(1) = GRIDS(ng)%mask_rho(Inode,Jnode)
			m(2) = GRIDS(ng)%mask_rho(Inode+1,Jnode)
			m(3) = GRIDS(ng)%mask_rho(Inode+1,Jnode+1)
			m(4) = GRIDS(ng)%mask_rho(Inode,Jnode+1)
#endif
			nwater=m(1)+m(2)+m(3)+m(4)
			if (nwater.LT.4.0) then
				call getcoeff(X,Y,m,cff)
			endif
		  CASE("temp")
				
			v(1,:) = HYDRODATA(ng)%temp(Inode,Jnode,:,t)
			v(2,:) = HYDRODATA(ng)%temp(Inode+1,Jnode,:,t)
			v(3,:) = HYDRODATA(ng)%temp(Inode+1,Jnode+1,:,t)	
			v(4,:) = HYDRODATA(ng)%temp(Inode,Jnode+1,:,t)	
#ifdef WETDRY	
			m(1) = HYDRODATA(ng)%wetdry_mask_rho(Inode,Jnode,t)
			m(2) = HYDRODATA(ng)%wetdry_mask_rho(Inode+1,Jnode,t)
			m(3) = HYDRODATA(ng)%wetdry_mask_rho(Inode+1,Jnode+1,t)
			m(4) = HYDRODATA(ng)%wetdry_mask_rho(Inode,Jnode+1,t)	
#else			
			m(1) = GRIDS(ng)%mask_rho(Inode,Jnode)
			m(2) = GRIDS(ng)%mask_rho(Inode+1,Jnode)
			m(3) = GRIDS(ng)%mask_rho(Inode+1,Jnode+1)
			m(4) = GRIDS(ng)%mask_rho(Inode,Jnode+1)
#endif
			nwater=m(1)+m(2)+m(3)+m(4)
			if (nwater.LT.4.0) then
				call getcoeff(X,Y,m,cff)
			endif
		  CASE("AKs")
				
			v(1,:) = HYDRODATA(ng)%AKs(Inode,Jnode,:,t)
			v(2,:) = HYDRODATA(ng)%AKs(Inode+1,Jnode,:,t)
			v(3,:) = HYDRODATA(ng)%AKs(Inode+1,Jnode+1,:,t)		
			v(4,:) = HYDRODATA(ng)%AKs(Inode,Jnode+1,:,t)
#ifdef WETDRY	
			m(1) = HYDRODATA(ng)%wetdry_mask_rho(Inode,Jnode,t)
			m(2) = HYDRODATA(ng)%wetdry_mask_rho(Inode+1,Jnode,t)
			m(3) = HYDRODATA(ng)%wetdry_mask_rho(Inode+1,Jnode+1,t)
			m(4) = HYDRODATA(ng)%wetdry_mask_rho(Inode,Jnode+1,t)	
#else			
			m(1) = GRIDS(ng)%mask_rho(Inode,Jnode)
			m(2) = GRIDS(ng)%mask_rho(Inode+1,Jnode)
			m(3) = GRIDS(ng)%mask_rho(Inode+1,Jnode+1)
			m(4) = GRIDS(ng)%mask_rho(Inode,Jnode+1)
#endif
			nwater=m(1)+m(2)+m(3)+m(4)
			if (nwater.LT.4.0) then
				call getcoeff(X,Y,m,cff)
			endif
		  CASE("turbaccel")
				
			v(1,:) = HYDRODATA(ng)%Accelstd_t(Inode,Jnode,:,t)
			v(2,:) = HYDRODATA(ng)%Accelstd_t(Inode+1,Jnode,:,t)
			v(3,:) = HYDRODATA(ng)%Accelstd_t(Inode+1,Jnode+1,:,t)		
			v(4,:) = HYDRODATA(ng)%Accelstd_t(Inode,Jnode+1,:,t)
#ifdef WETDRY	
			m(1) = HYDRODATA(ng)%wetdry_mask_rho(Inode,Jnode,t)
			m(2) = HYDRODATA(ng)%wetdry_mask_rho(Inode+1,Jnode,t)
			m(3) = HYDRODATA(ng)%wetdry_mask_rho(Inode+1,Jnode+1,t)
			m(4) = HYDRODATA(ng)%wetdry_mask_rho(Inode,Jnode+1,t)	
#else			
			m(1) = GRIDS(ng)%mask_rho(Inode,Jnode)
			m(2) = GRIDS(ng)%mask_rho(Inode+1,Jnode)
			m(3) = GRIDS(ng)%mask_rho(Inode+1,Jnode+1)
			m(4) = GRIDS(ng)%mask_rho(Inode,Jnode+1)
#endif
			nwater=m(1)+m(2)+m(3)+m(4)
			if (nwater.LT.4.0) then
				call getcoeff(X,Y,m,cff)
			endif
		  CASE("turbvort")
				
			v(1,:) = HYDRODATA(ng)%Vortstd_t(Inode,Jnode,:,t)
			v(2,:) = HYDRODATA(ng)%Vortstd_t(Inode+1,Jnode,:,t)
			v(3,:) = HYDRODATA(ng)%Vortstd_t(Inode+1,Jnode+1,:,t)		
			v(4,:) = HYDRODATA(ng)%Vortstd_t(Inode,Jnode+1,:,t)
#ifdef WETDRY	
			m(1) = HYDRODATA(ng)%wetdry_mask_rho(Inode,Jnode,t)
			m(2) = HYDRODATA(ng)%wetdry_mask_rho(Inode+1,Jnode,t)
			m(3) = HYDRODATA(ng)%wetdry_mask_rho(Inode+1,Jnode+1,t)
			m(4) = HYDRODATA(ng)%wetdry_mask_rho(Inode,Jnode+1,t)	
#else			
			m(1) = GRIDS(ng)%mask_rho(Inode,Jnode)
			m(2) = GRIDS(ng)%mask_rho(Inode+1,Jnode)
			m(3) = GRIDS(ng)%mask_rho(Inode+1,Jnode+1)
			m(4) = GRIDS(ng)%mask_rho(Inode,Jnode+1)
#endif
			nwater=m(1)+m(2)+m(3)+m(4)
			if (nwater.LT.4.0) then
				call getcoeff(X,Y,m,cff)
			endif
					 
		CASE("waveaccelu")
				
			v(1,:) = HYDRODATA(ng)%Accelustd_w(Inode,Jnode,:,t)
			v(2,:) = HYDRODATA(ng)%Accelustd_w(Inode+1,Jnode,:,t)
			v(3,:) = HYDRODATA(ng)%Accelustd_w(Inode+1,Jnode+1,:,t)		
			v(4,:) = HYDRODATA(ng)%Accelustd_w(Inode,Jnode+1,:,t)
#ifdef WETDRY	
			m(1) = HYDRODATA(ng)%wetdry_mask_rho(Inode,Jnode,t)
			m(2) = HYDRODATA(ng)%wetdry_mask_rho(Inode+1,Jnode,t)
			m(3) = HYDRODATA(ng)%wetdry_mask_rho(Inode+1,Jnode+1,t)
			m(4) = HYDRODATA(ng)%wetdry_mask_rho(Inode,Jnode+1,t)	
#else			
			m(1) = GRIDS(ng)%mask_rho(Inode,Jnode)
			m(2) = GRIDS(ng)%mask_rho(Inode+1,Jnode)
			m(3) = GRIDS(ng)%mask_rho(Inode+1,Jnode+1)
			m(4) = GRIDS(ng)%mask_rho(Inode,Jnode+1)
#endif
			nwater=m(1)+m(2)+m(3)+m(4)
			if (nwater.LT.4.0) then
				call getcoeff(X,Y,m,cff)
			endif
			
		CASE("waveaccelv")
				
			v(1,:) = HYDRODATA(ng)%Accelvstd_w(Inode,Jnode,:,t)
			v(2,:) = HYDRODATA(ng)%Accelvstd_w(Inode+1,Jnode,:,t)
			v(3,:) = HYDRODATA(ng)%Accelvstd_w(Inode+1,Jnode+1,:,t)		
			v(4,:) = HYDRODATA(ng)%Accelvstd_w(Inode,Jnode+1,:,t)
#ifdef WETDRY	
			m(1) = HYDRODATA(ng)%wetdry_mask_rho(Inode,Jnode,t)
			m(2) = HYDRODATA(ng)%wetdry_mask_rho(Inode+1,Jnode,t)
			m(3) = HYDRODATA(ng)%wetdry_mask_rho(Inode+1,Jnode+1,t)
			m(4) = HYDRODATA(ng)%wetdry_mask_rho(Inode,Jnode+1,t)	
#else			
			m(1) = GRIDS(ng)%mask_rho(Inode,Jnode)
			m(2) = GRIDS(ng)%mask_rho(Inode+1,Jnode)
			m(3) = GRIDS(ng)%mask_rho(Inode+1,Jnode+1)
			m(4) = GRIDS(ng)%mask_rho(Inode,Jnode+1)
#endif
			nwater=m(1)+m(2)+m(3)+m(4)
			if (nwater.LT.4.0) then
				call getcoeff(X,Y,m,cff)
			endif
		  CASE("w")
				
			v(1,:) = HYDRODATA(ng)%W(Inode,Jnode,:,t)
			v(2,:) = HYDRODATA(ng)%W(Inode+1,Jnode,:,t)
			v(3,:) = HYDRODATA(ng)%W(Inode+1,Jnode+1,:,t)		
			v(4,:) = HYDRODATA(ng)%W(Inode,Jnode+1,:,t)
#ifdef WETDRY	
			m(1) = HYDRODATA(ng)%wetdry_mask_rho(Inode,Jnode,t)
			m(2) = HYDRODATA(ng)%wetdry_mask_rho(Inode+1,Jnode,t)
			m(3) = HYDRODATA(ng)%wetdry_mask_rho(Inode+1,Jnode+1,t)
			m(4) = HYDRODATA(ng)%wetdry_mask_rho(Inode,Jnode+1,t)	
#else			
			m(1) = GRIDS(ng)%mask_rho(Inode,Jnode)
			m(2) = GRIDS(ng)%mask_rho(Inode+1,Jnode)
			m(3) = GRIDS(ng)%mask_rho(Inode+1,Jnode+1)
			m(4) = GRIDS(ng)%mask_rho(Inode,Jnode+1)
#endif
			nwater=m(1)+m(2)+m(3)+m(4)
		
		
			if (nwater.LT.4.0) then
				call getcoeff(X,Y,m,cff)
			endif
			
			
			
			
		
		  CASE DEFAULT
			write(*,*) 'Problem interpolating ',var
			write(*,*) ' '
			write(*,*) 'The Program Cannot Continue and Will Terminate'
			stop
		END SELECT


		do i=1,nz
		  zprof(i)= cff(1)*v(1,i)+cff(2)*v(2,i)+cff(3)*v(3,i)+cff(4)*v(4,i)
		enddo
		getInterp3D=(zcff(1)*zprof(Ibot)+zcff(2)*zprof(Itop))

		
		DEALLOCATE(v)
		DEALLOCATE(zprof)

	  END FUNCTION getInterp3D
	  
	  DOUBLE PRECISION FUNCTION getInterpStr(var,ng,Ipar,Jpar,t)
		!This Function returns the interpolated value at the particle's location 
		!  using the interpolation variables stored from function setInterp, and
		!  the hydrodynamic variables that have been read in

		use grid_mod, only: GRIDS
		use hydro_mod, only: HYDRODATA
		use param_mod, only: xi_rho,eta_rho,xi_u,eta_u,xi_v,eta_v
		
		IMPLICIT NONE
		CHARACTER(LEN=*), INTENT(IN) :: var
		INTEGER, INTENT(IN) :: ng,t
		DOUBLE PRECISION, INTENT(IN) :: Ipar,Jpar
		INTEGER :: Inode,Jnode,i
		DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: v
		DOUBLE PRECISION :: m(4),cff(4)
		DOUBLE PRECISION :: X,Y,uvIpar,uvJpar
		DOUBLE PRECISION :: Ival,nwater
	    
		
		ALLOCATE(v(4))
		! write(*,*) '-----'
		 !write(*,"(I,F10.4,F10.4,F10.4,F10.4,F10.4,F10.4)") ztype,zlev(1),ZPar,zlev(2),zeta,zcff(1),zcff(2)
		! write(*,*) Itop,Ibot
		Inode=floor(Ipar)
		Jnode=floor(Jpar)
		X=Ipar-Inode
		Y=Jpar-Jnode
		

	!	 f(x,y)=f(0,0)(1-x)(1-y)+f(1,0)x(1-y)+f(0,1)(1-x)y+f(1,1)xy,}
		cff(1)=(1-X)*(1-Y)
	    cff(2)=X*(1-Y)
	    cff(3)=X*Y
	    cff(4)=(1-X)*Y
		!Determine which data to interpolate from
		
		SELECT CASE(var)
		  CASE("bustr")
			uvIpar=Ipar-0.5d0
			Inode=floor(uvIpar)
			
			
			
			X=uvIpar-Inode
			cff(1)=(1-X)*(1-Y)
			cff(2)=X*(1-Y)
			cff(3)=X*Y
			cff(4)=(1-X)*Y
			
		  
			
			if (Inode.GE.xi_u(ng)) then
				
				v(1) = HYDRODATA(ng)%bustr(xi_u(ng),Jnode,t)
				v(2) = HYDRODATA(ng)%bustr(xi_u(ng),Jnode,t)
				v(3) = HYDRODATA(ng)%bustr(xi_u(ng),Jnode+1,t)	
				v(4) = HYDRODATA(ng)%bustr(xi_u(ng),Jnode+1,t)	
#ifdef WETDRY
				m(1) = HYDRODATA(ng)%wetdry_mask_u(xi_u(ng),Jnode,t)
				m(2) = 0
				m(3) = 0
				m(4) = HYDRODATA(ng)%wetdry_mask_u(xi_u(ng),Jnode+1,t)
#else			
				m(1) = GRIDS(ng)%mask_u(xi_u(ng),Jnode)
				m(2) = 0
				m(3) = 0
				m(4) = GRIDS(ng)%mask_u(xi_u(ng),Jnode+1)	
#endif					
			
			elseif  (Inode.LT.1) then
				v(1) = HYDRODATA(ng)%bustr(1,Jnode,t)
				v(2) = HYDRODATA(ng)%bustr(1,Jnode,t)
				v(3) = HYDRODATA(ng)%bustr(1,Jnode+1,t)	
				v(4) = HYDRODATA(ng)%bustr(1,Jnode+1,t)	
#ifdef WETDRY		
				m(1) = 0
				m(2) = HYDRODATA(ng)%wetdry_mask_u(1,Jnode,t)
				m(3) = HYDRODATA(ng)%wetdry_mask_u(1,Jnode+1,t)
				m(4) = 0
#else	
				m(1) = 0
				m(2) = GRIDS(ng)%mask_u(1,Jnode)
				m(3) = GRIDS(ng)%mask_u(1,Jnode+1)
				m(4) = 0	
#endif			
			else
				v(1) = HYDRODATA(ng)%bustr(Inode,Jnode,t)
				v(2) = HYDRODATA(ng)%bustr(Inode+1,Jnode,t)
				v(3) = HYDRODATA(ng)%bustr(Inode+1,Jnode+1,t)	
				v(4) = HYDRODATA(ng)%bustr(Inode,Jnode+1,t)
#ifdef WETDRY	
		        m(1) = HYDRODATA(ng)%wetdry_mask_u(Inode,Jnode,t)
				m(2) = HYDRODATA(ng)%wetdry_mask_u(Inode+1,Jnode,t)
				m(3) = HYDRODATA(ng)%wetdry_mask_u(Inode+1,Jnode+1,t)
				m(4) = HYDRODATA(ng)%wetdry_mask_u(Inode,Jnode+1,t)
#else			
				m(1) = GRIDS(ng)%mask_u(Inode,Jnode)
				m(2) = GRIDS(ng)%mask_u(Inode+1,Jnode)
				m(3) = GRIDS(ng)%mask_u(Inode+1,Jnode+1)
				m(4) = GRIDS(ng)%mask_u(Inode,Jnode+1)
#endif		
			endif
			
			nwater=m(1)+m(2)+m(3)+m(4)
			
			if ((nwater.LT.4.0).and.(nwater.GT.0.0)) then
				!call getcoeff(X,Y,m,cff)
				
					if (m(1).eq.0.0d0) v(1)=0.0
					if (m(2).eq.0.0d0) v(2)=0.0
					if (m(3).eq.0.0d0) v(3)=0.0
					if (m(4).eq.0.0d0) v(4)=0.0
				
			elseif (nwater.EQ.0.0) then
				cff=0
				v=0
			endif
		
				
				
		  CASE("bvstr")
			uvJpar=Jpar-0.5
			Jnode=floor(uvJpar)
			if (Jnode.lt.1) Jnode=1
			if (Jnode.gt.eta_v(ng)) Inode=eta_v(ng)
			Y=uvJpar-Jnode
			cff(1)=(1-X)*(1-Y)
			cff(2)=X*(1-Y)
			cff(3)=X*Y
			cff(4)=(1-X)*Y
			
			if (Jnode.GE.eta_v(ng)) then
				v(1) = HYDRODATA(ng)%bvstr(Inode,eta_v(ng),t)
				v(2) = HYDRODATA(ng)%bvstr(Inode+1,eta_v(ng),t)
				v(3) = HYDRODATA(ng)%bvstr(Inode+1,eta_v(ng),t)	
				v(4) = HYDRODATA(ng)%bvstr(Inode,eta_v(ng),t)
				
				
#ifdef WETDRY				
				m(1) = HYDRODATA(ng)%wetdry_mask_v(Inode,eta_v(ng),t)
				m(2) = HYDRODATA(ng)%wetdry_mask_v(Inode+1,eta_v(ng),t)
				m(3) = 0
				m(4) = 0
#else			
				m(1) = GRIDS(ng)%mask_v(Inode,eta_v(ng))
				m(2) = GRIDS(ng)%mask_v(Inode+1,eta_v(ng))
				m(3) = 0
				m(4) = 0
#endif			
			elseif (Jnode.LT.1) then
			
				v(1) = HYDRODATA(ng)%bvstr(Inode,1,t)
				v(2) = HYDRODATA(ng)%bvstr(Inode+1,1,t)
				v(3) = HYDRODATA(ng)%bvstr(Inode+1,1,t)	
				v(4) = HYDRODATA(ng)%bvstr(Inode,1,t)
#ifdef WETDRY				
				m(1) = 0
				m(2) = 0
				m(3) = HYDRODATA(ng)%wetdry_mask_v(Inode+1,1,t)
				m(4) = HYDRODATA(ng)%wetdry_mask_v(Inode,1,t)
#else				
				m(1) = 0
				m(2) = 0
				m(3) = GRIDS(ng)%mask_v(Inode+1,1)
				m(4) = GRIDS(ng)%mask_v(Inode,1)
#endif
			else
				v(1) = HYDRODATA(ng)%bvstr(Inode,Jnode,t)
				v(2) = HYDRODATA(ng)%bvstr(Inode+1,Jnode,t)
				v(3) = HYDRODATA(ng)%bvstr(Inode+1,Jnode+1,t)	
				v(4) = HYDRODATA(ng)%bvstr(Inode,Jnode+1,t)
#ifdef WETDRY	
				m(1) = HYDRODATA(ng)%wetdry_mask_v(Inode,Jnode,t)
				m(2) = HYDRODATA(ng)%wetdry_mask_v(Inode+1,Jnode,t)
				m(3) = HYDRODATA(ng)%wetdry_mask_v(Inode+1,Jnode+1,t)
				m(4) = HYDRODATA(ng)%wetdry_mask_v(Inode,Jnode+1,t)
#else
				m(1) = GRIDS(ng)%mask_v(Inode,Jnode)
				m(2) = GRIDS(ng)%mask_v(Inode+1,Jnode)
				m(3) = GRIDS(ng)%mask_v(Inode+1,Jnode+1)
				m(4) = GRIDS(ng)%mask_v(Inode,Jnode+1)
#endif
			endif
			nwater=m(1)+m(2)+m(3)+m(4)
			if ((nwater.LT.4.0).and.(nwater.GT.0.0)) then
				!call getcoeff(X,Y,m,cff)
					if (m(1).eq.0.0d0) v(1)=0.0
					if (m(2).eq.0.0d0) v(2)=0.0
					if (m(3).eq.0.0d0) v(3)=0.0
					if (m(4).eq.0.0d0) v(4)=0.0
			elseif (nwater.EQ.0.0) then
				cff=0
				v=0
			endif
			
		
		  CASE DEFAULT
			write(*,*) 'Problem interpolating ',var
			write(*,*) ' '
			write(*,*) 'The Program Cannot Continue and Will Terminate'
			stop
		END SELECT


		getInterpStr=cff(1)*v(1)+cff(2)*v(2)+cff(3)*v(3)+cff(4)*v(4)

		
		DEALLOCATE(v)

	  END FUNCTION getInterpStr
	subroutine getInterpAKs(ng,Ipar,Jpar,t,i,zprof)
		!This Function returns the interpolated value at the particle's location 
		!  using the interpolation variables stored from function setInterp, and
		!  the hydrodynamic variables that have been read in

		use grid_mod, only: GRIDS
		use hydro_mod, only: HYDRODATA
		use param_mod, only: s_rho,s_w,xi_rho,eta_rho,xi_u,eta_u,xi_v,eta_v
		
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: ng,t,i
		DOUBLE PRECISION, INTENT(IN) :: Ipar,Jpar
		INTEGER :: Inode,Jnode,nz
		DOUBLE PRECISION, INTENT(OUT) :: zprof
		DOUBLE PRECISION :: m(4),cff(4),v(4)
		DOUBLE PRECISION :: X,Y
		DOUBLE PRECISION :: Ival,zlev(2),nwater

			  

		nz=s_w(ng)


		Inode=floor(Ipar)
		Jnode=floor(Jpar)
		X=Ipar-Inode
		Y=Jpar-Jnode
		

	!	 f(x,y)=f(0,0)(1-x)(1-y)+f(1,0)x(1-y)+f(0,1)(1-x)y+f(1,1)xy,}
		cff(1)=(1-X)*(1-Y)
	    cff(2)=X*(1-Y)
	    cff(3)=X*Y
	    cff(4)=(1-X)*Y

		v(1) = HYDRODATA(ng)%AKs(Inode,Jnode,i,t)
		v(2) = HYDRODATA(ng)%AKs(Inode+1,Jnode,i,t)
		v(3) = HYDRODATA(ng)%AKs(Inode+1,Jnode+1,i,t)		
		v(4) = HYDRODATA(ng)%AKs(Inode,Jnode+1,i,t)
		m(1) = GRIDS(ng)%mask_rho(Inode,Jnode)
		m(2) = GRIDS(ng)%mask_rho(Inode+1,Jnode)
		m(3) = GRIDS(ng)%mask_rho(Inode+1,Jnode+1)
		m(4) = GRIDS(ng)%mask_rho(Inode,Jnode+1)
		nwater=m(1)+m(2)+m(3)+m(4)
		if (nwater.LT.4.0) then
				call getcoeff(X,Y,m,cff)
		endif

		zprof= cff(1)*v(1)+cff(2)*v(2)+cff(3)*v(3)+cff(4)*v(4)
	




	  END subroutine getInterpAKs
 SUBROUTINE getcoeff(X,Y,m,cff)
		!This Function GETS coefficients
		use grid_mod, only: GRIDS
		use hydro_mod, only: HYDRODATA
		IMPLICIT NONE
		DOUBLE PRECISION, INTENT(IN) :: X,Y
		INTEGER :: Inode,Jnode,t,i
		DOUBLE PRECISION, INTENT(IN) :: m(4)
		DOUBLE PRECISION, INTENT(INOUT) :: cff(4)
		DOUBLE PRECISION :: dist(4),W(4),WT
		
		
		cff=0.0
		dist(1)=sqrt(X**2+Y**2)
		W(1)=m(1)/dist(1)
		dist(2)=sqrt((1-X)**2+Y**2)
		W(2)=m(2)/dist(2)
		dist(3)=sqrt((1-X)**2+(1-Y)**2)
		W(3)=m(3)/dist(3)
		dist(4)=sqrt(X**2+(1-Y)**2)
		W(4)=m(4)/dist(4)
	    WT=W(1)+W(2)+W(3)+W(4)
		
		
		if ( WT.EQ.0.0D0) then
			do i=1,4 
				cff(i)=0.0D0
			enddo		
		else
			do i=1,4 
				cff(i)=W(i)/WT
			enddo
		endif
		

		



	  END SUBROUTINE getcoeff
END MODULE INT_MOD