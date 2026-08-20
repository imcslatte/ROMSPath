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
  PUBLIC :: bounds,zbounds,IsBeachPossible,WhetherToBeach,Unbeach

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


SUBROUTINE IsBeachPossible(newXpos,newYpos,ng,coord_close,node_dists,m,InBeachWindow)
use grid_mod,   only: GRIDS
use param_mod, only: xi_rho,eta_rho,beachDistWindow,nBeachCells
use int_mod, only: getInterp2D
DOUBLE PRECISION, INTENT(IN) :: newXpos,newYpos
INTEGER, INTENT(IN) :: ng
INTEGER :: I,J,icurr,jcurr,layerInd,centInd,ind_close(2,2),xi_ind,eta_ind,xind,yind,currsum
DOUBLE PRECISION :: nLon,nLat,bdist,s,currdist,tLat1,tLon1,tLat2,tLon2,dist2beach,parLon,parLat
DOUBLE PRECISION, INTENT(OUT) :: coord_close(2,2),node_dists(2) 
INTEGER, INTENT(OUT) :: m(2*nBeachCells,2*nBeachCells)
LOGICAL, INTENT(OUT) :: InBeachWindow

InBeachWindow = .FALSE.
node_dists = 9999.
coord_close = 99.

I = floor(newXpos)
J = floor(newYpos)
m = 99

do icurr = (I-nBeachCells+1),(I+nBeachCells)
        do jcurr = (J-nBeachCells+1),(J+nBeachCells)
                  if ((icurr.LT.2) .or. (jcurr.LT.2) .or. (icurr.GT.(xi_rho(ng)-1))) then !not a valid cell reference, but not land for current grid
                        m(icurr-(I-nBeachCells),jcurr-(J-nBeachCells)) = 1
                  elseif (jcurr.GT.eta_rho(ng)) then !land, with current grid
                  !ASSUMPTION: all particles that beach where there isn't a mask beach at the large end of the eta coordinates. 
                  !asummption repeated farther down
                  m(icurr-(I-nBeachCells),jcurr-(J-nBeachCells)) = 0
                  else !check mask to get m value
                        m(icurr-(I-nBeachCells),jcurr-(J-nBeachCells)) = INT(GRIDS(ng)%mask_rho(icurr,jcurr))  
                  endif    
        enddo
enddo 


if (SUM(m) < ((2*nBeachCells)**2)) then !near boundary, should check distance from beach

        !ID 2 closest boundary cells
        ! first, check 4 surrounding. if any of these are 0, then InBeachWindow
        ! (because beachDist must greater than any cell dimension)
        centInd = nBeachCells
        currsum = SUM(m(centInd:(centInd+1),centInd:(centInd+1)))
        if ((I.LT.2) .OR. (I.GT.(xi_rho(ng)-1))) then !might be able to reduce this now that I have I instead of icurr here (eg, I.LT.1)
                InBeachWindow = .FALSE. !essentially, do nothing because particle is in escape window
        elseif ((J+1).GT.eta_rho(ng)) then
                InBeachWindow = .TRUE.
                node_dists = 999.
        elseif (currsum<4) then
                InBeachWindow = .TRUE.
        else !now check cells that are farther out
                parLon = getInterp2D("lon",ng,newXpos,newYpos,1)
                parLat = getInterp2D("lat",ng,newXpos,newYpos,1)	
                CALL ID2Closest(ng,m,parLon,parLat,I,J,coord_close,node_dists)
                
                ! interpolate closest distance         
                if ((node_dists(1)<9998.) .AND. (node_dists(2)>9998.)) then
                        dist2beach = node_dists(1)
                        if (dist2beach<beachDistWindow) then
                                InBeachWindow = .TRUE.
                        endif
                elseif ((node_dists(1)>9998.) .AND. (node_dists(2)<9998.)) then
                        dist2beach = node_dists(2)
                        if (dist2beach<beachDistWindow) then
                                InBeachWindow = .TRUE.
                        endif
                elseif ((node_dists(1)<9998.) .AND. (node_dists(2)<9998.)) then
                        bdist = GCircle(coord_close(1,1),coord_close(1,2),coord_close(2,1),coord_close(2,2))
                        s = .5*(bdist+node_dists(1)+node_dists(2))
                        dist2beach = (2./bdist)*SQRT(s*(s-node_dists(1))*(s-node_dists(2))*(s-bdist))
                        if (dist2beach < beachDistWindow) then
                                InBeachWindow = .TRUE.
                        endif
                endif
        endif

endif


END SUBROUTINE IsBeachPossible


LOGICAL FUNCTION WhetherToBeach()
use param_mod, only: BeachTimescale,idt
DOUBLE PRECISION :: pb, rand_num

pb = 1.-EXP(-DBLE(idt)/(BeachTimescale*24.*60.*60.))
call RANDOM_NUMBER(rand_num)

if (rand_num<pb) then
        WhetherToBeach = .TRUE.
else
        WhetherToBeach = .FALSE.
endif


END FUNCTION WhetherToBeach


SUBROUTINE Unbeach(ng,oldXPos,oldYPos,newXPos,newYPos,coord_close,node_dists,m,ubXPos,ubYPos)
use int_mod, only: getInterp2D,LL2ij
use param_mod, only: xi_rho,eta_rho,nBeachCells,Ngrid
use grid_mod,   only: GRIDS
INTEGER, INTENT(IN) :: ng
DOUBLE PRECISION, INTENT(IN) :: oldXPos,oldYPos,NewXPos,NewYPos
DOUBLE PRECISION, INTENT(INOUT) :: coord_close(2,2),node_dists(2)
INTEGER, INTENT (INOUT)       :: m(2*nBeachCells,2*nBeachCells)
INTEGER                       :: I,J,icurr,jcurr,nmask,xi_ind(2),eta_ind(2),centInd
DOUBLE PRECISION              :: parLonNew,parLatNew,parLonOld,parLatOld,dLonL,dLatL,dLonP,dLatP
DOUBLE PRECISION              :: alpha,mL,mP,theta,dp,dProj,N1(2),N2(2),d1,d2,NLon(1),NLat(1),tempubXpos(1),tempubYpos(1)
DOUBLE PRECISION              :: currsum,km2londeg,km2latdeg,useLon,useLat,test
LOGICAL                       :: ingrid,obound
DOUBLE PRECISION, INTENT(OUT) :: ubXPos,ubYPos

if (.NOT.((node_dists(1)>998.).AND.(node_dists(1)<9998.).AND.(node_dists(2)>998.).AND.(node_dists(2)<9998.))) then
parLonOld = getInterp2D("lon",ng,oldXpos,oldYpos,1)
parLatOld = getInterp2D("lat",ng,oldXpos,oldYpos,1)	
parLonNew = getInterp2D("lon",ng,newXpos,newYpos,1)
parLatNew = getInterp2D("lat",ng,newXpos,newYpos,1)	
endif

!Useful if in beaching window and closest coordinates weren't already ID'ed
!Determine whether there are 2 close coordinates already. Get them if possible
!Leave at old position if not
if ((node_dists(1)>9998.).OR.(node_dists(2)>9998.)) then

        node_dists = 9999. !redoing!
        I = floor(newXpos)
        J = floor(newYpos)
        useLon = parLonNew
        useLat = parLatNew

        centInd = nBeachCells
        currsum = SUM(m(centInd:(centInd+1),centInd:(centInd+1)))
        
        if ((currsum.EQ.0).OR.(currsum.EQ.4)) then !retry from old position if fully within mask or ID2Closest already fully ran and didn't work

        !Try to find two closest from old position
                I = floor(oldXpos)
                J = floor(oldYpos)
                useLon = parLonOld
                useLat = parLatOld

                m = 99
        
                do icurr = (I-nBeachCells+1),(I+nBeachCells)
                         do jcurr = (J-nBeachCells+1),(J+nBeachCells)
                                if ((icurr.LT.2) .or. (jcurr.LT.2) .or. (icurr.GT.(xi_rho(ng)-1))) then !not a valid cell reference, but not land for current grid
                                        m(icurr-(I-nBeachCells),jcurr-(J-nBeachCells)) = 1
                                elseif (jcurr.GT.eta_rho(ng)) then !land, with current grid
                                !ASSUMPTION: all particles that beach where there isn't a mask beach at the large end of the eta coordinates. 
                                m(icurr-(I-nBeachCells),jcurr-(J-nBeachCells)) = 0
                                else !check mask to get m value
                                        m(icurr-(I-nBeachCells),jcurr-(J-nBeachCells)) = INT(GRIDS(ng)%mask_rho(icurr,jcurr))  
                                endif    
                        enddo
                enddo 

                currsum = SUM(m(centInd:(centInd+1),centInd:(centInd+1)))
        endif

        if (currsum.EQ.0) then !old position is still fully within mask
                WRITE(*,*) "Error: Old 'unbeached' position is fully within the mask"
        elseif (currsum.EQ.1) then !3 land nodes are in inner ring, take two adjacent to ocean node
                if ((m(centInd,centInd).EQ.1).OR.(m(centInd+1,centInd+1).EQ.1)) then
                        xi_ind(1) = I
                        eta_ind(1) = J+1
                        xi_ind(2) = I+1
                        eta_ind(2) = J 
                else 
                        xi_ind(1) = I
                        eta_ind(1) = J
                        xi_ind(2) = I+1
                        eta_ind(2) = J+1
                endif
                

                coord_close(1,1) = GRIDS(ng)%lon_rho(xi_ind(1),eta_ind(1))
                coord_close(1,2) = GRIDS(ng)%lat_rho(xi_ind(1),eta_ind(1))
                
                coord_close(2,1) = GRIDS(ng)%lon_rho(xi_ind(2),eta_ind(2))
                coord_close(2,2) = GRIDS(ng)%lat_rho(xi_ind(2),eta_ind(2))

                node_dists = 1. !just to record that they were found
        elseif (currsum.EQ.2) then !2 land nodes are in inner ring, these are the two closest
                xi_ind =  999999.
                eta_ind = 999999.
                if (m(centInd,centInd).EQ.0) then
                        xi_ind(1) = I
                        eta_ind(1) = J
                endif
                if (m(centInd+1,centInd).EQ.0) then
                        if (xi_ind(1) > 99999.) then
                                xi_ind(1) = I+1
                                eta_ind(1) = J
                        else
                                xi_ind(2) = I+1
                                eta_ind(2) = J
                        endif
                endif
                if (m(centInd,centInd+1).EQ.0) then
                        if (xi_ind(1) > 99999.) then
                                xi_ind(1) = I
                                eta_ind(1) = J+1
                        else
                                xi_ind(2) = I
                                eta_ind(2) = J+1
                        endif
                endif
                if(m(centInd+1,centInd+1).EQ.0) then
                        xi_ind(2) = I+1
                        eta_ind(2) = J+1
                endif
  
                coord_close(1,1) = GRIDS(ng)%lon_rho(xi_ind(1),eta_ind(1))
                coord_close(1,2) = GRIDS(ng)%lat_rho(xi_ind(1),eta_ind(1))
                
                coord_close(2,1) = GRIDS(ng)%lon_rho(xi_ind(2),eta_ind(2))
                coord_close(2,2) = GRIDS(ng)%lat_rho(xi_ind(2),eta_ind(2))

                node_dists = 1. !just to record they were found 
        elseif (currsum.EQ.3) then !1 "land" node is in inner ring, find two closest in next ring
                CALL ID2Closest(ng,m,useLon,useLat,I,J,coord_close,node_dists)
        elseif (SUM(m) < ((2*nBeachCells)**2)) then !near boundary, do check distance from beach
                !ID 2 closest boundary cells
                CALL ID2Closest(ng,m,useLon,useLat,I,J,coord_close,node_dists)
        else
                WRITE(*,*) "Error: not registering as near a boundary"
        endif
endif

!Again, check whether there are 2 close coordinates. If not, unbeach by leaving
!at old location
if ((node_dists(1)>9998.).OR.(node_dists(2)>9998.)) then
        ubXPos = oldXPos
        ubYPos = oldYPos
        !WRITE(*,*) "Particle left at old position because there were not two close coordinates"
elseif ((node_dists(1)>998.).OR.(node_dists(1)>998.)) then
        ubXPos = oldXPos
        ubYPos = oldYPos
        !WRITE(*,*) "Particle left at old position because new position was outside of grid"
else !move alongshore!

        !get slope of line connecting land points
        dLonL = GCircle(coord_close(1,1),coord_close(1,2),coord_close(2,1),coord_close(1,2))
        if (coord_close(2,1)<coord_close(1,1)) then
                dLonL = -dLonL
        endif
        dLatL = GCircle(coord_close(1,1),coord_close(1,2),coord_close(1,1),coord_close(2,2))
        if (coord_close(2,2)<coord_close(1,2)) then
                dLatL = -dLatL
        endif
        mL = dLatL/dLonL


        !get slope of line connecting particle positions
        dLonP = GCircle(parLonOld,parLatOld,parLonNew,parLatOld)
        if (parLonNew<parLonOld) then
                dLonP = -dLonP
        endif
        dLatP = GCircle(parLonOld,parLatOld,parLonOld,parLatNew)
        if (parLatNew<parLatOld) then
                dLatP = -dLatP
        endif
        mP = dLatP/dLonP

        !calculate angles
        alpha = ATAN(mL) !angle between beach and cartesian
        theta = ATAN( (mL-mP)/(1+mL*mP) ) !angle between particle trajectory and beach

        !calculate projected distance
        dp = GCircle(parLonOld,parLatOld,parLonNew,parLatNew)
        dproj = dp*COS(theta)

        !get km to degree conversions
        km2londeg = ABS((coord_close(2,1)-coord_close(1,1))/dLonL)
        km2latdeg = ABS((coord_close(2,2)-coord_close(1,2))/dLatL)

        !get two possible points
        N1(1) = parLonOld + km2londeg*dproj*COS(alpha)
        N1(2) = parLatOld + km2latdeg*dproj*SIN(alpha)  
        N2(1) = parLonOld - km2londeg*dproj*COS(alpha)
        N2(2) = parLatOld - km2latdeg*dproj*SIN(alpha)  

        !check which point is closer to the new particle location
        d1 = gcircle(N1(1),N1(2),parLonNew,parLatNew)
        d2 = gcircle(N2(1),N2(2),parLonNew,parLatNew)
        if (d1<d2) then !N1 is where the particle should unbeach to
                NLon(1) = N1(1)
                NLat(1) = N1(2)
        else !N2 is where the particle should unbeach to
                NLon(1) = N2(1)
                NLat(1) = N2(2)
        endif

                call ll2ij(grids(ng)%lon_rho,grids(ng)%lat_rho,grids(ng)%angle,NLon,NLat,	&
				1,xi_rho(ng),eta_rho(ng),tempubXPos,tempubYPos)
                ubXpos = tempubXPos(1)
                ubYpos = tempubYPos(1)
	!make sure new position is in bounds
        !note that this is not adapted for nested grids
	!do loopind=1,Ngrid
        !        tempX=GRIDS(int(par(n,pGID)))%scl(loopind,1)*ubXPos+GRIDS(int(par(n,pGID)))%off(loopind,1)
        !        tempY=GRIDS(int(par(n,pGID)))%scl(loopind,2)*ubYPos+GRIDS(int(par(n,pGID)))%off(loopind,2)
        	call bounds(1,ubXpos,ubYpos,nmask,ingrid,obound)
       !end do
        

        !if not in bounds, use old position instead
        if (.NOT.ingrid) then
                ubXPos = oldXPos
                ubYPos = oldYPos
        endif

endif


END SUBROUTINE Unbeach


SUBROUTINE ID2Closest(ng,m,parLon,parLat,I,J,coord_close,node_dists)
use grid_mod,   only: GRIDS
use param_mod, only: xi_rho,eta_rho,nBeachCells
use int_mod, only: getInterp2D
DOUBLE PRECISION, INTENT(IN) :: parLon,parLat
INTEGER, INTENT(IN) :: ng,m(2*nBeachCells,2*nBeachCells),I,J
INTEGER :: icurr,jcurr,layerInd,centInd,ind_close(2,2),xi_ind,eta_ind,currsum
DOUBLE PRECISION :: nLon,nLat,bdist,s,currdist,tLat1,tLon1,tLat2,tLon2,dist2beach
DOUBLE PRECISION, INTENT(INOUT) :: coord_close(2,2),node_dists(2) 

                centInd = nBeachCells

                do layerInd = 1,nBeachCells
                    currsum = SUM(m((centInd-layerInd+1):(centInd+layerInd),(centInd-layerInd+1):(centInd+layerInd)))    
                    if (currsum < (layerInd*2)*(layerInd*2)) then
                                                      
                        !compute distance from particle to any 0 in subarray
                        do icurr = centInd-layerInd+1,centInd+layerInd
                                do jcurr = centInd-layerInd+1,centInd+layerInd
                                        !only calculate for outer perimeter
                                        if ((icurr.EQ.(centInd-layerInd+1)).OR.(icurr.EQ.(centInd+layerInd)).OR.(jcurr.EQ.(centInd-layerInd+1)).OR.(jcurr.EQ.centInd+layerInd)) then
                                                if (m(icurr,jcurr)<1) then
                                                        xi_ind = I-centInd+icurr
                                                        eta_ind = J-centInd+jcurr
                                                        if (eta_ind.GT.eta_rho(ng)) then !check whether in grid. if not,extrapolate lat/lon
                                                               if ((eta_ind-1).GT.eta_rho(ng)) then 
                                                                        WRITE(*,*) "Error in IsBeachPossible: tried to index outside of mask"
                                                                endif
                                                                tLon1 = GRIDS(ng)%lon_rho(xi_ind,eta_ind-2)
                                                                tLat1 = GRIDS(ng)%lat_rho(xi_ind,eta_ind-2)
                                                                tLon2 = GRIDS(ng)%lon_rho(xi_ind,eta_ind-1)
                                                                tLat2 = GRIDS(ng)%lat_rho(xi_ind,eta_ind-1)
                                                                nLon = tLon2 + (tLon2-tLon1)
                                                                nLat = tLat2 + (tLat2-tLat1)
                                                        else
                                                                !get lat lon coordinates of point
                                                                nLon = GRIDS(ng)%lon_rho(xi_ind,eta_ind)
                                                                nLat = GRIDS(ng)%lat_rho(xi_ind,eta_ind)
                                                                !estimate distance between particle and those coordinates
                                                        endif
                                                                currdist = GCircle(parLon,parLat,nLon,nLat) 
                                                        !hold onto the two closest distances and associated i/j 
                                                        if (node_dists(1)>node_dists(2)) then
                                                                if (currdist<node_dists(1)) then
                                                                        node_dists(1) = currdist
                                                                        ind_close(1,1) = xi_ind
                                                                        ind_close(1,2) = eta_ind
                                                                        coord_close(1,1) = nLon
                                                                        coord_close(1,2) = nLat
                                                                endif
                                                        else
                                                                if (currdist<node_dists(2)) then
                                                                        node_dists(2) = currdist
                                                                        ind_close(2,1) = xi_ind
                                                                        ind_close(2,2) = eta_ind
                                                                        coord_close(2,1) = nLon
                                                                        coord_close(2,2) = nLat
                                                                endif
                                                        endif
                                                endif
                                        endif
                                enddo
                        enddo
                       
                        if ((node_dists(1) < 999.).AND.(node_dists(2)<999.)) then 
                                EXIT !exit do loop because we've found the two closest
                        endif
                    endif
                enddo

END SUBROUTINE ID2Closest

DOUBLE PRECISION FUNCTION GCircle(lon1,lat1,lon2,lat2)
DOUBLE PRECISION, INTENT(IN) :: lon1,lat1,lon2,lat2
DOUBLE PRECISION :: radius,PI, deg2rad,rad2deg,slon,slat,elon,elat,alpha
INTEGER :: ind

PI = 4*ATAN(1.d0)
radius = 6371.315
deg2rad = PI/180.
rad2deg = 180./PI

slon = lon1*deg2rad
slat = lat1*deg2rad
elon = lon2*deg2rad
elat = lat2*deg2rad

alpha = sin(slat)*sin(elat)+cos(slat)*cos(elat)*cos(elon-slon)

if (alpha>1.) then !this corresponds to 0 distance
        alpha = 1.
elseif (alpha<-1.) then !this corresponds to halfway around the world
        alpha = -1.
endif

alpha = ACOS(alpha)

GCircle = radius*alpha

END FUNCTION GCircle

!  ! This subroutine calculates the intersection between the particle
!  ! trajectory and the boundary line in a grid cell, and then calculates
!  ! the reflection, returning the new particle location
!   subroutine intersect_reflect(Xpos,Ypos,nXpos,nYpos,fintersectX,fintersectY,  &
!     freflectX,freflectY,intersectf,skipbound,isWater)
!     IMPLICIT NONE
!     INTEGER, INTENT(OUT) :: intersectf
!     INTEGER, INTENT(INOUT) :: skipbound
!     DOUBLE PRECISION, INTENT(IN) :: Xpos,Ypos,nXpos,nYpos
!     DOUBLE PRECISION, INTENT(OUT) :: fintersectX,fintersectY,freflectX,freflectY
!     LOGICAL, OPTIONAL, INTENT(OUT) :: isWater
!     INTEGER :: i,intersect,skipboundi
!     DOUBLE PRECISION :: crossk,dPBC,mBCperp,rx1,rx2,ry1,ry2,Bp,distBC,dist1,   &
!     dist2,intersctx,interscty,rPxyzX,rPxyzY,Mbc,Bbc,Mp,bcx1,bcy1,bcx2,bcy2,  &
!      bBCperp,xhigh,xlow,yhigh,ylow,d_Pinter,dtest,bxhigh,bxlow,byhigh,bylow
!
!     distBC=0.0
!     Mbc = 0.0
!     Bbc = 0.0
!     Mp = 0.0
!     Bp = 0.0
!     intersect=0
!     intersectf=0
!     skipboundi = skipbound
!     fintersectX = -999999.
!     fintersectY = -999999.
!     freflectX = -999999.
!     freflectY = -999999.
!     dtest = 999999.
!     isWater = .FALSE.
!
!     if (Xpos.GE.nXpos) then
!      xhigh = Xpos
!      xlow = nXpos
!     else
!      xhigh = nXpos
!      xlow = Xpos
!     endif  
!
!     if (Ypos.GE.nYpos) then
!      yhigh = Ypos
!      ylow = nYpos
!     else
!      yhigh = nYpos
!      ylow = Ypos
!     endif
!
!     do i=1,nbounds
!
!     if (i == skipbound) cycle
!
!        intersect = 0
!        bcx1=bnd_x(1,i)
!        bcy1=bnd_y(1,i)
!        bcx2=bnd_x(2,i)
!        bcy2=bnd_y(2,i)
!
!        !If the boundary segment end points are both east, west, north, or 
!        !  south of the particle's previous or new location, cycle to next 
!        !  boundary
!        if( ((bcx1 > xhigh) .AND. (bcx2 > xhigh)) .OR. &
!           ((bcx1 < xlow ) .AND. (bcx2 < xlow )) .OR. &
!           ((bcy1 > yhigh) .AND. (bcy2 > yhigh)) .OR. &
!           ((bcy1 < ylow ) .AND. (bcy2 < ylow ))      ) cycle
!       
!        if (bcx1.GE.bcx2) then
!         bxhigh = bcx1
!         bxlow = bcx2
!        else
!         bxhigh = bcx2
!         bxlow = bcx1
!        endif  
!
!        if (bcy1.GE.bcy2) then
!         byhigh = bcy1
!         bylow = bcy2
!        else
!         byhigh = bcy2
!         bylow = bcy1
!        endif
!
!        !First determine if an undefined denominator is possible
!        if (bcx1.EQ.bcx2 .OR. nXpos.EQ.Xpos ) then
!         !test if they both vertical, if so cycle because they cannot intersect
!         if (bcx1.EQ.bcx2 .AND. nXpos.EQ.Xpos ) cycle
!         !test if perpendicular and parrallel to coordinate axes
!         if (bcx1.EQ.bcx2 .AND. nYpos.EQ.Ypos ) then
!          !undefined denominator, perp. & || to axes
!          intersctx = bcx1
!          interscty = nYpos
!          if (intersctx.LE.xhigh  .AND. intersctx.GE.xlow .AND.              &
!           interscty.LE.yhigh  .AND. interscty.GE.ylow .AND.                &
!           intersctx.LE.bxhigh .AND. intersctx.GE.bxlow .AND.               &
!           interscty.LE.byhigh .AND. interscty.GE.bylow  ) then
!           dPBC=sqrt((intersctx-nXpos)**2+(interscty-nYpos)**2)
!           rx1=nXpos+(DBLE(2.0)*dPBC)
!           ry1=nYpos
!           rx2=nXpos-(DBLE(2.0)*dPBC)
!           ry2=nYpos
!           dist1 = sqrt( (intersctx-rx1)**2 + (interscty-ry1)**2 )
!           dist2 = sqrt( (intersctx-rx2)**2 + (interscty-ry2)**2 )
!           if(dist1.LT.dist2) then
!            rPxyzX= rx1
!            rPxyzY= ry1
!           elseif(dist1.GT.dist2) then
!            rPxyzX= rx2
!            rPxyzY= ry2
!           endif
!           intersect=1
!          endif
!         elseif (nXpos.EQ.Xpos .AND. bcy1.EQ.bcy2 ) then
!          !undefined denominator, perp. & || to axes
!          intersctx = nXpos
!          interscty = bcy1
!          if (intersctx.LE.xhigh .AND.  intersctx.GE.xlow .AND.              &
!           interscty.LE.yhigh  .AND. interscty.GE.ylow .AND.                &
!           intersctx.LE.bxhigh .AND. intersctx.GE.bxlow .AND.               &
!           interscty.LE.byhigh .AND. interscty.GE.bylow  ) then
!           dPBC=sqrt((intersctx-nXpos)**2+(interscty-nYpos)**2)
!           rx1=nXpos
!           ry1=nYpos+(DBLE(2.0)*dPBC)
!           rx2=nXpos
!           ry2=nYpos-(DBLE(2.0)*dPBC)
!           dist1 = sqrt( (intersctx-rx1)**2 + (interscty-ry1)**2 )
!           dist2 = sqrt( (intersctx-rx2)**2 + (interscty-ry2)**2 )
!           if(dist1.LT.dist2) then
!            rPxyzX= rx1
!            rPxyzY= ry1
!           elseif(dist1.GT.dist2) then
!            rPxyzX= rx2
!            rPxyzY= ry2
!           endif
!           intersect=1
!          endif
!         elseif (bcx1.EQ.bcx2 .AND. nYpos.NE.Ypos ) then
!          !undefined denominator, not perpendicular
!          Mp = (nYpos-Ypos)/(nXpos-Xpos)
!          Bp = Ypos - Mp*Xpos
!          intersctx = bcx1
!          interscty = Mp*intersctx + Bp
!          if (intersctx.LE.xhigh .AND.  intersctx.GE.xlow .AND.              &
!             interscty.LE.yhigh  .AND. interscty.GE.ylow .AND.              &
!             intersctx.LE.bxhigh .AND. intersctx.GE.bxlow .AND.             &
!             interscty.LE.byhigh .AND. interscty.GE.bylow  ) then
!           dPBC = nXpos-intersctx
!           rx1=nXpos+(DBLE(2.0)*dPBC)
!           ry1=nYpos
!           rx2=nXpos-(DBLE(2.0)*dPBC)
!           ry2=nYpos
!           dist1 = sqrt( (intersctx-rx1)**2 + (interscty-ry1)**2 )
!           dist2 = sqrt( (intersctx-rx2)**2 + (interscty-ry2)**2 )
!           if(dist1.LT.dist2) then
!            rPxyzX= rx1
!            rPxyzY= ry1
!           elseif(dist1.GT.dist2) then
!            rPxyzX= rx2
!            rPxyzY= ry2
!           endif
!           intersect=1
!          endif
!         elseif (nXpos.EQ.Xpos .AND. bcy1.NE.bcy2  ) then
!          !undefined denominator, not perpendicular
!          Mbc = (bcy2-bcy1)/(bcx2-bcx1)
!          Bbc = bcy2 - Mbc*bcx2
!          intersctx = nXpos
!          interscty = Mbc*intersctx + Bbc
!          if (intersctx.LE.xhigh .AND.  intersctx.GE.xlow .AND.              &
!             interscty.LE.yhigh  .AND. interscty.GE.ylow .AND.              &
!             intersctx.LE.bxhigh .AND. intersctx.GE.bxlow .AND.             &
!             interscty.LE.byhigh .AND. interscty.GE.bylow ) then
!           !Now use cross product to determine the distance of the particle 
!           !  from the boundary
!           distBC = sqrt((bcx1-bcx2)**2+(bcy1-bcy2)**2)
!           crossk= ((nXpos-bcx1)*(bcy2-bcy1)) - ((bcx2-bcx1)*(nYpos-bcy1))
!           dPBC = sqrt(crossk**2)/distBC
!           !find line perpendicular to boundary
!           mBCperp = DBLE(-1.0)/Mbc
!           bBCperp = nYpos - mBCperp*nXpos
!           !find two potential reflection points
!           rx1 = sqrt( ((DBLE(2.0)*dPBC)**2)/(DBLE(1.0)+mBCperp**2) ) +nXpos
!           ry1 = mBCperp*rx1 + bBCperp
!           rx2 = sqrt( ((DBLE(2.0)*dPBC)**2)/(DBLE(1.0)+mBCperp**2) )       &
!              * DBLE(-1.0) + nXpos
!           ry2 = mBCperp*rx2 + bBCperp
!           !point closest to intersection of boundary and particle trajectory
!           !  is the right one
!           dist1 = sqrt( (intersctx-rx1)**2 + (interscty-ry1)**2 )
!           dist2 = sqrt( (intersctx-rx2)**2 + (interscty-ry2)**2 )
!           if(dist1.LT.dist2) then
!            rPxyzX= rx1
!            rPxyzY= ry1
!           elseif(dist1.GT.dist2) then
!            rPxyzX= rx2
!            rPxyzY= ry2
!           endif
!           intersect=1
!          endif
!         endif
!        else
!
!         if(intersect == 0)then
!
!           Mbc = (bcy2-bcy1)/(bcx2-bcx1)
!           Bbc = bcy2 - Mbc*bcx2
!           Mp = (nYpos-Ypos)/(nXpos-Xpos)
!           Bp = Ypos - Mp*Xpos
!           intersctx = (Bbc - Bp)/(Mp - Mbc)
!           interscty = Mp*intersctx + Bp
!
!           !when bc parallel with x-axis, byhigh=bylow=intersecty
!           if (Mbc.EQ.0.0) interscty = byhigh
!       
!           if (intersctx.LE.xhigh .AND.  intersctx.GE.xlow .AND.              &
!              interscty.LE.yhigh  .AND. interscty.GE.ylow .AND.              &
!              intersctx.LE.bxhigh .AND. intersctx.GE.bxlow .AND.             &
!              interscty.LE.byhigh .AND. interscty.GE.bylow  ) then
!
!            if (Mbc.EQ.0.0) then  !inverse slope denominator not OK
!              dPBC = nYpos-bcy1
!              rx1=nXpos
!              ry1=nYpos+(DBLE(2.0)*dPBC)
!              rx2=nXpos
!              ry2=nYpos-(DBLE(2.0)*dPBC)
!              dist1 = sqrt( (intersctx-rx1)**2 + (interscty-ry1)**2 )
!              dist2 = sqrt( (intersctx-rx2)**2 + (interscty-ry2)**2 ) 
!              if(dist1.LT.dist2) then
!               rPxyzX= rx1
!               rPxyzY= ry1
!              elseif(dist1.GT.dist2) then
!               rPxyzX= rx2
!               rPxyzY= ry2
!              endif
!               intersect=1
!              endif  
!
!            if(intersect == 0)then
!
!              !Now use cross product to determine the distance of the
!              !  particle from the boundary
!              distBC = sqrt((bcx1-bcx2)**2+(bcy1-bcy2)**2)
!              crossk= ((nXpos-bcx1)*(bcy2-bcy1)) - ((bcx2-bcx1)*(nYpos-bcy1))
!              dPBC = sqrt(crossk**2)/distBC
!              !find line perpendicular to boundary
!              mBCperp = DBLE(-1.0)/Mbc
!              bBCperp = nYpos - mBCperp*nXpos
!              !find two potential reflection points
!              rx1 = sqrt(((DBLE(2.0)*dPBC)**2)/(DBLE(1.0)+mBCperp**2)) +nXpos
!              ry1 = mBCperp*rx1 + bBCperp
!              rx2 = sqrt(((DBLE(2.0)*dPBC)**2)/(DBLE(1.0)+mBCperp**2))       &
!                 * DBLE(-1.0) + nXpos
!              ry2 = mBCperp*rx2 + bBCperp
!              !point closest to intersection of boundary and particle 
!              !  trajectory is the right one
!              dist1 = sqrt( (intersctx-rx1)**2 + (interscty-ry1)**2 )
!              dist2 = sqrt( (intersctx-rx2)**2 + (interscty-ry2)**2 )
!              if(dist1.LT.dist2) then
!               rPxyzX= rx1
!               rPxyzY= ry1
!              elseif(dist1.GT.dist2) then
!               rPxyzX= rx2
!               rPxyzY= ry2
!              endif
!              intersect=1
!             
!            endif
!           endif  
!         endif
!        endif
!
!
!        d_Pinter = sqrt( (Xpos-intersctx)**2 + (Ypos-interscty)**2 )
!        if( (intersect .EQ. 1) .AND. (d_Pinter .LT. dtest) ) then
!         fintersectX = intersctx
!         fintersectY = interscty
!         freflectX = rPxyzX
!         freflectY = rPxyzY
!         intersectf = 1
!         dtest = d_Pinter
!         skipboundi = i
!         isWater = .NOT. land(i)
!        endif
!
!     enddo
!
!     skipbound = skipboundi
!   END SUBROUTINE intersect_reflect




END MODULE
