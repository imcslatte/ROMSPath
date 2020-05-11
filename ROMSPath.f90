! ROMSPath - Larval TRANSport Lagrangian model v.4                                
! Date: 5 March 2019
!
! Description: The Lagrangian TRANSport model (ROMSPath) is an 
! off-line particle-tracking model that runs with the stored predictions of
! a 3D hydrodynamic model, specifically the Regional Ocean Modeling System 
! (ROMS). Although ROMSPath was built to simulate oyster larvae, it can  
! be adapted to simulate passive particles and other planktonic organisms. 
! ROMSPath is written in Fortran 90 and is designed to track the trajectories 
! of particles in three dimensions. It includes a 4th order Runge-Kutta scheme 
! for particle advection and a random displacement model for vertical turbulent
! particle motion. Reflective boundary conditions, larval behavior, and 
! settlement routines are also included. Components of ROMSPath have been in 
! development since 2002 and are described in the following publications:
! North et al. 2004, North et al. 2006a, North et al. 2006b, 
! North et al. 2008, North et al. 2011, Schlag and North 2012.
!
! Developers:
!   Elizabeth North: enorth@umces.edu
!   Zachary Schlag: zschlag@umces.edu
!   Ian Mitchell: imitchell@umces.edu
!   Elias Hunter: hunter@marine.rutgers.edu
!
!   University of Maryland
!   Center for Envir. Science
!   Horn Point Laboratory
!   Cambridge, MD 21613 USA
!
! Funding was provided by the National Science Foundation Biological 
! and Physical Oceanography Programs, Maryland Department of Natural 
! Resources, NOAA Chesapeake Bay Office, NOAA Maryland Sea Grant College 
! Program, & NOAA-funded UMCP Advanced Study Institute for the Environment. 
! 
! **********************************************************************
! **********************************************************************
! **                      Copyright (c) 2013                          **
! **   The University of Maryland Center for Environmental Science    **
! **********************************************************************
! **                                                                  **
! ** This Software is open-source and licensed under the following    **
! ** conditions as stated by MIT/X License:                           **
! **                                                                  **
! **  (See http://www.opensource.org/licenses/mit-license.php ).      **
! **                                                                  **
! ** Permission is hereby granted, free of charge, to any person      **
! ** obtaining a copy of this Software and associated documentation   **
! ** files (the "Software"), to deal in the Software without          **
! ** restriction, including without limitation the rights to use,     **
! ** copy, modify, merge, publish, distribute, sublicense,            **
! ** and/or sell copies of the Software, and to permit persons        **
! ** to whom the Software is furnished to do so, subject to the       **
! ** following conditions:                                            **
! **                                                                  **
! ** The above copyright notice and this permission notice shall      **
! ** be included in all copies or substantial portions of the         **
! ** Software.                                                        **
! **                                                                  **
! ** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,  **
! ** EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE           **
! ** WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE  **
! ** AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT  **
! ** HOLDERS BE LIABLE FOR ANY CLAIMS, DAMAGES OR OTHER LIABILITIES,  **
! ** WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     **
! ** FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR    **
! ** OTHER DEALINGS IN THE SOFTWARE.                                  **
! **                                                                  **
! ** The most current official versions of this Software and          **
! ** associated tools and documentation are available at:             **
! **                                                                  **
! **  http://northweb.hpl.umces.edu/ROMSPath.htm                        **
! **                                                                  **
! ** We ask that users make appropriate acknowledgement of            **
! ** The University of Maryland Center for Environmental Science,     **
! ** individual developers, participating agencies and institutions,  **
! ** and funding agencies. One way to do this is to cite one or       **
! ** more of the relevant publications listed at:                     **
! **                                                                  **
! **  http://northweb.hpl.umces.edu/ROMSPath.htm#Description            **
! **                                                                  **
! **********************************************************************
! ********************************************************************** 

PROGRAM main

! ROMSPath.f90 contains the main structure of the particle-tracking program. 
! It executes the external time step, internal time step, and particle loops, 
! advects particles, and writes output. It calls modules that read in 
! hydrodynamic model information, move particles due to turbulence and 
! behavior, test if particles are in habitat polygons, and apply boundary 
! conditions to keep particles in the model domain. 
!
! Program created by:   Elizabeth North
! Modified by:          Elias Hunter
! Created on:           2004
! Last Modified on:     8 August 2019

IMPLICIT NONE
!   *************************************************************************
!   *                                                                       *
!   *                       Variable Declarations                           *
!   *                                                                       *
!   *************************************************************************

  INTEGER, PARAMETER :: nAttrib   = 20

  INTEGER, PARAMETER :: pX        =  1  ! Particle X-coordinate
  INTEGER, PARAMETER :: pY        =  2  ! Particle Y-coordinate
  INTEGER, PARAMETER :: pZ        =  3  ! Particle Z-coordinate
  INTEGER, PARAMETER :: pnX       =  4  ! Particle new X-coordinate
  INTEGER, PARAMETER :: pnY       =  5  ! Particle new Y-coordinate
  INTEGER, PARAMETER :: pnZ       =  6  ! Particle new Z-coordinate
  INTEGER, PARAMETER :: ppX       =  7  ! Particle previous X-coordinate
  INTEGER, PARAMETER :: ppY       =  8  ! Particle previous Y-coordinate
  INTEGER, PARAMETER :: ppZ       =  9  ! Particle previous Z-coordinate
  INTEGER, PARAMETER :: pStatus   = 10  ! Status of particle (previously Color)
  INTEGER, PARAMETER :: pDOB      = 11  ! Particle Date Of Birth
  INTEGER, PARAMETER :: pAge      = 12  ! Particle Age (s)
  INTEGER, PARAMETER :: pLifespan = 13  ! Age at which particle settled or died
  INTEGER, PARAMETER :: pGID 	  = 14  ! Current grid ID. 
  INTEGER, PARAMETER :: pSize 	  = 15  ! Current particle size
  INTEGER, PARAMETER :: pAcc 	  = 16  ! Modeled PArticle acceleration
  INTEGER, PARAMETER :: pVort 	  = 17  ! Modeled particle vorticity
  INTEGER, PARAMETER :: pbehaveW  = 18  ! MOdeled particle bahavioral velocity
  INTEGER, PARAMETER :: pSSF 	  = 19  ! Modeled Sink/Swim flag
  INTEGER, PARAMETER :: pWD 	  = 20  ! Water Depth

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: par
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: P_Salt,P_Temp,mean_salt,mean_temp
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: P_HOB,P_bustr,P_bvstr
  INTEGER, ALLOCATABLE, DIMENSION(:) :: startpoly,endpoly,hitBottom,hitLand
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: isIn
  DOUBLE PRECISION :: ex(3),ix(3)
  INTEGER :: printdt,ets,its,mI
  REAL :: timeCounts(8),times(9)
 
  INTEGER :: NCcount,NCstart,prcount
!   *************************************************************************
!   *                                                                       *
!   *                             Execution                                 *
!   *                                                                       *
!   *************************************************************************

  call run_ROMSPath()

contains



  subroutine run_ROMSPath()
    ! *************************************************************************
    ! *                                                                       *
    ! *                              Run Model                                *
    ! *                                                                       *
    ! *************************************************************************
    use param_mod, only: days,dt,tdim
    integer :: seconds,stepT
	DOUBLE PRECISION :: before,after,tdiff 	
    call ini_ROMSPath()

    write(*,'(/,A)') '****** BEGIN ITERATIONS *******'
	


      ! days*24*60*60 = total number of seconds to run the model
      ! divide that by dt to get the number of external time steps
      seconds = int(days*86400.0) !Total seconds to run model
      stepT   = seconds/dt        !number of external time steps

      do ets=1,stepT
			call run_External_Timestep()
      enddo

    call fin_ROMSPath()

  end subroutine run_ROMSPath



  subroutine ini_ROMSPath()
    ! *************************************************************************
    ! *                                                                       *
    ! *                           Initialize Model                            *
    ! *                                                                       *
    ! *************************************************************************
!    use behavior_mod, only: initBehave,setOut,die
    use boundary_mod, only: bounds,zbounds
    USE INT_MOD,    ONLY: getInterp2D,getInterp3D
    use random_mod,   only: init_genrand,init_random_seed
    use hydro_mod,    only: updatehydro

    use param_mod,    only: numpar,days,dt,idt,seed,parfile,settlementon,   &
                      Behavior,TrackCollisions,SaltTempOn,Ngrid,SaltTempMean,    &
					  xi_rho,eta_rho,t_b,t_c,t_f,tstep,initsize,WriteBottom,    &
                      WriteHeaders,WriteModelTiming,ErrorFlag,getParams,Behavior
	use grid_mod,   only: InitGrid,GRIDS
	use INT_MOD,   only: LL2ij,inside
    integer :: n,istat,ng,i,j,nmask
    logical :: ingrid,obound,inzgrid
    double precision, allocatable, dimension(:) :: pLon,pLat,Ipar,Jpar
	CHARACTER(len=200) :: filenm

    integer :: in_island,inbounds,test
    double precision:: tdepth,zeta,dpran
	

	
	t_b = 1    !Back step is 1st time step in arrays
	t_c = 2    !Center step is 2nd time step in arrays
	t_f = 3    !Forward step is 3rd time step in arrays
	tstep=1
  ! ***************************************************************************
  ! *                          Get Parameter Values                           *
  ! ***************************************************************************

    CALL getParams()
	
	CALL InitGrid()
	
    CALL writeModelInfo()

    write(*,*) ' '
    write(*,*) ' *************** ROMSPath INITIALIZATION ************** '
#ifdef WETDRY
       write(*,*) 'WETDRY ACTIVE'
#endif	
#ifdef GROWTH
       write(*,*) 'GROWTH ACTIVE'
#endif
#ifdef STOKES
       write(*,*) 'STOKES ACTIVE'
#endif
  ! ! ***************************************************************************
  ! ! *                       Allocate Dynamic Variables                        *
  ! ! ***************************************************************************

     ALLOCATE(par(numpar,nAttrib))

    IF(SettlementOn)THEN
      ALLOCATE(startpoly(numpar))
      ALLOCATE(endpoly(numpar))
      endpoly = 0               !initialize end polygon location to zero
    ENDIF

    IF(SaltTempOn)THEN
      ALLOCATE(P_Salt(numpar))
      ALLOCATE(P_Temp(numpar))
      P_Salt = 0.0
      P_Temp = 0.0
	   IF(SaltTempMean)THEN
		ALLOCATE(mean_salt(numpar))
		ALLOCATE(mean_temp(numpar))
		mean_salt = 0.0
		mean_temp = 0.0
	   ENDIF
    ENDIF

	IF(WriteBottom)THEN
      ALLOCATE(P_HOB(numpar))
      ALLOCATE(P_bustr(numpar))
      ALLOCATE(P_bvstr(numpar))
      P_HOB = 0.0
      P_bustr = 0.0
      P_bvstr = 0.0
    ENDIF
	
    IF(TrackCollisions)THEN
      ALLOCATE(hitBottom(numpar))
      ALLOCATE(hitLand(numpar))
      hitBottom = 0
      hitLand = 0
    ENDIF

    ! !Local variables for read-in of Latitude and Longitude
     ALLOCATE(pLon(numpar))
     ALLOCATE(pLat(numpar))
     ALLOCATE(isIn(numpar))
     ALLOCATE(Ipar(numpar))
     ALLOCATE(Jpar(numpar))

    ! ! *************************************************************************
    ! ! *         Initialize print counters and random number generator         *
    ! ! *************************************************************************

    ! ! THE FOLLOWING VARIABLE INITIALIZATIONS SHOULD NOT BE CHANGED:
     prcount=0                  !print counter; number of external time steps
     printdt=0                  !print counter
	 
    !set random random Seed Value (how inception is that)
	IF (seed .EQ. 0) THEN
	 call init_random_seed(seed)
	ENDIF

    CALL init_genrand(seed)!set random number generator Seed Value
	
    ! ! *************************************************************************
    ! ! *                                                                       *
    ! ! *                    Initialize Hydrodynamic data         		        *
    ! ! *                                                                       *
    ! ! *************************************************************************

    !! THIS IS TEMPORARY 
    IF (Behavior.eq.10) THEN
		write(*,*) "BEHAVIOR TYPE 10 IS UNAVAILABLE, STOPPING"
		STOP
    ENDIF
	
	call	updateHydro(.TRUE.,1,t_b)
	call	updateHydro(.FALSE.,2,t_c)
	call	updateHydro(.FALSE.,3,t_f)

    ! ! *************************************************************************
    ! ! *                   Initialize Particle Attributes                      *
    ! ! *************************************************************************

    ! ! Read-in lat/long of particles. If settlement module is on, read in    
    ! ! the habitat polygon on which the particle start                       
    write(*,*) 'read in particle locations', numpar

    OPEN(1,FILE=TRIM(parfile))

      do n=1,numpar
	  
        par(n,pDOB)      = -9.0
        if(settlementon)then
          read (1,*) pLon(n),pLat(n),par(n,pZ),par(n,pDOB),startpoly(n)
        else
          read (1,*) pLon(n),pLat(n),par(n,pZ),par(n,pDOB)
        endif
		
        par(n,pX)  = 1.0
        par(n,pY)  = 1.0
        par(n,pnX) = 1.0
        par(n,pnY) = 1.0
        par(n,pnZ) = 1.0
        par(n,ppX) = 1.0    !initialize to 0.0 to indicate no previous location
        par(n,ppY) = 1.0    !initialize to 0.0 to indicate no previous location
        par(n,ppZ) = 1.0    !initialize to 0.0 to indicate no previous location
        par(n,pStatus)   = 0.0
        par(n,pAge)      = 0.0
        par(n,pLifespan) = 0.0
        par(n,pGID)      = dble(Ngrid)
        par(n,pSize)      = initsize
        par(n,pAcc) = 0.0
        par(n,pVort) = 0.0
        par(n,pbehaveW) = 0.0
        par(n,pSSF) = 0.0
		
		Ipar(n)=0.0
		Jpar(n)=0.0
		IF(WriteBottom)THEN
		P_HOB(n)=9999.0
		ENDIF
		isIn(n)=.False.
      enddo

    CLOSE(1)

	
	
	
	
	write(*,*) '*********'
	
	do ng =1,Ngrid

		
		call LL2ij(GRIDS(ng)%lon_rho,GRIDS(ng)%lat_rho,GRIDS(ng)%angle,Plon,Plat,	&
				numpar,xi_rho(ng),eta_rho(ng),Ipar,Jpar)

	
		do n=1,numpar
				call bounds(ng,Ipar(n),Jpar(n),nmask,ingrid,obound)
	
				
				! par(n,pnX)=Ipar(n)
				! par(n,pnY)=Jpar(n)
				! par(n,ppX)=Ipar(n)
				! par(n,ppY)=Jpar(n)
				! if ((ingrid).and.(par(n,pGID).eq.dble(Ngrid))) then

	
				if (.not.(isIn(n))) then
						
					if (ingrid) then
						par(n,pX)=Ipar(n)
						par(n,pY)=Jpar(n)
						par(n,pGID)=dble(ng)
						par(n,pStatus)=0.0
						isIn(n)=ingrid	
						call zbounds(ng,Ipar(n),Jpar(n),par(n,pZ),inzgrid,t_b)
						if (inzgrid) then
						else
							par(n,pStatus)=9.0
						endif
					endif
				endif 
				

		enddo
		
		if (SaltTempOn) then
			do n=1,numpar
				if (isIn(n)) then
				
					tdepth = DBLE(-1.0)* getInterp2D("depth",int(par(n,pGID)),par(n,pX),par(n,pY),t_c)
					zeta =  getInterp2D("zeta",int(par(n,pGID)),par(n,pX),par(n,pY),t_c)
					 P_salt(n)=getInterp3d("salt",int(par(n,pGID)),par(n,pX),par(n,pY),par(n,pZ),t_c,1,zeta,tdepth)
					 P_temp(n)=getInterp3d("temp",int(par(n,pGID)),par(n,pX),par(n,pY),par(n,pZ),t_c,1,zeta,tdepth)
					 par(n,pWD)=(DBLE(-1.0)*tdepth)+zeta
					 if (SaltTempMean) then
						mean_salt(n)=P_salt(n)
						mean_temp(n)=P_temp(n)
						mI=1
					 endif
					 
				 endif
		
			enddo
		endif
	 
		  ! if (ng.eq.1) then
			! temp1=GRIDS(ng)%scl(ng+1,1)*Ipar(n)+GRIDS(ng)%off(ng+1,1)
			! temp2=GRIDS(ng)%scl(ng+1,2)*Jpar(n)+GRIDS(ng)%off(ng+1,2)
		  ! else
		  
			! temp1=GRIDS(ng)%scl(ng-1,1)*Ipar(n)+GRIDS(ng)%off(ng-1,1)
			! temp2=GRIDS(ng)%scl(ng-1,2)*Jpar(n)+GRIDS(ng)%off(ng-1,2)
		  ! endif
		  ! write(10,"(F20.10,F20.10,F20.10,F20.10,F20.10,F20.10)") Plon(n),Plat(n),Ipar(n),Jpar(n),temp1,temp2
		 ! enddo
		 ! close(10)
		 
		 ! WRITE(filenm,'(A,I3.3)') 'GRID_',ng		
		
		 ! OPEN(10,FILE=TRIM(filenm),STATUS='REPLACE')
		 ! do i=1,xi_rho(ng)
			 ! do j=1,eta_rho(ng)step
				! if (ng.eq.1) then
						! temp1=GRIDS(ng)%scl(ng+1,1)*GRIDS(ng)%xi(i,j)+GRIDS(ng)%off(ng+1,1)
						! temp2=GRIDS(ng)%scl(ng+1,2)*GRIDS(ng)%eta(i,j)+GRIDS(ng)%off(ng+1,2)
				! else
		  
						! temp1=GRIDS(ng)%scl(ng-1,1)*GRIDS(ng)%xi(i,j)+GRIDS(ng)%off(ng-1,1)
						! temp2=GRIDS(ng)%scl(ng-1,2)*GRIDS(ng)%eta(i,j)+GRIDS(ng)%off(ng-1,2)
				! endif
			 
				 ! write(10,"(F20.10,F20.10,I4,F20.10,F20.10)") GRIDS(ng)%xi(i,j),GRIDS(ng)%eta(i,j) ,GRIDS(ng)%mask_rho(i,j),temp1,temp2
			 ! enddo
	    ! 
		! close(10)
    enddo
	
	

	write(*,*) '*********'
	write(*,*) '  Particle n=1 Latitude=',pLat(1),'Longitude=',pLon(1)
    write(*,*) '  Particle n=1 Depth=',par(1,pZ)
    write(*,*) '  Particle n=1 X=',par(1,pX),'Y=',par(1,pY)
    ! if(settlementon) write(*,*) '  Particle n=5 Start Polygon=',startpoly(5)
  
    ! ! *******************************************************************
    ! ! *                    Initialize NetCDF Output                     *
    ! ! *******************************************************************

     !Create NetCDF Output File
 
     CALL initNetCDF()
     CALL createNetCDF(par(:,pDOB))
   

    prcount = 0
	call writeNetCDF(0,pLon,pLat)
      ! if (SaltTempOn) then
        ! if(TrackCollisions)then
          ! CALL writeNetCDF(0,par(:,pAge),pLon,pLat,par(:,pZ),par(:,pStatus),par(:,pGid),   &
                          ! SALT=P_Salt,TEMP=P_Temp,HITB=hitBottom,HITL=hitLand)
        ! else
          ! CALL writeNetCDF(0,par(:,pAge),pLon,pLat,par(:,pZ),par(:,pStatus),par(:,pGid),   &
                          ! SALT=P_Salt,TEMP=P_Temp)
        ! endif
      ! else
        ! if(TrackCollisions)then
          ! CALL writeNetCDF(0,par(:,pAge),pLon,pLat,par(:,pZ),par(:,pStatus),par(:,pGid),   &
                          ! HITB=hitBottom,HITL=hitLand)
        ! else
          ! CALL writeNetCDF(0,par(:,pAge),pLon,pLat,par(:,pZ),par(:,pStatus),par(:,pGid))
        ! endif
      ! endif





    ! !Initialize Behavior
    ! CALL initBehave()

    


    ! !Create files to output 'land hits' and 'bottom hits'
    ! IF(TrackCollisions) then
      ! OPEN(100,FILE='LandHits.csv',STATUS='REPLACE')
        ! write(100,*)'numpar,lon,lat,depth,age,time,hitLand'
      ! CLOSE(100)
      ! OPEN(101,FILE='BottomHits.csv',STATUS='REPLACE')
        ! write(101,*)'numpar,lon,lat,depth,age,time,hitBottom'
      ! CLOSE(101)
    ! ENDIF

  
    ! !Create file to track model timing
    ! IF(WriteModelTiming)then
      ! OPEN(300,FILE='Timing.csv',STATUS='REPLACE')
        ! write(300,*)'Daytime,Elapsed,Hydro,Hydro %,Set Up,Set Up %,',     &
        ! 'Advection,Advection %,HTurb,HTurb %,VTurb,VTurb %,Behavior,',    &
        ! 'Behavior %,Update,Update %'
      ! CLOSE(300)

    ! ENDIF

   

      ! IF(TrackCollisions)then
        ! OPEN(100,FILE='LandHits_Headers.txt',STATUS='REPLACE')
          ! write(100,*)'column 01: numpar  -Particle identification number ',   &
                      ! '(dimensionless)'
          ! write(100,*)'column 02: lon     -Longitude of particle at end of ',  &
                      ! 'time step (decimal °)'
          ! write(100,*)'column 03: lat     -Latitude  of particle at end of ',  &
                      ! 'time step (decimal °)'
          ! write(100,*)'column 04: depth   -Depth of particle at end of time ', &
                      ! 'step (meters)'
          ! write(100,*)'column 05: age     -Age of particle (in days since ',   &
                      ! 'released)'
          ! write(100,*)'column 06: time    -Model time (in days since the ',    &
                      ! 'start of the model)'
          ! write(100,*)'column 07: hitLand -Number of times the particle ',     &
                      ! 'struck land in the last print interval time step'
        ! CLOSE(100)

        ! OPEN(101,FILE='BottomHits_Headers.txt',STATUS='REPLACE')
          ! write(101,*)'column 01: numpar    -Particle identification number ', &
                      ! '(dimensionless)'
          ! write(101,*)'column 02: lon       -Longitude of particle at end ',   &
                      ! 'of time step (decimal °)'
          ! write(101,*)'column 03: lat       -Latitude  of particle at end ',   &
                      ! 'of time step (decimal °)'
          ! write(101,*)'column 04: depth     -Depth of particle at end of ',    &
                      ! 'time step (meters)'
          ! write(101,*)'column 05: age       -Age of particle (in days since ', &
                      ! 'released)'
          ! write(101,*)'column 06: time      -Model time (in days since the ',  &
                      ! 'start of the model)'
          ! write(101,*)'column 07: hitBottom -Number of times the particle ',   &
                      ! 'struck bottom in the last print interval time step'
        ! CLOSE(101)
      ! ENDIF
    ! ENDIF

    ! !Deallocate local variables
     DEALLOCATE(pLon,pLat)

    ! !Output time spent initializing model
    ! call CPU_TIME(times(1))
    ! write(*,'("Time to initialize model = ",f6.3," seconds.")') times(1)
    ! timeCounts = 0.0      !Initialize time counters to 0.0

  end subroutine ini_ROMSPath

  

  subroutine run_External_Timestep()
    use param_mod, only: dt,idt,WriteModelTiming,tdim,t_b,t_c,t_f,filenum,&
		numdigits,prefix,suffix,tstep,multifile
    use hydro_mod, only: updateHydro,HYDRODATA
	integer :: stepIT,ng
	
	real :: before,after,tdiff,ibefore,iafter
	ng=1
	if (tstep.gt.tdim(ng))	then
		tstep=1
		filenum=filenum+1.0		
		if (.NOT.multifile) then
			write(*,*) 'NOT ENOUGH TIMES STEP TO CONTINUE:'
			STOP
		endif
	endif
    stepIT  = int(dt/idt)                     !number of internal time steps

      !IF(WriteModelTiming) call CPU_TIME(before)
	  
      ! write(*,*) '---'
	  ! call CPU_TIME(before)
      !Read in hydrodynamic model data 
      IF(ets > 2) then			
	     t_b  = mod(t_b,3)+1  ! 1 -> 2 -> 3 -> 1
	     t_c  = mod(t_c,3)+1  ! 2 -> 3 -> 1 -> 2
	     t_f = mod(t_f,3)+1  ! 3 -> 1 -> 2 -> 3
		CALL updateHydro(.FALSE.,tstep,t_f)   !do not start updating until 3rd iteration
	  endif
	  ! call CPU_TIME(after)
	  ! tdiff=after-before
	  
	  ! write(*,*) tdiff
	  
	  
      !Prepare external time step values to be used for 
      !  calculating Advection and Turbulence
      ex=0.0
      ex(1) = (ets-2)*dt
      ex(2) = (ets-1)*dt
      ex(3) = ets*dt
	 ! call CPU_TIME(before)
	 call CPU_TIME(before)

		
       do its=1,stepIT
		
		! call CPU_TIME(ibefore)
        call run_Internal_Timestep()
	   ! call CPU_TIME(iafter)  
		! write(*,*) iafter-ibefore		
       enddo !ITloop
	  tstep=tstep+1
	 ! call CPU_TIME(after)  
	 ! tdiff=after-before
	 ! write(*,*) 'TIMES'
	 ! write(*,"(F10.2 ,F10.2 ,F10.2 ,F10.2 ,F10.2 ,F10.2 ,F10.2)") tdiff,timeCounts(1),timeCounts(2),timeCounts(3),timeCounts(4),timeCounts(5),timeCounts(6)
	 timeCounts=0.0
  end subroutine run_External_Timestep



  subroutine run_Internal_Timestep()
    use param_mod, only: idt,iPrint

	
    !  calculating Advection and Turbulence
    ix(1) = ex(2) + DBLE((its-2)*idt)
    ix(2) = ex(2) + DBLE((its-1)*idt)
    ix(3) = ex(2) + DBLE(its*idt)
    
    !********************************************************
    !*                    Particle Loop                     *
    !********************************************************
	
    call update_particles()
    !********************************************************
    !*                 PRINT OUTPUT TO FILE                 *
    !********************************************************
    mI=mI+1
    printdt=printdt+idt
	
    if(printdt.GE.iprint) then
      write(*,*) 'write output to file, day = ',(DBLE(ix(3))/DBLE(86400))
      !ix(3)/86400 = (current model time in seconds) /
      !              (# of seconds in a day)

      call dataOutput()

      printdt=0  !reset print counter
    endif

  end subroutine run_Internal_Timestep  



   subroutine fin_ROMSPath()
    ! use param_mod, only: numpar,outpathGiven,outpath,settlementon
    ! use behavior_mod, only: finBehave,getStatus
    ! use convert_mod, only: x2lon,y2lat
    ! use hydro_mod, only: finHydro

    ! !OUTPUT ENDFILE NAME CONSTRUCTION VARIABLE
    ! CHARACTER(LEN=100) :: efile

     integer :: n,d,h,m
     real :: fintime,s
    ! double precision :: pLon,pLat

    

    ! !DEALLOCATE LOCAL VARIABLES
    ! DEALLOCATE(par)
    ! IF(ALLOCATED(hitBottom)) DEALLOCATE(hitBottom)
    ! IF(ALLOCATED(startpoly)) DEALLOCATE(startpoly)
    ! IF(ALLOCATED(endpoly  )) DEALLOCATE(endpoly)
    ! IF(ALLOCATED(hitLand  )) DEALLOCATE(hitLand)
    ! IF(ALLOCATED(P_Salt   )) DEALLOCATE(P_Salt)
    ! IF(ALLOCATED(P_Temp   )) DEALLOCATE(P_Temp)

    ! !DEALLOCATE MODULE VARIABLES
    ! call finBehave()
    ! call finHydro()

    !Calculate model run time and output to screen before exiting
    call CPU_TIME(fintime)
    d = int(fintime/86400.0)             !# of full days that the model ran
    h = int((fintime - real(d*86400))/3600.0)       !# of hours   (minus days)
    m = int((fintime - real(d*86400 - h*3600))/60.0)!# of minutes (minus days and hours)
    s =  fintime - REAL(d*86400) - REAL(h*3600) - REAL(m*60) !# of seconds (- days, hrs and mins)

    11 format('Time to run model = ',i4,' days ',i4,' hours ',i4,              &
              ' minutes and ',f10.4,' seconds.')
    12 format('Time to run model = ',i4,' hours ',i4,' minutes and ',f10.4,     &
              ' seconds.')
    13 format('Time to run model = ',i4,' minutes and ',f10.4,' seconds.')
    14 format('Time to run model = ',f10.4,' seconds.')

    if(fintime > 86400.0)then
      write(*,11) d,h,m,s
    elseif(fintime > 3600.0)then
      write(*,12) h,m,s
    elseif(fintime > 60.0)then
      write(*,13) m,s
    else
      write(*,14) fintime
    endif

     write(*,'(/,A)') '****** END ROMSPath *******'

   end subroutine fin_ROMSPath

    

  subroutine update_particles()

    USE PARAM_MOD,      ONLY: numpar,xi_rho,eta_rho,s_rho,s_w,xi_u,eta_u,	   &
							  idt,HTurbOn,VTurbOn,settlementon,xi_v,eta_v,     &
                              Behavior,SaltTempOn,OpenOceanBoundary,Swimdepth, &
                              TrackCollisions,WriteModelTiming,mortality,      &
                              ErrorFlag,t_b,t_c,t_f,Ngrid,vertdist,scheme,nsb, &
							  SaltTempMean,WriteBottom,maxsize,Process_VA							  
    !USE SETTLEMENT_MOD, ONLY: isSettled,testSettlement
#ifdef GROWTH
	USE GROWTH_MOD,   ONLY:  growlarva
#endif

    USE BEHAVIOR_MOD,   ONLY: behave
    USE BOUNDARY_MOD,   ONLY: bounds
    USE GRID_MOD,    ONLY: getSlevel,getWlevel,GRIDS
    USE HTURB_MOD,      ONLY: HTurb
    USE VTURB_MOD,      ONLY: VTurb
    USE ADVECTION_MOD,  ONLY: RKAdvect
    USE INT_MOD,        ONLY: getinterp2d,getinterp3d,polintd,getInterpStr

    IMPLICIT NONE

    ! Iteration Variables
    INTEGER :: i,deplvl,n

    ! Particle tracking
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: Pwc_zb,Pwc_zc,Pwc_zf
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: Pwc_wzb,Pwc_wzc,Pwc_wzf
    DOUBLE PRECISION :: Xpar,Ypar,Zpar,newXpos,newYpos,newZpos,P_zb,P_zc,P_zf, &
      P_depth,P_zeta,ey(3)
    
    ! Behavior and Turbulence
    DOUBLE PRECISION :: TurbHx,TurbHy,TurbV,Behav,XBehav,YBehav,ZBehav
    LOGICAL :: bott   ! for Behavior 7 along with XBehav,YBehav,ZBehav

    ! Boundaries
    INTEGER :: intersectf,skipbound,inbounds,reflects,inpoly,nmask,Inode,Jnode &
	  ,ngid
    DOUBLE PRECISION :: reflect,fintersectX,fintersectY,freflectX,freflectY,   &
      Xpos,Ypos,nXpos,nYpos,pm,pn
    LOGICAL :: ingrid,obound,hbot,htop

    ! Advection
    DOUBLE PRECISION :: AdvectX,AdvectY,AdvectZ,maxpartdepth,minpartdepth,     &
      kn1_u,kn1_v,kn1_w,kn2_u,kn2_v,kn2_w,kn3_u,kn3_v,kn3_w,kn4_u,kn4_v,kn4_w, &
      P_V,P_U,P_W,UAD,VAD,WAD,x1,x2,x3,y1,y2,y3,z1,z2,z3,btemp,tbustr,tbvstr


	  
	DOUBLE PRECISION :: tempX,tempY,tdepth,zeta,zetab,zetac,zetaf,behout(4)
	  

   
    DO n=1,numpar

      ! *********************************************************
      ! *                                                       *
      ! *        Update Particle Age and Characteristics        *
      ! *                                                       *
      ! *********************************************************

      !If the particle is not yet released, set new location to 
      !  current location, and cycle to next particle
	  if(ix(3) <= par(n,pDOB))then
         cycle
	  endif
	  
	  if(par(n,pStatus).eq.9.0)then
         cycle
	  endif
#ifdef GROWTH	  
	  if(par(n,pSize).gt.maxsize)then
		 par(n,pStatus)=9.0
         cycle
	  endif
#endif
	  
      ! !If there are open ocean boundaries and the current
      ! !  particle has exited the model domain via them, skip it
	  if(OpenOceanBoundary)then
         if(.not.isIn(n)) cycle
       endif
		
       !Update particle age
       par(n,pAge) = par(n,pAge) + float(idt)
	   
	   

      ! !If particle settled or dead, skip tracking
      ! if(settlementon)then
        ! if ( isSettled(n) ) cycle
      ! endif


	  


      ! ! *********************************************************
      ! ! *                                                       *
      ! ! *          Find Element that Contains Particle          *
      ! ! *                                                       *
      ! ! *********************************************************

      ! !Get node Boundary ibformation
      Xpar = par(n,pX)
      Ypar = par(n,pY)
      Zpar = par(n,pZ)
	  Inode=floor(Xpar)
	  Jnode=floor(Ypar)



      ! ! *********************************************************
      ! ! *                                                       *
      ! ! *             Prepare for Particle Movement             *
      ! ! *                                                       *
      ! ! *********************************************************
		
	  pm = getInterp2D("pm",int(par(n,pGID)),Xpar,Ypar,1)
	  pn = getInterp2D("pn",int(par(n,pGID)),Xpar,Ypar,1)
      AdvectX = 0.0
      AdvectY = 0.0
      AdvectZ = 0.0
      TurbHx = 0.0
      TurbHy = 0.0
      TurbV = 0.0
      ! Behav = 0.0

      ! ! *********************************************************
      ! ! *                                                       *
      ! ! *                       ADVECTION                       *
      ! ! *                                                       *
      ! ! *********************************************************

	  
      ! !Find advection currents at original coordinates
	   call CPU_TIME(times(1))

	   SELECT CASE (scheme)
		CASE (1)
			call RKAdvect(Xpar,Ypar,Zpar,ex,ix,pm,pn,int(par(n,pGID)),ets,AdvectX,AdvectY,AdvectZ)
        CASE DEFAULT
			AdvectX = 0.0
			AdvectY = 0.0
			AdvectZ = 0.0
		END SELECT
		

		
      ! ! *********************************************************
      ! ! *                                                       *
      ! ! *                  Horizontal Turbulence                *
      ! ! *                                                       *
      ! ! *********************************************************

      ! IF (WriteModelTiming) call CPU_TIME(times(4))


	   call CPU_TIME(times(2))
       IF (HTurbOn) CALL HTurb(TurbHx,TurbHy,int(par(n,pGID)))
	   TurbHx=TurbHx*pm
	   TurbHy=TurbHy*pn	
		

      ! ! *********************************************************
      ! ! *                                                       *
      ! ! *                   Verticle Turbulence                 *
      ! ! *                                                       *
      ! ! ********************************************************* 


	  
	   call CPU_TIME(times(3))
	   
	  
	   
       IF (VTurbOn) CALL VTurb(Xpar,Ypar,Zpar,ets,ex,ix,int(par(n,pGID)),TurbV)

	   

  
      ! ! *********************************************************
      ! ! *                                                       *
      ! ! *                       Behavior                        *
      ! ! *                                                       *
      ! ! *********************************************************


	 
	   call CPU_TIME(times(4))
	   !write(*,*) times(4)-times(3)
	   
		 
		CALL behave(Xpar,Ypar,Zpar,XBehav,YBehav,ZBehav,par(n,pSize),ex,ix,int(par(n,pGID)),behout)
		par(n,pAcc)=behout(1)
		par(n,pVort)=behout(2)
		par(n,pbehaveW)=behout(3)
		par(n,pSSF)=behout(4)
          


	
	 
      ! ! *********************************************************
      ! ! *                                                       *
      ! ! *     Update Particle Locations and Check Boundaries    *
      ! ! *                                                       *
      ! ! *********************************************************

      ! IF(WriteModelTiming) call CPU_TIME(times(7))

	
	   call CPU_TIME(times(5))
      !Update due to Advection and Turbulence
      newXpos = par(n,pX) + AdvectX + TurbHx 
      newYpos = par(n,pY) + AdvectY + TurbHy
	 
       !

      ! !Assign new particle positions
	  
	  do i=1,Ngrid
	  	tempX=GRIDS(int(par(n,pGID)))%scl(i,1)*newXpos+GRIDS(int(par(n,pGID)))%off(i,1)
		tempY=GRIDS(int(par(n,pGID)))%scl(i,2)*newYpos+GRIDS(int(par(n,pGID)))%off(i,2)
		call bounds(i,tempX,tempY,nmask,ingrid,obound)
		if (ingrid) then  
			par(n,pGID)=i
			exit
		endif
	  enddo
	  tdepth = DBLE(-1.0)* getInterp2D("depth",int(par(n,pGID)),tempX,tempY,t_c)
	  
	  ey(1) =  DBLE(1.0)*getInterp2D("zeta",int(par(n,pGID)),tempX,tempY,t_b)
	  ey(2) =  DBLE(1.0)*getInterp2D("zeta",int(par(n,pGID)),tempX,tempY,t_c)
	  ey(3) =  DBLE(1.0)*getInterp2D("zeta",int(par(n,pGID)),tempX,tempY,t_f)
			  
	  P_zeta=polintd(ex,ey,3,ix(2))
				
				
	 par(n,pWD)=(DBLE(-1.0)*tdepth)+P_zeta
	 hbot=.FALSE.
	 htop=.FALSE.
	 
	 
	  SELECT CASE (nsb)
		CASE (0)
            newZpos = par(n,pZ) + AdvectZ + TurbV+Zbehav
			!newZpos = par(n,pZ) + AdvectZ +Zbehav
		     if (newZpos.LT.tdepth) 	then
				newZpos = tdepth +  ABS(newZpos-tdepth)
				hbot=.TRUE.
			endif
			 if (newZpos.GT.P_zeta) 	then
			  !  write(*,*) P_zeta,TurbV,ABS(newZpos-P_zeta)
				newZpos = P_zeta - ABS(newZpos-P_zeta)
				htop=.TRUE.
			endif
		CASE (1)  !Near-Surface
			newZpos = P_zeta-vertdist
		CASE (2) !Near-Bottom
			newZpos = P_depth+vertdist
		CASE DEFAULT
			WRITE(*,*) 'NO VALID BEHAVIOR SET' 
			EXIT
	  END SELECT

	  if (ingrid) then
		isIn(n)=ingrid
		par(n,pX) = tempX
		par(n,pY) = tempY
		par(n,pZ) = newZpos
		par(n,pStatus) = 1.0
		if (htop) par(n,pStatus) = 1.1
		if (hbot) par(n,pStatus) = 1.2
		
	  else
		if (obound) then
			isIn(n)=ingrid
			par(n,pGID)=Ngrid
			par(n,pStatus) = -8.0 
		else
			isIn(n)=.TRUE.
			par(n,pStatus) = -9.0 
		endif
	 endif	  

	 

    

  
#ifdef GROWTH
      ! ! *********************************************************
      ! ! *                                                       *
      ! ! *                      Growth                       *
      ! ! *                                                       *
      ! ! *********************************************************
	  
	  call growlarva(P_temp(n),P_salt(n),par(n,pAge),par(n,pSize),par(n,pStatus))
       
#endif
      ! ! *********************************************************
      ! ! *                                                       *
      ! ! *                      Settlement                       *
      ! ! *                                                       *
      ! ! *********************************************************

      ! if(settlementon) then

        ! CALL testSettlement(par(n,pAge),n,par(n,pX),par(n,pY),inpoly)
        ! if (inpoly .GT. 0) then
          ! par(n,pnZ) = P_depth
          ! endpoly(n) = inpoly
          ! par(n,pLifespan) = par(n,pAge)
        ! endif

      ! endif 

  	 ! tdepth = DBLE(-1.0)* getInterp2D("depth",int(par(n,pGID)),par(n,pX),par(n,pY),t_c)
	  
	  
  	 
	  ! ! *********************************************************
      ! ! *                                                       *
      ! ! *               Bottom Stuff              *
      ! ! *                                                       *
      ! ! *********************************************************
	  if (WriteBottom) then
		btemp=DBLE(-1.0)*tdepth+par(n,pZ)
		P_HOB(n)=min(btemp,P_HOB(n))
		
		if (btemp.EQ.P_HOB(n)) then 
		
			ey(1)=getInterpStr("bustr",int(par(n,pGID)),par(n,pX),par(n,pY),t_b)
			ey(2)=getInterpStr("bustr",int(par(n,pGID)),par(n,pX),par(n,pY),t_c)
			ey(3)=getInterpStr("bustr",int(par(n,pGID)),par(n,pX),par(n,pY),t_f)
			P_bustr(n)=polintd(ex,ey,3,ix(2))
			ey(1)=getInterpStr("bvstr",int(par(n,pGID)),par(n,pX),par(n,pY),t_b)
			ey(2)=getInterpStr("bvstr",int(par(n,pGID)),par(n,pX),par(n,pY),t_c)
			ey(3)=getInterpStr("bvstr",int(par(n,pGID)),par(n,pX),par(n,pY),t_f)
			P_bvstr(n)=polintd(ex,ey,3,ix(2))
		endif
	   
	  endif
	  
	  
  	 
	  ! ! *********************************************************
      ! ! *                                                       *
      ! ! *                Salinity and Temperature               *
      ! ! *                                                       *
      ! ! *********************************************************
	  
	   call CPU_TIME(times(6))
		if (SaltTempOn) then
				
				zetab =  getInterp2D("zeta",int(par(n,pGID)),par(n,pX),par(n,pY),t_b)
				zetac =  getInterp2D("zeta",int(par(n,pGID)),par(n,pX),par(n,pY),t_c)
				zetaf =  getInterp2D("zeta",int(par(n,pGID)),par(n,pX),par(n,pY),t_f)
			  
				ey(1)=getInterp3d("salt",int(par(n,pGID)),par(n,pX),par(n,pY),par(n,pZ),t_b,1,zetab,tdepth)
				ey(2)=getInterp3d("salt",int(par(n,pGID)),par(n,pX),par(n,pY),par(n,pZ),t_c,1,zetac,tdepth)
				ey(3)=getInterp3d("salt",int(par(n,pGID)),par(n,pX),par(n,pY),par(n,pZ),t_f,1,zetaf,tdepth)
				P_salt(n)=polintd(ex,ey,3,ix(2))
				
				
				ey(1)=getInterp3d("temp",int(par(n,pGID)),par(n,pX),par(n,pY),par(n,pZ),t_b,1,zetab,tdepth)
				ey(2)=getInterp3d("temp",int(par(n,pGID)),par(n,pX),par(n,pY),par(n,pZ),t_c,1,zetac,tdepth)
				ey(3)=getInterp3d("temp",int(par(n,pGID)),par(n,pX),par(n,pY),par(n,pZ),t_f,1,zetaf,tdepth)
				
				P_temp(n)=polintd(ex,ey,3,ix(2))	
				if (SaltTempMean) then
					mean_salt(n)=mean_salt(n)+P_salt(n)
					mean_temp(n)=mean_temp(n)+P_temp(n)
					
				endif
					 
			    
		endif
      ! ! *****************************************************************
      ! ! *                      End of Particle Loop                     *
      ! ! *****************************************************************

      ! IF(WriteModelTiming) then
         call CPU_TIME(times(7))

         timeCounts(1) = timeCounts(1) + (times(2)-times(1))
         timeCounts(2) = timeCounts(2) + (times(3)-times(2))
         timeCounts(3) = timeCounts(3) + (times(4)-times(3))
         timeCounts(4) = timeCounts(4) + (times(5)-times(4))
         timeCounts(5) = timeCounts(5) + (times(6)-times(5))
         timeCounts(6) = timeCounts(6) + (times(7)-times(6))
      ! ENDIF


    ENDDO !end loop for each particle

  


  end subroutine update_particles



  subroutine dataOutput()
    use param_mod,   only: numpar,SaltTempOn,TrackCollisions,stokesprefix,turbstd_v_a_prefix
	USE INT_MOD,        ONLY: getinterp2d
    integer :: n
    double precision, dimension(numpar) :: pLon,pLat

    ! increment file number
    ! prcount = prcount + 1
    do n=1,numpar
		  pLon(n) = getInterp2D("lon",int(par(n,pGID)),par(n,pX),par(n,pY),1)
		  pLat(n) = getInterp2D("lat",int(par(n,pGID)),par(n,pX),par(n,pY),1)	
		  
    enddo
   
   call writeNetCDF(int(ix(3)),pLon,pLat)
   
   
    !Based on user options, write specified data to output
    ! IF(SaltTempOn)THEN
      ! IF(TrackCollisions)THEN
        ! CALL writeOutput(par(:,pX),par(:,pY),par(:,pZ),par(:,pAge),par(:,pStatus),prcount,    &
             ! HITBOTTOM=hitBottom,HITLAND=hitLand,P_SALT=P_Salt,P_TEMP=P_Temp)
      ! ELSE
        ! CALL writeOutput(par(:,pX),par(:,pY),par(:,pZ),par(:,pAge),par(:,pStatus),prcount,    &
             ! P_SALT=P_Salt,P_TEMP=P_Temp)
      ! ENDIF
    ! ELSE
      ! IF(TrackCollisions)THEN
        ! CALL writeOutput(par(:,pX),par(:,pY),par(:,pZ),par(:,pAge),par(:,pStatus),prcount,    &
             ! HITBOTTOM=hitBottom,HITLAND=hitLand)
      ! ELSE
        ! CALL writeOutput(par(:,pX),par(:,pY),par(:,pZ),par(:,pAge),par(:,pStatus),prcount)
      ! ENDIF
    ! ENDIF

   
    !If Tracking Model Timing, write Time data to file
    ! IF(WriteModelTiming)then
      ! call CPU_TIME(times(9))

      ! timeCounts(8) = times(9)-times(1)

      ! OPEN(300,FILE='Timing.csv',POSITION='APPEND')

        ! write(300,"(15(F14.4,','),F14.4)") (DBLE(ix(3))/DBLE(86400)),timeCounts(8),    &
          ! timeCounts(1),(timeCounts(1)/timeCounts(8))*DBLE(100.00),            &
          ! timeCounts(2),(timeCounts(2)/timeCounts(8))*DBLE(100.00),            &
          ! timeCounts(3),(timeCounts(3)/timeCounts(8))*DBLE(100.00),            &
          ! timeCounts(4),(timeCounts(4)/timeCounts(8))*DBLE(100.00),            &
          ! timeCounts(5),(timeCounts(5)/timeCounts(8))*DBLE(100.00),            &
          ! timeCounts(6),(timeCounts(6)/timeCounts(8))*DBLE(100.00),            &
          ! timeCounts(7),(timeCounts(7)/timeCounts(8))*DBLE(100.00)

      ! CLOSE(300)

      ! timeCounts = 0
      ! call CPU_TIME(times(1))
    ! ENDIF

  
 end subroutine dataOutput

  ! SUBROUTINE writeOutput(Xpar,Ypar,Zpar,P_age,P_status,prcount,hitBottom,hitLand,P_Salt,P_Temp)
    ! USE PARAM_MOD, ONLY: numpar,SaltTempOn,outpathGiven,outpath,      &
                         ! TrackCollisions,Behavior
    ! ! USE BEHAVIOR_MOD, ONLY: getStatus
    ! USE INT_MOD,        ONLY: getinterp2d
    ! USE HYDRO_MOD, ONLY: writeNetCDF

    ! IMPLICIT NONE

    ! DOUBLE PRECISION, INTENT(IN) :: Xpar(numpar),Ypar(numpar),Zpar(numpar),P_age(numpar),P_status(numpar)
    ! INTEGER         , INTENT(IN) :: prcount
    ! INTEGER, DIMENSION(numpar), INTENT(IN), OPTIONAL :: hitBottom,hitLand
    ! DOUBLE PRECISION, DIMENSION(numpar), INTENT(IN), OPTIONAL :: P_Salt,P_Temp

    ! INTEGER :: n
    ! DOUBLE PRECISION :: statuses(numpar)
    ! double precision, dimension(numpar) :: pLon,pLat

    ! !INPUT/OUTPUT FILE NAME CONSTRUCTION VARIABLES
    ! CHARACTER(LEN=100) :: filenm2
    ! CHARACTER(LEN=4  ) :: prefix2,suffix2
    ! INTEGER :: counter2

    ! !Convert particle position (in meters) to latitude and longitude &
    ! !Find identification number that describes a particle's behavior 
    ! !  type or status for use in visualization routines
    ! do n=1,numpar
		  ! pLon(n) = getInterp2D("lon",int(par(n,pGID)),Xpar(n),Ypar(n),1)
		  ! pLat(n) = getInterp2D("lat",int(par(n,pGID)),Xpar(n),Ypar(n),1)	
		  ! statuses(n) =P_status(n)

   ! enddo

   



      ! !Based on user options, Write specified data to file
      ! if (SaltTempOn) then
        ! if(TrackCollisions)then
          ! CALL writeNetCDF(int(ix(3)),P_age(:),pLon,pLat,Zpar,statuses,par(:,pGid),      &
               ! SALT=P_Salt,TEMP=P_Temp,HITB=hitBottom,HITL=hitLand)
        ! else
          ! CALL writeNetCDF(int(ix(3)),P_age(:),pLon,pLat,Zpar,statuses,par(:,pGid),           &
               ! SALT=P_Salt,TEMP=P_Temp)
        ! endif
      ! else
        ! if(TrackCollisions)then
          ! CALL writeNetCDF(int(ix(3)),P_age(:),pLon,pLat,Zpar,statuses,par(:,pGid),           &
               ! HITB=hitBottom,HITL=hitLand)
        ! else
          ! CALL writeNetCDF(int(ix(3)),P_age(:),pLon,pLat,Zpar,statuses,par(:,pGid))
        ! endif
      ! endif


  ! END SUBROUTINE writeOutput  
  
	SUBROUTINE initNetCDF()
	  
		!Initialize NetCDF Counters
		NCcount = 0
		NCstart = 0

	  END SUBROUTINE initNetCDF

	  SUBROUTINE createNetCDF(dob)
		USE PARAM_MOD, ONLY: numpar,NCOutFile,outpath,outpathGiven,NCtime,         &
			RunName,ExeDir,OutDir,RunBy,Institution,StartedOn,SaltTempMean,         &
			TrackCollisions,SaltTempOn,Ngrid,days,idt,VTurbOn,HTurbOn,deltat,      &
			serr,smth,sub,AKSback,maxsize,tempcut,initsize,deadage,		&
			a0,a1,a2,a3,a4,a5,a6,a7,a8,TempOffset,WriteBottom,WriteWaterDepth,Behavior,		&
			vort_cr,vort_sat,b0pv,b1pv,b0wv,b1w,acc_cr,acc_sat,&
			b0pa,b1pa,b0wa,va_flag,OpenOceanBoundary,swimfast,Process_VA,	&
			WriteWaterDepth,seed
		USE GRID_MOD, ONLY: reftime, time_units
		USE netcdf
		IMPLICIT NONE

		DOUBLE PRECISION, DIMENSION(numpar), OPTIONAL, INTENT(IN) :: dob

		INCLUDE 'netcdf.inc'

		CHARACTER(LEN=200) :: ncFile
		
	    character(len=10) :: sdate,stime,zone
		
	    character(len=20) :: sdatetime
		INTEGER :: STATUS,NCID,numparID,timeID,pageID,modtimeID,lonID,latID,ngID,       &
				   depthID,statusID,hitBID,hitLID,dobID,saltID,tempID,date_time(8),		&
				   HOBID,bustrID,bvstrID,VORTID,ACCID,BWID,SSFID,WDID
#ifdef GROWTH		
		INTEGER ::		   sizeID
#endif

		!   NF90_CREATE           ! create netCDF dataset: enter define mode
		!        ...
		!      NF90_DEF_DIM       ! define dimensions: from name and length
		!        ...
		!      NF90_DEF_VAR       ! define variables: from name, type, dims
		!        ...
		!      NF90_PUT_ATT       ! assign attribute values
		!        ...
		!   NF90_ENDDEF           ! end definitions: leave define mode
		!        ...
		!      NF90_PUT_VAR       ! provide values for variable
		!        ...
		!   NF90_CLOSE            ! close: save new netCDF dataset


		!NF90_CREATE

		!Reset Print Counter to 0
		prcount = 0

		IF(outpathGiven)THEN
		  IF(NCtime == 0 ) THEN
			ncFile = TRIM(outpath) // TRIM(NCOutFile) // '.nc'
		  ELSE
			NCcount = NCcount + 1
			write(ncFile,"(A,A,A,I3.3,A)")TRIM(outpath),TRIM(NCOutFile),'_',       &
										  NCcount,'.nc'
		  ENDIF
		ELSE
		  IF(NCtime == 0 ) THEN
			ncFile = TRIM(NCOutFile) // '.nc'
		  ELSE
			NCcount = NCcount + 1
			write(ncFile,"(A,A,I3.3,A)")TRIM(NCOutFile),'_',NCcount,'.nc'
		  ENDIF
		ENDIF

		write(*,*)'Creating NetCDF Output File: ',TRIM(ncFile)

		STATUS = NF90_CREATE(TRIM(ncFile), NF90_NETCDF4, NCID)
		IF(STATUS /= NF90_NOERR) THEN
		  WRITE(*,*) 'Problem creating NetCDF output file'
		  WRITE(*,*) NF_STRERROR(STATUS)
		  STOP
		ENDIF

		!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		!NF90_DEF_DIM

			STATUS = NF90_DEF_DIM(NCID,'numpar',numpar,numparID)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: numpar dim'
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_DEF_DIM(NCID,'time',NF90_UNLIMITED,timeID)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: time dim'
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

		!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		!NF90_DEF_VAR

			STATUS = NF90_DEF_VAR(NCID,'model_time',NF_DOUBLE,(/timeID/),modtimeID,&
								  deflate_level=1,shuffle=.true.)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: time var'
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			IF( PRESENT(dob) )THEN
			  STATUS = NF90_DEF_VAR(NCID,'dob',NF_FLOAT,(/numparID/),dobID,       &
								  deflate_level=1,shuffle=.true.)
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: dob var'
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			ENDIF

			STATUS = NF90_DEF_VAR(NCID,'age',NF_DOUBLE,(/numparID,timeID/),pageID, &
								  deflate_level=1,shuffle=.true.)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: age var'
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_DEF_VAR(NCID,'lon',NF_FLOAT,(/numparID,timeID/),lonID,  &
								  deflate_level=1,shuffle=.true.)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: lon var'
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_DEF_VAR(NCID,'lat',NF_FLOAT,(/numparID,timeID/),latID,  &
								  deflate_level=1,shuffle=.true.)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: lat var'
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_DEF_VAR(NCID,'depth',NF_FLOAT,(/numparID,timeID/),      &
								  depthID,deflate_level=1,shuffle=.true.)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: depth var'
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_DEF_VAR(NCID,'status',NF_FLOAT,(/numparID,timeID/),      &
								  statusID,deflate_level=1,shuffle=.true.)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: status var'
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			IF (Ngrid .gt.1) THEN
				STATUS = NF90_DEF_VAR(NCID,'GID',NF_FLOAT,(/numparID,timeID/),      &
								  ngID,deflate_level=1,shuffle=.true.)
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: status var'
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


			ENDIF
			
			IF(TrackCollisions)THEN
			  STATUS =NF90_DEF_VAR(NCID,'hitBottom',NF_FLOAT,(/numparID,timeID/), &
								   hitBID,deflate_level=1,shuffle=.true.)
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: Bottom ', &
												  'Collision var'
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_DEF_VAR(NCID,'hitLand',NF_FLOAT,(/numparID,timeID/),  &
									hitLID,deflate_level=1,shuffle=.true.)
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: Land ',   &
												  'Collision var'
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			ENDIF

			IF(SaltTempOn)THEN
			  STATUS = NF90_DEF_VAR(NCID,'salinity',NF_FLOAT,(/numparID,timeID/), &
									saltID,deflate_level=1,shuffle=.true.)
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
												  'Salinity var'
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_DEF_VAR(NCID,'temperature',NF_FLOAT,                  &
									(/numparID,timeID/),tempID,deflate_level=1)
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
												  'Temperature var'
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			ENDIF
			
			
			IF(WriteBottom)THEN
			  STATUS = NF90_DEF_VAR(NCID,'HOB',NF_FLOAT,(/numparID,timeID/), &
									HOBID,deflate_level=1,shuffle=.true.)
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
												  'HIB var'
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_DEF_VAR(NCID,'bustr',NF_FLOAT,                  &
									(/numparID,timeID/),bustrID,deflate_level=1)
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
												  'Bottom Stress U var'
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			  
			  STATUS = NF90_DEF_VAR(NCID,'bvstr',NF_FLOAT,                  &
									(/numparID,timeID/),bvstrID,deflate_level=1)
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
												  'Bottom Stress V var'
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			ENDIF
			
			IF(WriteWaterDepth)THEN
			  STATUS = NF90_DEF_VAR(NCID,'WaterDepth',NF_FLOAT,(/numparID,timeID/), &
									WDID,deflate_level=1,shuffle=.true.)
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
												  'Waterdepth var'
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			ENDIF
			
			IF((Behavior .eq. 10) .OR. (Process_VA))THEN
				STATUS = NF90_DEF_VAR(NCID,'vorticity',NF_FLOAT,(/numparID,timeID/), &
									VORTID,deflate_level=1,shuffle=.true.)
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
												  'vorticity var'
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			  
				STATUS = NF90_DEF_VAR(NCID,'acceleration',NF_FLOAT,(/numparID,timeID/), &
									ACCID,deflate_level=1,shuffle=.true.)
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
												  'acceleration var'
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			  
				STATUS = NF90_DEF_VAR(NCID,'behave_w',NF_FLOAT,(/numparID,timeID/), &
									BWID,deflate_level=1,shuffle=.true.)
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
												  'behave_w var'
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			  
				STATUS = NF90_DEF_VAR(NCID,'swinksinkflag',NF_FLOAT,(/numparID,timeID/), &
									SSFID,deflate_level=1,shuffle=.true.)
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
												  'swinksinkflag var'
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			
			
			endif
			
#ifdef GROWTH
			
			STATUS = NF90_DEF_VAR(NCID,'size',NF_FLOAT,(/numparID,timeID/),sizeID, &
								  deflate_level=1,shuffle=.true.)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: size var'
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

#endif
			

		!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		!NF90_PUT_ATT

			!Particle Time
			STATUS = NF90_PUT_ATT(NCID, modtimeID, "long_name",                    &
								  "Model time")
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_PUT_ATT(NCID, modtimeID, "units", trim(time_units))
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_PUT_ATT(NCID, modtimeID, "field",                        &
								  "model_time, scalar, series")
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			IF( PRESENT(dob) )THEN
			  !Particle Date of Birth
			  STATUS = NF90_PUT_ATT(NCID, dobID, "long_name",                      &
					   "Date of Birth of particles in seconds from model start")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_PUT_ATT(NCID, dobID, "units", "seconds")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_PUT_ATT(NCID, dobID, "field", "age, scalar, series")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			ENDIF

			!Particle Age
			STATUS = NF90_PUT_ATT(NCID, pageID, "long_name", "age of particles")
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_PUT_ATT(NCID, pageID, "units", "seconds")
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_PUT_ATT(NCID, pageID, "field", "age, scalar, series")
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			!Longitude
			STATUS = NF90_PUT_ATT(NCID, lonID, "long_name",                        &
								  "longitude of particles")
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_PUT_ATT(NCID, lonID, "units", "decimal degrees E")
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_PUT_ATT(NCID, lonID, "field", "lon, scalar, series")
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


			!Latitude
			STATUS =NF90_PUT_ATT(NCID, latID, "long_name", "latitude of particles")
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_PUT_ATT(NCID, latID, "units", "decimal degrees N")
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_PUT_ATT(NCID, latID, "field", "lat, scalar, series")
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


			!Depth
			STATUS = NF90_PUT_ATT(NCID, depthID, "long_name", "depth of particles")
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_PUT_ATT(NCID, depthID, "units", "meters below ROMS z=0")
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_PUT_ATT(NCID, depthID, "field", "depth, scalar, series")
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


			!status
			STATUS = NF90_PUT_ATT(NCID, statusID, "long_name",                      &
					 "identification number for particle behavior or status")
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_PUT_ATT(NCID, statusID, "units",                          &
								  "nondimensional, see ROMSPath User Guide")
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_PUT_ATT(NCID, statusID, "field", "status, scalar, series")
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
#ifdef GROWTH


			STATUS = NF90_PUT_ATT(NCID, sizeID, "long_name", "size of particles")
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_PUT_ATT(NCID, sizeID, "units", "um")
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_PUT_ATT(NCID, sizeID, "field", "size, scalar, series")
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

#endif
			IF (Ngrid .gt.1) THEN
			
				STATUS = NF90_PUT_ATT(NCID, ngID, "long_name",                      &
						"Grid ID")
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

				STATUS = NF90_PUT_ATT(NCID, ngID, "units",                          &
								  "none")
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


			ENDIF

			IF(TrackCollisions)THEN
			  !hitBottom
			  STATUS = NF90_PUT_ATT(NCID, hitBID, "long_name",                     &
									"# of times Particle Collided with Bottom")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_PUT_ATT(NCID, hitBID, "units", "Number of Collisions")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_PUT_ATT(NCID, hitBID, "field",                         &
									"hitBottom, scalar, series")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


			  !hitLand
			  STATUS = NF90_PUT_ATT(NCID, hitLID, "long_name",                     &
									"# of times Particle Collided with Land")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_PUT_ATT(NCID, hitLID, "units", "Number of Collisions")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_PUT_ATT(NCID, hitLID, "field",                         &
									"hitLand, scalar, series")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			ENDIF

			IF(SaltTempOn)THEN
			  !salt
			  STATUS = NF90_PUT_ATT(NCID, saltID, "long_name",                     &
									"Salinity at the particle's location")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS=NF90_PUT_ATT(NCID,saltID, "field", "salinity, scalar, series")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


			  !temp
			  STATUS = NF90_PUT_ATT(NCID, tempID, "long_name",                     &
									"Temperature at the particle's location")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_PUT_ATT(NCID, tempID, "units", "° Celsius")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_PUT_ATT(NCID, tempID, "field",                         &
									"temperature, scalar, series")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			ENDIF

		   IF(WriteBottom)THEN
			  !HIB
			  STATUS = NF90_PUT_ATT(NCID, HOBID, "long_name",                     &
									"Particle Minimum height above bottom since last output time step")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS=NF90_PUT_ATT(NCID,HOBID, "field", "HOB, scalar, series")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
				
			  STATUS = NF90_PUT_ATT(NCID, HOBID, "units", "m")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  !BUSTR
			  STATUS = NF90_PUT_ATT(NCID, bustrID, "long_name",                     &
									"U bottom stress at minimum HOB")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_PUT_ATT(NCID, bustrID, "units",  "newton meter-2" )
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			  !BVSTR
			  STATUS = NF90_PUT_ATT(NCID, bvstrID, "long_name",                     &
									"V bottom stress at minimum HOB")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_PUT_ATT(NCID, bvstrID, "units",  "newton meter-2" )
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  
			ENDIF
			IF(WriteWaterDepth)THEN
			  !HOB
			  STATUS = NF90_PUT_ATT(NCID, WDID, "long_name",                     &
									"Local Water Depth at Particle Location")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS=NF90_PUT_ATT(NCID,WDID, "field", "WaterDepth, scalar, series")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
				
			  STATUS = NF90_PUT_ATT(NCID, WDID, "units", "m")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			 endif
			
			
			IF((Behavior .eq. 10) .OR. (Process_VA))THEN		  
			!VORTID,ACCID,BWID,SSFID
			!Vorticity
			  STATUS = NF90_PUT_ATT(NCID, VORTID, "long_name",                     &
									"Turbulent Vorticity Magnitude at particle location")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS=NF90_PUT_ATT(NCID,VORTID, "field", "vorticity, scalar, series")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
				
			  STATUS = NF90_PUT_ATT(NCID, VORTID, "units", "1/s")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			  
			  
			!Acceleration
			  STATUS = NF90_PUT_ATT(NCID, ACCID, "long_name",                     &
									"Acceleration Magnitude at particle location")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS=NF90_PUT_ATT(NCID,ACCID, "field", "acceleration, scalar, series")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
				
			  STATUS = NF90_PUT_ATT(NCID, ACCID, "units", "m/s^2")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			  
			  
			  !Behavior vertical velocity
			  STATUS = NF90_PUT_ATT(NCID, BWID, "long_name",                     &
									"Behavior vertical velocity at particle location")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS=NF90_PUT_ATT(NCID,BWID, "field", "behavior vertical velocity, scalar, series")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
				
			  STATUS = NF90_PUT_ATT(NCID, BWID, "units", "m/s")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
				
				
			  !Sinking/Swimming flag
			  STATUS = NF90_PUT_ATT(NCID, SSFID, "long_name",                     &
									"Sinking/Swimming flag at particle location")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS=NF90_PUT_ATT(NCID,SSFID, "field", "flag, scalar, series")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
				
			  STATUS = NF90_PUT_ATT(NCID, SSFID, "units", "none")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
				
			



		ENDIF			
			!Global
			

			STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "title", "Particle tracking output using ROMS data.")
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "author", "Elias Hunter (hunter@marine.rutgers.edu)")
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "adv_time_step", idt)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			
			STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "days_run", days)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			
			
			STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "Behavior", Behavior)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			
			If( OpenOceanBoundary) then
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "OpenOceanBoundary", 'TRUE')
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS) 
			endif
			
			STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "swimfast", swimfast)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			
			IF((Behavior .eq. 10) .OR. (Process_VA))THEN
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "vort_cr", vort_cr)
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "vort_sat", vort_sat)
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "b0pv", b0pv)
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "b1pv", b1pv)
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
				
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "acc_cr", acc_cr)
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "acc_sat", acc_sat)
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "b0pa", b0pa)
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "b1pa", b1pa)
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "b0wa", b0wa)
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "va_flag", va_flag)
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			
			
			endif
			
			if (VturbOn) THEN
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "Vertical_Turb", 'True')
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "seed", seed)
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "turb_time_step", deltat)
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "spline_error_cutoff", serr)
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "spline_smooth_param", smth)
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "subset_scale", sub)
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "AKs_bakcground",AKSback )
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			else
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "Vertical_Turb", 'False')
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			ENDIF
					
					
					
			if (HturbOn) THEN
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "Horiz_Turb", 'True')
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "seed", seed)
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			else
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "Horiz_Turb", 'False')
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			ENDIF
			
			 IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "run_name", RunName)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "executable", ExeDir)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "output_loc", OutDir)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "administrator", RunBy)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "institution", Institution)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			
			
            call date_and_time(sdate,stime,zone,date_time)
			write(sdatetime,"(I4'/'I0.2'/'I0.2' 'I0.2':'I0.2':'I0.2)") date_time(1),date_time(2),date_time(3),date_time(5),date_time(6),date_time(7)
			
			STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "Start_datetime", sdatetime)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			
#ifdef GROWTH
			STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "initsize", initsize)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			
		
			STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "maxsize", maxsize)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			
			STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "tempcut", tempcut)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			
			STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "a0", a0)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			
			STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "a1", a1)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			
			STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "a2", a2)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			
			STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "a3", a3)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			
			STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "a4", a4)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			
			STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "a5", a5)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			
			STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "a6", a6)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			
			STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "a7", a7)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			
			STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "a8", a8)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
#endif
			if (SaltTempMean) then
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "SaltTempMean", 'TRUE')
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			else
			
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "SaltTempMean", 'FALSE')
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			endif
			STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "TempOffset", TempOffset)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			! call date_and_time(sdate,stime,zone,date_time)
			! write(sdatetime,"(I4'/'I0.2'/'I0.2' 'I0.2':'I0.2':'I0.2)") date_time(1),date_time(2),date_time(3),date_time(5),date_time(6),date_time(7)
			! ! write(*,"(I4 I2 I2 I2 I2 I2)"),date_time(1),date_time(2),date_time(3),date_time(5),date_time(6),date_time(7)
			! write(*,*) sdatetime
			


		!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		!NF_ENDDEF

		  STATUS = NF90_ENDDEF(NCID)
		  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: EndDef'
		  IF(status /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

		!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		!NF_PUT_VAR

		  !Particle Date of Birth
		  IF( PRESENT(dob) )THEN
			STATUS = NF90_INQ_VARID(NCID, "dob", dobID)
			STATUS = NF90_PUT_VAR(NCID, dobID, dob)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put dob'
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)  
		  ENDIF

		!NF_CLOSE

		STATUS = NF_CLOSE(NCID)
		IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: Close'
		IF(status /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

	  END SUBROUTINE createNetCDF
  SUBROUTINE writeNetCDF(time,lon,lat)
		USE PARAM_MOD, ONLY: numpar,SaltTempOn,NCOutFile,outpath,outpathGiven,     &
			NCtime,TrackCollisions,Ngrid,SaltTempMean,WriteBottom,WriteWaterDepth,Behavior,		&
			Process_VA,WriteWaterDepth
		USE GRID_MOD, ONLY: reftime
		USE netcdf
		IMPLICIT NONE

		INTEGER, INTENT(IN) :: time
		DOUBLE PRECISION, INTENT(IN) :: lon(numpar),lat(numpar)
		! DOUBLE PRECISION, INTENT(IN), OPTIONAL :: Salt(numpar),Temp(numpar)
		! INTEGER, INTENT(IN), OPTIONAL :: hitB(numpar),hitL(numpar)

		INCLUDE 'netcdf.inc'

		CHARACTER(LEN=200) :: ncFile
		INTEGER :: STATUS,NCID,modtimeID,pageID,lonID,latID,depthID,hitBID,hitLID, &
				   statusID,saltID,tempID,ngID,n,HOBID,bustrID,bvstrID,BehaveID,WDID
#ifdef GROWTH		
		INTEGER ::		   sizeID
#endif
   
		INTEGER :: NCelapsed

		!If only one NetCDF output file is being written to:
		IF(NCtime == 0) THEN

		  IF(outpathGiven)THEN
			ncFile = TRIM(outpath) // TRIM(NCOutFile) // '.nc'
		  ELSE
			ncFile = TRIM(NCOutFile) // '.nc'
		  ENDIF

		!If sequentially numbered NetCDF output files are being written to:
		ELSE

		  NCelapsed = time - NCstart

		  !If specified time interval has been reached, create new NetCDF file
		  IF(NCelapsed >= NCtime) THEN
			NCstart = time
			call createNetCDF()
		  ENDIF

		  IF(outpathGiven)THEN
			write(ncFile,"(A,A,A,I3.3,A)")TRIM(outpath),TRIM(NCOutFile),'_',       &
										  NCcount,'.nc'
		  ELSE
			write(ncFile,"(A,A,I3.3,A)")TRIM(NCOutFile),'_',NCcount,'.nc'
		  ENDIF

		ENDIF

		prcount = prcount + 1

		STATUS = NF90_OPEN(TRIM(ncFile), NF90_WRITE, NCID)
		IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
		
		  !Particle Time
		  STATUS = NF90_INQ_VARID(NCID, "model_time", modtimeID)
		  STATUS = NF90_PUT_VAR(NCID, modtimeID, DBLE(time)+reftime, start = (/ prcount /))
		  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put model_time, time: ',time
		  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)  

		  !Particle Age
		  STATUS = NF90_INQ_VARID(NCID, "age", pageID)
		  STATUS = NF90_PUT_VAR(NCID, pageID, par(:,pAge),             &
								start = (/ 1, prcount /),     &
								count = (/ numpar,  1 /))
		  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put age, time: ',time
		  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

		  !Longitude
		  STATUS = NF90_INQ_VARID(NCID, "lon", lonID)
		  STATUS = NF90_PUT_VAR(NCID, lonID, lon,             &
								start = (/ 1, prcount /),     &
								count = (/ numpar,  1 /))
		  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put lon, time: ',time
		  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

		  !Latitude
		  STATUS = NF90_INQ_VARID(NCID, "lat", latID)
		  STATUS = NF90_PUT_VAR(NCID, latID, lat,             &
								start = (/ 1, prcount /),     &
								count = (/ numpar,  1 /))
		  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put lat, time: ',time
		  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

		  !Depth
		  STATUS = NF90_INQ_VARID(NCID, "depth", depthID)
		  STATUS = NF90_PUT_VAR(NCID, depthID, par(:,pZ),         &
								start = (/ 1, prcount /),     &
								count = (/ numpar,  1 /))
		  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put depth, time: ',time
		  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

		  !status
		  STATUS = NF90_INQ_VARID(NCID, "status", statusID)
		  STATUS = NF90_PUT_VAR(NCID, statusID, par(:,pStatus),        &
								start = (/ 1, prcount /),     &
								count = (/ numpar,  1 /))
		  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put status, time: ',time
		  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
		
		
		  IF (Ngrid .gt. 1) THEN 
			 STATUS = NF90_INQ_VARID(NCID, "GID", ngID)
			 
			STATUS = NF90_PUT_VAR(NCID, ngID,par(:,pGID),        &
								start = (/ 1, prcount /),     &
								count = (/ numpar,  1 /))
		  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put GID, time: ',time
		  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
		   ENDIF
			
		  !hitBottom
		  IF( TrackCollisions)THEN
			STATUS = NF90_INQ_VARID(NCID, "hitBottom", hitBID)
			STATUS = NF90_PUT_VAR(NCID, hitBID, hitBottom,       &
								  start = (/ 1, prcount /),   &
								  count = (/ numpar,  1 /))
			IF(STATUS /= NF90_NOERR) WRITE(*,*)                                    &
			  'Problem put # Bottom Collisions, time: ',time
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
		  ENDIF

		  !hitLand
		  IF( TrackCollisions) THEN
			STATUS = NF90_INQ_VARID(NCID, "hitLand", hitLID)
			STATUS = NF90_PUT_VAR(NCID, hitLID, hitLand,      &
								  start = (/ 1, prcount /),   &
								  count = (/ numpar,  1 /))
			IF(STATUS /= NF90_NOERR) WRITE(*,*)                                    &
			  'Problem put # Land Collisions, time: ',time
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
		  ENDIF

		  !SALTTEMP
		  if (SaltTempOn) then
		  	if (SaltTempMean) then
				do n=1,numpar
					mean_salt(n)=mean_salt(n)/mI
					mean_temp(n)=mean_temp(n)/mI
				enddo
					
				STATUS = NF90_INQ_VARID(NCID, "salinity", saltID)
				STATUS = NF90_PUT_VAR(NCID, saltID, mean_salt,       &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put salinity, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

				STATUS = NF90_INQ_VARID(NCID, "temperature", tempID)
				STATUS = NF90_PUT_VAR(NCID, tempID, mean_temp,      &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put temperature, time: ', &
													time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
				
			    do n=1,numpar
					mean_salt(n)=0.0
					mean_temp(n)=0.0
				enddo
				mI=0
					
			else
				STATUS = NF90_INQ_VARID(NCID, "salinity", saltID)
				STATUS = NF90_PUT_VAR(NCID, saltID, P_salt,       &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put salinity, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

				STATUS = NF90_INQ_VARID(NCID, "temperature", tempID)
				STATUS = NF90_PUT_VAR(NCID, tempID, P_temp,      &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put temperature, time: ', &
													time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
			ENDIF
		  ENDIF
		  
		  
		  if (WriteBottom) then
				STATUS = NF90_INQ_VARID(NCID, "HOB", HOBID)
				STATUS = NF90_PUT_VAR(NCID, HOBID, P_HOB,       &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put HOB, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
				do n=1,numpar
					P_HOB(n)=9999.0
				enddo
				STATUS = NF90_INQ_VARID(NCID, "bustr", bustrID)
				STATUS = NF90_PUT_VAR(NCID, bustrID, P_bustr,       &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put bustr, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
				
				STATUS = NF90_INQ_VARID(NCID, "bvstr", bvstrID)
				STATUS = NF90_PUT_VAR(NCID, bvstrID, P_bvstr,       &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put bvstr, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
		  endif
		  if (WriteWaterDepth) then
				STATUS = NF90_INQ_VARID(NCID, "WaterDepth", WDID)
				STATUS = NF90_PUT_VAR(NCID, WDID, par(:,pWD),       &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put Water Depth, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
				
			endif
#ifdef GROWTH
			
			STATUS = NF90_INQ_VARID(NCID, "size", sizeID)
			STATUS = NF90_PUT_VAR(NCID, sizeID,  par(:,pSize),       &
								  start = (/ 1, prcount /),   &
								  count = (/ numpar,  1 /))
			IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put size, time: ',time
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
#endif


		if ((Behavior.EQ.10) .OR. (Process_VA)) then
			if (prcount.GT.1) then
		!WRITING BEHAVIOR PARAMETERS INTO LAST TIME STEP, prcount-1
		
				STATUS = NF90_INQ_VARID(NCID, "acceleration", BehaveID)
				STATUS = NF90_PUT_VAR(NCID, BehaveID,  par(:,pAcc),       &
									start = (/ 1, prcount-1 /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put acceleration, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
			
				STATUS = NF90_INQ_VARID(NCID, "vorticity", BehaveID)
				STATUS = NF90_PUT_VAR(NCID, BehaveID,  par(:,pVort),       &
									start = (/ 1, prcount-1 /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put vorticity, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
			
			
				STATUS = NF90_INQ_VARID(NCID, "behave_w", BehaveID)
				STATUS = NF90_PUT_VAR(NCID, BehaveID,  par(:,pbehaveW),       &
								  start = (/ 1, prcount-1 /),   &
								  count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put behavior W, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
			
				STATUS = NF90_INQ_VARID(NCID, "swinksinkflag", BehaveID)
				STATUS = NF90_PUT_VAR(NCID, BehaveID,  par(:,pSSF),       &
								  start = (/ 1, prcount-1 /),   &
								  count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put Sink/Swim Flag, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
			endif 
		endif
		

		STATUS = NF_CLOSE(NCID)

	  END SUBROUTINE writeNetCDF
	  
	  
  SUBROUTINE writeModelInfo()
    !This subroutine simply writes model information to standard output
    USE PARAM_MOD
    IMPLICIT NONE

    CHARACTER(len=10) :: tmp !For Converting Integers to Characters
    CHARACTER(len=200) :: filenm
	INTEGER :: ng

    write(*,*) ' ******************** Model Info ******************** '
    write(*,*) ' '

    write(*,*) ' Run Name:              = ',TRIM(RunName)
    write(*,*) ' Executable Directory:  = ',TRIM(ExeDir)
    write(*,*) ' Output Directory:      = ',TRIM(OutDir)
    write(*,*) ' Run By:                = ',TRIM(RunBy)
    write(*,*) ' Institution:           = ',TRIM(Institution)
    write(*,*) ' Started On:            = ',TRIM(StartedOn)
    write(*,*) ' '

    write(tmp,'(F10.3)') days
    tmp = ADJUSTL(tmp)
    write(*,*) ' Days:                  = ',TRIM(tmp)
    write(tmp,'(I10)') numpar
    tmp = ADJUSTL(tmp)
    write(*,*) ' Particles:             = ',TRIM(tmp)
    write(*,*) ' Particle File:         = ',TRIM(parfile)
    write(*,*) ' '
    SELECT CASE(Behavior)
      CASE(0)
      write(*,*) ' Behavior:              = Passive'
      CASE(1)
      write(*,*) ' Behavior:              = Near-Surface'
      CASE(2)
      write(*,*) ' Behavior:              = Near-Bottom'
      CASE(3)
      write(*,*) ' Behavior:              = Diurnal Vertical Migration'
      CASE(4)
      write(*,*) ' Behavior:              = C.virginica oyster larvae'
      CASE(5)
      write(*,*) ' Behavior:              = C.ariakensis oyster larvae'
      CASE(6)
      write(*,*) ' Behavior:              = Constant sink/float'
      CASE(7)
      write(*,*) ' Behavior:              = Tidal Stream Transport'
    END SELECT

    if(mortality)then
      write(*,*) ' Particle Mortality:    = On'
    else
      write(*,*) ' Particle Mortality:    = Off'
    endif

    if(settlementon)then
      write(*,*) ' Settlement:            = On'
      write(*,*) ' Habitat File:          = ',TRIM(habitatfile)
      if(holesExist)write(*,*) ' Hole File:             = ',TRIM(holefile)
    else
      write(*,*) ' Settlement:            = Off'
    endif
    write(*,*) ' '

    if(HTurbOn)then
      write(*,*) ' Horizontal Turbulence: = On'
    else
      write(*,*) ' Horizontal Turbulence: = Off'
    endif
    if(VTurbOn)then
      write(*,*) ' Vertical Turbulence:   = On'
    else
      write(*,*) ' Vertical Turbulence:   = Off'
    endif
    if(OpenOceanBoundary)then
      write(*,*) ' Ocean Boundary:        = Open'
    else
      write(*,*) ' Ocean Boundary:        = Closed'
    endif
    if(SaltTempOn)then
      write(*,*) ' Salt & Temp Output:    = On'
    else
      write(*,*) ' Salt & Temp Output:    = Off'
    endif
    if(TrackCollisions)then
      write(*,*) ' Track Collisions:      = Yes'
    else
      write(*,*) ' Track Collisions:      = No'
    endif
    if(WriteModelTiming)then
      write(*,*) ' Track Model Timing:    = Yes'
    else
      write(*,*) ' Track Model Timing:    = No'
    endif
	do ng=1,Ngrid
		SELECT CASE(numdigits)
		  CASE(1)
			WRITE(filenm,'(A,I1.1,A)') TRIM(prefix(ng)),filenum,TRIM(suffix)
		  CASE(2)
			WRITE(filenm,'(A,I2.2,A)') TRIM(prefix(ng)),filenum,TRIM(suffix)
		  CASE(3)
			WRITE(filenm,'(A,I3.3,A)') TRIM(prefix(ng)),filenum,TRIM(suffix)
		  CASE(4)
			WRITE(filenm,'(A,I4.4,A)') TRIM(prefix(ng)),filenum,TRIM(suffix)
		  CASE(5)
			WRITE(filenm,'(A,I5.5,A)') TRIM(prefix(ng)),filenum,TRIM(suffix)
		  CASE(6)
			WRITE(filenm,'(A,I6.6,A)') TRIM(prefix(ng)),filenum,TRIM(suffix)
		  CASE(7)
			WRITE(filenm,'(A,I7.7,A)') TRIM(prefix(ng)),filenum,TRIM(suffix)
		  CASE(8)
			WRITE(filenm,'(A,I8.8,A)') TRIM(prefix(ng)),filenum,TRIM(suffix)
		  CASE DEFAULT
			WRITE(*,*) 'Model presently does not support numdigits of ',numdigits
			WRITE(*,*) 'Please use numdigit value from 1 to 8'
			WRITE(*,*) '  OR modify code in Hydrodynamic module'
			STOP
		END SELECT

		write(*,*) ' '
		write(*,*) ' First Hydro File:      = ',TRIM(filenm)

		write(*,*) ' '
		write(tmp,'(I10)') seed
		tmp = ADJUSTL(tmp)
		write(*,*) ' Seed:                  = ',TRIM(tmp)
		write(*,*) ' '
		if (Process_VA) then
			WRITE(*,*) 'Turbulence files: ',turbstd_v_a_prefix(ng)
			WRITE(*,*) 'Wave  files: ',wavestd_prefix(ng)
		endif
#ifdef STOKES
		WRITE(*,*) 'Stokes files: ',stokesprefix(ng)
#endif
		
	enddo 
  END SUBROUTINE writeModelInfo

 end program
