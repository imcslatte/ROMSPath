! ROMSPath - Offlinwe PArticle tracking model v1.1                               
! Date: 25 OCtober 2021
!
! Description: ROMSPath is an 
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
!   Rutgers The State University of New Jersey
!   Department of Marine and Coastal Sciences
!   New Brunswick, NJ 08901 USA
!
! Funding was provided by the National Science Foundation Biological 
! and Physical Oceanography Programs, Maryland Department of Natural 
! Resources, NOAA Chesapeake Bay Office, NOAA Maryland Sea Grant College 
! Program, & NOAA-funded UMCP Advanced Study Institute for the Environment. 
! 
! **********************************************************************
! **********************************************************************
! **                      Copyright (c) 2019                       **
! **   								    **
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
! **  	                     **
! **                                                                  **
! ** We ask that users make appropriate acknowledgement of            **
! ** The University of Maryland Center for Environmental Science,     **
! ** individual developers, participating agencies and institutions,  **
! ** and funding agencies. One way to do this is to cite one or       **
! ** more of the relevant publications listed at:                     **
! **                                                                  **
! ** 	          **
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
! Last Modified on:     25 October 2021
! ROMSPath Version: 1.0.1

IMPLICIT NONE
!   *************************************************************************
!   *                                                                       *
!   *                       Variable Declarations                           *
!   *                                                                       *
!   *************************************************************************

  INTEGER, PARAMETER :: nAttrib   = 22

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
  INTEGER, PARAMETER :: pZeta 	  = 21  ! ROMS Zeta
  INTEGER, PARAMETER :: pBath 	  = 22  ! ROMS Bathymetry (h)

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: par
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: P_Salt,P_Temp,mean_salt,mean_temp  
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: P_Light,mean_light,P_live_biofoul,P_dead_biofoul,P_total_biofoul, B_encounter,B_growth,B_grazing,B_mort,B_remin,P_settling_vel,mean_biofoul_live,mean_biofoul_dead,mean_encounter,mean_growth,mean_grazing,mean_mort,mean_remin
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
                      WriteHeaders,WriteModelTiming,ErrorFlag,getParams,Behavior,&
                      LightOn,BiofoulOn,LightMean,BiofoulMean, &
                      rhop
	use grid_mod,   only: InitGrid,GRIDS
	use INT_MOD,   only: LL2ij,inside
    use bf_mod,        only:biofoul_subr,settling_vel_func
    integer :: n,istat,ng,i,j,nmask
    logical :: ingrid,obound,inzgrid
    double precision, allocatable, dimension(:) :: pLon,pLat,Ipar,Jpar
	CHARACTER(len=200) :: filenm

    integer :: in_island,inbounds,test
    double precision:: tdepth,zeta,dpran
    double precision :: settling_vel
	

	
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

    IF(LightOn)THEN  
      ALLOCATE(P_Light(numpar))
      P_Light = 0.0
	   IF(LightMean)THEN
		ALLOCATE(mean_light(numpar))
		mean_light = 0.0
	   ENDIF
    ENDIF
    
    IF(BiofoulOn)THEN   
      ALLOCATE(P_live_biofoul(numpar))
      P_live_biofoul = 0.0
      ALLOCATE(P_dead_biofoul(numpar))
      P_dead_biofoul = 0.0
      ALLOCATE(P_total_biofoul(numpar))
      P_total_biofoul = 0.0
      ALLOCATE(B_encounter(numpar))
      B_encounter = 0.0
      ALLOCATE(B_growth(numpar))
      B_growth = 0.0
      ALLOCATE(B_grazing(numpar))
      B_grazing = 0.0
      ALLOCATE(B_mort(numpar))
      B_mort = 0.0
      ALLOCATE(B_remin(numpar))
      B_remin = 0.0
      ALLOCATE(P_settling_vel(numpar))
      P_settling_vel = 0.0
	   IF(BiofoulMean)THEN
		ALLOCATE(mean_biofoul_live(numpar))
		mean_biofoul_live = 0.0
		ALLOCATE(mean_biofoul_dead(numpar))
		mean_biofoul_dead = 0.0
		ALLOCATE(mean_encounter(numpar))
		mean_encounter = 0.0
		ALLOCATE(mean_growth(numpar))
		mean_growth = 0.0
		ALLOCATE(mean_grazing(numpar))
		mean_grazing = 0.0
		ALLOCATE(mean_mort(numpar))
		mean_mort = 0.0
		ALLOCATE(mean_remin(numpar))
		mean_remin = 0.0
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
        par(n,ppX) = 1.0    
        par(n,ppY) = 1.0    
        par(n,ppZ) = 1.0    
        par(n,pStatus)   = 9.0
        par(n,pAge)      = 0.0
        par(n,pLifespan) = 0.0
        par(n,pGID)      = dble(Ngrid)
        par(n,pSize)      = initsize
        par(n,pAcc) = 0.0
        par(n,pVort) = 0.0
        par(n,pbehaveW) = 0.0
        par(n,pSSF) = 0.0
        par(n,pWD) = 9999.0
        par(n,pBath) = 9999.0
        par(n,pZeta) = 9999.0
		
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
					else	
					
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
					 par(n,pZeta)=zeta
					 par(n,pBath)=(DBLE(-1.0)*tdepth)
					 if (SaltTempMean) then
						mean_salt(n)=P_salt(n)
						mean_temp(n)=P_temp(n)
						mI=1
					 endif
					 
				 endif
		
			enddo
		endif

	
		if (LightOn) then
			do n=1,numpar
				if (isIn(n)) then
				
					tdepth = DBLE(-1.0)* getInterp2D("depth",int(par(n,pGID)),par(n,pX),par(n,pY),t_c)
					zeta =  getInterp2D("zeta",int(par(n,pGID)),par(n,pX),par(n,pY),t_c)
					 P_light(n)=getInterp3d("light",int(par(n,pGID)),par(n,pX),par(n,pY),par(n,pZ),t_c,1,zeta,tdepth)
					 par(n,pWD)=(DBLE(-1.0)*tdepth)+zeta
					 par(n,pZeta)=zeta
					 par(n,pBath)=(DBLE(-1.0)*tdepth)
					 if (LightMean) then
						mean_light(n)=P_light(n)
						mI=1
					 endif
					 
				 endif
		
			enddo
		endif

                if (BiofoulOn) then
                        do n = 1,numpar
                                if (isIn(n)) then
                                        P_settling_vel(n) = settling_vel_func(par(n,pSize),P_total_biofoul(n),int(par(n,pgid)),par(n,px),par(n,py),par(n,pz),t_c)
                                endif
                        enddo
                endif



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
     

    ! !Deallocate local variables
     DEALLOCATE(pLon,pLat)

 

  end subroutine ini_ROMSPath

  

  subroutine run_External_Timestep()
    use param_mod, only: dt,idt,WriteModelTiming,tdim,t_b,t_c,t_f,filenum,&
		numdigits,prefix,suffix,tstep,multifile
    use hydro_mod, only: updateHydro,HYDRODATA
	USE GRID_MOD, ONLY: reftime
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

 

      !Read in hydrodynamic model data 
	 
      IF(ets > 2) then			
	     t_b  = mod(t_b,3)+1  ! 1 -> 2 -> 3 -> 1
	     t_c  = mod(t_c,3)+1  ! 2 -> 3 -> 1 -> 2
	     t_f = mod(t_f,3)+1  ! 3 -> 1 -> 2 -> 3
		CALL updateHydro(.FALSE.,tstep,t_f)   !do not start updating until 3rd iteration
	  endif

	  
      !Prepare external time step values to be used for 
      !  calculating Advection and Turbulence
      ex=0.0
      ex(1) = (ets-2)*dt
      ex(2) = (ets-1)*dt
      ex(3) = ets*dt
	 ! call CPU_TIME(before)
	 call CPU_TIME(before)

       do its=1,stepIT
		

        call run_Internal_Timestep()
       enddo !ITloop
	  tstep=tstep+1

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
							  idt,dt,HTurbOn,VTurbOn,settlementon,xi_v,eta_v,     &
                              Behavior,SaltTempOn,LightOn,BiofoulOn,OpenOceanBoundary,NoBounce,Swimdepth, &
                              TrackCollisions,WriteModelTiming,mortality,      &
                              ErrorFlag,t_b,t_c,t_f,Ngrid,vertdist,scheme,nsb, &
			      SaltTempMean,LightMean,BiofoulMean,WriteBottom,maxsize,Process_VA, &
                              rhop,nBeachCells,YesProbBeach      !rhop							  
    !USE SETTLEMENT_MOD, ONLY: isSettled,testSettlement
#ifdef GROWTH
	USE GROWTH_MOD,   ONLY:  growlarva
#endif

    USE BEHAVIOR_MOD,   ONLY: behave
    USE BOUNDARY_MOD,   ONLY: bounds,IsBeachPossible,WhetherToBeach,Unbeach
    USE GRID_MOD,       ONLY: getSlevel,getWlevel,GRIDS
    USE HTURB_MOD,      ONLY: HTurb
    USE VTURB_MOD,      ONLY: VTurb
    USE ADVECTION_MOD,  ONLY: RKAdvect
    USE INT_MOD,        ONLY: getinterp2d,getinterp3d,polintd,sinintd,getInterpStr
    USE BF_MOD,        ONLY: biofoul_subr,settling_vel_func
    IMPLICIT NONE

    ! Iteration Variables
    INTEGER :: i,deplvl,n

    ! Particle tracking
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: Pwc_zb,Pwc_zc,Pwc_zf
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: Pwc_wzb,Pwc_wzc,Pwc_wzf
    DOUBLE PRECISION :: Xpar,Ypar,Zpar,newXpos,newYpos,newZpos,P_zb,P_zc,P_zf, &
      P_zeta,ey(3)
    
    ! Behavior and Turbulence
    DOUBLE PRECISION :: TurbHx,TurbHy,TurbV,Behav,XBehav,YBehav,ZBehav
    LOGICAL :: bott   ! for Behavior 7 along with XBehav,YBehav,ZBehav

    ! Boundaries
    INTEGER :: intersectf,skipbound,inbounds,reflects,inpoly,nmask,Inode,Jnode &
	  ,ngid,m(2*nBeachCells,2*nBeachCells)
    DOUBLE PRECISION :: reflect,fintersectX,fintersectY,freflectX,freflectY,   &
      Xpos,Ypos,nXpos,nYpos,pm,pn,parLon,parLat,coord_close(2,2),node_dists(2)
    LOGICAL :: ingrid,obound,hbot,htop,InBeachWindow,yesBeach

    ! Advection
    DOUBLE PRECISION :: AdvectX,AdvectY,AdvectZ,maxpartdepth,minpartdepth,     &
      kn1_u,kn1_v,kn1_w,kn2_u,kn2_v,kn2_w,kn3_u,kn3_v,kn3_w,kn4_u,kn4_v,kn4_w, &
      P_V,P_U,P_W,UAD,VAD,WAD,x1,x2,x3,y1,y2,y3,z1,z2,z3,btemp,tbustr,tbvstr


	  
	DOUBLE PRECISION :: tempX,tempY,tdepth,zeta,zetab,zetac,zetaf,behout(4)
	 
     ! Biofouling
     DOUBLE PRECISION :: settling_vel,daylength_val  
     INTEGER :: ng
   
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
        if (NoBounce) then
                if (ABS(par(n,pStatus)- 1.2) < 0.01) then
                        cycle
                elseif (ABS(par(n,pStatus)+9.) < 0.01) then
                        cycle
                endif
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

	   

        !moved below beaching algorithm 
      !! ! *********************************************************
      !! ! *                                                       *
      !! ! *                       Behavior                        *
      !! ! *                                                       *
      !! ! *********************************************************


      !   
      !     call CPU_TIME(times(4))
      !     !write(*,*) times(4)-times(3)
      !     
      !  	 
      !  	CALL behave(Xpar,Ypar,Zpar,XBehav,YBehav,ZBehav,par(n,pSize),ex,ix,int(par(n,pGID)),P_total_biofoul(n),behout)
      !  	par(n,pAcc)=behout(1)
      !  	par(n,pVort)=behout(2)
      !  	par(n,pbehaveW)=behout(3)
      !  	par(n,pSSF)=behout(4)
      !          WRITE(*,*) 'ZBehav: ', ZBehav 
      !          WRITE(*,*) 'P_total_biofoul: ', P_total_biofoul(n)

	
	 
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

if (YesProbBeach) then
        !Probabilistic beaching
        !Note that beaching is not fully integrated with nested grids
        !check if beached/in beaching window
                CALL IsBeachPossible(newXPos,newYPos,int(par(n,pGID)),coord_close,node_dists,m,InBeachWindow)
    

        !determine whether to beach
        yesBeach = .FALSE.
        if (InBeachWindow) then
                yesBeach = WhetherToBeach()
        endif

        !if beached but shouldn't be, unbeach
        !if not beached but should be, beach
        if (InBeachWindow .AND. .NOT.ingrid .AND. .NOT.yesBeach) then
                CALL Unbeach(int(par(n,pGID)),par(n,pX),par(n,pY),newXPos,newYPos,coord_close,node_dists,m,tempX,tempY)
                ingrid = .TRUE.
        elseif (ingrid .AND. yesBeach) then
                !not necessary to find exact location since normal beaching doesn't update final location?
                !update status
                ingrid = .FALSE.
                obound = .TRUE.
        endif
endif
      ! ! *********************************************************
      ! ! *                                                       *
      ! ! *                       Behavior                        *
      ! ! *                                                       *
      ! ! *********************************************************


	 
	   call CPU_TIME(times(4))
	   !write(*,*) times(4)-times(3)
	   
		 
		CALL behave(Xpar,Ypar,Zpar,XBehav,YBehav,ZBehav,par(n,pSize),ex,ix,int(par(n,pGID)),P_total_biofoul(n),behout)
		par(n,pAcc)=behout(1)
		par(n,pVort)=behout(2)
		par(n,pbehaveW)=behout(3)
		par(n,pSSF)=behout(4)

      ! ! *********************************************************
      ! ! *                                                       *
      ! ! *     Update Particle Locations and Check Boundaries in Z   *
      ! ! *                                                       *
      ! ! *********************************************************

	  tdepth = DBLE(-1.0)* getInterp2D("depth",int(par(n,pGID)),tempX,tempY,t_c)
	  
	  ey(1) =  DBLE(1.0)*getInterp2D("zeta",int(par(n,pGID)),tempX,tempY,t_b)
	  ey(2) =  DBLE(1.0)*getInterp2D("zeta",int(par(n,pGID)),tempX,tempY,t_c)
	  ey(3) =  DBLE(1.0)*getInterp2D("zeta",int(par(n,pGID)),tempX,tempY,t_f)
			  
	  P_zeta=polintd(ex,ey,3,ix(2))

	 par(n,pWD)=(DBLE(-1.0)*tdepth)+P_zeta
	 par(n,pZeta)=P_zeta
	 par(n,pBath)=(DBLE(-1.0)*tdepth)
	 
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
			newZpos = tdepth+vertdist
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


	  ! ! *********************************************************
      ! ! *                                                       *
      ! ! *                    Light 	                          *
      ! ! *                                                       *
      ! ! *********************************************************
	  
	   call CPU_TIME(times(6))
		if (LightOn) then
				
				zetac =  getInterp2D("zeta",int(par(n,pGID)),par(n,pX),par(n,pY),t_c)
				zetaf =  getInterp2D("zeta",int(par(n,pGID)),par(n,pX),par(n,pY),t_c)
			  
				ey(2)=getInterp3d("light",int(par(n,pGID)),par(n,pX),par(n,pY),par(n,pZ),t_c,1,zetac,tdepth)	
				ey(3)=getInterp3d("light",int(par(n,pGID)),par(n,pX),par(n,pY),par(n,pZ),t_f,1,zetaf,tdepth)	
                                ng = int(par(n,pGID)) 
                                !daylength_val = HYDRODATA(int(ng)%daylength(t_c)
                                !daylength_val = daylength_val*60.D0*60.D0       !convert units

                                if (ey(2) > 0) then !.and. (ix(2) < (ex(2)+daylength_val/2.D0))) then ! afternoon                                 
                                        P_light(n)=sinintd(int(par(n,pGID)),ex(2),ey(2),ix(2),t_c)
                                        !if (i .EQ. 1) then
                                        !if (ex(2) > 0) then 
                                        !       WRITE(*,*) 'in if 1'
                                        !       WRITE(*,*) 'P_light(n):', P_light(n)
                                        ! endif
                                else if (ey(3) > 0) then !.and. (ix(2) >= (ex(2)+daylength_val/2.D0))) then ! morning
                                        P_light(n)=sinintd(int(par(n,pGID)),ex(3),ey(3),ix(2),t_c)
                                        !if (i .EQ. 1) then
                                         !     WRITE(*,*) 'ex(3)', ex(3)
                                          !    WRITE(*,*) 'ey(3)', ey(3)
                                           !   WRITE(*,*) 'ix(2)', ix(2)
                                            !  WRITE(*,*) 'P_light(n): ',P_light(n)
                                       ! endif
                                else
                                        P_light(n) = 0.
                                endif
			        
	
				if (LightMean) then
					mean_light(n)=mean_light(n)+P_light(n)
				endif
					 
			    
		endif

	  ! ! *********************************************************
      ! ! *                                                         *
      ! ! *                    Biofoul		                        *
      ! ! *                                                         *
      ! ! *********************************************************
	  
	   call CPU_TIME(times(6))
		if (BiofoulOn) then
			
                                ! get PP (primary productivity, used as growth
                                ! rate)	
				!zetac =  getInterp2D("zeta",int(par(n,pGID)),par(n,pX),par(n,pY),t_c)
			  
				!ey(2)=getInterp3d("PP",int(par(n,pGID)),par(n,pX),par(n,pY),par(n,pZ),t_c,1,zetac,tdepth)	
				!ey(3)=getInterp3d("PP",int(par(n,pGID)),par(n,pX),par(n,pY),par(n,pZ),t_f,1,zetac,tdepth)	
                                

                                !if (isnan(ey(2)) .and. isnan(ey(3))) then
                                !        PP = 0
                                !if (ey(2) > 0) then !.and. (ix(2) < (ex(2)+daylength_val/2.d0))) then ! afternoon  was dt/2 instead of daylength               
                                  !      WRITE(*,*) "ey(2)", ey(2)
                                   !     WRITE(*,*) "ey(2) > 0"
                                !        PP=sinintd(int(par(n,pGID)),ex(2),ey(2),ix(2),t_c)
                                !else if (ey(3) > 0) then !.and. (ix(2) >= (ex(2)+daylength_val/2.D0))) then ! morning
                                !        PP=sinintd(int(par(n,pGID)),ex(3),ey(3),ix(2),t_c)
                                    !    WRITE(*,*) "ey(3) > 0"
                               ! else
                                !        PP=0.d0
                                !endif

                                P_settling_vel(n) = settling_vel_func(par(n,pSize),P_total_biofoul(n),int(par(n,pgid)),par(n,px),par(n,py),par(n,pz),t_c)
                                !calculate biofouling! (cells/m^2 plastic)
                                CALL biofoul_subr(P_live_biofoul(n),P_dead_biofoul(n),P_total_biofoul(n),par(n,pSize),P_settling_vel(n),idt,int(par(n,pgid)),par(n,px),par(n,py),par(n,pz),ix,ex,P_light(n),B_encounter(n),B_growth(n),B_grazing(n),B_mort(n),B_remin(n)) 

				if (BiofoulMean) then
					mean_Biofoul_live(n)=mean_biofoul_live(n)+P_live_biofoul(n)
					mean_Biofoul_dead(n)=mean_biofoul_dead(n)+P_dead_biofoul(n)
                                        mean_encounter(n) = mean_encounter(n)+B_encounter(n)
                                        mean_growth(n) = mean_growth(n)+B_growth(n)
                                        mean_grazing(n) = mean_grazing(n)+B_grazing(n)
                                        mean_mort(n) = mean_mort(n)+B_mort(n)
                                        mean_remin(n) = mean_remin(n)+B_remin(n)
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
    use param_mod,   only: numpar,SaltTempOn,LightOn,BiofoulOn,TrackCollisions,stokesprefix,turbstd_v_a_prefix
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
			RunName,ExeDir,OutDir,RunBy,Institution,StartedOn,SaltTempMean,LightMean,BiofoulMean,         &
			TrackCollisions,SaltTempOn,LightOn,BiofoulOn,Ngrid,days,idt,VTurbOn,HTurbOn,deltat,      &
			serr,smth,sub,AKSback,maxsize,tempcut,initsize,deadage,		&
			a0,a1,a2,a3,a4,a5,a6,a7,a8,TempOffset,WriteBottom,WriteWaterDepth,Behavior,		&
			vort_cr,vort_sat,b0pv,b1pv,b0wv,b1w,acc_cr,acc_sat,&
			b0pa,b1pa,b0wa,va_flag,OpenOceanBoundary,NoBounce,swimfast,Process_VA,	&
			WriteWaterDepth,seed,WriteZeta,WriteBath
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
				   HOBID,bustrID,bvstrID,VORTID,ACCID,BWID,SSFID,WDID,BATHID,ZETAID, &
                                   lightID,biofoulLiveID,biofoulDeadID,encounterID,growthID,grazingID,mortID,reminID,settling_velID,sizeID 
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

			STATUS = NF90_DEF_VAR(NCID,'zp',NF_FLOAT,(/numparID,timeID/),      &
								  depthID,deflate_level=1,shuffle=.true.)
			IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: zp var'
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
			
			IF(LightOn)THEN
			  STATUS = NF90_DEF_VAR(NCID,'light',NF_FLOAT,(/numparID,timeID/), &
									lightID,deflate_level=1,shuffle=.true.)
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
												  'Light var'
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			ENDIF

			IF(BiofoulOn)THEN
			  STATUS = NF90_DEF_VAR(NCID,'biofoul_live',NF_DOUBLE,(/numparID,timeID/), &
									biofoulLiveID,deflate_level=1,shuffle=.true.)
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
												  'Living biofoul var'
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_DEF_VAR(NCID,'biofoul_dead',NF_DOUBLE,(/numparID,timeID/), &
									biofoulDeadID,deflate_level=1,shuffle=.true.)
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
												  'Dead biofoul var'
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_DEF_VAR(NCID,'encounter',NF_DOUBLE,(/numparID,timeID/), &
									encounterID,deflate_level=1,shuffle=.true.)
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
												  'Encounter var'
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_DEF_VAR(NCID,'growth',NF_DOUBLE,(/numparID,timeID/), &
									growthID,deflate_level=1,shuffle=.true.)
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
												  'Growth var'
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_DEF_VAR(NCID,'grazing',NF_DOUBLE,(/numparID,timeID/), &
									grazingID,deflate_level=1,shuffle=.true.)
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
												  'Grazing var'
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_DEF_VAR(NCID,'mort',NF_DOUBLE,(/numparID,timeID/), &
									mortID,deflate_level=1,shuffle=.true.)
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
												  'Mort var'
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_DEF_VAR(NCID,'remineralization',NF_DOUBLE,(/numparID,timeID/), &
									reminID,deflate_level=1,shuffle=.true.)
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
												  'Remineralization var'
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


			  STATUS = NF90_DEF_VAR(NCID,'settling_vel',NF_DOUBLE,(/numparID,timeID/), &
									settling_velID,deflate_level=1,shuffle=.true.)
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
												  'Settling_vel var'
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
			
			
			IF(WriteZeta)THEN
			  STATUS = NF90_DEF_VAR(NCID,'zeta',NF_FLOAT,(/numparID,timeID/), &
									ZETAID,deflate_level=1,shuffle=.true.)
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
												  'Zeta var'
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			ENDIF
			
			IF(WriteBath)THEN
			  STATUS = NF90_DEF_VAR(NCID,'h',NF_FLOAT,(/numparID,timeID/), &
									BATHID,deflate_level=1,shuffle=.true.)
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
												  'H var'
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
			STATUS = NF90_PUT_ATT(NCID, depthID, "long_name", "vertical position of the particle in ROMS z-coordinates, zp=0 is mean sea level, positive is up.")
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_PUT_ATT(NCID, depthID, "units", "meters")
			IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			STATUS = NF90_PUT_ATT(NCID, depthID, "field", "zp, scalar, series")
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

			  STATUS = NF90_PUT_ATT(NCID, tempID, "units", "� Celsius")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_PUT_ATT(NCID, tempID, "field",                         &
									"temperature, scalar, series")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			ENDIF
						
                        IF(LightOn)THEN
			  !light
			  STATUS = NF90_PUT_ATT(NCID, lightID, "long_name",                     &
									"Light at the particle's location")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS=NF90_PUT_ATT(NCID,lightID, "field", "light, scalar, series")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_PUT_ATT(NCID, lightID, "units", "watt m-2")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			ENDIF
                        
                        IF(BiofoulOn)THEN
			  !biofouling
			  STATUS = NF90_PUT_ATT(NCID, biofoulLiveID, "long_name",                     &
									"Living biofouling on the particle")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS=NF90_PUT_ATT(NCID,biofoulLiveID, "field", "live biofouling, scalar, series")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_PUT_ATT(NCID, biofoulLiveID, "units", "number/m^2")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_PUT_ATT(NCID, biofoulDeadID, "long_name",                     &
									"Dead biofouling on the particle")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS=NF90_PUT_ATT(NCID,biofoulDeadID, "field", "dead biofouling, scalar, series")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_PUT_ATT(NCID, biofoulDeadID, "units", "number/m^2")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_PUT_ATT(NCID, encounterID, "long_name",                     &
									"Encounter term in Kooi model")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS=NF90_PUT_ATT(NCID,encounterID, "field", "encounter, scalar, series")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_PUT_ATT(NCID, encounterID, "units", "number/m^2")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


			  STATUS = NF90_PUT_ATT(NCID, growthID, "long_name",                     &
									"Growth term in Kooi model")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS=NF90_PUT_ATT(NCID,growthID, "field", "growth, scalar, series")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_PUT_ATT(NCID, growthID, "units", "number/m^2")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


			  STATUS = NF90_PUT_ATT(NCID, grazingID, "long_name",                     &
									"Grazing term in Kooi model")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS=NF90_PUT_ATT(NCID,grazingID, "field", "grazing, scalar, series")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_PUT_ATT(NCID, grazingID, "units", "number/m^2")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


			  STATUS = NF90_PUT_ATT(NCID, mortID, "long_name",                     &
									"Mort term in Kooi model")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS=NF90_PUT_ATT(NCID,mortID, "field", "mort, scalar, series")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_PUT_ATT(NCID, mortID, "units", "number/m^2")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


			  STATUS = NF90_PUT_ATT(NCID, reminID, "long_name",                     &
									"Remineralization term in Kooi model")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS=NF90_PUT_ATT(NCID,reminID, "field", "remin, scalar, series")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_PUT_ATT(NCID, reminID, "units", "number/m^2")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)



			  STATUS = NF90_PUT_ATT(NCID, settling_velID, "long_name",                     &
									"Settling velocity of the particle")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS=NF90_PUT_ATT(NCID,settling_velID, "field", "settling velocity, scalar, series")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS = NF90_PUT_ATT(NCID, settling_velID, "units", "m/s")
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
			
			  STATUS = NF90_PUT_ATT(NCID, WDID, "long_name",                     &
									"Local Water Depth at Particle Location")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS=NF90_PUT_ATT(NCID,WDID, "field", "WaterDepth, scalar, series")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
				
			  STATUS = NF90_PUT_ATT(NCID, WDID, "units", "m")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			ENDIF
			
			IF(WriteZeta)THEN
			
			  STATUS = NF90_PUT_ATT(NCID, ZETAID, "long_name",                     &
									"free surface displacement at particle locations")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS=NF90_PUT_ATT(NCID,ZETAID, "field", "zeta, scalar, series")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
				
			  STATUS = NF90_PUT_ATT(NCID, ZETAID, "units", "m")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			ENDIF
			
			
			IF(WriteBath)THEN
			
			  STATUS = NF90_PUT_ATT(NCID, BATHID, "long_name",                     &
									"bathymetry at particle locations")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

			  STATUS=NF90_PUT_ATT(NCID,BATHID, "field", "H, scalar, series")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
				
			  STATUS = NF90_PUT_ATT(NCID, BATHID, "units", "m")
			  IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			ENDIF
			
			
			
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
			
			If( NoBounce) then
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "NoBounce", 'TRUE')
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


			if (LightMean) then
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "LightMean", 'TRUE')
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			else
			
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "LightMean", 'FALSE')
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			endif
			
            if (BiofoulMean) then
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "BiofoulMean", 'TRUE')
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			else
			
				STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "BiofoulMean", 'FALSE')
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
			endif
                        
                        
			
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
		USE PARAM_MOD, ONLY: numpar,SaltTempOn,LightOn,BiofoulOn,NCOutFile,outpath,outpathGiven,     &
			NCtime,TrackCollisions,Ngrid,SaltTempMean,LightMean,BiofoulMean,WriteBottom,WriteWaterDepth,Behavior,		&
			Process_VA,WriteWaterDepth,WriteZeta,WriteBath
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
				   statusID,saltID,tempID,ngID,n,HOBID,bustrID,bvstrID,BehaveID,WDID, &
				   ZETAID,BATHID,lightID, &
                           biofoulLiveID,biofoulDeadID,encounterID,growthID,grazingID,mortID,reminID,settling_velID,sizeID
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
		  STATUS = NF90_INQ_VARID(NCID, "zp", depthID)
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
		  
		  !Light
		  if (LightOn) then
		  	if (LightMean) then
				do n=1,numpar
					mean_light(n)=mean_light(n)/mI
				enddo
					
				STATUS = NF90_INQ_VARID(NCID, "light", lightID)
				STATUS = NF90_PUT_VAR(NCID, lightID, mean_light,       &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put light, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

			    do n=1,numpar
					mean_light(n)=0.0
				enddo

			else
				STATUS = NF90_INQ_VARID(NCID, "light", lightID)
				STATUS = NF90_PUT_VAR(NCID, lightID, P_light,       &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put light, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

			ENDIF
		  ENDIF
		  
                  !Biofouling
		  if (BiofoulOn) then
		  	if (BiofoulMean) then
				do n=1,numpar
					mean_biofoul_live(n)=mean_biofoul_live(n)/mI
					mean_biofoul_dead(n)=mean_biofoul_dead(n)/mI
                                        mean_encounter(n)=mean_encounter(n)/mI
                                        mean_growth(n)=mean_growth(n)/mI
                                        mean_grazing(n)=mean_grazing(n)/mI
                                        mean_mort(n)=mean_mort(n)/mI
                                        mean_remin(n)=mean_remin(n)/mI
				enddo
					
				STATUS = NF90_INQ_VARID(NCID, "biofoul_live", biofoulLiveID)
				STATUS = NF90_PUT_VAR(NCID, biofoulLiveID, mean_biofoul_live,       &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put live biofouling, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

				STATUS = NF90_INQ_VARID(NCID, "biofoul_dead", biofoulDeadID)
				STATUS = NF90_PUT_VAR(NCID, biofoulDeadID, mean_biofoul_dead,       &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put biofouling, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

				STATUS = NF90_INQ_VARID(NCID, "encounter", encounterID)
				STATUS = NF90_PUT_VAR(NCID, encounterID, mean_encounter,       &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put encounter, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

				STATUS = NF90_INQ_VARID(NCID, "growth", growthID)
				STATUS = NF90_PUT_VAR(NCID, growthID, mean_growth,       &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put growth, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

				STATUS = NF90_INQ_VARID(NCID, "grazing", grazingID)
				STATUS = NF90_PUT_VAR(NCID, grazingID, mean_grazing,       &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put grazing, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

				STATUS = NF90_INQ_VARID(NCID, "mort", mortID)
				STATUS = NF90_PUT_VAR(NCID, mortID, mean_mort,       &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put mort, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

				STATUS = NF90_INQ_VARID(NCID, "remineralization", reminID)
				STATUS = NF90_PUT_VAR(NCID, reminID, mean_remin,       &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put remineralization, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

			    do n=1,numpar
					mean_biofoul_live(n)=0.0
					mean_biofoul_dead(n)=0.0
                                        mean_encounter(n)=0.0
                                        mean_growth(n)=0.0
                                        mean_grazing(n)=0.0
                                        mean_mort(n)=0.0
                                        mean_remin(n)=0.0
				enddo
				

			else
				STATUS = NF90_INQ_VARID(NCID, "biofoul_live", biofoulLiveID)
				STATUS = NF90_PUT_VAR(NCID, biofoulLiveID, P_live_biofoul,       &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put live biofoul, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

				STATUS = NF90_INQ_VARID(NCID, "biofoul_dead", biofoulDeadID)
				STATUS = NF90_PUT_VAR(NCID, biofoulLiveID, P_dead_biofoul,       &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put dead biofoul, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

				STATUS = NF90_INQ_VARID(NCID, "encounter", encounterID)
				STATUS = NF90_PUT_VAR(NCID, encounterID, B_encounter,       &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put encounter, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

				STATUS = NF90_INQ_VARID(NCID, "growth", growthID)
				STATUS = NF90_PUT_VAR(NCID, growthID, B_growth,       &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put growth, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

				STATUS = NF90_INQ_VARID(NCID, "grazing", grazingID)
				STATUS = NF90_PUT_VAR(NCID, grazingID, B_grazing,       &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put grazing, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

				STATUS = NF90_INQ_VARID(NCID, "mort", mortID)
				STATUS = NF90_PUT_VAR(NCID, mortID, B_mort,       &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put mort, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

				STATUS = NF90_INQ_VARID(NCID, "remineralization", reminID)
				STATUS = NF90_PUT_VAR(NCID, reminID, B_remin,       &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put remineralization, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

			ENDIF
				STATUS = NF90_INQ_VARID(NCID, "settling_vel", settling_velID)
				STATUS = NF90_PUT_VAR(NCID, settling_velID, P_settling_vel,       &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put settling_vel, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
		  ENDIF

		  mI=0        		  
		  
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
		  if (WriteZeta) then
				STATUS = NF90_INQ_VARID(NCID, "zeta", ZETAID)
				STATUS = NF90_PUT_VAR(NCID, ZETAID, par(:,pZeta),       &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put Zeta, time: ',time
				IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
				
		  endif
		  if (WriteBath) then
				STATUS = NF90_INQ_VARID(NCID, "h", BATHID)
				STATUS = NF90_PUT_VAR(NCID, BATHID, par(:,pBath),       &
									start = (/ 1, prcount /),   &
									count = (/ numpar,  1 /))
				IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put h, time: ',time
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
	USE HYDRO_MOD, ONLY: getFileNames
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
    if(NoBounce)then
      write(*,*) ' Bounces Allowed:        = No'
    else
      write(*,*) ' Bounces Allowed:        = Yes'
    endif
    if(SaltTempOn)then
      write(*,*) ' Salt & Temp Output:    = On'
    else
      write(*,*) ' Salt & Temp Output:    = Off'
    endif
    if(LightOn)then
      write(*,*) ' Light Output:    = On'
    else
      write(*,*) ' Light Output:    = Off'
    endif
    if(BiofoulOn)then
      write(*,*) ' Biofoul Output:    = On'
    else
      write(*,*) ' Biofoul Output:    = Off'
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
		call getFileNames(filenm,prefix(ng),filenum)
		! SELECT CASE(numdigits)
		  ! CASE(1)
			! WRITE(filenm,'(A,I1.1,A)') TRIM(prefix(ng)),filenum,TRIM(suffix)
		  ! CASE(2)
			! WRITE(filenm,'(A,I2.2,A)') TRIM(prefix(ng)),filenum,TRIM(suffix)
		  ! CASE(3)
			! WRITE(filenm,'(A,I3.3,A)') TRIM(prefix(ng)),filenum,TRIM(suffix)
		  ! CASE(4)
			! WRITE(filenm,'(A,I4.4,A)') TRIM(prefix(ng)),filenum,TRIM(suffix)
		  ! CASE(5)
			! WRITE(filenm,'(A,I5.5,A)') TRIM(prefix(ng)),filenum,TRIM(suffix)
		  ! CASE(6)
			! WRITE(filenm,'(A,I6.6,A)') TRIM(prefix(ng)),filenum,TRIM(suffix)
		  ! CASE(7)
			! WRITE(filenm,'(A,I7.7,A)') TRIM(prefix(ng)),filenum,TRIM(suffix)
		  ! CASE(8)
			! WRITE(filenm,'(A,I8.8,A)') TRIM(prefix(ng)),filenum,TRIM(suffix)
		  ! CASE DEFAULT
			! WRITE(*,*) 'Model presently does not support numdigits of ',numdigits
			! WRITE(*,*) 'Please use numdigit value from 1 to 8'
			! WRITE(*,*) '  OR modify code in Hydrodynamic module'
			! STOP
		! END SELECT

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
