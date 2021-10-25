!
!This file defines many variables that are read from the GRID.data and ROMSPath.data files
!And also groups the input parameters/variables into namelists
!
!NOTE: variables in a namelist can NOT be dynamic variables!!!!! 
!      dynamic namelists are NOT yet supported in FORTRAN90/95 standard
!

  INTEGER,parameter :: MAXNgrid=5
  DOUBLE PRECISION,parameter :: Eradius = 6371315.0D0
  DOUBLE PRECISION, parameter :: pi = 3.14159265358979323846D0
  DOUBLE PRECISION, parameter :: deg2rad = pi / 180.0D0
  DOUBLE PRECISION, parameter :: rad2deg =180.0D0 / pi
  DOUBLE PRECISION, parameter :: g =980D0 ! cm s^-2 Gravity
  DOUBLE PRECISION, parameter :: rhof =1.026D0 ! (g cm^-3) Density
  DOUBLE PRECISION, parameter :: nu =0.01D0 !  (cm^2 s^-1) kinematic viscosity of seawater
  DOUBLE PRECISION, parameter :: mu =0.010260D0 ! (g cm^-1 s^-1) dynamic viscosity of seawater(nu*rhof)
  DOUBLE PRECISION, parameter :: rhop =1.06D0 ! (g cm^-3) Density of larva
  
  
!--these are the grid dimension/paramters read from the history/averages files. 

  INTEGER :: xi_rho(MAXNgrid)                 ! 
  INTEGER :: eta_rho(MAXNgrid)                 ! 
  INTEGER :: xi_u(MAXNgrid)                 !       
  INTEGER :: eta_u(MAXNgrid)                 !      
  INTEGER :: xi_v(MAXNgrid)                 !       
  INTEGER :: eta_v(MAXNgrid)                 !         
  INTEGER :: s_rho(MAXNgrid)                 !         
  INTEGER :: s_w(MAXNgrid)                 !      
  INTEGER :: Vtransform(MAXNgrid)            !       
  INTEGER :: Vstretching(MAXNgrid)            !       
  INTEGER :: theta_s(MAXNgrid)                 !       
  INTEGER :: theta_b(MAXNgrid)                 !       
  INTEGER :: tline(MAXNgrid)                 !      
  INTEGER :: zob(MAXNgrid)                 !    
  INTEGER :: tdim(MAXNgrid)		           !        
  INTEGER :: t_b
  INTEGER :: t_c 
  INTEGER :: t_f
  INTEGER :: tstep
  DOUBLE PRECISION :: hc(MAXNgrid)	
  
!  INTEGER :: rho_nodes(MAXNgrid)          ! number rho nodes
!  INTEGER :: u_nodes(MAXNgrid)            ! number of u nodes
!  INTEGER :: v_nodes(MAXNgrid)            ! number of v nodes!

!  INTEGER :: max_rho_elements(MAXNgrid)   ! total number of rho elements 
!  INTEGER :: max_u_elements(MAXNgrid)     ! total number of u elements
!  INTEGER :: max_v_elements(MAXNgrid)     ! total           v

!  INTEGER :: rho_elements(MAXNgrid)       ! number of rho elements with at least 1 vertex is water
!  INTEGER :: u_elements(MAXNgrid)         !            u 
!  INTEGER :: v_elements(MAXNgrid)         !            v

!group the grid info section in a namelist:

  namelist/gridinfo/ xi_rho,eta_rho,xi_u,eta_u,xi_v,eta_v,t_b,t_c,t_f,tstep,   &
     &				s_rho,s_w,Vtransform,Vstretching,theta_s,theta_b,tline,zob        



!The following used to be in ROMSPath.inc:


!*** NUMBER OF PARTICLES ***

  INTEGER :: numpar

  namelist/numparticles/numpar      ! Number of particles



!*** TIME PARAMETERS ***

  REAL             :: days          ! Number of days to run the model
  INTEGER          :: iprint        ! Print interval for ROMSPath output (s); 3600 = every hour
  INTEGER          :: dt            ! External time step (duration between hydro model predictions) (s) 
  INTEGER          :: idt           ! Internal (particle tracking) time step (s)
!  DOUBLE PRECISION :: dstart        ! Start time in data relative to ROMS model initilization 
  
!  namelist/timeparam/ days,iprint,dt,idt,dstart
  namelist/timeparam/ days,iprint,dt,idt


!!*** ROMS HYDRODYNAMIC MODULE PARAMETERS ***
  LOGICAL          :: readZeta      ! If .TRUE. read in sea-surface height   (zeta) from NetCDF file, else use constZeta
  DOUBLE PRECISION :: constZeta     ! Constant value for Zeta if readZeta is .FALSE.
  LOGICAL          :: readSalt      ! If .TRUE. read in salinity             (salt) from NetCDF file, else use constSalt
  DOUBLE PRECISION :: constSalt     ! Constant value for Salt if readSalt is .FALSE.
  LOGICAL          :: readTemp      ! If .TRUE. read in temperature          (temp) from NetCDF file, else use constTemp
  DOUBLE PRECISION :: constTemp     ! Constant value for Temp if readTemp is .FALSE.
  LOGICAL          :: readU         ! If .TRUE. read in u-momentum component (U   ) from NetCDF file, else use constU
  DOUBLE PRECISION :: constU        ! Constant value for U if readU is .FALSE.
  LOGICAL          :: readV         ! If .TRUE. read in v-momentum component (V   ) from NetCDF file, else use constV
  DOUBLE PRECISION :: constV        ! Constant value for V if readV is .FALSE.
  LOGICAL          :: readW         ! If .TRUE. read in w-momentum component (W   ) from NetCDF file, else use constW
  DOUBLE PRECISION :: constW        ! Constant value for W if readW is .FALSE.
  LOGICAL          :: readAks       ! If .TRUE. read in salinity vertical diffusion coefficient (Aks ) from NetCDF file, else use constAks
  DOUBLE PRECISION :: constAks      ! Constant value for Aks if readAks is .FALSE.
  LOGICAL          :: readDens
  DOUBLE PRECISION :: constDens
  CHARACTER(LEN=200),Dimension(MAXNgrid) :: stokesprefix
  CHARACTER(LEN=200),Dimension(MAXNgrid) :: turbstd_v_a_prefix
  CHARACTER(LEN=200),Dimension(MAXNgrid) :: wavestd_prefix
  
  LOGICAL :: Process_VA           ! PROCESS Vort./accel and write to netcdf file.
  LOGICAL :: Process_WA           ! PROCESS Wave Accel and write to netcdf file.
!
  namelist/hydroparam/readZeta,constZeta,readSalt,   &
                    & constSalt,readTemp,constTemp,readU,readU,constU,readV,     &
                    & constV,readW,constW,readAks,constAks,readDens,constDens,   &
					& stokesprefix,turbstd_v_a_prefix,wavestd_prefix,Process_VA,Process_WA


!*** TURBULENCE MODULE PARAMETERS ***

  LOGICAL          :: VTurbOn       ! Vertical   Turbulence on (.TRUE.) or off (.FALSE.)
  DOUBLE PRECISION :: serr 			!Cubic spline error Cutoff
  DOUBLE PRECISION :: smth			!Cubic spline smoothing
  DOUBLE PRECISION :: sub		    !Resolution multiplier for cubic spline smoothing
  DOUBLE PRECISION :: deltat		!Vertical turbulence time step
  DOUBLE PRECISION :: AKSback		!Vertical turbulence AKS background
  LOGICAL          :: HTurbOn       ! Horizontal Turbulence on (.TRUE.) or off (.FALSE.)
  DOUBLE PRECISION :: ConstantHTurb(MAXNgrid) ! Constant value of horizontal turbulence (m2/s)

  namelist/turbparam/VTurbOn,serr,smth,sub,deltat,AKSback,HTurbOn,ConstantHTurb

!*** Advection MODULE PARAMETERS ***
	INTEGER :: scheme    !Advection Scheme
	INTEGER :: nsb         ! Neutral, Surface or bottom 
	DOUBLE PRECISION :: vertdist  ! Vertical distance from surface or bottom for behavior types 1 or 2.   
	
  namelist/advectparam/scheme,nsb,vertdist
!*** BEHAVIOR MODULE PARAMETERS ***
  INTEGER :: Behavior               ! Behavior type (specify a number)
                                    !   Note: The behavior types numbers are: 
                                    !     0 Passive, 1 near-surface, 2 near-bottom, 
                                    ! 
  LOGICAL :: OpenOceanBoundary      ! Note: If you want to allow particles to "escape" via open ocean 
                                    !   boundaries, set this to TRUE; Escape means that the particle 
                                    !   will stick to the boundary and stop moving
  DOUBLE PRECISION :: pediage       ! Age when particle reaches max swim speed and can settle (s)
                                    !   Note: for oyster larvae behavior (types 4 & 5):
                                    !     pediage = age at which a particle becomes a pediveliger
                                    !   Note: pediage does not cause particles to settle if the Settlement module is not on
  DOUBLE PRECISION :: swimstart     ! Age that swimming or sinking begins (s) 1 day = 1.*24.*3600.
  DOUBLE PRECISION :: swimslow      ! Swimming speed when particle begins to swim (m/s)
  DOUBLE PRECISION :: swimfast      ! Maximum swimming speed (m/s)  0.05 m/s for 5 mm/s
                                    !   Note: for constant swimming speed for behavior types 1,2 & 3: 
                                    !     set swimslow = swimfast = constant speed
  DOUBLE PRECISION :: Sgradient     ! Salinity gradient threshold that cues larval behavior (psu/m)
                                    !   Note: This parameter is only used if Behavior = 4 or 5. 
  DOUBLE PRECISION :: sink          ! Sinking velocity for behavior type 6
                                    !   Note: This parameter is only used if Behavior = 6.
! Tidal Stream Transport behavior type:
  DOUBLE PRECISION :: Hswimspeed    ! Horizontal swimming speed (m/s)
  DOUBLE PRECISION :: Swimdepth     ! Depth at which fish swims during flood time 
                                    ! in meters above bottom (this should be a positive value)


  namelist/behavparam/Behavior,OpenOceanBoundary,pediage,swimstart,swimslow,swimfast,Sgradient,sink,Hswimspeed,Swimdepth
  
  
  
!*** BEHAVIOR MODULE PARAMETERS , FUCHS PARAMATERIZATION***
  
!%%%%%%%%%%%%%%%%% VORTICITY RESPONSES
  DOUBLE PRECISION :: vort_cr   !%(s^-1) critical vorticity for inducing response
  DOUBLE PRECISION :: vort_sat    !%(s^-1) vorticity where response saturates

  DOUBLE PRECISION :: b0pv    !% min probability of swimming vs. vorticity
  DOUBLE PRECISION :: b1pv    !% max probability of swimming vs. vorticity
  DOUBLE PRECISION :: b0wv    !% (cm s^-1) max swimming velocity vs vorticity -- keep it for flexibility
  DOUBLE PRECISION :: b1w     !% (cm s^-1) neutral buoyancy (no response)

!%%%%%%%%%%%%%%%%% ACCELERATION RESPONSES
  DOUBLE PRECISION :: acc_cr     !%(cm s^-2) critical acceleration for inducing response
  DOUBLE PRECISION :: acc_sat	!%(cm s^-2) acceleration where response saturates

  DOUBLE PRECISION :: b0pa   !% min probability of swimming vs. acceleration
  DOUBLE PRECISION :: b1pa   !% max probability of swimming vs. acceleration
  DOUBLE PRECISION :: b0wa   !% (cm s^-1) max swimming velocity vs acceleration


  INTEGER          :: va_flag ! 0=Both, 1=Vorticity Only, 1=Acceleration Only
   
  namelist/fuchsparam/vort_cr,vort_sat,b0pv,b1pv,b0wv,b1w,acc_cr,acc_sat,b0pa,b1pa,b0wa,va_flag
  
!*** Growth MODULE PARAMETERS ***

  INTEGER :: Growth               ! Griwth/Age type (specify a number)
                                    !   Note: The aging types numbers are: 
                                    !     0 none, 1 Use deadage, 2 Use Growth, 
  
  LOGICAL :: mortality              ! TRUE if particles can die; else FALSE
  DOUBLE PRECISION :: deadage       ! Age at which a particle stops moving (i.e., dies) (s)
  DOUBLE PRECISION :: initsize       ! Initial size of Larva(Egg size?)
  DOUBLE PRECISION :: maxsize       ! Maximum size of larva. (Stop moving after this) 
  DOUBLE PRECISION :: tempcut       ! Temperature cutoff for growth
  DOUBLE PRECISION :: a0       !  Growth Coefficient 0
  DOUBLE PRECISION :: a1       !  Growth Coefficient 1
  DOUBLE PRECISION :: a2       !  Growth Coefficient 2
  DOUBLE PRECISION :: a3       !  Growth Coefficient 3
  DOUBLE PRECISION :: a4       !  Growth Coefficient 4
  DOUBLE PRECISION :: a5       !  Growth Coefficient 5
  DOUBLE PRECISION :: a6       !  Growth Coefficient 6
  DOUBLE PRECISION :: a7       !  Growth Coefficient 7
  DOUBLE PRECISION :: a8       !  Growth Coefficient 8
  
  
  
  namelist/growparam/Growth,mortality,deadage,initsize,maxsize,tempcut,a0,a1,a2,a3,a4,a5,a6,a7,a8
  
  
!*** DVM. The following are parameters for the Diurnal Vertical Migration (DVM) behavior type:
  DOUBLE PRECISION :: twistart      ! Time of twilight start (hr) **
  DOUBLE PRECISION :: twiend        ! Time of twilight end (hr) **
  DOUBLE PRECISION :: daylength     ! Length of day (hr) **
  DOUBLE PRECISION :: Em            ! Irradiance at solar noon (microE m^-2 s^-1) **
  DOUBLE PRECISION :: Kd            ! Vertical attenuation coefficient
  DOUBLE PRECISION :: thresh        ! Light threshold that cues behavior (microE m^-2 s^-1)
  !  Note: These values were calculated for September 1 at the latitude of 37.0 (Chesapeake Bay mouth)
  !  Note: Variables marked with ** were calculated with light_v2BlueCrab.f (not included in ROMSPath yet)
  !  Note: These parameters are only used if Behavior = 3 

  namelist/dvmparam/twistart,twiend,daylength,Em,Kd,thresh



!*** SETTLEMENT MODULE PARAMETERS ***
  LOGICAL :: settlementon           ! settlement module on (.TRUE.) or off (.FALSE.)
  !  Note: If settlement is off: set minholeid, maxholeid, minpolyid, maxpolyid, pedges, & hedges to 1
  !        to avoid both wasted variable space and errors due to arrays of size 0.
  !        If settlement is on and there are no holes: set minholeid, maxholeid, & hedges to 1
  LOGICAL :: holesExist             ! Are there holes in habitat? yes(TRUE) no(FALSE)
  INTEGER :: minpolyid              ! Lowest habitat polygon id number
  INTEGER :: maxpolyid              ! Highest habitat polygon id number
  INTEGER :: minholeid              ! Lowest hole id number
  INTEGER :: maxholeid              ! Highest hole id number
  INTEGER :: pedges                 ! Number of habitat polygon edge points (# of rows in habitat polygon file)
  INTEGER :: hedges                 ! Number of hole edge points (number of rows in holes file)

  namelist/settleparam/settlementon,holesExist,minpolyid,maxpolyid,minholeid,maxholeid,pedges,hedges





!*** INPUT FILE NAME AND LOCATION PARAMETERS ***; 

!  ** ROMS NetCDF Model Grid file **
  
  INTEGER :: Ngrid                 ! Number of grids
  DOUBLE PRECISION :: refine(MAXNgrid)                ! Number of grids

  namelist/romsgrid/Ngrid,refine



!  ** ROMS Predictions NetCDF Input File **
!  Filename = prefix + filenum + suffix
  CHARACTER(LEN=200),Dimension(MAXNgrid) :: prefix      ! NetCDF Input Filename prefix
  CHARACTER(LEN=200) :: suffix      ! NetCDF Input Filename suffix
  CHARACTER(LEN=100) :: time_vname      ! NetCDF Input Time vairable name (usually ocean_time)
  CHARACTER(LEN=100) :: time_dname      ! NetCDF Input Time dimension name (usually ocean_time)
  INTEGER :: filenum                ! Number in First NetCDF Input Filename
  INTEGER :: numdigits              ! Number of digits in number portion of file name (with leading zeros)
!  LOGICAL :: startfile              ! .TRUE. means the first file has an additional time step
  LOGICAL :: multifile              ! .TRUE. means multiple files are used. .False. means only a single file/url
  !Note: the path to the file is necessary if the file is not in the same folder as the code
  !Note: if .nc file in separate folder in Windows, then include path in prefix. For example:
  !      CHARACTER(LEN=15), PARAMETER :: prefix='D:\ROMS\y95hdr_'   
  !      if .nc file in separate folder in Linux, then include path in prefix. For example:
  !      CHARACTER(LEN=26), PARAMETER :: prefix='/share/lzhong/1995/y95hdr_'   

  namelist/romsoutput/prefix,suffix,filenum,numdigits,multifile,time_vname,time_dname




!  ** Particle Location Input File **
  CHARACTER(LEN=200) :: parfile     ! Particle locations file
  !Note: the path to the file is necessary if the file is not in the same folder as the code

  namelist/parloc/parfile



!  ** Habitat Polygon Location Input Files **
  CHARACTER(LEN=200) :: habitatfile ! Habitat polygon file
  CHARACTER(LEN=200) :: holefile    ! Holes in habitat file
  !Note: the path to the file is necessary if the file is not in the same folder as the code

  namelist/HabPolyLoc/habitatfile,holefile



!  ** Output Related Variables **
  CHARACTER(LEN=200) :: outpath     ! Location to write output files
  CHARACTER(LEN=100) :: NCOutFile   ! Name of the NetCDF output file if outputting to NetCDF
  LOGICAL :: outpathGiven           ! If TRUE files are written to the path given in outpath
  INTEGER :: NCtime                 ! Time interval between creation of new NetCDF output files

  !NetCDF Model Metadata:

  CHARACTER(LEN=100) :: RunName     ! Unique Identifier for this particular model run
  CHARACTER(LEN=200) :: ExeDir      ! Location of the model run executable
  CHARACTER(LEN=200) :: OutDir      ! Location of the model run output files
  CHARACTER(LEN=100) :: RunBy       ! Name of person who setup/run the model
  CHARACTER(LEN=100) :: Institution ! Place the model is run
  CHARACTER(LEN=200) :: StartedOn   ! Date the model run began

  namelist/output/outpath,NCOutFile,outpathGiven,NCtime, &
                  RunName,ExeDir,OutDir,RunBy,Institution,StartedOn



!*** OTHER PARAMETERS *** 
  INTEGER :: seed                   ! Seed value for random number generator (Mersenne Twister)
  INTEGER :: ErrorFlag              ! What to do if an error is encountered: 0=stop, 1=return particle to previous location
                                    ! 2=kill particle & stop tracking, 3=set particle out of bounds & stop tracking
									! Note: Options 1-3 will output information to ErrorLog.txt
  LOGICAL :: SaltTempOn             ! Calculate salinity and temperature at particle 
                                    ! location: yes (.TRUE.) or no (.FALSE.)
  LOGICAL :: SaltTempMean           ! Average Salinity and temperature
  DOUBLE PRECISION :: TempOffset    ! Temperature offset applied to input
  
  LOGICAL :: WriteBottom            ! Write out bottom stress and height above bottom
  LOGICAL :: WriteWaterDepth        ! Write Total water depth
  LOGICAL :: WriteZeta              ! Write zeta
  LOGICAL :: WriteBath              ! Write ROMS bathymetry(H)
  LOGICAL :: TrackCollisions        ! Write Bottom and Land Hit Files? .TRUE.=yes, .FALSE.=no
  LOGICAL :: WriteHeaders           ! Write .txt files with column headers? .TRUE.=yes, .FALSE.=no
  LOGICAL :: WriteModelTiming       ! Write .csv file with model timing data? .TRUE.=yes, .FALSE.=no
  LOGICAL :: WriteProblemFile       ! Write a file with problem particles for debugging? .TRUE.=yes, .FALSE.=no

  INTEGER :: ijbuff                 ! number of extra elements to read in on every side of the particles

  LOGICAL :: FreeSlip               ! use free slip condition?

  namelist/other/seed,SaltTempOn,TrackCollisions,WriteHeaders,SaltTempMean,WriteBottom, &
                 WriteModelTiming,WriteProblemFile,ijbuff,ErrorFlag,FreeSlip,TempOffset,&
				 WriteWaterDepth,WriteZeta,WriteBath
