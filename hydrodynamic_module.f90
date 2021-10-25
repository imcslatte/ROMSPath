	MODULE HYDRO_MOD

	!  This module handles all the input from the hydrodynamic NetCDF input files.
	!  It is the only module that interacts with NetCDF input files.  It contains
	!  all the variables read in from the NetCDF files.  It also contains all the
	!  information and variables related to the grid elements.
	!
	!  Created by:            Zachary Schlag        
	!  Created on:            07 Aug 2008
	!  Last Modified on:         Feb 2013

	  IMPLICIT NONE
	  PRIVATE
	  PUBLIC ::getFileNames

	  SAVE

	 

	  !Used for reading in NetCDF variables one time step at a time
	  INTEGER :: STARTr(4),COUNTr(4),STARTz(3),COUNTz(3)

	  !These variables keep track of the interpolation method and weights
	  ! INTEGER :: tOK
	   DOUBLE PRECISION :: t,u,Wgt1,Wgt2,Wgt3,Wgt4
	  ! !The Rho, U, and V nodes that make up the Rho, U, and V element that 
	  ! !  the particle is in
	  ! INTEGER, ALLOCATABLE,DIMENSION(:) :: rnode1,rnode2,rnode3,rnode4,unode1,unode2,unode3,unode4,vnode1,   &
				 ! vnode2,vnode3,vnode4
		
	  !read in zeta,salinity,temperature,vertical diffusivity, and U,V,W velocities 
	  !  at hydrodynamic back, center, and forward time
	 
	  INTEGER, ALLOCATABLE,DIMENSION(:) :: stepf   !Keeps track of the forward time step

	 
	  
	  ! !S-Level location variables
	  ! DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SC,CS,SCW,CSW
	  ! !Depth at each rho node location
	  ! DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: depth 
	  ! !Rho, U, and V grid wet elements(four node numbers that make up the element)
	  ! !  (wet means at least 1 node is masked as water)
	  ! INTEGER, ALLOCATABLE, DIMENSION(:,:) :: RE,UE,VE
	  ! !Keeps track of the Rho, U, and V element that each particle is in
	  ! INTEGER, ALLOCATABLE, DIMENSION(:) :: P_r_element,P_u_element,P_v_element
	  ! !For each element, a list containing itself and all the elements that share a 
	  ! !  node with that element;  used to speed up determining which element the 
	  ! !  particle has moved to, if it has moved at all
	  

	  TYPE HDATA
		DOUBLE PRECISION,pointer :: zeta(:,:,:)
		DOUBLE PRECISION,pointer :: bustr(:,:,:)
		DOUBLE PRECISION,pointer :: bvstr(:,:,:)
		DOUBLE PRECISION,pointer :: salt(:,:,:,:)
		DOUBLE PRECISION,pointer :: temp(:,:,:,:)
		DOUBLE PRECISION,pointer :: AKs(:,:,:,:)
		DOUBLE PRECISION,pointer :: U(:,:,:,:)
		DOUBLE PRECISION,pointer :: V(:,:,:,:)
		DOUBLE PRECISION,pointer :: W(:,:,:,:)
		DOUBLE PRECISION,pointer :: time(:)
		DOUBLE PRECISION,pointer :: Accelstd_t(:,:,:,:)
		DOUBLE PRECISION,pointer :: Vortstd_t(:,:,:,:)
		DOUBLE PRECISION,pointer :: Accelustd_w(:,:,:,:)
		DOUBLE PRECISION,pointer :: Accelvstd_w(:,:,:,:)
		DOUBLE PRECISION,pointer :: Accelwstd_w(:,:,:,:)
#ifdef STOKES
		DOUBLE PRECISION,pointer :: SU(:,:,:,:)
		DOUBLE PRECISION,pointer :: SV(:,:,:,:)
#endif

#ifdef WETDRY
		DOUBLE PRECISION,pointer :: wetdry_mask_u(:,:,:)
		DOUBLE PRECISION,pointer :: wetdry_mask_v(:,:,:)
		DOUBLE PRECISION,pointer :: wetdry_mask_rho(:,:,:)
#endif
	  END TYPE HDATA		
	  TYPE (HDATA), allocatable :: HYDRODATA(:) 	


	  !Keeps track if the grid has been read in yet or not
	  !  If the grid hasn't been read in, the boundaries can't be made
	  ! LOGICAL :: GRD_SET = .FALSE.

	  !The concatenated hydrodynamic input file name
	  CHARACTER(len=200) :: filenm,turbfilenm,wavefilenm
#ifdef STOKES
	  CHARACTER(len=200) :: stokesfilenm
#endif
	  character(len=256) :: Iname
	  !Counters for NetCDF files
	  
	  
	  !The following procedures have been made public:
	 ! PUBLIC :: initGrid,initHydro,updateHydro,setEle,setEle_all,setInterp,        &
	 !   getInterp,interp,WCTS_ITPI,getSlevel,getWlevel,getMask_Rho,getUVxy,        &
	 !   getR_ele,getP_r_element,finHydro,initNetCDF,createNetCDF,writeNetCDF
		
	  PUBLIC :: updateHydro,HYDRODATA        
	! &
		! getMask_Rho,getUVxy,finHydro,setEle,       &
		! initNetCDF,createNetCDF,writeNetCDF,MODGRID, &
		! getInterp,interp,setInterp,WCTS_ITPI,getSlevel,getWlevel

	CONTAINS

	 


	  SUBROUTINE updateHydro(FIRST,tstep,tind)
		!This Subroutine reads in the hydrodynamic information for the first 
		!  iteration
		USE PARAM_MOD, ONLY: numpar,xi_rho,eta_rho,s_rho,s_w,suffix,&
			prefix,filenum,numdigits,readZeta,constZeta,readSalt,constSalt, &
			readTemp,constTemp,readDens,constDens,readU,constU,readV,constV,readW, &
			constW,readAks,constAks,Ngrid,xi_u,eta_u,xi_v,eta_v,tdim,t_b,t_c,t_f,&
			stokesprefix,TempOffset,WriteBottom,turbstd_v_a_prefix,wavestd_prefix,&
			Behavior,Process_VA,Process_WA,time_vname,time_dname 	
		USE netcdf
		IMPLICIT NONE

		INCLUDE 'netcdf.inc'

		INTEGER :: STATUS,NCID,VID,DIMID,dimcount
		LOGICAL,INTENT(IN) :: FIRST

		INTEGER ,INTENT(IN) :: tstep,tind
		INTEGER :: i,j,k,ng
		real :: before,after,tdiff
		

			
		do ng=1,Ngrid
			if (FIRST) then 
				if (ng.eq.1) allocate(HYDRODATA(Ngrid))
					!ALLOCATE MODULE VARIABLES
					ALLOCATE(HYDRODATA(ng)%zeta(xi_rho(ng),eta_rho(ng),3))
					ALLOCATE(HYDRODATA(ng)%bustr(xi_u(ng),eta_u(ng),3))
					ALLOCATE(HYDRODATA(ng)%bvstr(xi_v(ng),eta_v(ng),3))
					ALLOCATE(HYDRODATA(ng)%salt(xi_rho(ng),eta_rho(ng),s_rho(ng),3))
					ALLOCATE(HYDRODATA(ng)%temp(xi_rho(ng),eta_rho(ng),s_rho(ng),3))
					ALLOCATE(HYDRODATA(ng)%W(xi_rho(ng),eta_rho(ng),s_w(ng),3))
					ALLOCATE(HYDRODATA(ng)%AKs(xi_rho(ng),eta_rho(ng),s_w(ng),3))
					ALLOCATE(HYDRODATA(ng)%U(xi_u(ng),eta_u(ng),s_rho(ng),3))
					ALLOCATE(HYDRODATA(ng)%V(xi_v(ng),eta_v(ng),s_rho(ng),3))
					ALLOCATE(HYDRODATA(ng)%AKs(xi_rho(ng),eta_rho(ng),s_w(ng),3))
					ALLOCATE(HYDRODATA(ng)%Accelstd_t(xi_rho(ng),eta_rho(ng),s_w(ng),3))
					ALLOCATE(HYDRODATA(ng)%Vortstd_t(xi_rho(ng),eta_rho(ng),s_w(ng),3))
					ALLOCATE(HYDRODATA(ng)%Accelustd_w(xi_rho(ng),eta_rho(ng),s_w(ng),3))
					ALLOCATE(HYDRODATA(ng)%Accelvstd_w(xi_rho(ng),eta_rho(ng),s_w(ng),3))
					ALLOCATE(HYDRODATA(ng)%Accelwstd_w(xi_rho(ng),eta_rho(ng),s_w(ng),3))
#ifdef STOKES
					ALLOCATE(HYDRODATA(ng)%SU(xi_u(ng),eta_u(ng),s_rho(ng),3))
					ALLOCATE(HYDRODATA(ng)%SV(xi_v(ng),eta_v(ng),s_rho(ng),3))
#endif
#ifdef WETDRY
					ALLOCATE(HYDRODATA(ng)%wetdry_mask_u(xi_u(ng),eta_u(ng),3))
					ALLOCATE(HYDRODATA(ng)%wetdry_mask_v(xi_v(ng),eta_v(ng),3))
					ALLOCATE(HYDRODATA(ng)%wetdry_mask_rho(xi_rho(ng),eta_rho(ng),3))
					HYDRODATA(ng)%wetdry_mask_u 	= 0
					HYDRODATA(ng)%wetdry_mask_v 	= 0
					HYDRODATA(ng)%wetdry_mask_rho 	= 0
#endif
					HYDRODATA(ng)%zeta 	= 0
					HYDRODATA(ng)%salt 	= 0
					HYDRODATA(ng)%temp 	= 0
					HYDRODATA(ng)%AKs 	= 0
					HYDRODATA(ng)%U 	= 0
					HYDRODATA(ng)%V 	= 0
					HYDRODATA(ng)%W 	= 0
					HYDRODATA(ng)%Accelstd_t 	= 0
					HYDRODATA(ng)%Vortstd_t 	= 0
					HYDRODATA(ng)%Accelustd_w 	= 0
					HYDRODATA(ng)%Accelvstd_w 	= 0
					HYDRODATA(ng)%Accelwstd_w 	= 0
			
					
				
			endif
			
		 
		  ! !Open netCDF file
		   call getFileNames(filenm,prefix(ng),filenum)
		   call getFileNames(turbfilenm,turbstd_v_a_prefix(ng),filenum)
		   call getFileNames(wavefilenm,wavestd_prefix(ng),filenum)
#ifdef STOKES
		   call getFileNames(stokesfilenm,stokesprefix(ng),filenum)
#endif


			if (tstep.eq.1)then
				write(*,*) "New ROMS File:"
				write(*,*) filenm
		    endif
			
	
		
		  ! Read in data for first three external time steps
		 
		  STATUS = NF90_OPEN(TRIM(filenm), NF90_NOWRITE, NCID)
		  if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem NF90_OPEN HYDROFILE'
		  if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)




			STATUS = NF90_INQ_DIMID(NCID,trim(time_dname),DIMID)
			STATUS = NF90_INQUIRE_DIMENSION(NCID,DIMID,len=dimcount)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem with dimid:'
			  write(*,*) time_dname
			  write(*,*) NF90_STRERROR(STATUS)
			endif
			tdim(ng) = dimcount
			
			startz(1)=tstep
			countz(1)=1
		  STATUS = NF90_INQ_VARID(NCID,trim(time_vname),VID)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem finding time variable:'
			  write(*,*) trim(time_vname)
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif
		!  STATUS = NF90_GET_VAR(NCID,VID,ttime,startz,countz)
		  

		
		  
		  if(readZeta)then  
			! **** Zeta ****
			startz(1)=1
			startz(2)=1
			startz(3)=tstep

			countz(1)=xi_rho(ng)
			countz(2)=eta_rho(ng)
			countz(3)=1

			STATUS = NF90_INQ_VARID(NCID,'zeta',VID)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem find zeta'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif

			STATUS = NF90_GET_VAR(NCID,VID,HYDRODATA(ng)%zeta(:,:,tind),STARTz,COUNTz)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem read zeta array 1'	
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif
		  else
			HYDRODATA(ng)%zeta= constZeta
		  endif
		
		  if(readSalt)then
			! **** Salt ****
			startr(1)=1
			startr(2)=1
			startr(3)=1
			startr(4)=tstep

			countr(1)=xi_rho(ng)
			countr(2)=eta_rho(ng)
			countr(3)=s_rho(ng)
			countr(4)=1
			STATUS = NF90_INQ_VARID(NCID,'salt',VID)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem find salt'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif

			STATUS = NF90_GET_VAR(NCID,VID,HYDRODATA(ng)%salt(:,:,:,tind),STARTr,COUNTr)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem read salt array'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif
		  else
			HYDRODATA(ng)%salt = constSalt
		  endif
		  
		  if(readTemp)then  
			! **** Temp ****
			startr(1)=1
			startr(2)=1
			startr(3)=1
			startr(4)=tstep

			countr(1)=xi_rho(ng)
			countr(2)=eta_rho(ng)
			countr(3)=s_rho(ng)
			countr(4)=1
			STATUS = NF90_INQ_VARID(NCID,'temp',VID)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem find temp'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif

			STATUS = NF90_GET_VAR(NCID,VID,HYDRODATA(ng)%temp(:,:,:,tind),STARTr,COUNTr)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem read temp array'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif

			do i=1,xi_rho(ng)
				do j=1,eta_rho(ng)
					do k=1,s_rho(ng)
							HYDRODATA(ng)%temp(i,j,k,tind)=HYDRODATA(ng)%temp(i,j,k,tind)+TempOffset
					enddo
				enddo
			enddo	
		  else
			HYDRODATA(ng)%temp = constTemp
		  endif
		
		
		 
		  
		  ! call CPU_TIME(before)
		  if(readU)then  
			! **** U velocity ****
			startr(1)=1
			startr(2)=1
			startr(3)=1
			startr(4)=tstep

			countr(1)=xi_u(ng)
			countr(2)=eta_u(ng)
			countr(3)=s_rho(ng)
			countr(4)=1
			
			

			
			STATUS = NF90_INQ_VARID(NCID,'u',VID)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem find u'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif

			STATUS = NF90_GET_VAR(NCID,VID,HYDRODATA(ng)%U(:,:,:,tind),STARTr,COUNTr)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem read u array'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif

			
		  else
			HYDRODATA(ng)%U = constU
		  endif

		  
		  
		  
		  
		  if(readV)then  
			! **** V velocity ****
			startr(1)=1
			startr(2)=1
			startr(3)=1
			startr(4)=tstep
			
			countr(1)=xi_v(ng)
			countr(2)=eta_v(ng)
			countr(3)=s_rho(ng)
			countr(4)=1
			
			STATUS = NF90_INQ_VARID(NCID,'v',VID)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem find v'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif

			STATUS = NF90_GET_VAR(NCID,VID,HYDRODATA(ng)%V(:,:,:,tind),STARTr,COUNTr)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem read v array'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif
		
			
		  else
			HYDRODATA(ng)%V = constV
		  endif

		  ! call CPU_TIME(after)
		  ! tdiff=after-before
		  ! write(*,*) '****'
		  ! write(*,*) tdiff
		  
		 if(WriteBottom)then  
			! **** U stress ****
			startr(1)=1
			startr(2)=1
			startr(3)=tstep

			countr(1)=xi_u(ng)
			countr(2)=eta_u(ng)
			countr(3)=1
			STATUS = NF90_INQ_VARID(NCID,'bustr',VID)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem find bustr'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif

			STATUS = NF90_GET_VAR(NCID,VID,HYDRODATA(ng)%bustr(:,:,tind),STARTr,COUNTr)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem read bustr array'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif
			
			! **** v stress ****
			startr(1)=1
			startr(2)=1
			startr(3)=tstep

			countr(1)=xi_v(ng)
			countr(2)=eta_v(ng)
			countr(3)=1
			STATUS = NF90_INQ_VARID(NCID,'bvstr',VID)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem find bvstr'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif

			STATUS = NF90_GET_VAR(NCID,VID,HYDRODATA(ng)%bvstr(:,:,tind),STARTr,COUNTr)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem read bvstr array'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif


			
		  else
			HYDRODATA(ng)%bustr = 0.0
			HYDRODATA(ng)%bvstr = 0.0
		  endif

			
		  if(readW)then  
			! **** W velocity ****
			startr(1)=1
			startr(2)=1
			startr(3)=1
			startr(4)=tstep

			countr(1)=xi_rho(ng)
			countr(2)=eta_rho(ng)
			countr(3)=s_w(ng)
			countr(4)=1
			STATUS = NF90_INQ_VARID(NCID,'w',VID)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem find w'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif

			STATUS = NF90_GET_VAR(NCID,VID,HYDRODATA(ng)%W(:,:,:,tind),STARTr,COUNTr)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem read w array'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif
		  else
			HYDRODATA(ng)%W = constW
		  endif

		  if(readAks)then  
			! **** Vertical diffusivity for salt (Aks) ****
			startr(1)=1
			startr(2)=1
			startr(3)=1
			startr(4)=tstep

			countr(1)=xi_rho(ng)
			countr(2)=eta_rho(ng)
			countr(3)=s_w(ng)
			countr(4)=1
			STATUS = NF90_INQ_VARID(NCID,'AKs',VID)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem find AKs'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif

			STATUS = NF90_GET_VAR(NCID,VID,HYDRODATA(ng)%AKs(:,:,:,tind),STARTr,COUNTr)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem read AKs array'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif
		  else
			HYDRODATA(ng)%AKs = constAks
		  endif

#ifdef WETDRY
		
		   startz(1)=1
		   startz(2)=1
		   startz(3)=tstep
			
			countz(1)=xi_v(ng)
			countz(2)=eta_v(ng)
			countz(3)=1
			STATUS = NF90_INQ_VARID(NCID,'wetdry_mask_v',VID)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem find wetdry_mask_v'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif


			STATUS = NF90_GET_VAR(NCID,VID,HYDRODATA(ng)%wetdry_mask_v(:,:,tind),STARTz,COUNTz)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem read wetdry_mask_v array'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif
			
			countz(1)=xi_u(ng)
			countz(2)=eta_u(ng)
			countz(3)=1
			STATUS = NF90_INQ_VARID(NCID,'wetdry_mask_u',VID)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem find wetdry_mask_u'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif

			STATUS = NF90_GET_VAR(NCID,VID,HYDRODATA(ng)%wetdry_mask_u(:,:,tind),STARTz,COUNTz)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem read wetdry_mask_u array'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif
			
			
			countz(1)=xi_rho(ng)
			countz(2)=eta_rho(ng)
			countz(3)=1
			STATUS = NF90_INQ_VARID(NCID,'wetdry_mask_rho',VID)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem find wetdry_mask_rho'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif

			STATUS = NF90_GET_VAR(NCID,VID,HYDRODATA(ng)%wetdry_mask_rho(:,:,tind),STARTz,COUNTz)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem read wetdry_mask_rho array'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif
			
			
#endif
		  ! !close the dataset and reassign the NCID
		  STATUS = NF90_CLOSE(NCID)
#ifdef STOKES
          
		! call CPU_TIME(before)
		  STATUS = NF90_OPEN(TRIM(stokesfilenm), NF90_NOWRITE, NCID)
		  if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem NF90_OPEN STOKESFILE'
		  if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
		 
		 

			! **** STOKES U velocity ****
			startr(1)=1
			startr(2)=1
			startr(3)=1
			startr(4)=tstep

			countr(1)=xi_u(ng)
			countr(2)=eta_u(ng)
			countr(3)=s_rho(ng)
			countr(4)=1
			
			
			STATUS = NF90_INQ_VARID(NCID,'ustokes',VID)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem find u'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif

			STATUS = NF90_GET_VAR(NCID,VID,HYDRODATA(ng)%SU(:,:,:,tind),STARTr,COUNTr)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem read stokes u array'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif

		  	do i=1,xi_u(ng)
				do j=1,eta_u(ng)
					do k=1,s_rho(ng)
							HYDRODATA(ng)%U(i,j,k,tind)=HYDRODATA(ng)%U(i,j,k,tind)+HYDRODATA(ng)%SU(i,j,k,tind)
					enddo
				enddo
			enddo
			
			
			! **** STOKES V velocity ****
			startr(1)=1
			startr(2)=1
			startr(3)=1
			startr(4)=tstep
			
			countr(1)=xi_v(ng)
			countr(2)=eta_v(ng)
			countr(3)=s_rho(ng)
			countr(4)=1
			
			STATUS = NF90_INQ_VARID(NCID,'vstokes',VID)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem find v'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif

			STATUS = NF90_GET_VAR(NCID,VID,HYDRODATA(ng)%SV(:,:,:,tind),STARTr,COUNTr)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem read stokes v array'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif
			
			

		  
		  STATUS = NF90_CLOSE(NCID)	
		  	  ! call CPU_TIME(after)
		  ! tdiff=after-before
		  ! write(*,*) '^^^'
		  ! write(*,*) tdiff

			do i=1,xi_v(ng)
				do j=1,eta_v(ng)
					do k=1,s_rho(ng)
							HYDRODATA(ng)%V(i,j,k,tind)=HYDRODATA(ng)%V(i,j,k,tind)+HYDRODATA(ng)%SV(i,j,k,tind)
					enddo
				enddo
			enddo		  
	

#endif	

		if ((Behavior.EQ.10) .OR. (Process_VA)) then
! call CPU_TIME(before)
		 
			startr(1)=1
			startr(2)=1
			startr(3)=1
			startr(4)=tstep

			countr(1)=xi_rho(ng)
			countr(2)=eta_rho(ng)
			countr(3)=s_w(ng)
			countr(4)=1
			
		  STATUS = NF90_OPEN(TRIM(turbfilenm), NF90_NOWRITE, NCID)
		  if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem NF90_OPEN turbfilenm'
		  if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

			
			
			STATUS = NF90_INQ_VARID(NCID,'vortstd_cmpnt_turb',VID)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem find vortstd_cmpnt_turb'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif

			STATUS = NF90_GET_VAR(NCID,VID,HYDRODATA(ng)%Vortstd_t(:,:,:,tind),STARTr,COUNTr)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem read vortstd_cmpnt_turb array'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif

		
			
			STATUS = NF90_INQ_VARID(NCID,'accelstd_cmpnt_turb',VID)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem find accelstd_cmpnt_turb'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif

			STATUS = NF90_GET_VAR(NCID,VID,HYDRODATA(ng)%Accelstd_t(:,:,:,tind),STARTr,COUNTr)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem read accelstd_cmpnt_turb array'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif
			
		  STATUS = NF90_CLOSE(NCID)	
			  
		endif
		
		
		
		if (Process_WA) then
	
		  STATUS = NF90_OPEN(TRIM(wavefilenm), NF90_NOWRITE, NCID)
		  if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem NF90_OPEN wavefilenm'
		  if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
	
			STATUS = NF90_INQ_VARID(NCID,'accelstd_u_wave',VID)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem find accelstd_u_wave'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif
			
			STATUS = NF90_GET_VAR(NCID,VID,HYDRODATA(ng)%Accelustd_w(:,:,:,tind),STARTr,COUNTr)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem read accelstd_u_waves array'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif
			
	
			STATUS = NF90_INQ_VARID(NCID,'accelstd_v_wave',VID)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem find accelstd_v_wave'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif
			
			STATUS = NF90_GET_VAR(NCID,VID,HYDRODATA(ng)%Accelvstd_w(:,:,:,tind),STARTr,COUNTr)
			if (STATUS .NE. NF90_NOERR) then
			  write(*,*) 'Problem read accelstd_u_waves array'
			  write(*,*) NF90_STRERROR(STATUS)
			  stop
			endif
		  
		  STATUS = NF90_CLOSE(NCID)	
		endif
		
		
		
	enddo
		

	  END SUBROUTINE updateHydro


	 

	  SUBROUTINE finHydro()
		!This subroutine closes all the module's allocatable variables
		IMPLICIT NONE


		DEALLOCATE(HYDRODATA)

	  END SUBROUTINE finHydro

	   SUBROUTINE getFileNames(filenm,prefix,filenum)
		USE PARAM_MOD, ONLY: numdigits,suffix,multifile
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: filenum 
		CHARACTER(len=200), INTENT(INOUT) :: filenm,prefix
		
		if (multifile) then
			SELECT CASE(numdigits)
			CASE(1)
				WRITE(filenm,'(A,I1.1,A)') TRIM(prefix),filenum,TRIM(suffix)
			CASE(2)
				WRITE(filenm,'(A,I2.2,A)') TRIM(prefix),filenum,TRIM(suffix)
			CASE(3)
				WRITE(filenm,'(A,I3.3,A)') TRIM(prefix),filenum,TRIM(suffix)
			CASE(4)
				WRITE(filenm,'(A,I4.4,A)') TRIM(prefix),filenum,TRIM(suffix)
			CASE(5)
				WRITE(filenm,'(A,I5.5,A)') TRIM(prefix),filenum,TRIM(suffix)
			CASE(6)
				WRITE(filenm,'(A,I6.6,A)') TRIM(prefix),filenum,TRIM(suffix)
			CASE(7)
				WRITE(filenm,'(A,I7.7,A)') TRIM(prefix),filenum,TRIM(suffix)
			CASE(8)
				WRITE(filenm,'(A,I8.8,A)') TRIM(prefix),filenum,TRIM(suffix)
			CASE DEFAULT
				WRITE(*,*) 'Model presently does not support numdigits of ',numdigits
				WRITE(*,*) 'Please use numdigit value from 1 to 8'
				WRITE(*,*) '  OR modify code in Hydrodynamic module'
				STOP
			END SELECT
		else
			filenm=TRIM(prefix)
		endif
	
	  END SUBROUTINE

	END MODULE HYDRO_MOD
