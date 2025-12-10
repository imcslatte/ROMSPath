MODULE GRID_MOD 
 
! Grid MOdule

! ROMSPath Version: 1.0.1
 
IMPLICIT NONE 
PUBLIC 
SAVE 


DOUBLE PRECISION ::  reftime
CHARACTER(len=100) :: time_units
TYPE GRIDDATA
	DOUBLE PRECISION,pointer :: s_rho(:)
	DOUBLE PRECISION,pointer :: s_w(:)
	DOUBLE PRECISION,pointer :: cs_r(:)
	DOUBLE PRECISION,pointer :: cs_w(:)
	DOUBLE PRECISION,pointer :: H(:,:) 
	DOUBLE PRECISION,pointer :: lon_rho(:,:) 
	DOUBLE PRECISION,pointer :: lat_rho(:,:) 	
	DOUBLE PRECISION,pointer :: xi(:,:) 	
	DOUBLE PRECISION,pointer :: eta(:,:) 
	DOUBLE PRECISION,pointer :: pm(:,:)
	DOUBLE PRECISION,pointer :: pn(:,:) 	
	DOUBLE PRECISION,pointer :: angle(:,:) 	 
	DOUBLE PRECISION,pointer :: scl(:,:)
	DOUBLE PRECISION,pointer :: off(:,:)
	DOUBLE PRECISION,pointer :: z_rho(:,:,:)
	DOUBLE PRECISION,pointer :: z_w(:,:,:) 	  	 		
	INTEGER,pointer :: mask_rho(:,:)
	INTEGER,pointer :: mask_u(:,:)
	INTEGER,pointer :: mask_v(:,:)
END TYPE GRIDDATA		
TYPE (GRIDDATA), allocatable :: GRIDS(:) 	
PUBLIC ::GRIDS
CONTAINS


  SUBROUTINE InitGrid()

  
    USE netcdf
	USE PARAM_MOD, ONLY: xi_rho,eta_rho,xi_u,eta_u,xi_v,eta_v,    &
		s_rho,s_w,Vtransform,Vstretching,theta_s,theta_b,tline,zob,	&
        prefix,suffix,filenum,numdigits,Ngrid,refine,hc,time_vname 
	USE HYDRO_MOD, ONLY: getFileNames
    IMPLICIT NONE

!    INTEGER, INTENT(OUT), OPTIONAL :: IOSTAT

    INCLUDE 'netcdf.inc'

    !NetCDF Variables
    INTEGER :: NCID,STATUS,VID,dimid,dimcount,ng,ng2,nref,dng

	CHARACTER(len=200) :: filenm,header
	CHARACTER(len=100) :: strtmp
	DOUBLE PRECISION :: grefine(Ngrid,Ngrid)
    !Grid File Output Variables
 !   INTEGER :: nR,nU,nV,maxR,maxU,maxV,wetR,wetU,wetV

    !Iteration Variables
    INTEGER :: i,j,err,ii,jj

    err = 0
	header="PROBLEM READING GRID INFORMATION"
    ! *********************** GET GRID INFO ***********************
	
    ! OPEN NETCDF FILE - GET NCID VALUE
	do ng=1,Ngrid
		if (ng.eq.1) allocate(GRIDS(Ngrid))
		
		call getFileNames(filenm,prefix(ng),filenum)
	
	
		
		STATUS = NF90_OPEN(filenm,NF90_NOWRITE,NCID)
		if (STATUS .NE. NF90_NOERR) then
		write(*,*) 'Problem with NF90_OPEN:'
		write(*,*) 'File not found:'
		write(*,*) filenm
		err = 10
		call errorHandler(header,-1)
		endif

    ! GET VALUES FOR xi_rho,xi_u,xi_v,eta_rho,eta_u,eta_v

		STATUS = NF90_INQ_DIMID(NCID,'xi_rho',dimid)
		STATUS = NF90_INQUIRE_DIMENSION(NCID,dimid,len=dimcount)
		if (STATUS .NE. NF90_NOERR) then
			write(*,*) 'Problem dimid xi_rho'
			err = 20 
		call errorHandler(header,-1)
		endif
		xi_rho(ng) = dimcount

		STATUS = NF90_INQ_DIMID(NCID,'eta_rho',dimid)
		STATUS = NF90_INQUIRE_DIMENSION(NCID,dimid,len=dimcount)
		if (STATUS .NE. NF90_NOERR) then
			write(*,*) 'Problem dimid eta_rho'
			err = 20 
			call errorHandler(header,-1)
		endif
		eta_rho(ng) = dimcount
		!print *, eta_rho(ng)	!ELI
		!print *, dimcount	!ELI
		
		STATUS = NF90_INQ_DIMID(NCID,'xi_u',dimid)
		STATUS = NF90_INQUIRE_DIMENSION(NCID,dimid,len=dimcount)
		if (STATUS .NE. NF90_NOERR) then
			write(*,*) 'Problem dimid xi_u'
			err = 20 
			call errorHandler(header,-1)
		endif
		xi_u(ng) = dimcount

		!STATUS = NF90_INQ_DIMID(NCID,'eta_u',dimid)
		STATUS = NF90_INQ_DIMID(NCID,'eta_rho',dimid)	!ELI   
		STATUS = NF90_INQUIRE_DIMENSION(NCID,dimid,len=dimcount)
		if (STATUS .NE. NF90_NOERR) then
			write(*,*) 'Problem dimid eta_u'
			err = 20 
			call errorHandler(header,-1)
		endif
		eta_u(ng) = dimcount

		!STATUS = NF90_INQ_DIMID(NCID,'xi_v',dimid)
		STATUS = NF90_INQ_DIMID(NCID,'xi_rho',dimid)	!ELI   
		STATUS = NF90_INQUIRE_DIMENSION(NCID,dimid,len=dimcount)
		if (STATUS .NE. NF90_NOERR) then
			write(*,*) 'Problem dimid xi_v'
			call errorHandler(header,-1)
			err = 20 
		endif
		xi_v(ng) = dimcount

		STATUS = NF90_INQ_DIMID(NCID,'eta_v',dimid)
		STATUS = NF90_INQUIRE_DIMENSION(NCID,dimid,len=dimcount)
		if (STATUS .NE. NF90_NOERR) then
			write(*,*) 'Problem dimid eta_v'
			call errorHandler(header,-1)
			err = 20 
		endif
		eta_v(ng) = dimcount
		
		STATUS = NF90_INQ_DIMID(NCID,'s_rho',dimid)
		STATUS = NF90_INQUIRE_DIMENSION(NCID,dimid,len=dimcount)
		if (STATUS .NE. NF90_NOERR) then
			write(*,*) 'Problem dimid s_rho'
			call errorHandler(header,-1)
			err = 20 
		endif
		s_rho(ng) = dimcount
		
		STATUS = NF90_INQ_DIMID(NCID,'s_w',dimid)
		STATUS = NF90_INQUIRE_DIMENSION(NCID,dimid,len=dimcount)
		if (STATUS .NE. NF90_NOERR) then
			write(*,*) 'Problem dimid s_w'
			call errorHandler(header,-1)
			err = 20 
		endif
		s_w(ng) = dimcount
		

    ! READ IN Grid vertical trasnformation paramters
	
		STATUS = NF90_INQ_VARID(NCID,'Vtransform',VID)
		STATUS = NF90_GET_VAR(NCID,VID,Vtransform(ng))
		if (STATUS .NE. NF90_NOERR) then
			write(*,*) 'Problem read Vtransform'
			err = 40 
			call errorHandler(header,-1)
		endif
		
		STATUS = NF90_INQ_VARID(NCID,'Vstretching',VID)
		STATUS = NF90_GET_VAR(NCID,VID,Vstretching(ng))
		if (STATUS .NE. NF90_NOERR) then
			write(*,*) 'Problem read Vstretching'
			err = 40 
			call errorHandler(header,-1)
		endif
		
		STATUS = NF90_INQ_VARID(NCID,'theta_s',VID)
		STATUS = NF90_GET_VAR(NCID,VID,theta_s(ng))
		if (STATUS .NE. NF90_NOERR) then
			write(*,*) 'Problem read theta_s'
			err = 40 
			call errorHandler(header,-1)
		endif
		
		STATUS = NF90_INQ_VARID(NCID,'theta_b',VID)
		STATUS = NF90_GET_VAR(NCID,VID,theta_b(ng))
		if (STATUS .NE. NF90_NOERR) then
			write(*,*) 'Problem read theta_b'
			err = 40 
			call errorHandler(header,-1)
		endif
		
		STATUS = NF90_INQ_VARID(NCID,'hc',VID)
		STATUS = NF90_GET_VAR(NCID,VID,hc(ng))
		if (STATUS .NE. NF90_NOERR) then
			write(*,*) 'Problem read hc'
			err = 40 
			call errorHandler(header,-1)
		endif

    ! ! ALLOCATE VARIABLE ARRAY DIMENSIONS	
		ALLOCATE(GRIDS(ng)%mask_rho(xi_rho(ng),eta_rho(ng)))
		ALLOCATE(GRIDS(ng)%H(xi_rho(ng),eta_rho(ng)))
		ALLOCATE(GRIDS(ng)%angle(xi_rho(ng),eta_rho(ng)))
		ALLOCATE(GRIDS(ng)%lon_rho(xi_rho(ng),eta_rho(ng)))
		ALLOCATE(GRIDS(ng)%lat_rho(xi_rho(ng),eta_rho(ng)))
		ALLOCATE(GRIDS(ng)%pm(xi_rho(ng),eta_rho(ng)))
		ALLOCATE(GRIDS(ng)%pn(xi_rho(ng),eta_rho(ng)))
		ALLOCATE(GRIDS(ng)%mask_u(xi_u(ng),eta_u(ng)))
		ALLOCATE(GRIDS(ng)%mask_v(xi_v(ng),eta_v(ng)))
		ALLOCATE(GRIDS(ng)%s_rho(s_rho(ng)))
		ALLOCATE(GRIDS(ng)%s_w(s_w(ng)))
		ALLOCATE(GRIDS(ng)%cs_r(s_rho(ng)))
		ALLOCATE(GRIDS(ng)%cs_w(s_w(ng)))
		ALLOCATE(GRIDS(ng)%z_w(s_w(ng),Ngrid,3))
		ALLOCATE(GRIDS(ng)%z_rho(s_rho(ng),Ngrid,3))
		ALLOCATE(GRIDS(ng)%xi(xi_rho(ng),eta_rho(ng)))
		ALLOCATE(GRIDS(ng)%eta(xi_rho(ng),eta_rho(ng)))
		ALLOCATE(GRIDS(ng)%scl(Ngrid,2))
		ALLOCATE(GRIDS(ng)%off(Ngrid,2))
        
     
		
    ! READ IN DATA FROM NETCDF FILE TO VARIABLES
	
	
      ! rho grid mask
		STATUS = NF90_INQ_VARID(NCID,'mask_rho',VID)
		STATUS = NF90_GET_VAR(NCID,VID,GRIDS(ng)%mask_rho)
		if (STATUS .NE. NF90_NOERR) then
			write(*,*) 'Problem read mask_rho'
			err = 40 
			call errorHandler(header,-1)
		endif

      ! u grid mask
		STATUS = NF90_INQ_VARID(NCID,'mask_u',VID)
		STATUS = NF90_GET_VAR(NCID,VID,GRIDS(ng)%mask_u)
		if(STATUS .NE. NF90_NOERR) then
			write(*,*)'Problem read mask_u'
			err = 40 
			call errorHandler(header,-1)
		endif

      ! v grid mask
		STATUS = NF90_INQ_VARID(NCID,'mask_v',VID)
		STATUS = NF90_GET_VAR(NCID,VID,GRIDS(ng)%mask_v)
		if(STATUS .NE. NF90_NOERR) then
			write(*,*)'Problem read mask_v'
			err = 40 
			call errorHandler(header,-1)
		endif

		 ! Longitude
		STATUS = NF90_INQ_VARID(NCID,'lon_rho',VID)
		STATUS = NF90_GET_VAR(NCID,VID,GRIDS(ng)%lon_rho)
		if(STATUS .NE. NF90_NOERR) then
			write(*,*)'Problem read lon_rho'
			err = 40 
			call errorHandler(header,-1)
		endif
		
			 ! Latitude
		STATUS = NF90_INQ_VARID(NCID,'lat_rho',VID)
		STATUS = NF90_GET_VAR(NCID,VID,GRIDS(ng)%lat_rho)
		if(STATUS .NE. NF90_NOERR) then
			write(*,*)'Problem read lat_rho'
			err = 40 
			call errorHandler(header,-1)
		endif	
		
			 ! Latitude
		STATUS = NF90_INQ_VARID(NCID,'h',VID)
		STATUS = NF90_GET_VAR(NCID,VID,GRIDS(ng)%H)
		if(STATUS .NE. NF90_NOERR) then
			write(*,*)'Problem read H'
			err = 40 
			call errorHandler(header,-1)
		endif	
			 ! PM
		STATUS = NF90_INQ_VARID(NCID,'pm',VID)
		STATUS = NF90_GET_VAR(NCID,VID,GRIDS(ng)%pm)
		if(STATUS .NE. NF90_NOERR) then
			write(*,*)'Problem read pm'
			err = 40 
			call errorHandler(header,-1)
		endif
			 ! PN
		STATUS = NF90_INQ_VARID(NCID,'pn',VID)
		STATUS = NF90_GET_VAR(NCID,VID,GRIDS(ng)%pn)
		if(STATUS .NE. NF90_NOERR) then
			write(*,*)'Problem read pn'
			err = 40 
			call errorHandler(header,-1)
		endif
			 ! Angle
		STATUS = NF90_INQ_VARID(NCID,'angle',VID)
		STATUS = NF90_GET_VAR(NCID,VID,GRIDS(ng)%angle)
		if(STATUS .NE. NF90_NOERR) then
			write(*,*)'Problem read pn'
			err = 40 
			call errorHandler(header,-1)
		endif
		
		STATUS = NF90_INQ_VARID(NCID,'s_rho',VID)
		STATUS = NF90_GET_VAR(NCID,VID,GRIDS(ng)%s_rho)
		if (STATUS .NE. NF90_NOERR) then
			write(*,*) 'Problem read s_rho'
			write(*,*) NF90_STRERROR(STATUS)
			call errorHandler(header,-1)
		endif

      ! Cs value on rho grid (Cs_r)
		STATUS = NF90_INQ_VARID(NCID,'Cs_r',VID)
		
		STATUS = NF90_GET_VAR(NCID,VID,GRIDS(ng)%cs_r)
		if (STATUS .NE. NF90_NOERR) then
			write(*,*) 'Problem read CS_r'
			write(*,*) NF90_STRERROR(STATUS)
			call errorHandler(header,-1)
		endif
      ! s-coordinate on w grid (sc_w)
		STATUS = NF90_INQ_VARID(NCID,'s_w',VID)
		STATUS = NF90_GET_VAR(NCID,VID,GRIDS(ng)%s_w)
		if (STATUS .NE. NF90_NOERR) then
			write(*,*) 'Problem read s_w'
			write(*,*) NF90_STRERROR(STATUS)
			call errorHandler(header,-1)
		endif
      ! Cs value on w grid (Cs_w)
		STATUS = NF90_INQ_VARID(NCID,'Cs_w',VID)
		STATUS = NF90_GET_VAR(NCID,VID,GRIDS(ng)%cs_w)
		if (STATUS .NE. NF90_NOERR) then
			write(*,*) 'Problem read cs_w'
			write(*,*) NF90_STRERROR(STATUS)
			call errorHandler(header,-1)
		endif
		 ! Ocean_time
		STATUS = NF90_INQ_VARID(NCID,trim(time_vname),VID)
		STATUS = NF90_GET_VAR(NCID,VID,reftime)
		if (STATUS .NE. NF90_NOERR) then
			write(*,*) 'Problem reading Time variable:'
			write(*,*) trim(time_vname)
			write(*,*) NF90_STRERROR(STATUS)
			call errorHandler(header,-1)
	    endif
		
		
		
		STATUS = nf90_get_att(NCID, VID,'units', strtmp) 
		if (STATUS .NE. NF90_NOERR) then
			write(*,*) 'NO time units'
			write(*,*) NF90_STRERROR(STATUS)
			call errorHandler(header,-1)
		endif

		
		if(index(strtmp,'days') .gt. 0)then
			reftime=reftime*dble(86400.0)
			write(*,*) 'time units in input file are days'
			time_units(1:7)='seconds'
			time_units(8:100)=strtmp(5:100-3)
		elseif(index(strtmp,'hours') .gt. 0)then
			write(*,*) 'time units in input file are hours'
			reftime=reftime*dble(3600.0)
			time_units(1:7)='seconds'
			time_units(8:100)=strtmp(6:100-2)
		elseif(index(strtmp,'minutes') .gt. 0)then
			reftime=reftime*dble(60.0)
			write(*,*) 'time units in input file are minutes'
			time_units(1:7)='seconds'
			time_units(8:100)=strtmp(8:100)
		!elseif(index(strtmp,'seconds') .gt. 0)then	
		elseif(index(strtmp,'second') .gt. 0)then	!ELI
			write(*,*) 'time units in input file are seconds'
			time_units=strtmp
		else
			write(*,*) 'Not a recognized time unit'
			call errorHandler(header,-1)
		endif
		STATUS = NF90_CLOSE(NCID)
		if(STATUS /= NF90_NOERR) then
		write(*,*)'Problem closing NCID'
		err = 50
	    call errorHandler(header,-1)
		endif
		

	    write(*,*) time_units
		write(*,*) reftime

  ! ********************** MAKE Xi,eta griod **********************
  
  		 do i=1,xi_rho(ng)
			 do j=1,eta_rho(ng)
				GRIDS(ng)%xi(i,j)=i
				GRIDS(ng)%eta(i,j)=j
			 enddo
	    enddo

  
		GRIDS(ng)%scl=1.0d0
		GRIDS(ng)%off=0.0d0
		grefine=1.0d0
	enddo
	do ng=1,Ngrid
		do ng2=1,Ngrid
			dng=abs(ng-ng2)
			if (dng.GT.0.0) then
				do nref=1,dng
					grefine(ng,ng2)=refine(dng)*grefine(ng,ng2)
				enddo
			endif
		enddo
	enddo
		
	!!!!!!!!!!!THIS IS a bit kludgy. 
	!!!!! setting scales and offsets for coordinate transform. 
	do ng2=1,Ngrid
		do ng=1,Ngrid
			loop1:do i=1,xi_rho(ng)
					do j=1,eta_rho(ng)
						do ii=1,xi_rho(ng2)
							do jj=1,eta_rho(ng2)
							if ((GRIDS(ng)%lon_rho(i,j).EQ.GRIDS(ng2)%lon_rho(ii,jj)).and. &
							(GRIDS(ng)%lat_rho(i,j).EQ.GRIDS(ng2)%lat_rho(ii,jj))) then
							
								exit loop1
							endif
						enddo
					enddo
				enddo
			enddo loop1
			GRIDS(ng)%scl(ng2,:)=1.0d0/grefine(ng,ng2)
			GRIDS(ng2)%scl(ng,:)=grefine(ng2,ng)
		
			GRIDS(ng)%off(ng2,1)=dble(ii)-GRIDS(ng)%scl(ng2,1)*dble(i)
			GRIDS(ng)%off(ng2,2)=dble(jj)-GRIDS(ng)%scl(ng2,2)*dble(j)
			GRIDS(ng2)%off(ng,1)=dble(i)-GRIDS(ng2)%scl(ng,1)*dble(ii)
			GRIDS(ng2)%off(ng,2)=dble(j)-GRIDS(ng2)%scl(ng,2)*dble(jj)
		
		
		
		
		enddo 
	enddo
	
	! do ng=1,Ngrid
		! write(*,*) '1--------------'
		! do ng2=1,Ngrid
			! write(*,*) ng,ng2
			! write(*,*) grefine(ng,ng2)
			! write(*,*) GRIDS(ng)%off(ng2,1),GRIDS(ng2)%off(ng,2)
			! write(*,*) GRIDS(ng)%scl(ng2,1),GRIDS(ng2)%scl(ng,2)
		
		! enddo
	! enddo
    !If IOSTAT is present, set return value to error code
    !IF(PRESENT(IOSTAT)) IOSTAT = err
    !  0=No Errors                 30=Error allocating arrays
    ! 10=Error Opening NCgridfile  40=Error getting variables
    ! 20=Error getting dimensions  50=Error Closing NCgridfile

  END SUBROUTINE InitGrid
  
  SUBROUTINE errorHandler(header, flag)
    IMPLICIT NONE
    CHARACTER(LEN=120), INTENT(IN) :: header
    INTEGER, INTENT(IN) :: flag
	
    WRITE(*,"(A120)")header               !print error message in report.txt
    STOP

   
  END SUBROUTINE errorHandler

  
   DOUBLE PRECISION FUNCTION getSlevel(zeta,depth,ng,i)
		!This function returns the depth of the current s-level
		USE PARAM_MOD, ONLY: s_rho,s_w,Vtransform,Vstretching,hc 
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: ng,i 
		DOUBLE PRECISION, INTENT(IN) :: zeta,depth
		DOUBLE PRECISION :: h,S


		! convert negative depth to positive depth
		h = DBLE(-1.0) * depth
		

		SELECT CASE(Vtransform(ng))

		  CASE(1)  !Rutgers-ROMS formulation, eqn (1) of 
			!https://www.myroms.org/wiki/index.php/Vertical_S-coordinate

			S = hc(ng)*GRIDS(ng)%s_rho(i)+(h-hc(ng))*GRIDS(ng)%cs_r(i)
			getSlevel = S+zeta*(DBLE(1.0)+S/h)

		  CASE(2)  !UCLA-formulation, eqn(2) of 
			!https://www.myroms.org/wiki/index.php/Vertical_S-coordinate

			S = (hc(ng)*GRIDS(ng)%s_rho(i)+h*GRIDS(ng)%cs_r(i)) / (hc(ng)+h)
			getSlevel = zeta+(zeta+h)*S

		  CASE(3)  !Song, Y. and D. B. Haidvogel, 1994: A semi-implicit
			!ocean circulation model using a generalized topography-following 
			!coordinate system, J. Comp. Phys., 115 (1), 228-244.

			getSlevel = zeta*(DBLE(1.0)+GRIDS(ng)%s_rho(i))+hc(ng)*GRIDS(ng)%s_rho(i)+(h-hc(ng))*GRIDS(ng)%cs_r(i)

		  CASE DEFAULT
			write(*,*) 'ERROR: Illegal Vtransform number'
			write(*,*) ' '
			write(*,*) 'The Program Cannot Continue and Will Terminate'
			stop

		END SELECT

	  END FUNCTION getSlevel


  DOUBLE PRECISION FUNCTION getWlevel(zeta,depth,ng,i)
		!This function returns the depth of the current s-level
		USE PARAM_MOD, ONLY: s_rho,s_w,Vtransform,Vstretching,hc 
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: ng,i 
		DOUBLE PRECISION, INTENT(IN) :: zeta,depth
		DOUBLE PRECISION :: h,S


		! convert negative depth to positive depth
		h = DBLE(-1.0) * depth
		

		SELECT CASE(Vtransform(ng))

		  CASE(1)  !Rutgers-ROMS formulation, eqn (1) of 
			!https://www.myroms.org/wiki/index.php/Vertical_S-coordinate

			S = hc(ng)*GRIDS(ng)%s_w(i)+(h-hc(ng))*GRIDS(ng)%cs_w(i)
			getWlevel = S+zeta*(DBLE(1.0)+S/h)

		  CASE(2)  !UCLA-formulation, eqn(2) of 
			!https://www.myroms.org/wiki/index.php/Vertical_S-coordinate

			S = (hc(ng)*GRIDS(ng)%s_w(i)+h*GRIDS(ng)%cs_w(i)) / (hc(ng)+h)
			getWlevel = zeta+(zeta+h)*S

		  CASE(3)  !Song, Y. and D. B. Haidvogel, 1994: A semi-implicit
			!ocean circulation model using a generalized topography-following 
			!coordinate system, J. Comp. Phys., 115 (1), 228-244.

			getWlevel = zeta*(DBLE(1.0)+GRIDS(ng)%s_w(i))+hc(ng)*GRIDS(ng)%s_w(i)+(h-hc(ng))*GRIDS(ng)%cs_w(i)

		  CASE DEFAULT
			write(*,*) 'ERROR: Illegal Vtransform number'
			write(*,*) ' '
			write(*,*) 'The Program Cannot Continue and Will Terminate'
			stop

		END SELECT

	  END FUNCTION getWlevel
	

END MODULE GRID_MOD 
