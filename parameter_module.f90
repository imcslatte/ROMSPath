MODULE PARAM_MOD 
 
!  The Parameter Module reads in the include file, ROMSPath.data, making the 
!  parameters declared within available to all the other modules. It also 
!  reads in information from the NetCDF grid file and calculates values of
!  grid specific parameters, making them available to all the other modules.
! 
!  Created by:            Zachary Schlag 
!  Created on:            28 Jul 2008 
!  Last Modified on:         Feb 2013
 
IMPLICIT NONE 
PUBLIC 
SAVE 

include 'ROMSPath.h'

CONTAINS


  SUBROUTINE getParams()
  !Subroutine to read all input parameters from ROMSPath.data 

    character(len=120) :: header
	character(len=256) :: Iname
    integer :: istat,err

    err = 0
	
	CALL getarg(1,Iname)
	
	WRITE(6,*) '---------'
	WRITE(6,*) trim(Iname)
	WRITE(6,*) '----------'
    OPEN(1,file=Iname)                  !--- read control variables:
!    OPEN(1,file='ROMSPath.data')                  !--- read control variables:


	
      IF(err == 0) THEN
        READ(1,nml=numparticles ,IOSTAT=istat)  !--- number of particles
        IF(istat/=0)err = 10
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=timeparam    ,IOSTAT=istat)  !--- time info
        IF(istat/=0)err = 20
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=hydroparam   ,IOSTAT=istat)  !--- hydrodynamics info
        IF(istat/=0)err = 30
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=turbparam    ,IOSTAT=istat)  !--- turbulence info
        IF(istat/=0)err = 40
      ENDIF      
	  IF(err == 0) THEN
        READ(1,nml=advectparam    ,IOSTAT=istat)  !--- Advection info
      IF(istat/=0)err = 45
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=behavparam   ,IOSTAT=istat)  !--- behavior info
        IF(istat/=0)err = 50
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=fuchsparam   ,IOSTAT=istat)  !--- FUCHS  info
        IF(istat/=0)err = 55
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=growparam   ,IOSTAT=istat)  !--- growth info
        IF(istat/=0)err = 60
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=dvmparam     ,IOSTAT=istat)  !--- diurnal vertical migration 
        IF(istat/=0)err = 60
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=settleparam  ,IOSTAT=istat)  !--- settlement info
        IF(istat/=0)err = 70
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=romsgrid     ,IOSTAT=istat)  !--- roms grid
        IF(istat/=0)err = 90
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=romsoutput   ,IOSTAT=istat)  !--- roms history output file
        IF(istat/=0)err = 100
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=parloc       ,IOSTAT=istat)  !--- particle locations
        IF(istat/=0)err = 110
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=HabPolyLoc   ,IOSTAT=istat)  !--- habitat polygon info
        IF(istat/=0)err = 120
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=output       ,IOSTAT=istat)  !--- output related info
        IF(istat/=0)err = 130
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=other        ,IOSTAT=istat)  !--- other misc 
        IF(istat/=0)err = 140
      ENDIF
    CLOSE(1)




    SELECT CASE(err)
      CASE(0)
        header='No Errors'
      CASE(10)
        header='Error when reading numparticles, pls check ROMSPath.data'
      CASE(20)
        header='Error when reading timeparam, pls check ROMSPath.data'
      CASE(30)
        header='Error when reading hydroparam, pls check ROMSPath.data'
      CASE(40)
        header='Error when reading turbparam, pls check ROMSPath.data'
      CASE(45)
        header='Error when reading advectparam, pls check ROMSPath.data'
      CASE(50)
        header='Error when reading behavparam, pls check ROMSPath.data'
      CASE(55)
        header='Error when reading fuchsparam, pls check ROMSPath.data'
      CASE(60)
        header='Error when reading growthparam, pls check ROMSPath.data'
      CASE(70)
        header='Error when reading settleparam, pls check ROMSPath.data'
      CASE(90)
        header='Error when reading romsgrid, pls check ROMSPath.data'
      CASE(100)
        header='Error when reading romsoutput, pls check ROMSPath.data'
      CASE(110)
        header='Error when reading parloc, pls check ROMSPath.data'
      CASE(120)
        header='Error when reading HabPolLoc, pls check ROMSPath.data'
      CASE(130)
        header='Error when reading output, pls check ROMSPath.data'
      CASE(140)
        header='Error when reading other, pls check ROMSPath.data'
      CASE(150)
        header='Error when reading gridinfo, pls check GRID.data'
      CASE DEFAULT
        header='Error: unexpected err number'
    END SELECT

    IF(err/=0) CALL errorHandler(Header,-1)  !print the error message and stop


  END SUBROUTINE getParams


  SUBROUTINE errorHandler(header, flag)
    IMPLICIT NONE
    CHARACTER(LEN=120), INTENT(IN) :: header
    INTEGER, INTENT(IN) :: flag
	
    IF (flag .eq. -1) THEN
      WRITE(*,"(A120)")header               !print error message in report.txt
      STOP
    ELSE
      WRITE(*,"('***** WARNING *****')")    !print warning message to screen
      WRITE(*,"(A120)")header
    ENDIF
   
  END SUBROUTINE errorHandler


END MODULE PARAM_MOD 
