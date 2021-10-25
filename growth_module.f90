MODULE GROWTH_MOD

! GRowTH MODULE
!
! Growth algorithms and code created by: Heidi FUchs, implemented by Elias Hunter,
! Created on:                              2019
! Last Modified on:                        22 March 2019
! ROMSPath Version: 1.0.1

  IMPLICIT NONE
  PRIVATE
  SAVE

  
  !The following procedures have been made public:
  PUBLIC :: growlarva

CONTAINS



  SUBROUTINE growlarva(P_temp,P_salt,P_age,P_size,P_status)  !Update particle status

    USE PARAM_MOD, ONLY: Growth,mortality,deadage,initsize,maxsize, &
				a0,a1,a2,a3,a4,a5,a6,a7,a8,idt,tempcut
    IMPLICIT NONE

  
    DOUBLE PRECISION, INTENT(IN) ::P_salt,P_temp,P_age 
	DOUBLE PRECISION, INTENT(INOUT) ::P_size,P_status
	DOUBLE PRECISION :: GR,DS

	if (Growth.eq. 1) then
		if (P_age.ge.deadage) then
			P_status=9.0
		endif	
	elseif (Growth.eq. 2) then
		
		if (P_temp.gt.tempcut) then
			GR = a0 + a1*P_temp + a2*P_salt + a3*P_temp*P_salt &
				+ a4*(P_temp**2) + a5*(P_salt**2) + a6*(P_temp**2)*P_salt + &
				a7*P_temp*(P_salt**2) + a8*(P_temp**2)*(P_salt**2)
			if (GR.LT.0.0) then
				GR=0.0
			endif
		else
			GR=0.0
		endif
		
		
		
		DS=GR*DBLE(idt)/86400.0D0
		P_size=P_size+DS
	endif
 

  END SUBROUTINE growlarva

 


END MODULE GROWTH_MOD
