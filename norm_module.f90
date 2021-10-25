MODULE PDF_MOD

!  The Norm Module contains the function Norm, which returns a random number 
!    (deviate) drawn from a normal distribution with zero mean and unit 
!    variance (i.e., standard deviation = 1).
!
!  Created by:            Elizabeth North
!  Modified by:           Zachary Schlag
!  Created on:            22 Aug 2008
!  Last Modified on:         Feb 2011
! ROMSPath Version: 1.0.1

IMPLICIT NONE
PUBLIC

CONTAINS

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~                                                               ~~
! ~~                       FUNCTION norm                           ~~
! ~~                                                               ~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  DOUBLE PRECISION FUNCTION norm()
    ! This function returns a normally distributed deviate with 
    !   zero mean and unit variance (i.e., standard deviation = 1). 
    !   By E. North, 8/22/08
    USE PARAM_MOD,  ONLY: PI
    USE RANDOM_MOD, ONLY: genrand_real3
    IMPLICIT NONE

    DOUBLE PRECISION :: dev1,dev2

    dev1 = genrand_real3()
    dev2 = genrand_real3()
    norm = sqrt(DBLE(-2.)*log(dev1)) * cos(DBLE(2.)*PI*dev2)

  END FUNCTION norm

END MODULE NORM_MOD