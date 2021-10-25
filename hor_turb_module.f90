MODULE HTURB_MOD

!   A random walk model is used to simulate turbulent particle motion in the 
!   horizontal direction (x- and y- directions).
!
! Created by:             Elizabeth North
! Modified by:            Zachary Schlag
! Created on:             2003
! Last Modified on:       18 Aug 2008
! ROMSPath Version: 1.0.1

IMPLICIT NONE
PUBLIC

CONTAINS

!       ***************** Horizontal Turbulence (RWM) *********************
!     ***********************************************************************
!    ** Random Walk Model (Visser 1997 MEPS 158:275-281) for simulating     **
!    ** turbulent diffusion, applied in the horizontal direction            **
!    **    z(t+1)= z + R[2/r K(z)dt]**0.5                                   **
!    **    where z = particle vertical location at time t                   **
!    **      K  = horizontal diffusivity (KM from QUODDY)                   **
!    **      dt = time step of RWM (deltat)                                 **
!    **      R  = random process with mean = 0 and standard deviation = r.  **
!    **                                                                     **
!    ** Programmed by EW North February 2003 UMCES HPL enorth@hpl.umces.edu **
!    *************************************************************************

  SUBROUTINE HTurb(TurbHx,TurbHy,ng)
    USE PARAM_MOD, ONLY: ConstantHTurb,idt
    USE PDF_MOD,  ONLY: norm
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(OUT) :: TurbHx,TurbHy	
	INTEGER, INTENT(in) ::ng
    DOUBLE PRECISION :: devX,devY,r
    DOUBLE PRECISION :: KM

    r=1.0             ! standard deviation of the random deviate
    KM=ConstantHTurb(ng)  ! constant horizontal diffusivity
	
    devX=norm()       ! the random deviate in the X direction
    devY=norm()       ! the random deviate in the Y direction

    !Apply random walk model to calculate horizontal turbulent 
    !  particle displacement
    TurbHx= devX*(DBLE(2.0)/r * KM *idt)**0.5 
    TurbHy= devY*(DBLE(2.0)/r * KM *idt)**0.5

  END SUBROUTINE HTurb

END MODULE HTURB_MOD
