MODULE RANDOM_MOD
IMPLICIT NONE
PUBLIC
SAVE

! ROMSPath Version: 1.0.1
! The following Mersenne Twister program, mt19937ar.f, is used to generate  
! random numbers between 0 and 1 from a uniform distribution. The Mersenne
! Twister is a fast random number generator with a period of 2^19937-1. 
! It was downloaded from the following website:
!   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/FORTRAN/mt19937ar.f
! See the Mersenne Twister Home Page for more information:
!   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
!
! Zachary Schlag converted the program to F90 for use in ROMSPath (8/29/08). 

! ************* Mersenne Twister ****************
!
!  A C-program for MT19937, with initialization improved 2002/1/26.
!  Coded by Takuji Nishimura and Makoto Matsumoto.
!
!  Before using, initialize the state by using init_genrand(seed)  
!  or init_by_array(init_key, key_length).
!
!  Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
!  All rights reserved.                          
!  Copyright (C) 2005, Mutsuo Saito,
!  All rights reserved.                          
!
!  Redistribution and use in source and binary forms, with or without
!  modification, are permitted provided that the following conditions
!  are met:
!
!    1. Redistributions of source code must retain the above copyright
!       notice, this list of conditions and the following disclaimer.
!
!    2. Redistributions in binary form must reproduce the above copyright
!       notice, this list of conditions and the following disclaimer in the
!       documentation and/or other materials provided with the distribution.
!
!    3. The names of its contributors may not be used to ENDorse or promote 
!       products derived from this software without specific prior written 
!       permission.
!
!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
!  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER
!  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!
!  Any feedback is very welcome.
!  http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
!  email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
!
!-----------------------------------------------------------------------
!  FORTRAN77 translation by Tsuyoshi TADA. (2005/12/19)
!
!  FORTRAN90 translation by Zachary Schlag (2008/08/29)
!
!     ---------- initialize routines ----------
!  SUBROUTINE init_genrand(seed): initialize with a seed
!  SUBROUTINE init_by_array(init_key,key_length): initialize by an array
!
!     ---------- generate FUNCTIONs ----------
!  INTEGER FUNCTION genrand_int32(): signed 32-bit INTEGER
!  INTEGER FUNCTION genrand_int31(): unsigned 31-bit INTEGER
!  DOUBLE PRECISION FUNCTION genrand_real1(): [0,1] with 32-bit resolution
!  DOUBLE PRECISION FUNCTION genrand_real2(): [0,1) with 32-bit resolution
!  DOUBLE PRECISION FUNCTION genrand_real3(): (0,1) with 32-bit resolution
!  DOUBLE PRECISION FUNCTION genrand_res53(): (0,1) with 53-bit resolution
!
!  This program uses the following non-standard intrinsics.
!    ishft(i,n): If n>0, shifts bits in i by n positions to left.
!                If n<0, shifts bits in i by n positions to right.
!    iand (i,j): Performs logical AND on corresponding bits of i and j.
!    ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
!    ieor (i,j): Performs exclusive OR on corresponding bits of i and j.
!

  INTEGER, PRIVATE, PARAMETER :: N = 624
  INTEGER, PRIVATE, PARAMETER :: M = 397
  INTEGER, PRIVATE, PARAMETER :: DONE = 123456789

  INTEGER, PRIVATE :: ALLBIT_MASK
  INTEGER, PRIVATE :: TOPBIT_MASK
  INTEGER, PRIVATE :: UPPER_MASK
  INTEGER, PRIVATE :: LOWER_MASK
  INTEGER, PRIVATE :: MATRIX_A
  INTEGER, PRIVATE :: T1_MASK
  INTEGER, PRIVATE :: T2_MASK
  INTEGER, PRIVATE :: mag01(0:1)
  INTEGER, PRIVATE :: mt(0:N-1)
  INTEGER, PRIVATE :: mti
  INTEGER, PRIVATE :: initialized

CONTAINS


!-----------------------------------------------------------------------
!     initialize mt(0:N-1) with a seed
!-----------------------------------------------------------------------
  SUBROUTINE init_genrand(s)
    INTEGER, INTENT(IN) :: s

    CALL mt_initln
    mt(0)=iand(s,ALLBIT_MASK)
    do mti=1,N-1
      mt(mti)=1812433253*ieor(mt(mti-1),ishft(mt(mti-1),-30))+mti
      mt(mti)=iand(mt(mti),ALLBIT_MASK)
    enddo
    initialized=DONE

  END SUBROUTINE init_genrand


!-----------------------------------------------------------------------
!     initialize by an array with array-length
!     init_key is the array for initializing keys
!     key_length is its length
!-----------------------------------------------------------------------
  SUBROUTINE init_by_array(init_key,key_length)
    INTEGER, INTENT(IN) :: init_key(0:*)
    INTEGER, INTENT(IN) :: key_length

    INTEGER :: i,j,k

    CALL init_genrand(19650218)
    i=1
    j=0
    do k=max(N,key_length),1,-1
      mt(i)=ieor(mt(i),ieor(mt(i-1),ishft(mt(i-1),-30))*1664525)+init_key(j)+j
      mt(i)=iand(mt(i),ALLBIT_MASK)
      i=i+1
      j=j+1
      if(i.ge.N)then
        mt(0)=mt(N-1)
        i=1
      endif
      if(j.ge.key_length)then
        j=0
      endif
    enddo
    do k=N-1,1,-1
      mt(i)=ieor(mt(i),ieor(mt(i-1),ishft(mt(i-1),-30))*1566083941)-i
      mt(i)=iand(mt(i),ALLBIT_MASK)
      i=i+1
      if(i.ge.N)then
        mt(0)=mt(N-1)
        i=1
      endif
    enddo
    mt(0)=TOPBIT_MASK

  END SUBROUTINE init_by_array


!-----------------------------------------------------------------------
!     generates a random number on [0,0xffffffff]-interval
!-----------------------------------------------------------------------
  INTEGER FUNCTION genrand_int32()
    INTEGER :: y,kk

    if(initialized.ne.DONE)then
      CALL init_genrand(21641)
    endif

    if(mti.ge.N)then
      do kk=0,N-M-1
        y=ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
        mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
      enddo
      do kk=N-M,N-1-1
        y=ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
        mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
      enddo
      y=ior(iand(mt(N-1),UPPER_MASK),iand(mt(0),LOWER_MASK))
      mt(kk)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
      mti=0
    endif

    y=mt(mti)
    mti=mti+1

    y=ieor(y,ishft(y,-11))
    y=ieor(y,iand(ishft(y,7),T1_MASK))
    y=ieor(y,iand(ishft(y,15),T2_MASK))
    y=ieor(y,ishft(y,-18))

    genrand_int32=y

  END FUNCTION genrand_int32


!-----------------------------------------------------------------------
!     generates a random number on [0,0x7fffffff]-interval
!-----------------------------------------------------------------------
  INTEGER FUNCTION genrand_int31()

    genrand_int31=int(ishft(genrand_int32(),-1))

  END FUNCTION genrand_int31


!-----------------------------------------------------------------------
!     generates a random number on [0,1]-real-interval
!-----------------------------------------------------------------------
  DOUBLE PRECISION FUNCTION genrand_real1()
    DOUBLE PRECISION :: r

    r=dble(genrand_int32())
    if(r.lt.0.d0)r=r+2.d0**32
    genrand_real1=r/4294967295.d0

  END FUNCTION genrand_real1


!-----------------------------------------------------------------------
!     generates a random number on [0,1)-real-interval
!-----------------------------------------------------------------------
  DOUBLE PRECISION FUNCTION genrand_real2()
    DOUBLE PRECISION :: r

    r=dble(genrand_int32())
    if(r.lt.0.d0)r=r+2.d0**32
    genrand_real2=r/4294967296.d0

  END FUNCTION genrand_real2


!-----------------------------------------------------------------------
!     generates a random number on (0,1)-real-interval
!-----------------------------------------------------------------------
  DOUBLE PRECISION FUNCTION genrand_real3()
    DOUBLE PRECISION :: r

    r=dble(genrand_int32())
    if(r.lt.0.d0)r=r+2.d0**32
    genrand_real3=(r+0.5d0)/4294967296.d0

  END FUNCTION genrand_real3


!-----------------------------------------------------------------------
!     generates a random number on [0,1) with 53-bit resolution
!-----------------------------------------------------------------------
  DOUBLE PRECISION FUNCTION genrand_res53()
    DOUBLE PRECISION :: a,b

    a=dble(ishft(genrand_int32(),-5))
    b=dble(ishft(genrand_int32(),-6))
    if(a.lt.0.d0)a=a+2.d0**32
    if(b.lt.0.d0)b=b+2.d0**32
    genrand_res53=(a*67108864.d0+b)/9007199254740992.d0

  END FUNCTION genrand_res53


!-----------------------------------------------------------------------
!     initialize large number (over 32-bit constant number)
!-----------------------------------------------------------------------
  SUBROUTINE mt_initln()

    TOPBIT_MASK=1073741824
    TOPBIT_MASK=ishft(TOPBIT_MASK,1)
    ALLBIT_MASK=2147483647
    ALLBIT_MASK=ior(ALLBIT_MASK,TOPBIT_MASK)
    UPPER_MASK=TOPBIT_MASK
    LOWER_MASK=2147483647
    MATRIX_A=419999967
    MATRIX_A=ior(MATRIX_A,TOPBIT_MASK)
    T1_MASK=489444992
    T1_MASK=ior(T1_MASK,TOPBIT_MASK)
    T2_MASK=1875247104
    T2_MASK=ior(T2_MASK,TOPBIT_MASK)
    mag01(0)=0
    mag01(1)=MATRIX_A

  END SUBROUTINE mt_initln
  
  SUBROUTINE init_random_seed(seed)
			
            implicit none
#ifdef IFORT
			integer, external :: GETPID
#endif	
			
			INTEGER, INTENT(OUT) :: seed

            integer :: i, un, istat, dt(8), pid, t(2), s
            integer(8) :: count, tms
          
            ! First try if the OS provides a random number generator
            open(newunit=un, file="/dev/urandom", access="stream", &
                 form="unformatted", action="read", status="old", iostat=istat)
			
            if (istat == 0) then
				write(*,*) 'USING SYSTEM RANDOM SEED'	
               read(un) seed
               close(un)
            else
			
				write(*,*) 'USING  TIME/PID SEED'
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
               call system_clock(count)
               if (count /= 0) then
                  t = transfer(count, t)
               else
                  call date_and_time(values=dt)
                  tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
                  t = transfer(tms, t)
               end if
               s = ieor(t(1), t(2))
               pid = getpid() + 1099279 ! Add a prime
               s = ieor(s, pid)
               ! if (n >= 3) then
                  ! seed(1) = t(1) + 36269
                  ! seed(2) = t(2) + 72551
                  ! seed(3) = pid
                  ! if (n > 3) then
                     ! seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
                  ! end if
               ! else
               seed = s 
               ! end if
            end if
  END SUBROUTINE init_random_seed
  
  

END MODULE RANDOM_MOD