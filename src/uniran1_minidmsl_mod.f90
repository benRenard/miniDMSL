!******************************************************************
! (C) Copyright 2002-2020  ---  Dmitri Kavetski  ---  All rights reserved
!******************************************************************************
!           UNIFORM RANDOM NUMBER GENERATOR: SERIAL IMPLEMENTATION
!******************************************************************************
module uniran1_dmsl_mod
! Purpose: Mersenne Twister random number generator, wrapped into DMSL uniran interface
! see module mt19937-64.f95 for details
!-----
!=======================================================================
! TABLE OF CONTENTS
!
! seed_uniran(put,get):sub=seeds the sequence with user (or default) seed
! uniran_s(harvest):sub=harvest->scalar uniform random numbers
!
!=======================================================================
use kinds_dmsl_kit    ! numeric kind definitions
! NB: do not declare additional modules for global use
implicit none
!**********************************
! Procedure availability
private ! to protect low level engines and specific procedure names
!-----
! List of public generic procedures (see Table of Contents for description)
public::seed_uniran,uniran_s
!***
! private internal variables (current seed)
! The module is designed to provide maximum security for the seed, but
! also allow (controlled) access to its value.
!integer(mik),parameter::uniran1_idum_0=-460053  ! default initial seed
!integer(mik),SAVE::idum=uniran1_idum_0          ! current seed
!logical(mlk),SAVE::restart=.true.
!----------------------------------------------------
contains
!----------------------------------------------------
subroutine seed_uniran(put,get,CPUtime)
! Purpose: Processes the serial uniform random number seed.
! NB: Implements the only secure way to access/modify the random seed.
! ---
! Usage
!   * When put is present will reset the seed to "put"
!   * When get is present will return current seed in "get"
!   * If nothing present, do nothing
! ---
! Programmer: Dmitri Kavetski
! Kreated: 1 April 2002 AD
! Hisstory: various
!           18 Feb 2010 AD, Hammo
! ---
! Comments:
! 1. Differs from the F-90 intrinsic random_seed in that
!    no error is signalled when both put and get are present.
! 2. [MT] CPU_time - intrinsic, which allows user to set the random seed the CPU_time
! ---
use mt19937_64,only: init_genrand64
implicit none
! dummies
integer(mik),intent(in), optional::put  ! user-specified seed
integer(mik),intent(out),optional::get  ! current seed
logical(mlk),intent(in), optional::CPUtime  ! current seed
! locals
logical(mlk)::presentPut,presentGet
integer(8)::idum
! Start procedure here
presentPut=present(put); presentGet=present(get)
if(present(CPUtime))then;if(CPUtime)then !MT additions - 3/06/2008
  CALL SYSTEM_CLOCK(COUNT=IDUM)
  call init_genrand64(idum)
  if(presentGet)get=idum
  return
endif;endif
if(presentPut.and.presentGet) then      ! reset to user-provided seed
  idum=put
  call init_genrand64(idum)
  get=idum
  !restart=.true.    ! and return its value.
  !if(put==0)idum=uniran1_idum_0
elseif(presentPut)then  ! reset to user provided seed and restart
  idum=put
  call init_genrand64(idum)
  !restart=.true.
  !if(put==0)idum=uniran1_idum_0
elseif(presentGet)then  ! return current seed value, but do not restart
  get=undefRN
  write(*,*) 'Warning: get seed alone is not implemented'
  !get=idum
else                    ! reset to default seed and restart
  !idum=uniran1_idum_0; restart=.true.
endif
! End procedure here
endsubroutine seed_uniran
!----------------------------------------------------
subroutine uniran_s(harvest)
! Purpose: Packs uniform random numbers into a scalar "harvest"
implicit none
! dummies
real(mrk),intent(out)::harvest
! Start procedure
harvest=rangen()
! End procedure here
endsubroutine uniran_s
!----------------------------------------------------
function rangen()
use mt19937_64,only: genrand64_real3
implicit none
! dummies
real(mrk)::rangen
! locals
rangen=genrand64_real3()
! Start procedure
! End procedure here
endfunction rangen
!----------------------------------------------------
endmodule uniran1_dmsl_mod
!******************************************************************************
!        END UNIFORM RANDOM NUMBER GENERATOR: SERIAL IMPLEMENTATION
!******************************************************************************
