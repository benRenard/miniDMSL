!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! WARNING !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!
! Different from original DMSL file: all compiler directives (!DEC$) have been
! commented out to allow compilation on gfortran 


!*** pure subroutine quicksort_rv(arr,arr2,ascnd,err)
!*** ! Sorts arr into ascending numerical order using Quicksort.
!*** ! err/=EXIT_SUCCESS indicates scratch memory allocation failure.
!*** ! Settings: nn is the size of subarrays sorted by straight insertion
!*** !           nstack is the required auxiliary storage
!*** ! Adapted from Press et al. (1992) Numerical Recipes in Fortran-77.
!*** use utilities_dmsl_kit,only:quickif,swap
!*** implicit none
!*** ! dummies
!*** real(mrk),intent(inout)::arr(:)
!*** real(mrk),intent(inout)::arr2(:) ! only if !DEC$ DEFINED (HAVE_2)
!*** logical(mlk),intent(in),optional::ascnd
integer(mik),intent(out)::err
! local settings
integer(mik),parameter::nn=15,nstack=64
logical(mlk),parameter::ascndDef=.true.
! locals
!*** real(mrk)::a,(a2 - only if HAVE_2)
integer(mik)::istack(nstack)
integer(mik)::n,k,i,j,jstack,L,r
logical(mlk)::ascnd0,doOp
! start procedure here
ascnd0=quickif(ascnd,ascndDef)
err=EXIT_SUCCESS; n=size(arr); jstack=0; L=1; r=n
outerLoop:do
  if(r-L<nn)then      ! subarray small enough:
    do j=L+1,r        ! use straight insertion
      a =arr (j)
!!DEC$ IF DEFINED(HAVE_2)
!      a2=arr2(j)
!!DEC$ END IF
      do i=j-1,L,-1
        if(ascnd0)then; doOp=arr(i)<=a
        else;           doOp=arr(i)>=a; endif
        if(doOp)exit
        arr (i+1)=arr (i)
!!DEC$ IF DEFINED(HAVE_2)
!        arr2(i+1)=arr2(i)
!!DEC$ END IF
      enddo
      arr (i+1)=a
!!DEC$ IF DEFINED(HAVE_2)
!      arr2(i+1)=a2
!!DEC$ END IF
    enddo
    if(jstack==0)return
    r=istack(jstack)   ! pop stack and begin new
    L=istack(jstack-1) ! round of partitioning
    jstack=jstack-2
  else        ! use quicksort: choose median of left, center and right
    k=(L+r)/2 ! elements as partitioning element a. Also rearrange
    call swap(arr (k),arr (L+1))   ! so that a(L) <= a(L+1) <= a(r)
!!DEC$ IF DEFINED(HAVE_2)
!    call swap(arr2(k),arr2(L+1))
!!DEC$ END IF
    if(icomp(L,  r  ))then
      call swap(arr (L  ),arr (r  ))
!!DEC$ IF DEFINED(HAVE_2)
!      call swap(arr2(L  ),arr2(r  ))
!!DEC$ END IF
    endif
    if(icomp(L+1,r  ))then
      call swap(arr (L+1),arr (r  ))
!!DEC$ IF DEFINED(HAVE_2)
!      call swap(arr2(L+1),arr2(r  ))
!!DEC$ END IF
    endif
    if(icomp(L  ,L+1))then
      call swap(arr (L  ),arr (L+1))
!!DEC$ IF DEFINED(HAVE_2)
!      call swap(arr2(L  ),arr2(L+1))
!!DEC$ END IF
    endif
    i=L+1; j=r      ! initialise pointers for partitioning
    a =arr (L+1)    ! partitioning element
!!DEC$ IF DEFINED(HAVE_2)
!    a2=arr2(L+1)
!!DEC$ END IF
    do      ! here is the meat :-)
      do    ! scan up to find element >= a
        i=i+1
        if(ascnd0)then; doOp=arr(i)>=a
        else;           doOp=arr(i)<=a; endif
        if(doOp)exit
      enddo
      do    ! scan down to find element <=a
        j=j-1
        if(ascnd0)then; doOp=arr(j)<=a
        else;           doOp=arr(j)>=a; endif
        if(doOp)exit
      enddo
      if(j<i)exit   ! pointers crossed. partitioning complete
      call swap(arr (i),arr (j)) ! exchange elements
!!DEC$ IF DEFINED(HAVE_2)
!      call swap(arr2(i),arr2(j))
!!DEC$ END IF
    enddo
    arr (L+1)=arr (j); arr (j)=a ! insert partitioning element
!!DEC$ IF DEFINED(HAVE_2)
!    arr2(L+1)=arr2(j); arr2(j)=a2
!!DEC$ END IF
    jstack=jstack+2 ! push pointers on stack; process smaller subarray immediately
    if(jstack>nstack)then   ! stack size exceeded
      err=20; return
    endif
    if(r-i+1>=j-L)then      ! famous "is it an l or a 1" bug? wth ...
      istack(jstack)=r;   istack(jstack-1)=i; r=j-1
    else
      istack(jstack)=j-1; istack(jstack-1)=L; L=i
    endif
  endif
enddo outerLoop
! end main body here
contains
!-------
pure function icomp(i,j)result(doOp)
implicit none
integer(mik),intent(in)::i,j
integer(mik)::swp
logical(mlk)::doOp
if(ascnd0)then; doOp=arr(i)>arr(j)
else;           doOp=arr(i)<arr(j); endif
endfunction icomp
!-------
!*** ! end procedure here
!*** endsubroutine quicksort_rv
