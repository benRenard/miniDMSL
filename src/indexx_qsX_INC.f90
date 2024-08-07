!*** pure subroutine indexx_qsr(arr,indx,ascnd,err,message)
!*** ! Purpose: Creates an index of arr, such that arr(indx) is sorted
!*** ! Algorithm: Quicksort sorting.
!*** ! Adapted from Press et al Numerical Recipes in F-77 and F-90.
!*** use utilities_dmsl_kit,only:assertEq,arthsi,quickif,swap
!*** implicit none
!*** ! dummies
!*** real(mrk),intent(in)::arr(:)
!*** logical(mlk),intent(in),optional::ascnd
integer(mik),intent(out)::indx(:)
integer(mik),intent(out)::err
character(*),intent(out)::message
! local settings
!*** character(*),parameter::procnam="indexx_qsr"
integer(mik),parameter::nn=15,nstack=64
logical(mlk),parameter::ascndDef=.true.
! locals
!*** real(mrk)::a
integer(mik)::istack(nstack)
integer(mik)::n,k,i,j,jstack,L,r,indxt
logical(mlk)::ascnd0,doOp
! start procedure here
err=EXIT_SUCCESS
ascnd0=quickif(ascnd,ascndDef)
call assertEq(size(arr),size(indx),doOp,n)
if(.not.doOp)then  ! dimension error
  message="f-"//procnam//"/arrayDimError"
  err=10; indx=-10; return
endif
indx=arthsi(n)
jstack=0; L=1; r=n
outerLoop:do
  if(r-L<nn)then
    do j=L+1,r
      indxt=indx(j); a=arr(indxt)
      do i=j-1,L,-1
        if(ascnd0)then; doOp=arr(indx(i))<=a
        else;           doOp=arr(indx(i))>=a; endif
        if(doOp)exit
        indx(i+1)=indx(i)
      enddo
      indx(i+1)=indxt
    enddo
    if(jstack==0)return
    r=istack(jstack); L=istack(jstack-1); jstack=jstack-2
  else
    k=(L+r)/2
    call swap(indx(k),indx(L+1))
    call icomp_xchg(indx(L  ),indx(r  ))
    call icomp_xchg(indx(L+1),indx(r  ))
    call icomp_xchg(indx(L  ),indx(L+1))
    i=L+1; j=r; indxt=indx(L+1); a=arr(indxt)
    do
      do
        i=i+1
        if(ascnd0)then; doOp=arr(indx(i))>=a
        else;           doOp=arr(indx(i))<=a; endif
        if(doOp)exit
      enddo
      do
        j=j-1
        if(ascnd0)then; doOp=arr(indx(j))<=a
        else;           doOp=arr(indx(j))>=a; endif
        if(doOp)exit
      enddo
      if(j<i)exit
      call swap(indx(i),indx(j))
    enddo
    indx(L+1)=indx(j); indx(j)=indxt; jstack=jstack+2
    if(jstack>nstack)then
      write(message,'(i0,a,i0)')"f-"//procnam//"/exceededStack[nstack=",nstack,"]"
      err=20; indx=-20; return
    endif
    if(r-i+1>=j-L)then
      istack(jstack)=r;   istack(jstack-1)=i; r=j-1
    else
      istack(jstack)=j-1; istack(jstack-1)=L; L=i
    endif
  endif
enddo outerLoop
! end main body here
contains
!-------
pure subroutine icomp_xchg(i,j)
implicit none
integer(mik),intent(inout)::i,j
integer(mik)::swp
logical(mlk)::doOp
if(ascnd0)then; doOp=arr(i)>arr(j)
else;           doOp=arr(i)<arr(j); endif
if(doOp)then
  swp=i; i=j; j=swp
endif
endsubroutine icomp_xchg
!-------
!*** ! end procedure here
!*** endsubroutine indexx_qsr
