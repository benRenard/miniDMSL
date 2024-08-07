module linalg_dmsl_kit

use kinds_dmsl_kit    ! numeric kind definitions
implicit none
private

! Publicly-available stuff
public::&
  choles_dcmp,choles_invrt,&
  lu_dcmp,lu_invrt

! Overloads
!--
interface choles_dcmp
  module procedure choles_dcmp_engine1,&
    !choles_dcmp_engine2,choles_dcmp_engine3,&
    !choles_dcmp_gershgorin,&
    choles_dcmp_var1,choles_dcmp_var2
endinterface choles_dcmp
!--
interface choles_invrt
  module procedure choles_invrt_dcmp_engine,choles_invrt1,choles_invrt2, &
    choles_dcmp_invrt2
endinterface choles_invrt

!--
interface lu_dcmp
  module procedure lu_dcmp1
endinterface lu_dcmp
!--
interface lu_invrt
  module procedure lu_invrt1
endinterface lu_invrt
!--
interface lu_fwbw
  module procedure lu_fwbw1
endinterface lu_fwbw
!--
interface getCondest
  module procedure getCondest1,getCondest2,getCondest3
endinterface getCondest

contains

!----------------------------------------------------
! 3. CHOLESKI METHODS
!----------------------------------------------------
pure subroutine choles_dcmp_var1(a,posDefinite,iBad,logDet,condest,err,message)
! Purpose: overloaded for in-place decomposition of a
use utilities_dmsl_kit,only:zero,assertEq
implicit none
! dummies
real(mrk),intent(inout)::a(:,:)
logical(mlk),intent(out)::posDefinite
integer(mik),intent(out),optional::iBad
real(mrk),intent(out),optional::logDet
real(mrk),intent(out),optional::condest
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
character(*),parameter::procnam="choles_dcmp_var1"
real(mrk)::Ld(size(a,dim=1))
integer(mik)::i,n
! Start procedure here
call choles_dcmp_engine1(a,Ld,posDefinite,iBad,logDet,condest,err,message)
if(err/=EXIT_SUCCESS)then
  message="f-"//procnam//"/&"//message
endif
if(.not.posDefinite)return  ! decomposition failed
n=size(Ld)
forall(i=1:n)
  a(i,i)    =Ld(i) ! insert Ld into diagonal of a
  a(i,i+1:n)=zero  ! zero upper triangle (remainder of row)
endforall
! End procedure here
endsubroutine choles_dcmp_var1
!----------------------------------------------------
pure subroutine choles_dcmp_var2(a,Lmat,posDefinite,iBad,logDet,condest,err,message)
! Purpose: overloaded to decompose a and store in Lmat
use utilities_dmsl_kit,only:assertEq
implicit none
! dummies
real(mrk),intent(in)::a(:,:)
real(mrk),intent(out)::Lmat(:,:)
logical(mlk),intent(out)::posDefinite
integer(mik),intent(out),optional::iBad
real(mrk),intent(out),optional::logDet
real(mrk),intent(out),optional::condest
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
character(*),parameter::procnam="choles_dcmp_var2"
logical(mlk)::conform
integer(mik)::n
! Start procedure here
call assertEq(size(a,1),size(a,2),size(Lmat,1),size(Lmat,2),conform,n)
if(.not.conform) then
  err=10; message="f-"//procnam//"/nonConform(a,Lmat)"; return
elseif(n<=0) then
  err=20; message="f-"//procnam//"/bug?dim<=0"; return
endif
Lmat(:,:)=a(:,:)  ! copy and overwrite
call choles_dcmp_var1(Lmat,posDefinite,iBad,logDet,condest,err,message)
if(err/=EXIT_SUCCESS)then
  message="f-"//procnam//"/&"//message
endif
! End procedure here
endsubroutine choles_dcmp_var2
!----------------------------------------------------
pure subroutine choles_dcmp_engine1(A,Ld,posDefinite,iBad,logDet,condest,err,message)
! Purpose: Implements the classic Choleski decomposition algorithm of a
!          positive definite matrix A.
! Programmer: Dmitri Kavetski.
! NB: This routine accesses the upper triangle of a only and
!     overwrites the lower triangle of A. Hence, if a non-symmetric
!     matrix is supplied, an error will not be detected.
! IN: upper triangle of A
! OUT: Choleski factor L in compressed storage:
!         lower triangle of L stored in lower triangle of A
!         main diagonal of L stored in Ld
!         original matrix still in upper triangle of A
!      iBad         = row at which (first) non-positive pivot encountered (n+1 if ok)
!      posDefinite  = true if A was numerically positive definite
!      logdet       = logarithm of determinant of "A" (optional)
!      condest      = fast estimate of condition number (using max/min diag**2)
!      err and message = status description.
! Algorithm:
! * Standard Choleski decomposition without pivoting. Lack of pivoting
!   does not induce any instabilities, yet pivoting can reduce roundoff
!   errors and make the solution invariant wrt variable indexing
! * The row-by-row matrix access pattern (following NR-90 and Dennis and Schnabel)
!   is actually slightly more complex than in the perturbed Choleski-Gershgorin
!   (revised modified robust Choleski) since it works with the entire matrix.
!   This makes pivoting _harder_ to code than in variants that acess only a
!   tringular portion of A.
! For more powerful modified Choleski decomposition see choles_dcmp_engine2/3/4.
! Ref: - Press et al, NR-90.
!      - Dennis and Schanbel,1996,"Numerical methods for unconstrained
!        optimization and nonlinear equations",p.319 (simplified RC algorithm)
use utilities_dmsl_kit,only:zero,two,assertEq,ns=>number_string
implicit none
! dummies
real(mrk),intent(inout)::A(:,:)
real(mrk),intent(out)::Ld(:)
logical(mlk),intent(out)::posDefinite
integer(mik),intent(out),optional::iBad
real(mrk),intent(out),optional::logDet,condest
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
character(*),parameter::procnam="choles_dcmp_engine1"
integer(mik)::n,j
real(mrk)::summ
logical(mlk)::conform
! Start procedure here
err=EXIT_SUCCESS; message=procnam//"/ok"; posDefinite=.true. ! presume innocence
call assertEq(size(A,1),size(A,2),size(Ld),conform,n)
if(.not.conform)then
  err=10; message="f-"//procnam//"/nonConform(A,Ld)"; return
elseif(n==0)then
  err=-10; message="w-"//procnam//"/dim[A]==0"; return
elseif(n<0)then
  err=20; message="f-"//procnam//"/dim[A]<0"; return
endif
if(present(iBad))iBad=n+1
do j=1,n
  summ=A(j,j)-dot_product(A(j,1:j-1),A(j,1:j-1))
  if(summ<=zero)then  ! standard Choleski iteration. In the robust Choleski scheme
    posDefinite=.false.; err=-20  ! zero is replaced with tolerance
    message="f-"//procnam//"/nonPosDefinite:"//ns(j)
    if(present(iBad))iBad=j
    if(present(logDet))logDet=undefRN
    if(present(condEst))condEst=undefRN
    return
  endif
  Ld(j)=sqrt(summ)
  A(j+1:n,j)=(A(j,j+1:n)-matmul(A(j+1:n,1:j-1),A(j,1:j-1)))/Ld(j)
enddo
call setLogDetCondEstCholes(pd=posDefinite,d=Ld,&
  logDet=logDet,condEst=condEst,err=err,message=message)
if(err/=EXIT_SUCCESS)then
  message="f-"//procnam//"/&"//message
endif
! End procedure here
endsubroutine choles_dcmp_engine1
!----------------------------------------------------
pure subroutine setLogDetCondEstCholes(pd,d,minD,logDet,condEst,err,message)
! Purpose: Sets the log(det) and condEst from the Choleski diagonal.
! Programmer: Dmitri Kavetski
! Kreated: The Golden Age circa 2003-2004. Or during P-time.
! History: 10 July 2013 AD, Hua Hin White Sands - make minD optional.
!          10 July 2013 AD, Hua Hin White Sands - implement 'pd'
use utilities_dmsl_kit,only:zero,two
implicit none
! dummies
logical(mlk),intent(in)::pd
real(mrk),intent(in)::d(:)
real(mrk),intent(in),optional::minD
real(mrk),intent(out),optional::logDet,condEst
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
character(*),parameter::procnam="setLogDetCondEstCholes"
logical(mlk)::haveDet,haveCond
real(mrk)::minDuse
! start procedure here
err=EXIT_SUCCESS; haveDet=present(logDet); haveCond=present(condEst)
if(haveDet.or.haveCond)then
  if(pd)then  ! check to avoid crashes
    if(present(minD))then; minDuse=minD
    else;                  minDuse=minval(d); endif
    if(minDuse<=zero)then
      message="f-"//procnam//"/unexpected[pd.but.minDuse<=0]"
      err=EXIT_FAILURE; return
    endif
  endif
endif
if(haveDet)then
  if(pd)then
    logDet=two*sum(log(d))
  else
    logDet=undefRN
  endif
endif
if(haveCond)then
  if(pd)then
    condEst=(maxval(d)/minval(d))**2
  else
    condEst=undefRN
  endif
endif
! end procedure here
endsubroutine setLogDetCondEstCholes
!----------------------------------------------------
pure subroutine choles_invrt1(A,Ld,err,message)
! Purpose: Inverts a positive definite matrix A given its Choleski
! factors stored in A and diagonal in Ld as follows:
! INPUT:  Choleski factors L:
!         Lower triangle of L stored in lower half of array A
!         Diagonal of L stored in array Ld.
! OUTPUT: Inverse matrix Ainv:
!         Upper triangle of Ainv returned in upper triangle of array A.
! COMMENT:
! 1. Returned matrix is not symmetric because it contains both the
!    (original) Choleski factor of A and its (computed) inverse Ainv.
use utilities_dmsl_kit,only:one,assertEq
implicit none
! dummy variables definitions
real(mrk),dimension(:,:),intent(inout)::A
real(mrk),dimension(:),intent(in)::Ld
integer(mik),intent(out)::err
character(*),intent(out)::message
! local variables definitions
integer(mik)::n,i,j
logical(mlk)::compat
! Start procedure here
call assertEq(size(A,1),size(A,2),size(Ld),compat,n)
if(.not.compat.or.n<1) then
  err=10; message="choles_invrt1/dimError"; return
endif ! Do customised forward/backward substitution
forall(j=1:n)A(j,j)=one/Ld(j) ! diag of inverse
do j=1, n ! loop over columns in the inverse
  do i=j+1,n  ! forward substitution adapted for matrix inversion
    A(j,i)=-dot_product(A(i,j:i-1),A(j,j:i-1))/Ld(i)
  enddo   ! note accounting for leading zero's in rhs (unit matrix)
  do i=n,j,-1 ! back sub, exploiting symmetry of A and its inverse
    A(j,i)=(A(j,i)-dot_product(A(i+1:n,i),A(j,i+1:n)))/Ld(i)
  enddo
enddo
! End procedure here
endsubroutine choles_invrt1
!----------------------------------------------------
pure subroutine choles_invrt2(Lmat,ainv,err,message)
! Usage: overloaded for different interface.
! This routine operates on a matrix already decomposed into its Choleski factor Lmat.
! The returned matrix Ainv is the complete inverse (symmetric).
! ---
! Comments:
! 1. If Lmat not present, then assumes ainv contains the
!    Choleski factor of the matrix to be inverted.
!    (routine accesses lower triangle only);
! 2. See argument list of 'choles_invrt1'.
use utilities_dmsl_kit,only:assertEq,getDiag
implicit none
! dummy variables definitions
real(mrk),dimension(:,:),optional,intent(in)::Lmat
real(mrk),dimension(:,:),intent(inout)::ainv
integer(mik),intent(out)::err
character(*),intent(out)::message
! local variables definitions
integer(mik)::i,j,n
logical(mlk)::compat
! Start procedure here
if(present(Lmat))then
  call assertEq(size(ainv,1),size(ainv,2),size(Lmat,1),size(Lmat,2),compat,n)
  if(.not.compat.or.n<1) then
    err=10; message="choles_invrt2/dimError"; return
  endif
  ainv=Lmat ! a bit wasteful, but very simple code!
  call choles_invrt1(ainv,getDiag(Lmat),err,message)
else
  call assertEq(size(ainv,1),size(ainv,2),compat,n)
  if(.not.compat.or.n<1) then
    err=20; message="choles_invrt2/dimError"; return
  endif ! again, very simple code!
  call choles_invrt1(ainv,getDiag(ainv),err,message)
endif
forall(i=1:n,j=1:n,i>j)ainv(i,j)=ainv(j,i)  ! make symmetric
! End procedure here
endsubroutine choles_invrt2
!----------------------------------------------------
pure subroutine choles_invrt_dcmp_engine(A,Ld,posDefinite,logDet,iBad,condest,&
  err,message)
! Purpose: Inverts a positive definite matrix A using Choleski scheme
! Programmer: Dmitri Kavetski
! Last modified: 30 August 2000 AD
! Performance
! IN: matrix A to be inverted (completely destroyed)
! OUT: 1. Inverse matrix stored in upper triangle of A
!      2. Choleski factors of A stored in lower half of A and diag Ld
!      3. (optional) logarithm of the determinant
!      4. (optional) estimated condition number of A
implicit none
! dummy variables definitions
real(mrk),dimension(:,:),intent(inout)::A
real(mrk),dimension(:),intent(out)::Ld
logical(mlk),intent(out)::posDefinite
real(mrk),intent(out),optional::logDet,condest
integer(mik),intent(out),optional::iBad
integer(mik),intent(out)::err
character(*),intent(out)::message
! Start procedure here
call choles_dcmp_engine1(a=A,Ld=Ld,posDefinite=posDefinite,&
  logDet=logDet,iBad=iBad,condest=condest,err=err,message=message) ! decomposition
if(.not.posDefinite.or.err/=0) return
call choles_invrt1(A,Ld,err,message) ! customised forward/backward substitution
! End procedure here
endsubroutine choles_invrt_dcmp_engine
!----------------------------------------------------
pure subroutine choles_dcmp_invrt2(a,ainv,Lmat,posDefinite,logDet,iBad,condest,&
  err,message)
! Purpose: Improved jacket for Choleski decomposition / inverter.
! Programmer: Dmitri Kavetski.
! Interface: (Use of keywords strongly recommended):
!   a or ainv: matrix to be inverted/decomposed using Choleski algorithm
! NB: a is optional, ainv is NOT optional
!   if a provided, inverts a and stores inverse in ainv;
!   if a not provided, inverts and overwrites ainv unto itself
!   Lmat: lower diagonal factor of A (optional by-product)
!   logDet: logarithm of determinant (optional by-product)
!   condest = condition number of the matrix (fast estimate)
!   err: status indicator, 0=ok, -1=a is not not positive definite
!   message: performance/result message
! Comments:
! Employs economic methods for array storage: regardless of run mode,
! only needs O(n) scratch space for the diagonal of the L-factor.
use utilities_dmsl_kit,only:zero,assertEq
implicit none
! dummy variables definitions
real(mrk),dimension(:,:),intent(in),optional::a
real(mrk),dimension(:,:),intent(inout)::ainv
real(mrk),dimension(:,:),intent(out),optional::Lmat
logical(mlk),intent(out)::posDefinite
real(mrk),intent(out),optional::logDet,condest
integer(mik),intent(out),optional::iBad
integer(mik),intent(out)::err
character(*),intent(out)::message
! local variables definitions
logical(mlk)::ok
real(mrk)::diag(size(ainv,dim=1))
integer(mik)::i,j,n
! Start procedure here
if(present(a))then  ! keep original matrix a, put inverse in ainv.
  call assertEq(size(a,1),size(a,2),size(ainv,1),size(ainv,2),ok,n)
  if(.not.ok.or.n<1) then
    err=10; message="f-choles_dcmp_invrt2/dimError"; return
  endif               ! get inverse and L factor in compressed form
  ainv(:,:)=a(:,:)    ! assume ainv has junk in it
else    ! do not keep original matrix, overwrite ainv
  call assertEq(size(ainv,dim=1),size(ainv,dim=2),ok,n)
  if(.not.ok.or.n<1) then  ! dimension error
    err=20; message="f-choles_dcmp_invrt2/dimError"; return
  endif               ! get inverse and L factor in compressed form
endif
call choles_invrt_dcmp_engine(ainv,diag,posDefinite,logDet,iBad,condest,err,message)
if(.not.posDefinite.or.err/=0) return
if(present(Lmat))then ! keep L factor in separate matrix
  call assertEq(n,size(Lmat,1),size(Lmat,2),ok,n)
  if(.not.ok.or.n<1) then
    err=30; message="f-choles_dcmp_invrt2/dimError"; return
  endif
  forall(i=1:n)           Lmat(i,i)=diag(i)   ! set diagonal
  forall(i=1:n,j=1:n,i>j) Lmat(i,j)=ainv(i,j) ! and upper triangle
  forall(i=1:n,j=1:n,i<j) Lmat(i,j)=zero      ! zero upper triangle
endif
forall(i=1:n,j=1:n,i>j)   ainv(i,j)=ainv(j,i) ! make Ainv symmetric
err=0; message="choles_dcmp_invrt2/ok" ! let's hope for the best
! End procedure here
endsubroutine choles_dcmp_invrt2
!----------------------------------------------------

!----------------------------------------------------
! 2. LU METHODS
!----------------------------------------------------
pure subroutine lu_dcmp1(a,indx,d,sing,logDet,signDet,condK,condest,err,message)
! Purpose: Given a matrix "a", this routine computes the LU decomposition
! of a row-wise permutation of "a"; uses partial pivoting and implicit scaling,
! "outer product" Gaussian reduction of Golub and van Loan (Press et al, 1996)
! IN:   a       = matrix to be decomposed
!       condK   = (optional) array index of norms for condition number (typically,1)
! OUT:  a       = LU factors overwritten into "a"
!       indx    = records the row permutation
!       d       = +/-1 if even/odd number of row permutations (for determinant)
!       sing    = true if matrix singular
!       logDet/signDet  = (optional) logarithm and sign of the determinant
!       condest = (optional) condition number estimate (corresp to condK)
!       err     = status: err=0->ok, +1->small pivot used
!       message = optional message
! Ref: adapted from Press et al. (1996), NR-90, NR-77.
! Comments
! * Uses implicit scaling in the pivoting strategy. This is not always optimal,
!   but does make the algorithm more or less invariant under rescaling of the
!   equations. Set impScaling=.false. to disable implicit scaling.
! * If matrix contains a row of zeroes the routine will quit with err=20 without
!   touching "a" and return huge condition number. If a perturbed LU decomposition
!   is still necessary (see next comment), the the calling unit should
!   perturb at least one of the row elements and call this routine again.
! * If matrix singular, replaces zero pivot with "tiny" scaled pivot, which
!   can be useful in inverse iteration polishing of eigenvalues. These cases
!   are flagged by sing=.true. but err=0. The condition number will be based
!   on the modified matrix. This strategy is not unlike modified Choleski factorizations.
use utilities_dmsl_kit,only:assertEq,imaxloc,swap,outerProd,zero,one,getKnorm
implicit none
! dummies
real(mrk),intent(inout)::a(:,:)
integer(mik),intent(out)::indx(:)
real(mrk),intent(out)::d
integer(mik),intent(in),optional::condK(:)
real(mrk),optional,intent(out)::signDet,logDet
real(mrk),intent(out),optional::condEst(:)
logical(mlk),intent(out)::sing
integer(mik),intent(out)::err
character(*),intent(out)::message
! local parameters
real(mrk),parameter::tiny=1.e-20_mrk ! small for singular pivot replacing
logical(mlk),parameter::impScaling=.true. ! requests implicit scaling
! locals
real(mrk),dimension(size(indx))::vv  ! implicit scaling of each row
logical(mlk)::conform,jsing
integer(mik)::j,n,imax
character(1)::jmsg
real(mrk),allocatable::normA(:)
! Start procedure here
err=0; message="lu_dcmp1/ok"; sing=.false. ! presume innocence
d=one  ! no row interchanges yet
call assertEq(size(a,dim=1),size(a,dim=2),size(indx),conform,n)
if(.not.conform.or.n<1) then
  err=-10; message="f-lu_dcmp1/dimError(a,indx)"; return
endif
if(present(condest))then
  allocate(normA(size(condest)))  ! compute matrix norms for condition estimation
  forall(j=1:size(condEst))normA(j)=getKnorm(a,condK(j))
endif
vv(:)=maxval(abs(a),dim=2) ! get implicit scaling in each row
if(any(vv==zero)) then     ! there is a row of zeros in the matrix
  err=-20; sing=.true.; message="w-lu_dcmp1/singular:row=0"
  if(present(condest))condest=hugeRe
  return
endif
if(impScaling)then    ! - enable implicit scaling
  vv(:)=one/vv(:)     !   save the scaling to pivot based on "implicit scaling"
else                  ! - disable implicit scaling
  vv=one
endif
do j=1, n
  imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j))) ! find the pivot row
  if(j/=imax) then  ! yes, interchange rows
    call swap(a(imax,:),a(j,:))
    d=-d            ! ... and change the parity of "d"
    vv(imax)=vv(j)  ! also interchange the scale factor
  endif
  indx(j)=imax
  if(a(j,j)==zero) then ! matrix is (algorithmically) singular,however,
! set pivot to a small number (useful in some singular matrix applications)
! (Reference: Press et al.1992). Modified for improved scale-invariance by DK.
! this strategy is not unlike modified Choleski factorizations
    sing=.true.
    a(j,j)=tiny*maxval(abs(a(j+1:n,j)))         ! first scale based on column values
    if(a(j,j)==zero)a(j,j)=tiny*maxval(abs(a))  ! emergency scale based on entire matrix.
    err=1; message="w-lu_dcmp1/singular/usedSmallPivot"
  endif
  a(j+1:n,j)=a(j+1:n,j)/a(j,j) ! normalise column and reduce remaining submatrix
  a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
enddo
if(present(logDet).and.present(signDet)) then ! get log(determinant)
  call getlogDet_trimat(a,logDet,signDet,jsing,err,jmsg)
  signDet=d*signDet   ! accounting for row permutations
endif
! Estimate condition number if requested
if(present(condest))then
  call getCondest(normA=normA,lu=a,indx=indx,knorm=condK,&
                  condest=condest,err=err,message=message)
  if(err/=0)then
    err=30;message="f-lu_dcmp1/&"//message
  else
    message="lu_dcmp1/ok"
  endif
  deallocate(normA)
endif
! End procedure here
endsubroutine lu_dcmp1
!----------------------------------------------------
pure subroutine lu_fwbw1(a,indx,b,ynorm1,err,message)
! Purpose: solves the linear system a.x=b, with solution stored in "b"
! "a" is input as returned by lu_dcmp. "indx" is the permutation vector
! "a" and "indx" are not modified and can be reused with another "b"
! IN:  a      = LU decomposed matrix
!      b      = rhs vector
! OUT: b      = solution vector (overwrites b)
!      ynorm1 = ||y||1 where y is the intermediate vector.
!      err    = status: err=0->ok
!      message= optional message
! Comment: the algorithm anticipates that rhs may begin with zero's.
!          it is hence efficient for use as a matrix inverter.
use utilities_dmsl_kit,only:zero,assertEq
implicit none
! dummy variables definitions
real(mrk),intent(in)::a(:,:)
integer(mik),intent(in)::indx(:)
real(mrk),intent(inout)::b(:)
real(mrk),intent(out),optional::ynorm1
integer(mik),intent(out)::err
character(*),intent(out)::message
! local variables definitions
real(mrk)::summ
integer(mik)::i,ii,n,ll
logical(mlk)::conform
! Start procedure here
err=0; message="lu_fwbw1/ok"! presume innocence
call assertEq(size(a,dim=1),size(a,dim=2),size(indx),size(b),conform,n)
if(.not.conform.or.n<1) then ! dimension error
  err=-10; message="f-lu_fwbw1/dimError(a,indx)"
  b(:)=-hugeRE; return
endif
ii=0  ! Forward substitution, get first non-vanishing element of the rhs
do i=1, n     ! when ii is set to a +ve number, it will be the index
  ll=indx(i)  ! of the first non-vanishing element of the rhs
  summ=b(ll)  ! unscramble permutation
  b(ll)=b(i)
  if(ii/=0) then ! do full dot product
    summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
  elseif(summ/=zero) then ! a nonzero entry in b detected
    ii=i  ! will have to do the full dot product next time
  endif
  b(i)=summ
enddo
if(present(ynorm1))ynorm1=sum(abs(b))
do i=n, 1, -1 ! Back-substitution
  b(i)=(b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
enddo
! End procedure here
endsubroutine lu_fwbw1
!----------------------------------------------------
pure subroutine lu_invrt1(a,ainv,lu,indx,sing,logDet,signDet,condK,condest,err,message)
! Purpose: Inverts matrix "ainv" or "a" using LU decomposition
! STRONGLY RECOMMENDED: use argument keywords!!!
! NB:   a is optional, ainv is NOT optional
! USE:  provide "a" and "lu"=>orig. "a", invrs "ainv", LU's "lu"
!       provide "a",not"lu"=>orig "a", invrs "ainv". LU's discarded
!       not"a",provide"lu"=>inverts/overwrites "ainv", LU's "lu"
!       not"a",not"lu"=>inverts/overwrites "ainv", LU's discarded
!       indx=row permutation vector (must get it if LU factors needed)
!       logDet=logarithm of the determinant (optional)
!       signDet=sign of the determinant (optional)
!       condest = condition number estimate (condK=array index of norms)
!       err=status: err=0->ok
!       message=optional message
! NB: must provide array "indx" if LU factors requested (present "lu")
!     if determinant requested, must provide both logDet and signDet.
! Comments: an advantage of this code is that a is intent(in).
use utilities_dmsl_kit,only:assertEq,assertEqLog,unitMatrix
implicit none
! dummy variables definitions
real(mrk),intent(in),optional::a(:,:)
real(mrk),intent(inout)::ainv(:,:)
real(mrk),intent(out),optional::lu(:,:)
logical(mlk),intent(out)::sing
integer(mik),optional,intent(out)::indx(:)
integer(mik),intent(out)::err
character(*),intent(out)::message
real(mrk),optional,intent(out)::logDet,signDet
integer(mik),intent(in),optional::condK(:)
real(mrk),intent(out),optional::condEst(:)
! local variables definitions
integer(mik)::i,n
real(mrk)::d
integer(mik),allocatable::indx_temp(:)
real(mrk),allocatable::a_temp(:,:)
logical(mlk)::conform
character(100)::msg
! Start procedure here
call assertEq(size(ainv,dim=1),size(ainv,dim=2),conform,n)
if(.not.conform) then
  err=-10; message="f-lu_invrt1/dimError(ainv)"; return
elseif(n<=0) then
  err=-11; message="f-lu_invrt1/bug/dim<=0"; return
endif
! branches depending on requested tasks
if(.not.present(a).and..not.present(lu)) then  ! overwrite ainv with inverse
  allocate(a_temp(n,n),indx_temp(n),stat=err)  ! need temporary LU storage
  if(err/=0) then
    message="f-lu_invrt1/allocError(a_temp,indx_temp)";return
  endif
  a_temp(:,:)=ainv(:,:) ! copy ainv->a_temp and overwrite a_temp with LU factors
  call lu_dcmp(a=a_temp,indx=indx_temp,d=d,sing=sing,logdet=logDet,signdet=signDet,&
    condK=condK,condest=condest,err=err,message=msg)
  if(sing) then
    err=0;message="w-lu_invrt1/singularMatrix"; return
  elseif(err/=0) then  ! oops, error return
    message="f-lu_invrt1/:"//msg; return
  endif
  call unitmatrix(ainv(:,:)) ! initialise ainv as a unit matrix
  do i=1, n
    call lu_fwbw(a_temp,indx_temp,ainv(:,i),err=err,message=msg) ! inverse by columns
  enddo
  deallocate(a_temp,indx_temp,stat=err)  ! garbage collection
  if(err/=0) then
    message="f-lu_invrt1/deallocError(a_temp,indx_temp)";return
  endif
elseif(.not.present(a).and.present(lu)) then ! inverse/LUfactors
  if(.not.present(indx)) then
    err=-10; message="f-lu_invrt1/need indx"; return
  elseif(.not.assertEqLog(size(lu,1),size(lu,2),n)) then
    err=-10; message="f-lu_invrt1/size(lu)error"; return
  endif
  lu(:,:)=ainv(:,:)
  call lu_dcmp(a=lu,indx=indx,d=d,sing=sing,logdet=logDet,signdet=signDet,&
    condK=condK,condest=condest,err=err,message=msg)
  if(sing) then
    err=0;message="w-lu_invrt1/singularMatrix"; return
  elseif(err/=0) then  ! oops, error return
    message="f-lu_invrt1/:"//msg; return
  endif
  call unitmatrix(ainv(:,:)) ! initialise ainv as a unit matrix
  do i=1, n
    call lu_fwbw(lu,indx,ainv(:,i),err=err,message=msg) ! inverse by columns
  enddo
elseif(present(a).and.present(lu)) then  ! keep "a", get LU
  if(.not.present(indx)) then
    err=-10; message="f-lu_invrt1/missingIndx"; return
  elseif(.not.assertEqLog(size(lu,1),size(lu,2),size(a,1),size(a,2),n)) then
    err=-10; message="f-lu_invrt1/nonConform(lu,a)"; return
  endif
  lu(:,:)=a(:,:)  ! lu will contain the LU factors of a
  call lu_dcmp(a=lu,indx=indx,d=d,sing=sing,logdet=logDet,signdet=signDet,&
    condK=condK,condest=condest,err=err,message=msg)
  if(sing) then
    err=0;message="w-lu_invrt1/singularMatrix"; return
  elseif(err/=0) then  ! oops, error return
    message="f-lu_invrt1/:"//msg; return
  endif
  call unitmatrix(ainv(:,:)) ! initialise ainv as a unit matrix
  do i=1, n
    call lu_fwbw(lu,indx,ainv(:,i),err=err,message=msg) ! construct inverse by columns
  enddo
elseif(present(a).and..not.present(lu) ) then  ! keep "a", no LU
  if(.not.assertEqLog(size(a,1),size(a,2),n)) then
    err=-10; message="f-lu_invrt1/nonConform(a,ainv)"; return
  endif
  allocate(a_temp(n,n),indx_temp(n),stat=err) ! need temporary storage for LU's
  if(err/=0) then
    err=-10; message="f-lu_invrt1/allocateError(a_temp,idx_temp)"; return
  endif
  a_temp(:,:)=a(:,:)  ! initialise a_temp.
  call lu_dcmp(a=a_temp,indx=indx_temp,d=d,sing=sing,logdet=logDet,signdet=signDet,&
    condK=condK,condest=condest,err=err,message=msg)
  if(sing) then
    err=0;message="w-lu_invrt1/singularMatrix"; return
  elseif(err/=0) then  ! oops, error return
    message="f-lu_invrt1/"//msg; return
  endif
  call unitmatrix(ainv(:,:)) ! initialise ainv as a unit matrix
  do i=1, n   ! construct inverse by columns
    call lu_fwbw(a_temp,indx_temp,ainv(:,i),err=err,message=msg)
  enddo
  deallocate(a_temp,indx_temp,stat=err)  ! garbage collection
  if(err/=0) then
    message="f-lu_invrt1/deallocError(a_temp,indx_temp)";return
  endif
endif
err=0; message="lu_invrt1/ok" ! ughh...
! End procedure here
endsubroutine lu_invrt1
!----------------------------------------------------
pure subroutine getCondest1(a,knorm,sing,condest,logdet,signdet,err,message)
! Purpose: Estimates the condition number of matrix "a" given
! its LU decomposition (stored in lu). Uses k-norm of matrices.
! NB: the resulting estimate is an underestimate of the actual condition number.
! Programmer: Dmitri Kavetski
! Last modified:8March2003
! Ref: Golub and van Loan(1983) "Matrix computations",p.78.
use utilities_dmsl_kit,only:getKnorm
implicit none
! dummies
real(mrk),intent(in)::a(:,:)
integer(mik),intent(in)::knorm(:)
logical(mlk),intent(out)::sing
real(mrk),intent(out)::condest(:)
real(mrk),intent(out),optional::logdet,signdet
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
integer(mik)::i,indx(size(a,1))
real(mrk)::d,normA(size(knorm)),lu(size(a,1),size(a,1))
! Start procedure here
forall(i=1:size(knorm))normA(i)=getKnorm(a,knorm(i))
LU=a  ! make a scratch copy of a to LU-destroy
call lu_dcmp(a=LU,indx=indx,d=d,sing=sing,logDet=logdet,signdet=signDet,&
             err=err,message=message)
if(sing)then
  err=-10;message="f-getCondest1/&"//message
  condest=hugeRe; return
elseif(err/=0)then
  err=10;message="f-getCondest1/&"//message
  condest=hugeRe; return
endif
call getCondest(normA=normA,lu=LU,indx=indx,knorm=knorm,&
                condest=condest,err=err,message=message)
if(err/=0)then
  err=20;message="f-getCondest1/&"//message
  condest=hugeRe; return
endif
! End procedure here
endsubroutine getCondest1
!----------------------------------------------------
pure subroutine getCondest2(normA,lu,indx,knorm,condest,err,message)
! Purpose: Estimates the condition number of matrix "a" given
! its LU decomposition (stored in lu) and its matrix norms normA.
! NB: the resulting estimate is an underestimate of the actual condition number,
! computed from an estimated maximum amplification of a vector induced
! by the matrix. Uses k-norm of matrices.
! Algorithm:
! (i)   Solve a(t).r=d with special d to induce large vector growth
! (ii)  Solve a.y=r  (in the code y and r stored in y)
! (iii) Estimate ||ainv||=||v||/||y||
! Programmer: Dmitri Kavetski
! Last modified: 1 April 2004 (fixed up pivoting in step ib)
! Ref:
! - Golub and van Loan (1983) "Matrix computations",p.78.
! - Kahaner, D., C. Moler, and S. Nash (1989) "Numerical Methods and Software".
! - Dennis and Schnabel (1996) "Numerical methods for unconstrained optimization
!   and nonlinear equations". algorithm 3.3.1 (condest)
use utilities_dmsl_kit,only:zero,one,assertEq,getKnorm,lower_rsolv,swap
implicit none
! dummies
real(mrk),intent(in)::normA(:),lu(:,:)
integer(mik),intent(in)::indx(:),knorm(:)
real(mrk),intent(out)::condest(:)
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
integer(mik)::i,n,nk
real(mrk)::y(size(lu,1)),normR(size(knorm)),normY(size(knorm)) !,normT(size(knorm))
logical(mlk)::ok
logical(mlk),parameter::useGenericLNS=.false.
! Start procedure here
call assertEq(size(lu,1),size(lu,2),size(indx),ok,n)
if(.not.ok)then
  err=10;message="f-getCondest2/dimError(lu,indx)"
  condest=hugeRe; return
endif
call assertEq(size(normA),size(knorm),size(condest),ok,nk)
if(ok)then
  err=0;message="getCondest2/ok"
else
  err=20;message="f-getCondest2/inError/dimKnorm/=dimCondest"
  condest=hugeRe; return
endif
! * Obtain a large norm solution to U(t).w=d, with special "d" (step ia)
call largeNormSolve(T=LU,y=y,triang="Ut",err=err,message=message)
if(err/=0)then
  err=30;message="f-getCondest2/st1/&"//message
  condest=hugeRe; return
endif
!forall(i=1:nk)normT(i)=getKnorm(y,knorm(i))   ! norm of large-norm solution (not needed)
! * Solve L(t).r=w (store r in y) (step ib)
call lower_rsolv(a=LU,unitD=.true.,x=y,transp=.true.,err=err)
if(err/=0)then
  err=40;message="f-getCondest2/st2/&"//message
  condest=hugeRe; return
endif
do i=n,1,-1 ! unscramble pivoting permutations
  if(indx(i)/=i)call swap(y(indx(i)),y(i))
enddo
forall(i=1:nk)normR(i)=getKnorm(y,knorm(i))   ! compute and save vector norms
! * Solve L.U.z = P.r using standard LU forward/backward sub (step ii).
call lu_fwbw(a=LU,indx=indx,b=y,err=err,message=message)
if(err/=0)then
  err=50;message="f-getCondest2/st3/&"//message
  condest=hugeRe; return
endif
! * Estimate condition number using requested norms (step iii)
forall(i=1:nk)normY(i)=getKnorm(y,knorm(i))
forall(i=1:nk,normR(i)/=zero)condest(i)=normA(i)*normY(i)/normR(i)
forall(i=1:nk,normR(i)==zero)condest(i)=zero    ! safeguard zero division
! End procedure here
endsubroutine getCondest2
!----------------------------------------------------
pure subroutine getCondest3(a,adiag,triang,transpStore,condest,err,message)
! Purpose: Estimates the L1 condition number of a triangular matrix "a"
! triang denotes the type of matrix: "U","Ut".
! condest is an underestimate (ie, lower bound) of cond[A], computed
! from the estimated maximum amplification induced by a on a vector.
! The lower triangle is not accessed.
! Optionally, transpStore specifies that matrix transposed in lower triangle.
! Programmer: Dmitri Kavetski
! Last modified: 1 April 2004
! Comments:
! - Assumes a is not permuted. Hence the routine may yield incorrect
!   results if a is a triangular factor of a permuted general matrix.
!   (eg, pivoted Choleski factorization).
! Ref: Dennis and Schnabel (1996) Numerical methods for unconstrained
! optimization and nonlinear equations. Algorithm A.3.3.1. (p.309-311)
use utilities_dmsl_kit,only:zero,quickif,assertEq,upper_rsolv,lower_rsolv
implicit none
! dummies
real(mrk),intent(in)::a(:,:)
real(mrk),intent(in),optional::adiag(:)
character(*),intent(in)::triang
logical(mlk),intent(in),optional::transpStore
real(mrk),intent(out)::condest
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
integer(mik)::j,n
real(mrk)::normA,normR,normY,y(size(a,1))
character(2)::triang0
logical(mlk)::ok,transpStore0
logical(mlk),parameter::transpStoreDef=.false.
! Start procedure here
call assertEq(size(a,1),size(a,2),ok,n)
if(ok)then
  err=0;message="getCondest3/ok"
elseif(.not.ok)then
  err=10;message="f-getCondest3/dimError(a)"
  condest=hugeRe; return
endif
if(present(adiag))then
  if(size(adiag)/=n)then
    err=20;message="f-getCondest3/dimError(adiag)"
    condest=hugeRe; return
  endif
endif
transpStore0=quickif(transpStore,transpStoreDef)
! * Compute 1-norm of a, ||a||1
if(transpStore0)then
  if(present(adiag))then
    forall(j=1:n)y(j)=sum(abs(a(j,1:j-1)))+abs(adiag(j))
  else
    forall(j=1:n)y(j)=sum(abs(a(j,1:j)))
  endif
else
  if(present(adiag))then
    forall(j=1:n)y(j)=sum(abs(a(1:j-1,j)))+abs(adiag(j))
  else
    forall(j=1:n)y(j)=sum(abs(a(1:j,j)))
  endif
endif
normA=maxval(y) ! norm of matrix A, ||A||1
! * Obtain a large-norm solution to a(t).r=d, with special "d"
selectcase(triang)
case("u","U")               ! upper triangular
  triang0="Ut"
case("Ut","uT","ut","UT")   ! transposed U
  triang0="U"
case default
  err=100;message="f-getCondest3/unknownTriang"
  condest=hugeRe; return
endselect
call largeNormSolve(T=a,Tdiag=adiag,y=y,triang=triang0,&
                    transpStore=transpStore,err=err,message=message)
if(err/=0)then
  err=30;message="f-getCondest3/st1/&"//message
  condest=hugeRe; return
endif
normR=sum(abs(y))           ! compute ||r||1
! * Solve a.r=y (store r in y)
if(transpStore0)then
  call lower_rsolv(a=a,d=adiag,x=y,transp=.not.(triang0=="U"),err=err)
else
  call upper_rsolv(a=a,d=adiag,x=y,transp=.not.(triang0=="U"),err=err)
endif
if(err/=0)then
  err=40;message="f-getCondest3/st2/&"//message
  condest=hugeRe; return
endif
! * Estimate condition number from vector amplification factor
normY=sum(abs(y))           ! compute ||y||1
if(normR==zero)then ! prevent division by zero
  condest=zero
else
  condest=normA*normY/normR
endif
! End procedure here
endsubroutine getCondest3
!----------------------------------------------------
pure subroutine getlogDet_trimat(a,logDet,signDet,sing,err,message)
! auxiliary subroutine to compute the determinant of a triangular matrix
! in: matrix "a", out: log and sign of determinant
use utilities_dmsl_kit,only:assertEq,getDiag,zero,one
implicit none
! Dummies
real(mrk),intent(in)::a(:,:)
real(mrk),intent(out)::logDet,signDet
logical(mlk),intent(out)::sing
integer(mik),intent(out)::err
character(*),intent(out)::message
! Locals
integer(mik)::ich,n
logical(mlk)::ok
real(mrk),dimension(size(a,dim=1))::diag
! Start procedure here
call assertEq(size(a,dim=1),size(a,dim=2),ok,n)
if(.not.ok) then
  err=-10; message="f-getlogDet_trimat/nonSquareMat"; return
else
  err=0; message="getlogDet_trimat/ok"
endif                     ! get number of negative diagonal elements
diag=getDiag(a); ich=count(diag<zero)
diag=abs(diag); sing=any(diag==zero) ! singulars will have zero diagonal(s)
if(.not.sing) then        ! non singular matrix
  logDet=sum(log(diag))
  signdet=merge(one,-one,modulo(ich,2)==0) ! 1 if even, -1 if odd
endif
! End procedure here
endsubroutine getlogDet_trimat
!----------------------------------------------------
pure subroutine largeNormSolve(T,Tdiag,w,y,triang,transpStore,err,message)
! Purpose: given T(n,n) non-singular and upper-triangular, and
! set of weights, computes vector y so that ||y||inf~||Tinv||inf.
! Another words, computes large norm solution to U.y=d, choosing
! d(i)={+1,-1} so that y is large.
! triang ("U" or "Ut") specifies whether to solve "U.y=d" or "U(t).y=d"
! Optionally, diagonal of T may be stored in Tdiag
! Optionally, transpStore allows storing matrices in lower triangle of T.
! If weight "w" not provided, assumes w=1/|diag[T]| (Golub & van Loan,1983)
! Programmer: Dmitri Kavetski, 1 April 2004.
! Typical use:
! * Estimate condition number of triangular matrices
! * Estimate condition number of matrices from LU factorisations.
! Ref:
! * Golub and van Loan(1983)"Matrix Computations", Algorithm 4.5-1, p.77.
! * Dennis and Schnabel (1996) Numerical methods for unconstrained
!   optimization and nonlinear equations. Algorithm A.3.3.1., p.309-311.
use utilities_dmsl_kit,only:quickif
implicit none
! dummies
real(mrk),intent(in)::T(:,:)
real(mrk),intent(in),optional::w(:),Tdiag(:)
character(*),intent(in)::triang
logical(mlk),intent(in),optional::transpStore
real(mrk),intent(out)::y(:)
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
logical(mlk),parameter::transpStoreDef=.false.
logical(mlk)::transpStore0
! Start procedure here
transpStore0=quickif(transpStore,transpStoreDef)
selectcase(triang)
case("U","u")
  call largeNormSolve_u (T,Tdiag,transpStore0,w,y,err,message)
case("Ut","uT","ut","UT")
  call largeNormSolve_ut(T,Tdiag,transpStore0,w,y,err,message)
case default
  err=+10;message="f-largeNormSolve/unknownTriang:"//triang
endselect
! End procedure here
endsubroutine largeNormSolve
!----------------------------------------------------
pure subroutine largeNormSolve_u(T,Tdiag,transpStore,w,y,err,message)
! Purpose: Large-norm solution of U.y=d, with d={+/-1}.
! Comments:
! * No effort taken to rescale vectors, hence overflow may occur in some cases.
! See Kahaner et al. routine "SGECO" for implementation with frequent
! rescaling. For an input matrix of reasonable size, overflow would
! indicate extreme ill-conditioning.
use utilities_dmsl_kit,only:assertEq,one,zero,getdiag
implicit none
! dummies
real(mrk),intent(in)::T(:,:)
real(mrk),intent(in),optional::w(:),Tdiag(:)
logical(mlk),intent(in)::transpStore
real(mrk),intent(out)::y(:)
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
integer(mik)::k,n
real(mrk)::p(size(T,1)),pm(size(T,1)),ykp,ykm,sp,sm
real(mrk),allocatable,dimension(:)::lw
logical(mlk)::ok,useW,useDiag
! Start procedure here
useW=present(w); n=size(T,1)
useDiag=present(Tdiag)
if(useDiag)then;if(size(Tdiag)/=n)then
  err=20;message="f-largeNormSolve_u/dimError(Tdiag)"; return
endif;endif
if(useW)then
  call assertEq(size(T,1),size(T,2),size(w),size(y),ok,n)
else
  call assertEq(size(T,1),size(T,2),size(y),ok,n)
  allocate(lw(size(T,1)))
  if(present(Tdiag))then; lw(:)=one/abs(Tdiag)
  else;                   lw(:)=one/abs(getdiag(T)); endif
endif
if(ok)then; err=0;message="largeNormSolve_u/ok"
else;       err=10;message="f-largeNormSolve_u/dimError"; endif
if(useDiag)then; y(n)=one/Tdiag(n)
else;            y(n)=one/T(n,n); endif
k=n
if(transpStore)then         ! first iteration k=1 is fairly predictable
  p(1:k-1)=T(k,1:k-1)*y(k)
else
  p(1:k-1)=T(1:k-1,k)*y(k)  ! choose e(j)=+1 (arbitrary initial condition)
endif
do k = n-1, 1, -1
!p(:)=zero  ! this is the first iteration written in standard form
!do k = n, 1, -1
  if(useDiag)then
    ykp=(+one-p(k))/Tdiag(k); ykm=(-one-p(k))/Tdiag(k)
  else
    ykp=(+one-p(k))/T(k,k);   ykm=(-one-p(k))/T(k,k)
  endif
  if(transpStore)then         ! matrix stored in lower triangle
    pm(1:k-1)=p(1:k-1)+T(k,1:k-1)*ykm; p(1:k-1)=p(1:k-1)+T(k,1:k-1)*ykp
  else
    pm(1:k-1)=p(1:k-1)+T(1:k-1,k)*ykm; p(1:k-1)=p(1:k-1)+T(1:k-1,k)*ykp
  endif
  if(useW)then                ! use user-supplied weights
    sp=abs(ykp)+dot_product( w(1:k-1),abs(p (1:k-1)))
    sm=abs(ykm)+dot_product( w(1:k-1),abs(pm(1:k-1)))
  else
    sp=abs(ykp)+dot_product(lw(1:k-1),abs(p (1:k-1)))
    sm=abs(ykm)+dot_product(lw(1:k-1),abs(pm(1:k-1)))
  endif
  y(k)=merge(ykp,ykm,sp>=sm)
  if(sp<sm)p(1:k-1)=pm(1:k-1) ! e(j)=-1, else e(j)=+1
enddo
if(.not.useW)deallocate(lw)
! End procedure here
endsubroutine largeNormSolve_u
!----------------------------------------------------
pure subroutine largeNormSolve_ut(T,Tdiag,transpStore,w,y,err,message)
! Purpose: Large-norm solution of U(t).y=d, with d={+/-1}
use utilities_dmsl_kit,only:assertEq,one,zero,getdiag
implicit none
! dummies
real(mrk),intent(in)::T(:,:)
real(mrk),intent(in),optional::w(:),Tdiag(:)
logical(mlk),intent(in)::transpStore
real(mrk),intent(out)::y(:)
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
integer(mik)::k,n
real(mrk)::p(size(T,1)),pm(size(T,1)),ykp,ykm,sp,sm
real(mrk),allocatable,dimension(:)::lw
logical(mlk)::ok,useW,useDiag
! Start procedure here
useW=present(w); n=size(T,1)
useDiag=present(Tdiag)
if(useDiag)then;if(size(Tdiag)/=n)then
  err=20;message="f-largeNormSolve_ut/dimErrorTdiag"; return
endif;endif
if(useW)then
  call assertEq(size(T,1),size(T,2),size(w),size(y),ok,n)
else
  call assertEq(size(T,1),size(T,2),size(y),ok,n)
  allocate(lw(size(T,1)))
  if(present(Tdiag))then; lw(:)=one/abs(Tdiag)
  else;                   lw(:)=one/abs(getdiag(T)); endif
endif
if(ok)then; err=0;message="largeNormSolve_ut/ok"
else;       err=10;message="f-largeNormSolve_ut/dimError";endif
if(useDiag)then;    y(1)=one/Tdiag(1)
else;               y(1)=one/T(1,1);     endif
if(transpStore)then;p(2:n)=T(2:n,1)*y(1)
else;               p(2:n)=T(1,2:n)*y(1);endif
do k=2,n
!p(:)=zero
!do k = 1, n
  if(useDiag)then
    ykp=(+one-p(k))/Tdiag(k); ykm=(-one-p(k))/Tdiag(k)
  else
    ykp=(+one-p(k))/T(k,k);   ykm=(-one-p(k))/T(k,k)
  endif
  if(transpStore)then
    pm(k+1:n)=p(k+1:n)+T(k+1:n,k)*ykm; p(k+1:n)=p(k+1:n)+T(k+1:n,k)*ykp
  else
    pm(k+1:n)=p(k+1:n)+T(k,k+1:n)*ykm; p(k+1:n)=p(k+1:n)+T(k,k+1:n)*ykp
  endif
  if(useW)then
    sp=abs(ykp)+dot_product( w(k+1:n),abs(p (k+1:n)))
    sm=abs(ykm)+dot_product( w(k+1:n),abs(pm(k+1:n)))
  else
    sp=abs(ykp)+dot_product(lw(k+1:n),abs(p (k+1:n)))
    sm=abs(ykm)+dot_product(lw(k+1:n),abs(pm(k+1:n)))
  endif
  y(k)=merge(ykp,ykm,sp>=sm)
  if(sp<sm)p(k+1:n)=pm(k+1:n)
enddo
if(.not.useW)deallocate(lw)
! End procedure here
endsubroutine largeNormSolve_ut
!----------------------------------------------------

end module linalg_dmsl_kit
