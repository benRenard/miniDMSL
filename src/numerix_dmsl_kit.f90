module numerix_dmsl_kit

!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! WARNING !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!
! Different from original DMSL file: see changes in subroutine uniran_s


use kinds_dmsl_kit    ! numeric kind definitions
implicit none
private

! Publicly-available stuff
public::&
    ! stats
    getmean,getvar,getCV,getmoments,getmeanvar,&
    ! probability distributions
    normal_logp,exp_logp,invGamma_logp,invChi2_logp,beta_logp,&
    poisson_logPmf,binomial_Pmf,negBinomial_pmf,&
    normal_cdf,exp_cdf,&
    poisson_cdf,binomial_cdf,&
    normal_icdf,exp_icdf,&
    uniran,normaldev,invGamdev,gamdev,invChi2dev,expdev,&
    poidev,bnldev,negbindev,&
    ! sorting
    quicksort,indexx_qsort,&
    ! other
    znewton

! Overloads
!--
interface getmean
  module procedure getmean_s
endinterface getmean
!--
interface getvar
  module procedure getvar_s,getcovar_v2
endinterface getvar
!--
interface getCV
  module procedure getCV_s
endinterface getCV
!--
interface getmoments
  module procedure getmoments_s1,getmoments_s2
endinterface getmoments
!--
interface getmeanvar
  module procedure getmeanvar_s,getmeancovar_v,getmeanvar_v
endinterface getmeanvar
!--
interface quicksort
  module procedure quicksort1_rv,quicksort1_iv
endinterface quicksort
!--
interface indexx_qsort
  module procedure indexx_qsr,indexx_qsi
endinterface indexx_qsort
!--
interface normal_cdf
  module procedure normal_cdf_2
endinterface normal_cdf
!--
interface normal_icdf
  module procedure ppnd16_func,normal_icdf_v
endinterface normal_icdf
!--
interface uniran
  module procedure uniran_s,uniran_v,uniranHyper_s,uniranHyper_v,&
    uniran_is,uniran_iv,uniranHyper_is,uniranHyper_iv
endinterface uniran
!--
interface normaldev
  module procedure normaldev1_s,normaldev1_1m1v_vd,normaldev1_NmNv_vd,&
    normaldev_v,normaldev2_v
endinterface normaldev
!--
interface gamdev
  module procedure gamdevint_s,gamdevr_s
endinterface gamdev
!--
interface invGamdev
  module procedure invGamdev_s
endinterface invGamdev
!--
interface invChi2dev
  module procedure invChi2dev_s,invscChi2dev_s
endinterface invChi2dev
!--
interface expdev
  module procedure expdev_s,expdev_v,expdevb_s,expdevb_v
endinterface expdev
!--
interface poidev
  module procedure poidev2_s
endinterface poidev
!--
interface bnldev
  module procedure bnldev2_s
endinterface bnldev
!--
interface negbindev
  module procedure negbindev_s
endinterface negbindev
!--

contains

!----------------------------------------------------------------------
pure function getmean_s(x)
! Purpose: numerical approx. of the expectation of a sample of a
! a random scalar variable
! IN: vector x(1:ns) - note packaging
! OUT: numexpectation
implicit none
! Dummies
real(mrk),intent(in)::x(:)
real(mrk)::getmean_s
! local variables definitions
real(mrk)::nsr
! Start procedure here
nsr=real(size(x),mrk); getmean_s=sum(x)/nsr ! arithmetic average
! End procedure here
endfunction getmean_s
!----------------------------------------------------------------------
pure function getvar_s(x,method)
! Purpose: reformats output from the library function getmeanvar to
! deliver the variance
implicit none
real(mrk),intent(in)::x(:)
character(*),intent(in)::method
real(mrk)::getvar_s
! Locals
real(mrk)::junkmean
integer(mik)::jerr
character(1)::jmsg
! Start procedure here
call getmeanvar(x,junkmean,getvar_s,method,jerr,jmsg)
!End procedure here
endfunction getvar_s
!----------------------------------------------------
pure function getCV_s(x,method,setZero,smallMean)
! Purpose: returns the coefficient of variation of x.
! setZero allows to set the CV to zero if mean<smallmean
use utilities_dmsl_kit,only:zero,quickif
implicit none
real(mrk),intent(in)::x(:)
character(*),intent(in)::method
logical(mlk),intent(in),optional::setZero
real(mrk),intent(in),optional::smallMean
real(mrk)::getCV_s
! Locals
real(mrk)::junkmean,junkvar,smallmeana
integer(mik)::jerr
character(1)::jmsg
real(mrk),parameter::smallmeandef=1.e-8_mrk
! Start procedure here
call getmeanvar(x,junkmean,junkvar,method,jerr,jmsg)
junkmean=abs(junkmean)
smallmeana=quickif(smallmean,smallmeandef)
if(junkmean>=smallmeana)then
  getCV_s=sqrt(junkvar)/junkmean
else
  if(.not.present(setZero))then
    getCV_s=hugeRE
  elseif(setZero)then
    getCV_s=zero
  else
    getCV_s=hugeRE
  endif
endif
!End procedure here
endfunction getCV_s
!----------------------------------------------------
pure subroutine getmoments_s1(x,nonmiss,mean,adev,sdev,var,skew,kurt,err,message)
! Purpose: Computes the moments of random variable sample "x"
! IN:  real array "x",optional "nonmiss"
! OUT: mean "mean", average deviation "adev", standard deviation "sdev"
!      variance "var", skew "skew" and kurtosis "kurt"
! ---
! Kreated: DK, circa 2001-2003
! History:
!   2001 - 2012: various changes
!   16 June 2012 AD: Sonnenthal: implemented 'nonmiss' capability
!                                implemented eng link
! ---
! Comments:
use utilities_dmsl_kit,only:zero,three,setVars,getDiffFromMean,getXnonMiss
implicit none
! Dummies
real(mrk),intent(in)::x(:)
logical(mlk),intent(in),optional::nonmiss(:)
real(mrk),intent(out)::mean,adev,sdev,var,skew,kurt
integer(mik),intent(out)::err
character(*),intent(out)::message
! Locals
character(*),parameter::procnam="getmoments_s1"
logical(mlk)::haveNmiss
integer(mik)::n,ntot
real(mrk)::sumx,nr
real(mrk),allocatable::p(:),s(:)
! Start procedure here
call getXnonMiss(x=x,nonmiss=nonmiss,haveNmiss=haveNmiss,&
  nx=ntot,nOk=n,mean=mean,err=err,message=message)
if(err/=EXIT_SUCCESS)then
  message="f-"//procnam//"/&"//message
  call setVars(undefRN,adev,sdev,var,skew,kurt)
  return
endif  ! First pass to get mean, further passes for higher moments
allocate(p(n),s(n))
if(haveNmiss)then; sumx=sum(x,mask=nonmiss)
else;              sumx=sum(x); endif
nr=real(n,mrk); mean=sumx/real(n,mrk)
if(haveNmiss)then
  call getDiffFromMean(x,nonmiss,mean,s)
else
  s=x-mean
endif
call getmoments_eng(s,p,adev,sdev,var,skew,kurt,err,message)
if(err/=EXIT_SUCCESS)then
  message="f-"//procnam//"/&"//message
endif
deallocate(p,s)
! End procedure here
endsubroutine getmoments_s1
!----------------------------------------------------
pure subroutine getmoments_s2(x,nonmiss,mean,adev,sdev,var,skew,err,message)
! Purpose: overload for getmoments with optional arguments
use utilities_dmsl_kit,only:zero,getDiffFromMean,getXnonMiss
implicit none
! Dummies
real(mrk),intent(in)::x(:)
logical(mlk),intent(in),optional::nonmiss(:)
real(mrk),intent(out)::skew
real(mrk),intent(out),optional::mean,adev,var,sdev
integer(mik),intent(out)::err
character(*),intent(out)::message
! Locals
character(*),parameter::procnam="getmoments_s2"
logical(mlk)::haveNmiss
integer(mik)::n,ntot
real(mrk)::sumx,nr,meanL,sdevL,varL
real(mrk),allocatable::p(:),s(:)
! Start procedure here
call getXnonMiss(x=x,nonmiss=nonmiss,haveNmiss=haveNmiss,&
  nx=ntot,nOk=n,mean=mean,err=err,message=message)
selectcase(n)
case(1)       ! single element
  if(present(adev))adev=undefRN
  if(present(var))  var=undefRN; if(present(sdev))sdev=undefRN
endselect
if(err/=EXIT_SUCCESS)then
  message="f-"//procnam//"/&"//message
  return
endif  ! First pass to get mean, further passes for higher moments
allocate(p(n),s(n))
if(haveNmiss)then; sumx=sum(x,mask=nonmiss)
else;              sumx=sum(x); endif
nr=real(n,mrk); meanL=sumx/nr
if(present(mean))mean=meanL
if(haveNmiss)then
  call getDiffFromMean(x,nonmiss,meanL,s)
else
  s=x-meanL
endif
call getmoments_eng(s,p,adev,sdevL,varL,skew,err=err,message=message)
if(err/=EXIT_SUCCESS)then
  message="f-"//procnam//"/&"//message
endif
if(present(sdev))sdev=sdevL; if(present(var))var=varL
! End procedure here
endsubroutine getmoments_s2
!----------------------------------------------------
pure subroutine getmoments_eng(s,p,adev,sdev,var,skew,kurt,err,message)
! Purpose: Engine room of getmoments.
! Programmer: Dmitri Kavetski
! History: 16 June 2012 AD, Sonnenthal
!          8 March 2013 AD, Adelaide N153 - convention that skew=kurt=0 for detz
! Comments
! 1. Uses a tolerance to check for s being identical ("detz")
! ---
! Adapted from Press et al. (1996)
! Modified to implemented population, rather than sample, estimates.
! Additional sources: Wikipedia and Excel :-)
! ---
use utilities_dmsl_kit,only:zero,one,two,three
implicit none
! dummies
real(mrk),intent(in)::s(:)
real(mrk),intent(out)::p(:),var,sdev
real(mrk),intent(out),optional::adev,skew,kurt
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
character(*),parameter::procnam="getmoments_eng"
real(mrk),parameter::sdevTol=epsRe
logical(mlk),parameter::popVar=.true.,popSkew=.true.,popKurt=.true.
real(mrk)::nr,nr1,scal,ep,sumS2,sumS3,sumS4,skewFac,kurtFac1,kurtFac2
logical(mlk)::haveSkew,haveKurt
integer(mik)::n
! Start procedure here
err=0; message=procnam//"/ok"; haveSkew=present(skew); haveKurt=present(kurt)
n=size(s); nr=real(n,mrk); scal=sum(abs(s))/nr  ! "typical" scale of s (for tol checks)
if(present(adev))adev=scal
nr1=real(n-1,mrk); ep=sum(s); p=s**2; sumS2=sum(p)-ep**2/nr
if(popVar)then; var=sumS2/nr1
else;           var=sumS2/nr; endif
sdev=sqrt(var)
if(haveSkew.or.haveKurt)p=p*s;
if(haveSkew)sumS3=sum(p)
if(haveKurt)then
  p=p*s; sumS4=sum(p)
endif
if(sdev>sdevTol*scal)then ! can complete the calculations of skew und kurt
  if(haveSkew)then
!     skew=skew/(nr*sdev**3)      ! NR-based
    skew=sumS3/sqrt(sumS2**3/nr)  ! sample skewness
    if(popSkew)then               ! population skewness (Wiki, Excel)
      skewFac=sqrt(nr*nr1)/(nr-two); skew=skewFac*skew
    endif
  endif
  if(haveKurt)then
    if(popKurt)then               ! population skurtosis
      kurtFac1=nr1*nr*(nr+one)/((nr-two)*(nr-three))
      kurtFac2=nr1**2/((nr-two)*(nr-three))
      kurt=kurtFac1*sumS4/sumS2**2-kurtFac2*three
    else                          ! sample kurtosis
      kurt=sumS4/(nr*var**2)-three
    endif
  endif
else                    ! yo bro we r looking @ a deterministic var
  skew=zero; kurt=zero  ! convention: deterministics have zero var, skew, kurto, etc
endif
! End procedure here
endsubroutine getmoments_eng
!----------------------------------------------------
pure subroutine getmeanvar_s(x,mean,var,method,err,message)
! Purpose: Computes the mean and variance of a sample of a scalar random variable.
! IN:  vector x(1:ns), method="f"-fast scheme,"c"-classic,else 2-pass
! OUT: mean "mean" and variance "var"
! Comments:
! 1. When size(x) is very large, it is possible to overflow the
!    stack for intermediate computations of deviations, since on-the-fly
!    arrays go to the stack. In this case, need to introduce allocatable
!    arrays to store intermediate results on the heap
implicit none
! Dummies
real(mrk),intent(in)::x(:)
real(mrk),intent(out)::mean
real(mrk),intent(out),optional::var
character(*),intent(in)::method
integer(mik),intent(out)::err
character(*),intent(out)::message
! Locals
character(*),parameter::procnam="getmeanvar_s"
real(mrk),allocatable::s(:)
integer(mik)::ns
real(mrk)::nsr
! Start procedure here
err=0; message=procnam//"/ok"; ns=size(x); nsr=real(ns,mrk)
if(ns==0)then
  mean=hugeRe; if(present(var))var=hugeRe
  err=-100; message="f-"//procnam//"/emptySample"; return
endif
mean=sum(x)/nsr ! get approximate expectation of x first
if(present(var))then ! get approximate variance if requested
  if(ns==1)then
    err=-200; var=hugeRe; message="w-"//procnam//"/singleSample"; return
  endif
  selectcase(method)
  case("f")     ! fast scheme:    var[x]=E[x2]-E[x].E[x]
    var=sum(x**2)/real(ns-1,mrk)-mean**2
  case("c")     ! classic scheme: var[x]=E[(x-mean)^2]
    var=sum((x-mean)**2)/real(ns-1,mrk)
  case default  ! 2-pass formula to minimise round-off errors
    allocate(s(ns),stat=err)
    if(err/=0)then
      err=10;message="f-"//procnam//"/allocError";return
    endif
    s=x-mean ! compute deviations for each point
    var=(sum(s**2)-sum(s)**2/nsr)/real(ns-1,mrk)
    deallocate(s,stat=err)
    if(err/=0)then
      err=20;message="f-"//procnam//"/deAllocError";return
    endif
  endselect
endif
! End procedure here
endsubroutine getmeanvar_s
!----------------------------------------------------------------------
pure subroutine getmeancovar_v(x,package,mean,covar,method,err,message)
! Purpose: approximates the first 2 moments of a vector random variable
! x based on a sample
! IN: vector x(1:nv,1:ns)-note packaging (package=1 by default)
!     method="f"-fast scheme,"c"-classic,else 2-pass method
!     package (optional) if package==2, then x(1:ns,1:nv) storage (transposed)
! OUT: mean "mean(1:nv)" and covariance "covar(1:nv,1:nv)"
!      err/=0 -> dimension error
use utilities_dmsl_kit,only:assertEq,zero,quickif
implicit none
! Dummies
real(mrk),intent(in)::x(:,:)
real(mrk),intent(out)::mean(:)
real(mrk),intent(out),optional::covar(:,:)
integer(mik),intent(in),optional::package
character(*),intent(in)::method
integer(mik),intent(out)::err
character(*),intent(out)::message
! Locals
character(*),parameter::procnam="getmeancovar_v"
real(mrk),allocatable::s(:,:)
real(mrk)::nsr,residj
integer(mik)::i,j,k,nv,ns
logical(mlk)::ok
integer(mik)::storage
integer(mik),parameter::var_sam=1,sam_var=2,packageDef=var_sam
! Start procedure here
storage=quickif(package,packageDef)
selectcase(storage)
case(var_sam)
  call assertEq(size(x,dim=1),size(covar,1),size(covar,2),size(mean),ok,nv)
  if(.not.ok)then
    err=-1; message="f-"//procnam//"/arraysDimError"; return
  endif
  ns=size(x,dim=2); nsr=real(ns,mrk)
  mean=sum(x,dim=2)/real(ns,mrk) ! get expectation of x first
  if(present(covar))then ! if requested, approximate the covariance matrix
    selectcase(method)
    case("f","F","l","L"); continue
    case default ! compute deviations for each point
      allocate(s(nv,ns),stat=err)
      if(err/=0)then
        err=-10;message="f-"//procnam//"/allocError";return
      endif
      forall(i=1:nv) s(i,:)=x(i,:)-mean(i)
    endselect
    selectcase(method)
    case("f","F") ! fast scheme: var[x(i,j)]=E[x(i)x(j)]-E[x(i)].E[x(j)]
      forall(i=1:nv)covar(i,i)=sum(x(i,:)**2)/nsr-mean(i)**2
      forall(i=1:nv,j=1:nv,i<j)covar(i,j)=sum(x(i,:)*x(j,:))/nsr-mean(i)*mean(j)
    case("c","C") ! classic scheme: var[x(i,j)]=E[(x(i)-mean)*(x(j)-mean)]
      forall(i=1:nv)covar(i,i)=sum(s(i,:)**2)/nsr
      forall(i=1:nv,j=1:nv,i<j)covar(i,j)=sum(s(i,:)*s(j,:))/nsr
    case("l","L") ! low memory scheme: var[x(i,j)]=E[(x(i)-mean)*(x(j)-mean)]
      covar=zero
      do k=1,ns
        do j=1,nv
          residj=x(j,k)-mean(j)
          covar(j,j) = covar(j,j)+residj**2
          do i=j+1,nv
            covar(i,j) = covar(i,j)+residj*(x(i,k)-mean(i))
          enddo
        enddo
      enddo
      forall(i=1:nv,j=1:nv,i<j)covar(i,j)=covar(i,j)/nsr
    case default  ! 2-pass formula to minimise round-off errors
      forall(i=1:nv)covar(i,i)=(sum(s(i,:)**2)-sum(s(i,:))**2/nsr)/real(ns-1,mrk)
      forall(i=1:nv,j=1:nv,i<j)covar(i,j)=&
        (sum(s(i,:)*s(j,:))-sum(s(i,:))*sum(s(j,:))/nsr)/real(ns-1,mrk)
    endselect
    forall(i=1:nv,j=1:nv,i<j)covar(j,i)=covar(i,j) ! reflect onto upper half
    selectcase(method)
    case("f","F","l","L"); continue
    case default ! deallocate extra storage
      deallocate(s,stat=err)
      if(err/=0)then
        err=-10;message="f-"//procnam//"/deAllocError";return
      endif
    endselect
  endif
case(sam_var)
  call assertEq(size(x,dim=2),size(covar,1),size(covar,2),size(mean),ok,nv)
  if(.not.ok)then
    err=-1; message="f-"//procnam//"/arraysDimError"; return
  endif
  ns=size(x,dim=1); nsr=real(ns,mrk)
  mean=sum(x,dim=1)/real(ns,mrk) ! get expectation of x first
  if(present(covar))then ! if requested, approximate the covariance matrix
    selectcase(method)
    case("f","F","l","L"); continue
    case default ! compute deviations for each point
      allocate(s(ns,nv),stat=err)
      if(err/=0)then
        err=-10;message="f-"//procnam//"/allocError";return
      endif
      forall(i=1:nv) s(:,i)=x(:,i)-mean(i)
    endselect
    selectcase(method)
    case("f")     ! fast scheme: var[x(i,j)]=E[x(i)x(j)]-E[x(i)].E[x(j)]
      forall(i=1:nv)covar(i,i)=sum(x(:,i)**2)/nsr-mean(i)**2
      forall(i=1:nv,j=1:nv,i<j)covar(i,j)=sum(x(:,i)*x(:,j))/nsr-mean(i)*mean(j)
    case("c")     ! classic scheme: var[x(i,j)]=E[(x(i)-mean)*(x(j)-mean)]
      forall(i=1:nv)covar(i,i)=sum(s(:,i)**2)/nsr
      forall(i=1:nv,j=1:nv,i<j)covar(i,j)=sum(s(:,i)*s(:,j))/nsr
    case("l","L") ! low memory scheme: var[x(i,j)]=E[(x(i)-mean)*(x(j)-mean)]
      covar=zero
      do k=1,ns
        do j=1,nv
          residj=x(k,j)-mean(j)
          covar(j,j) = covar(j,j)+residj**2
          do i=j+1,nv
            covar(i,j) = covar(i,j)+residj*(x(k,i)-mean(i))
          enddo
        enddo
      enddo
      forall(i=1:nv,j=1:nv,i<j)covar(i,j)=covar(i,j)/nsr
    case default  ! 2-pass formula to minimise round-off errors
      forall(i=1:nv)covar(i,i)=(sum(s(:,i)**2)-sum(s(:,i))**2/nsr)/real(ns-1,mrk)
      forall(i=1:nv,j=1:nv,i<j)covar(i,j)=&
        (sum(s(:,i)*s(:,j))-sum(s(:,i))*sum(s(:,j))/nsr)/real(ns-1,mrk)
    endselect
    forall(i=1:nv,j=1:nv,i<j)covar(j,i)=covar(i,j) ! reflect onto upper half
    selectcase(method)
    case("f","F","l","L"); continue
    case default ! deallocate extra storage
      deallocate(s,stat=err)
      if(err/=0)then
        err=-10;message="f-"//procnam//"/deAllocError";return
      endif
    endselect
  endif
endselect
err=0; message=procnam//"/ok"
! End procedure here
endsubroutine getmeancovar_v
!----------------------------------------------------------------------
pure subroutine getmeanvar_v(x,package,mean,var,method,err,message)
! Purpose: similar to getmeancovar_v, but if only vars are required
! (does not calculate covariances).
! Efficient for vectorised evaluation of many var's
use utilities_dmsl_kit,only:assertEq,quickif
implicit none
! Dummies
real(mrk),intent(in)::x(:,:)
real(mrk),intent(out)::mean(:)
integer(mik),intent(in),optional::package
real(mrk),intent(out)::var(:)
character(*),intent(in)::method
integer(mik),intent(out)::err
character(*),intent(out)::message
! Locals
character(*),parameter::procnam="getmeanvar_v"
real(mrk),allocatable::s(:,:)
real(mrk)::nsr
integer(mik)::i,nv,ns
logical(mlk)::ok
integer(mik)::storage
integer(mik),parameter::var_sam=1,sam_var=2,packageDef=var_sam
! Start procedure here
storage=quickif(package,packageDef)
selectcase(storage)
case(var_sam)
  call assertEq(size(x,dim=1),size(var),size(mean),ok,nv)
  if(.not.ok)then
    err=1; message="f-"//procnam//"/arraysDimError"; return
  endif
  ns=size(x,dim=2); nsr=real(ns,mrk)
  mean=sum(x,dim=2)/real(ns,mrk) ! get expectation of x first
! approximate the variance vector (diagonal of covariance)
  selectcase(method)
  case("f")     ! fast scheme: var[x(i,j)]=E[x(i)x(j)]-E[x(i)].E[x(j)]
    forall(i=1:nv)var(i)=sum(x(i,:)**2)/nsr-mean(i)**2
  case("c")     ! classic scheme: var[x(i,j)]=E[(x(i)-mean)*(x(j)-mean)]
    forall(i=1:nv)var(i)=sum((x(i,:)-mean(i))**2)/nsr
  case default  ! 2-pass formula to minimise round-off errors
    allocate(s(nv,ns),stat=err)
    if(err/=0)then
      err=10;message="f-"//procnam//"/allocError";return
    endif
    forall(i=1:nv) s(i,:)=x(i,:)-mean(i) ! compute deviations for each point
    forall(i=1:nv)var(i)=(sum(s(i,:)**2)-sum(s(i,:))**2/nsr)/real(ns-1,mrk)
    deallocate(s,stat=err)
    if(err/=0)then
      err=20;message="f-"//procnam//"/deAllocError";return
    endif
  endselect
case(sam_var)
  call assertEq(size(x,dim=2),size(var),size(mean),ok,nv)
  if(.not.ok)then
    err=2; message="f-"//procnam//"/arraysDimError"; return
  endif
  ns=size(x,dim=1); nsr=real(ns,mrk)
  mean=sum(x,dim=1)/real(ns,mrk) ! get expectation of x first
  ! approximate the variance vector (diagonal of covariance)
  selectcase(method)
  case("f")     ! fast scheme: var[x(i,j)]=E[x(i)x(j)]-E[x(i)].E[x(j)]
    forall(i=1:nv)var(i)=sum(x(:,i)**2)/nsr-mean(i)**2
  case("c")     ! classic scheme: var[x(i,j)]=E[(x(i)-mean)*(x(j)-mean)]
    forall(i=1:nv)var(i)=sum((x(:,i)-mean(i))**2)/nsr
  case default  ! 2-pass formula to minimise round-off errors
    allocate(s(ns,nv),stat=err)
    if(err/=0)then
      err=-10;message="f-"//procnam//"/allocError";return
    endif
    forall(i=1:nv) s(:,i)=x(:,i)-mean(i) ! compute deviations for each point
    forall(i=1:nv)var(i)=(sum(s(:,i)**2)-sum(s(:,i))**2/nsr)/real(ns-1,mrk)
    deallocate(s,stat=err)
    if(err/=0)then
      err=-10;message="f-"//procnam//"/deAllocError";return
    endif
  endselect
endselect
err=0; message=procnam//"/ok"
! End procedure here
endsubroutine getmeanvar_v
!----------------------------------------------------------------------
pure function getcovar_v2(x,package,method)
! Purpose: reformats output from the library function getmeanvar to
! deliver the covariance matrix only.
! Allows use of package; 1=nv is leading dimension, x(nv,ns)
implicit none
real(mrk),intent(in)::x(:,:)
character(*),intent(in)::method
integer(mik),intent(in)::package
real(mrk)::getcovar_v2(size(x,package),size(x,package))
! Locals
real(mrk)::junkmean(size(x,dim=package))
character(1)::jmsg
integer(mik)::junkerr
! Start procedure here
call getmeanvar(x,package,junkmean,getcovar_v2,method,junkerr,jmsg)
!End procedure here
endfunction getcovar_v2
!----------------------------------------------------
pure subroutine quicksort1_rv(arr,ascnd,err)
! Sorts arr into ascending numerical order using Quicksort.
! err/=EXIT_SUCCESS indicates scratch memory allocation failure.
use utilities_dmsl_kit,only:quickif,swap
implicit none
! dummies
real(mrk),intent(inout)::arr(:)
logical(mlk),intent(in),optional::ascnd
! locals
real(mrk)::a
INCLUDE "quicksort_X_INC.f90"
! end procedure here
endsubroutine quicksort1_rv
!----------------------------------------------------
pure subroutine quicksort1_iv(arr,ascnd,err)
! Purpose: integer overload
use utilities_dmsl_kit,only:quickif,swap
implicit none
! dummies
integer(mik),intent(inout)::arr(:)
logical(mlk),intent(in),optional::ascnd
! locals
integer(mik)::a
INCLUDE "quicksort_X_INC.f90"
! end procedure here
endsubroutine quicksort1_iv
!----------------------------------------------------
pure subroutine indexx_qsr(arr,indx,ascnd,err,message)
! Purpose: Creates an index of real array arr,
! such that arr(indx) is sorted (ascending or descending)
! Algorithm: employs Quicksort sorting algorithm
use utilities_dmsl_kit,only:assertEq,arthsi,quickif,swap
implicit none
! dummies
real(mrk),intent(in)::arr(:)
logical(mlk),intent(in),optional::ascnd
! locals
character(*),parameter::procnam="indexx_qsr"
real(mrk)::a
INCLUDE "indexx_qsX_INC.f90"
! end procedure here
endsubroutine indexx_qsr
!----------------------------------------------------
pure subroutine indexx_qsi(arr,indx,ascnd,err,message)
! Purpose: overloaded for integer arrays
use utilities_dmsl_kit,only:assertEq,arthsi,quickif,swap
implicit none
! dummies
integer(mik),intent(in)::arr(:)
logical(mlk),intent(in),optional::ascnd
! locals
character(*),parameter::procnam="indexx_qsi"
integer(mik)::a
INCLUDE "indexx_qsX_INC.f90"
! end procedure here
endsubroutine indexx_qsi
!----------------------------------------------------
elemental function normal_logp(x,mean,var)
! Purpose: evaluates the log-probability of a normal variate x
! with expectation mean and variance var.
use utilities_dmsl_kit,only:half,zero,lntwopi
implicit none
! Dummies
real(mrk),intent(in)::x
real(mrk),intent(in),optional::mean,var
real(mrk)::normal_logp
! Locals
real(mrk)::dev
real(mrk),parameter::mhalf_lntwopi=-half*lntwopi
! Start procedure
if(present(mean))then; dev=x-mean
else;                  dev=x; endif
if(present(var))then
  if(dev==zero)then
    normal_logp=-half*(lntwopi+log(var))
  else
    normal_logp=-half*(dev**2/var+lntwopi+log(var))
  endif
else ! assume unit variance
  if(dev==zero)then
    normal_logp=mhalf_lntwopi
  else
    normal_logp=-half*(dev**2+lntwopi)
  endif
endif
! End procedure
endfunction normal_logp
!----------------------------------------------------
elemental function invGamma_logp(x,alpha,beta)
! Purpose: evaluates the log-probability of an inverse gamma variate x
! from inv-G(alpha,beta), alpha=shape, beta=scale
! Ref: Gelman (1995) "Bayesian Data Analysis" pg. 474-5.
use utilities_dmsl_kit,only:gammaln,zero,one
implicit none
! dummies
real(mrk),intent(in)::x,alpha,beta
real(mrk)::invGamma_logp
! Start procedure
if(alpha<=zero.or.beta<=zero.or.x<=zero)then
  invGamma_logp=-hugeRE; return
endif
invGamma_logp=alpha*log(beta)-(alpha+one)*log(x)-gammaln(alpha)-beta/x
! End procedure
endfunction invGamma_logp
!----------------------------------------------------
elemental function invChi2_logp(x,v,s2,getConst)
! Purpose: Evaluates the log-probability of sample x from
! the inverse (scaled) Chi^2 distribution with v degrees of freedom
! and optional scale s^2.
! Ref: Gelman (1995) "Bayesian Data Analysis" pg. 474-5.
! p(x|v,s2) = (v/2)^(v/2) s^v x^(-(v/2+1)) exp(-v s2/(2x)) / GAMMA(v/2)
! same as invGamma(alpha=v/2,beta=v/2 s2)
use utilities_dmsl_kit,only:zero,half,one,quickif
implicit none
! dummies
real(mrk),intent(in)::x,v
real(mrk),optional,intent(in)::s2
logical(mlk),intent(in),optional::getConst
real(mrk)::invChi2_logp
! locals
real(mrk)::hv
logical(mlk)::getConst0
logical(mlk),parameter::getConstDef=.true.
! Start procedure
if(x<=zero.or.v<=zero)then
  invChi2_logp=-hugeRE; return
endif
hv=half*v
getConst0=quickif(getConst,getConstDef)
if(getConst0)then         ! * include constant into computation
  if(present(s2))then     !   - scaled inv-Chi2 distribution
    if(s2>zero)then
      invChi2_logp=invGamma_logp(x=x,alpha=hv,beta=hv*s2)
    else
      invChi2_logp=-hugeRE
    endif
  else                    !   - standard inv-Chi2
      invChi2_logp=invGamma_logp(x=x,alpha=hv,beta=hv)
  endif
else                      ! * do not bother with constant
  if(present(s2))then
    if(s2>zero)then
      invChi2_logp=-(hv+one)*log(x)-hv*s2/x
    else
      invChi2_logp=-hugeRE
    endif
  else
      invChi2_logp=-(hv+one)*log(x)-hv/x
  endif
endif
! End procedure
endfunction invChi2_logp
!----------------------------------------------------
elemental function beta_logp(x,alpha,beta)
! Purpose: returns log-probability of sample x from Beta pdf B(alpha,beta)
! where the prior sample sizes alpha and beta > 0
! Ref: Gelman (1995) "Bayesian Data Analysis" pg. 476-477.
use utilities_dmsl_kit,only:gammaln,zero,one
implicit none
! Dummies
real(mrk),intent(in)::x,alpha,beta
real(mrk)::beta_logp
! Start procedure
if(alpha>zero.and.beta>zero.and.x>=zero.and.x<=one)then
  beta_logp=gammaln(alpha+beta)-gammaln(alpha)-gammaln(beta)+ &
            (alpha-one)*log(x)+(beta-one)*log(one-x)
else
  beta_logp=-hugeRe
endif
! End procedure
endfunction beta_logp
!----------------------------------------------------
elemental function exp_logp(x,beta)
! Purpose: Evaluates the log-probability of sample x from
! the exponential distribution with inverse scale beta.
! Ref: Gelman (1995) "Bayesian Data Analysis" pg. 474-5.
! p(x|beta) = beta exp(-beta x)
! same as gamma(alpha=1,beta=beta)
use utilities_dmsl_kit,only:zero,one
implicit none
! dummies
real(mrk),intent(in)::x,beta
real(mrk)::exp_logp
! Start procedure
if(x>=zero.and.beta>zero)then
  exp_logp=-beta*x+log(beta)
else
  exp_logp=-hugeRE
endif
! End procedure
endfunction exp_logp
!----------------------------------------------------
elemental function poisson_logPmf(theta,rate)
! Purpose: Evaluates the log-probability mass function of
! sample theta from the Poisson distribution with given rate.
! logarithm implemented to avoid overflow on factorial calculation
! for large theta.
use utilities_dmsl_kit,only:factln,zero
implicit none
! dummies
integer(mik),intent(in)::theta
real(mrk),intent(in)::rate
real(mrk)::poisson_logPmf
! Start procedure
if(theta>=0.and.rate>zero)then
  poisson_logPmf=-factln(theta)+theta*log(rate)-rate
else
  poisson_logPmf=-hugeRE
endif
! End procedure
endfunction poisson_logPmf
!----------------------------------------------------
elemental function binomial_pmf(theta,n,p)
! Purpose: Evaluates the probability mass function of
! sample theta from the Binomial distribution with given
! number of trials n and individual probability p.
use utilities_dmsl_kit,only:bico,zero,one
implicit none
! dummies
integer(mik),intent(in)::theta,n
real(mrk),intent(in)::p
real(mrk)::binomial_Pmf
! locals
real(mrk)::thetar,nr
! Start procedure
if(theta>=0.and.theta<=n.and.n>0.and.p>=zero.and.p<=one)then
  thetar=real(theta,mrk); nr=real(n,mrk)
  binomial_Pmf=bico(n=n,k=theta)*p**thetar*(one-p)**(nr-thetar)
else
  binomial_Pmf=-hugeRE
endif
! End procedure
endfunction binomial_Pmf
!----------------------------------------------------
elemental function negBinomial_pmf(theta,alpha,beta)
! Purpose: Evaluates the probability mass function of
! sample theta from the Negative-Binomial distribution with given
! shape alpha and inverse scale beta
use utilities_dmsl_kit,only:bico,zero,one
implicit none
! dummies
integer(mik),intent(in)::theta,alpha
real(mrk),intent(in)::beta
real(mrk)::negBinomial_pmf
! locals
real(mrk)::thetar,alphar,b1
! Start procedure
if(theta>=0.and.alpha>0.and.beta>zero)then
  thetar=real(theta,mrk); alphar=real(alpha,mrk); b1=one/(beta+one)
  negBinomial_pmf=bico(n=theta+alpha-1,k=alpha-1)*&
                  (beta*b1)**alpha*b1**theta
else
  negBinomial_pmf=-hugeRE
endif
! End procedure
endfunction negBinomial_pmf
!----------------------------------------------------
elemental function normal_cdf_1(x)
! Purpose: Evaluates the Gaussian cdf.
! Warning: definitely not full double precision
! Source: Abramowitz and Stegun 26.2.17, absolute error<7.5e-8
use utilities_dmsl_kit,only:zero,one,half,sqrt2pi
implicit none
! Dummies
real(mrk), intent(in) :: x
real(mrk) :: normal_cdf_1
! Locals
real(mrk),parameter::b1= 0.319381530_mrk,&  ! polynomial coefficients
                     b2=-0.356563782_mrk,&
                     b3= 1.781477937_mrk,&
                     b4=-1.821255978_mrk,&
                     b5= 1.330274429_mrk,&
                     p=0.2316419_mrk,c=one/sqrt2pi
real(mrk) :: xabs, t, z, b
! Start procedure
if(x==zero)then
  normal_cdf_1=half; return
endif
xabs=abs(x); t=one/(one+p*xabs); z=c*exp(-xabs**2*half)
b=((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t
if(x>zero)then
  normal_cdf_1=one-z*b
else
  normal_cdf_1=z*b
endif
! End procedure
endfunction normal_cdf_1
!----------------------------------------------------
elemental function normal_cdf_2(x,mean,sdev)
! Purpose: wrapper for non-normalised normal cdf
implicit none
! Dummies
real(mrk),intent(in)::x
real(mrk),intent(in),optional::mean,sdev
real(mrk)::normal_cdf_2
! Locals
real(mrk)::dev
! Start procedure
if(present(mean))then; dev=x-mean
else;                  dev=x; endif
if(present(sdev))then
  normal_cdf_2=normal_cdf_1(dev/sdev)
else ! assume unit variance
  normal_cdf_2=normal_cdf_1(dev)
endif
! End procedure
endfunction normal_cdf_2
!----------------------------------------------------
elemental function exp_cdf(x,beta)
! Purpose: Evaluates the cdf of an exponential deviate
use utilities_dmsl_kit,only:zero,one
implicit none
! dummies
real(mrk),intent(in)::x,beta
real(mrk)::exp_cdf
! Start procedure
if(x>=zero.and.beta>zero)then
  exp_cdf=one-exp(-beta*x)
else
  exp_cdf=-hugeRE
endif
! End procedure
endfunction exp_cdf
!----------------------------------------------------
elemental function poisson_cdf(x,Ex)
! Purpose: evaluates the cumulative Poisson distribution P(X<=x|Ex) with expected
! value Ex. Statistical meaning x="number of events"
! NB: in NR-90, P(X<x)
use utilities_dmsl_kit,only:gammq
implicit none
! Dummies
integer(mik),intent(in)::x
real(mrk),intent(in)::Ex
real(mrk)::poisson_cdf
! Start procedure
if(x>=0)then  ! non-negative x
  poisson_cdf=gammq(real(x+1,mrk),Ex) ! note x+1 as P(X<=x), not P(X<x)
else          ! negative x
  poisson_cdf=-hugeRE
endif
! End procedure
endfunction poisson_cdf
!----------------------------------------------------------------------
elemental function binomial_cdf(x,n,p)
! Purpose: evaluates the cumulative Binomial distribution P(X<=x|n,p)
! with trial probability p. x=number of events in n trials.
! NB: In NR-90, P(X>=x|n,p) (complement Q of the cumulative function)
use utilities_dmsl_kit,only:one,betai
implicit none
! Dummies
integer(mik),intent(in)::x,n
real(mrk),intent(in)::p
real(mrk)::binomial_cdf
! Start procedure
if(x>=0)then  ! non-negative x
  binomial_cdf=one-betai(real(x,mrk),real(n-x+2,mrk),p)
else          ! negative x
  binomial_cdf=-hugeRE
endif
! End procedure
endfunction binomial_cdf
!----------------------------------------------------------------------
elemental function PPND16_func(P)result(PPND16)
! Purpose: ALGORITHM AS241 APPL. STATIST. (1988) VOL. 37, NO. 3.
! Produces the normal deviate Z corresponding to a given lower
! tail area of P; Z is accurate to about 1 part in 10**16.
! The hash sums below are the sums of the mantissas of the coefficients.
! They are included for use in checking transcription.
! Comments:
! 1. Full 64-bit version of ppnd routine.
! 2. In order to maintain PUREity, no strong error traps for 0<=p>=1;
use utilities_dmsl_kit,only:zero,one,half
implicit none
! dummies
real(mrk),intent(in)::P
real(mrk)::PPND16
! locals
real(mrk)::Q,R
integer(mik)::IFAULT
real(mrk),parameter::SPLIT1=0.425_mrk,SPLIT2=5._mrk,&
                     CONST1=0.180625_mrk,CONST2=1.6_mrk
real(mrk):: & ! parameters of the approximation
  A0, A1, A2, A3, A4, A5, A6, A7, B1, B2, B3, B4, B5, B6, B7,&
  C0, C1, C2, C3, C4, C5, C6, C7, D1, D2, D3, D4, D5, D6, D7,&
  E0, E1, E2, E3, E4, E5, E6, E7, F1, F2, F3, F4, F5, F6, F7
! Coefficients for P close to 0.5
parameter(A0 = 3.3871328727963666080e+0_mrk,&
          A1 = 1.3314166789178437745e+2_mrk,&
          A2 = 1.9715909503065514427e+3_mrk,&
          A3 = 1.3731693765509461125e+4_mrk,&
          A4 = 4.5921953931549871457e+4_mrk,&
          A5 = 6.7265770927008700853e+4_mrk,&
          A6 = 3.3430575583588128105e+4_mrk,&
          A7 = 2.5090809287301226727e+3_mrk,&
          B1 = 4.2313330701600911252e+1_mrk,&
          B2 = 6.8718700749205790830e+2_mrk,&
          B3 = 5.3941960214247511077e+3_mrk,&
          B4 = 2.1213794301586595867e+4_mrk,&
          B5 = 3.9307895800092710610e+4_mrk,&
          B6 = 2.8729085735721942674e+4_mrk,&
          B7 = 5.2264952788528545610e+3_mrk)
! HASH sum AB 55.8831928806149014439
! Coefficients for P not close to 0, 0.5 or 1.
parameter(C0 = 1.42343711074968357734e+0_mrk,&
          C1 = 4.63033784615654529590e+0_mrk,&
          C2 = 5.76949722146069140550e+0_mrk,&
          C3 = 3.64784832476320460504e+0_mrk,&
          C4 = 1.27045825245236838258e+0_mrk,&
          C5 = 2.41780725177450611770e-1_mrk,&
          C6 = 2.27238449892691845833e-2_mrk,&
          C7 = 7.74545014278341407640e-4_mrk,&
          D1 = 2.05319162663775882187e+0_mrk,&
          D2 = 1.67638483018380384940e+0_mrk,&
          D3 = 6.89767334985100004550e-1_mrk,&
          D4 = 1.48103976427480074590e-1_mrk,&
          D5 = 1.51986665636164571966e-2_mrk,&
          D6 = 5.47593808499534494600e-4_mrk,&
          D7 = 1.05075007164441684324e-9_mrk)
! HASH sum CD 49.33206503301610289036
! Coefficients for P near 0 or 1.
parameter(E0 = 6.65790464350110377720e+0_mrk,&
          E1 = 5.46378491116411436990e+0_mrk,&
          E2 = 1.78482653991729133580e+0_mrk,&
          E3 = 2.96560571828504891230e-1_mrk,&
          E4 = 2.65321895265761230930e-2_mrk,&
          E5 = 1.24266094738807843860e-3_mrk,&
          E6 = 2.71155556874348757815e-5_mrk,&
          E7 = 2.01033439929228813265e-7_mrk,&
          F1 = 5.99832206555887937690e-1_mrk,&
          F2 = 1.36929880922735805310e-1_mrk,&
          F3 = 1.48753612908506148525e-2_mrk,&
          F4 = 7.86869131145613259100e-4_mrk,&
          F5 = 1.84631831751005468180e-5_mrk,&
          F6 = 1.42151175831644588870e-7_mrk,&
          F7 = 2.04426310338993978564e-15_mrk)
! HASH sum EF 47.52583317549289671629
! Start procedure here
IFAULT=0
if(P==half)then ! trivial case
  PPND16=zero; return
endif
Q=P-HALF
if(abs(Q)<=SPLIT1)then
  R=CONST1-Q**2
  PPND16=Q*(((((((A7*R+A6)*R+A5)*R+A4)*R+A3)*R+A2)*R+A1)*R+A0)/&
           (((((((B7*R+B6)*R+B5)*R+B4)*R+B3)*R+B2)*R+B1)*R+ONE)
else
  if(Q<ZERO)then; R=P
  else;           R=ONE-P; endif
  if(R<=ZERO)then
    PPND16=hugeRe;IFAULT=1;return
  endif
  R=sqrt(-log(R))
  if(R<=SPLIT2)then
    R=R-CONST2
    PPND16=(((((((C7*R+C6)*R+C5)*R+C4)*R+C3)*R+C2)*R+C1)*R+C0)/&
           (((((((D7*R+D6)*R+D5)*R+D4)*R+D3)*R+D2)*R+D1)*R+ONE)
  else
    R=R-SPLIT2
    PPND16=(((((((E7*R+E6)*R+E5)*R+E4)*R+E3)*R+E2)*R+E1)*R+E0)/&
           (((((((F7*R+F6)*R+F5)*R+F4)*R+F3)*R+F2)*R+F1)*R+ONE)
  endif
  if(Q<ZERO)PPND16=-PPND16
endif
! End procedure here
endfunction PPND16_func
!----------------------------------------------------
elemental function normal_icdf_v(p,mean,sdev)
! Generic inverse normal cdf with mean and sdev
implicit none
! dummies
real(mrk),intent(in)::mean,sdev,p
real(mrk)::normal_icdf_v
! Start procedure here
normal_icdf_v=normal_icdf(p)*sdev+mean
! End procedure here
endfunction normal_icdf_v
!----------------------------------------------------
elemental function exp_icdf(p,beta)
! Purpose: Evaluates the icdf of an exponential deviate
use utilities_dmsl_kit,only:zero,one
implicit none
! dummies
real(mrk),intent(in)::p,beta
real(mrk)::exp_icdf
! Start procedure
if(p>=zero.and.p<one.and.beta>zero)then
  exp_icdf=-log(one-p)/beta
else
  exp_icdf=-hugeRE
endif
! End procedure
endfunction exp_icdf
!----------------------------------------------------
subroutine uniran_s(harvest,ixs)
! WARNING: Different from original DMSL file
! Disabled multiple random streams option to avoid using uniran1_parallel_dmsl_mod
! Original uniran_s from DMSL is commented out below
use uniran1_dmsl_mod,only:generator=>uniran_s
implicit none
! dummies
real(mrk),intent(out)::harvest
integer(mik),intent(in),optional::ixs
! Start procedure here
if(present(ixs))then
  write(*,*) 'WARNING: multiple-streams option has been disabled.'
  write(*,*) 'Running single-stream option instead, provided [ixs] not used.'
endif
call generator(harvest)
! End procedure here
endsubroutine uniran_s
!----------------------------------------------------
!subroutine uniran_s(harvest,ixs)
!! Purpose: Enhanced to work with multiple random streams
!use uniran1_dmsl_mod,only:generator=>uniran_s
!use uniran1_parallel_dmsl_mod,only:generator_parallel=>uniran_parallel_s
!implicit none
!! dummies
!real(mrk),intent(out)::harvest
!integer(mik),intent(in),optional::ixs
!! Start procedure here
!if(present(ixs))then
!  call generator_parallel(harvest,ixs)
!else
!  call generator(harvest)
!endif
!! End procedure here
!endsubroutine uniran_s
!----------------------------------------------------
subroutine uniran_v(harvest)
! Purpose: Generate a standard uniform random vector
implicit none
real(mrk),intent(out)::harvest(:)
integer(mik)::i
! Start procedure here
do i=1,size(harvest)
  call uniran_s(harvest(i))
enddo
! End procedure here
endsubroutine uniran_v
!----------------------------------------------------
subroutine uniranHyper_s(harvest,low,high)
! Purpose: Generate a random point in hyperspace bound by low->high
implicit none
real(mrk),intent(in)::low,high
real(mrk),intent(out)::harvest
! Start procedure here
call uniran(harvest)
harvest=low+harvest*(high-low)
! End procedure here
endsubroutine uniranHyper_s
!----------------------------------------------------
subroutine uniranHyper_v(harvest,low,high)
! Purpose: Vector version of uniranHyper_s
implicit none
real(mrk),intent(in)::low(:),high(:)
real(mrk),intent(out)::harvest(:)
! Start procedure here
call uniran(harvest)
harvest=low+harvest*(high-low)
! End procedure here
endsubroutine uniranHyper_v
!----------------------------------------------------
subroutine uniran_is(harvest)
! Purpose: random integer with uniform probability on 0 and 1
use utilities_dmsl_kit,only:half
implicit none
! Dummies
integer(mik),intent(out)::harvest
! Locals
real(mrk)::temp
! Start procedure here
call uniran(temp)
harvest=merge(0,1,temp<=half)
! End procedure here
endsubroutine uniran_is
!----------------------------------------------------
subroutine uniran_iv(harvest)
! Purpose: overloaded for vector results
use utilities_dmsl_kit,only:half
implicit none
! Dummies
integer(mik),intent(out)::harvest(:)
! Locals
real(mrk),dimension(size(harvest))::temp
! Start procedure here
call uniran(temp)
harvest=merge(0,1,temp<=half)
! End procedure here
endsubroutine uniran_iv
!----------------------------------------------------
subroutine uniranHyper_is(harvest,low,high)
! Purpose: scalar random integer uniform between lower and upper inclusive
implicit none
! Dummies
integer(mik),intent(in)::low,high
integer(mik),intent(out)::harvest
! Locals
real(mrk)::temp
! Start procedure here
call uniran(temp)
harvest=int(low+temp*(high+1-low),mik)
! End procedure here
endsubroutine uniranHyper_is
!----------------------------------------------------
subroutine uniranHyper_iv(harvest,low,high)
! Purpose: overloaded for vector results
implicit none
! Dummies
integer(mik),intent(in)::low(:),high(:)
integer(mik),intent(out)::harvest(:)
! Locals
real(mrk)::temp(size(harvest))
! Start procedure here
call uniran(temp)
harvest=int(low+temp*(high+1-low),mik)
! End procedure here
endsubroutine uniranHyper_iv
!----------------------------------------------------
subroutine normaldev1_s(mean,sdev,var,ixs,gdev,err,message)
! Purpose: Returns a scalar normal deviate "gdev" with mean "mean"
! and variance "var", i.e., gdev ~ N(mean,sdev^2) or ~ N(mean,var)
! if mean/var(sdev) omitted assume standard values
! both sdev/var can not be specified.
! sdev=0 or var=0 yield deterministic result (mean), ie, Dirac function.
! err=status indicator, message=return message
! Comments:
! 1. ixs allows working with multiple random streams
use utilities_dmsl_kit,only:zero
implicit none
! dummies
real(mrk),intent(in),optional::mean
real(mrk),intent(in),optional::sdev,var
integer(mik),intent(in),optional::ixs
real(mrk),intent(out)::gdev
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
real(mrk)::udev,gStDev
! Start procedure
call uniran(udev,ixs); gStDev=normal_icdf(udev); gdev=undefRN
err=0; message="normaldev1_s/ok"
if(present(sdev).and.present(var))then
    err=-10; message="f-normaldev1_s/fatal/both(sdev,var)specified"
elseif(present(mean).and.present(sdev))then
  if(sdev>=zero)then  ! use sdev
    gdev=mean+sdev*gStDev
  else
    err=-1;message="f-normaldev1_s/fatal/sdev<0"
  endif
elseif(present(mean).and.present(var))then
  if(var>=zero)then   ! use var
    gdev=mean+sqrt(var)*gStDev
  else
    err=-2;message="f-normaldev1_s/fatal/var<0"
  endif
elseif(present(mean))then  ! assume unit variance
    gdev=mean+gStDev
elseif(present(sdev))then
  if(sdev>=zero)then  ! zero mean
    gdev=sdev*gStDev
  else
    err=-3;message="f-normaldev1_s/fatal/sdev<0"
  endif
elseif(present(var))then
  if(var>=zero)then   ! use var
    gdev=sqrt(var)*gStDev
  else
    err=-4;message="f-normaldev1_s/fatal/var<0"
  endif
else  ! nothing supplied: return standard normal deviate
    gdev=gStDev
endif
! End procedure
endsubroutine normaldev1_s
!----------------------------------------------------
subroutine normaldev1_1m1v_vd(mean,sdev,var,gdev,err,message)
! Purpose: Overload: vector normal deviate with scalar mean and scalar variance.
use utilities_dmsl_kit,only:zero,assertEq
implicit none
! Dummies
real(mrk),intent(in)::mean
real(mrk),optional,intent(in)::sdev,var
real(mrk),intent(out)::gdev(:)
integer(mik),intent(out)::err
character(*),intent(out)::message
! Locals
! Start procedure
call uniran(gdev); gdev=normal_icdf(gdev)
err=0; message="normaldev1_1m1v_vd/ok"
if(present(sdev).and.present(var))then
    err=-10;message="f-normaldev1_1m1v_vd/fatal/sdev&var";gdev=-hugeRE
elseif(present(sdev))then   ! use sdev
  if(sdev<zero)then
    err=-2; message="f-normaldev1_1m1v_vd/fatal/sdev<0"; gdev=-hugeRE
  else
    gdev=mean+sdev*gdev
  endif
elseif(present(var))then    ! use variance
  if(var<zero)then
    err=-4; message="f-normaldev1_1m1v_vd/fatal/var<0"; gdev=-hugeRE
  else
    gdev=mean+sqrt(var)*gdev
  endif
else  ! assume unit variance
  gdev=mean+gdev
endif
! End procedure
endsubroutine normaldev1_1m1v_vd
!----------------------------------------------------
subroutine normaldev1_NmNv_vd(mean,sdev,var,gdev,err,message)
! Purpose: vector normal deviate with vector mean and vector variance
! Mean must be supplied, then choice - EITHER var or sdev
use utilities_dmsl_kit,only:zero,assertEqLog
implicit none
! Dummies
real(mrk),dimension(:),optional,intent(in)::mean
real(mrk),dimension(:),optional,intent(in)::sdev,var
real(mrk),dimension(:),intent(out)::gdev
integer(mik),intent(out)::err
character(len=*),intent(out)::message
! Locals
! Start procedure
call uniran(gdev); gdev=normal_icdf(gdev)
err=0; message="normaldev1_NmNv_vd/ok"
if(present(sdev).and.present(var))then ! error
    err=-10;message="f-normaldev1_NmNv_vd/fatal/sdev&var";gdev=-hugeRE
elseif(present(sdev).and.present(mean))then   ! * user sdev and user mean
  if(.not.assertEqLog(size(mean),size(sdev),size(gdev)))then
    err=-1; message="f-normaldev1_NmNv_vd/dimError";    gdev=-hugeRE
  elseif(any(sdev<zero))then
    err=-2; message="f-normaldev1_NmNv_vd/fatal/sdev<0";gdev=-hugeRE
  else
    gdev=mean+sdev*gdev
  endif
elseif(present(sdev))then                     ! * user sdev and zero mean
  if(.not.assertEqLog(size(sdev),size(gdev)))then
    err=-3; message="f-normaldev1_NmNv_vd/dimError";    gdev=-hugeRE
  elseif(any(sdev<zero))then
    err=-4; message="f-normaldev1_NmNv_vd/fatal/sdev<0";gdev=-hugeRE
  else
    gdev=sdev*gdev
  endif
elseif(present(var).and.present(mean))then    ! * user variance and user mean
  if(.not.assertEqLog(size(mean),size(var),size(gdev)))then
    err=-5; message="f-normaldev1_NmNv_vd/dimError";    gdev=-hugeRE
  elseif(any(var<zero))then
    err=-4; message="f-normaldev1_NmNv_vd/fatal/var<0"; gdev=-hugeRE
  else
    gdev=mean+sqrt(var)*gdev
  endif
elseif(present(var))then                      ! * user variance and zero mean
  if(.not.assertEqLog(size(var),size(gdev)))then
    err=-6; message="f-normaldev1_NmNv_vd/dimError";    gdev=-hugeRE
  elseif(any(var<zero))then
    err=-7; message="normaldev1_NmNv_vd/fatal/var<0"; gdev=-hugeRE
  else
    gdev=sqrt(var)*gdev
  endif
elseif(present(mean))then                     ! * user mean and unit variance
  if(size(mean)/=size(gdev))then
    err=-8; message="f-normaldev1_NmNv_vd/dimError";    gdev=-hugeRE
  else
    gdev=mean+gdev
  endif
!else                                         ! * zero mean and unit variance
endif
! End procedure
endsubroutine normaldev1_NmNv_vd
!----------------------------------------------------
subroutine normaldev_v(mean,covar,scal,do_dcmp,covar_dcmp,gdev,err,message)
! Purpose: Returns a vector normal deviate "gdev"~ N(mean,covar)
! with mean "mean" and covariance matrix "covar"
! err=status indicator, message=return message
! if do_dcmp then need to Choleski the covariance matrix
!    if covar provided, it will store its Choleski factor in covar_dcmp
!    else it will overwrite covar_dcmp with its Choleski factor
! if .not.do_dcmp then assumes the decomposed form is in covar_dcmp.
! NB: covar is optional, covar_dcmp is NOT optional.
! Comments:
! - thanx to Frosty for pointing out the "LV" vs "L" bug in this routines
use utilities_dmsl_kit,only:fmatmul_mv
use linalg_dmsl_kit,only:choles_dcmp
implicit none
! Dummies
real(mrk),intent(in)::mean(:)
real(mrk),optional,intent(in)::covar(:,:),scal
logical(mlk),intent(in)::do_dcmp
real(mrk),intent(inout)::covar_dcmp(:,:)
real(mrk),intent(out)::gdev(:)
integer(mik),intent(out)::err
character(*),intent(out)::message
! Locals
logical(mlk)::pd
character(100)::msg
! Start procedure
err=0; message="normaldev_v/ok"
if(do_dcmp)then   ! do Choleski decomposition of supplied covar matrix
  if(present(covar))then ! store Chol factor of "covar" in "covar_dcmp"
    call choles_dcmp(a=covar,Lmat=covar_dcmp,posDefinite=pd,&
                     err=err,message=msg)
  else            ! overwrite "covar_dcmp" with its Choleski factor
    call choles_dcmp(a=covar_dcmp,posDefinite=pd,err=err,message=msg)
  endif
  if(.not.pd.or.err/=0)then   ! negative semi-definite matrix
    message="f-normaldev_v/"//msg; return
  endif
endif                                 ! generate uncorrelated vector
call normaldev1_NmNv_vd(gdev=gdev,err=err,message=msg)
if(err/=0)then
  message="f-normaldev_v/"//msg; return
endif
gdev=fmatmul_mv(covar_dcmp,gdev,"LV") ! introduce correlation
if(present(scal))gdev=scal*gdev       ! scale
gdev=mean+gdev                        ! and shift
! End procedure
endsubroutine normaldev_v
!----------------------------------------------------
subroutine normaldev2_v(mean,covar,scal,gdev,err,message)
! Purpose: Compact version of normaldev_v, with Choleski decomposition
! at each sample. Should only be used when covar changes all the time.
implicit none
! Dummies
real(mrk),intent(in)::mean(:),covar(:,:)
real(mrk),intent(out)::gdev(:)
real(mrk),intent(in),optional::scal
integer(mik),intent(out)::err
character(*),intent(out)::message
! Locals
real(mrk)::covar_dcmp(size(covar,dim=1),size(covar,dim=1))
logical(mlk),parameter::do_dcmp=.true.
! Start procedure
call normaldev_v(mean,covar,scal,do_dcmp,covar_dcmp,gdev,err,message)
! End procedure
endsubroutine normaldev2_v
!----------------------------------------------------
subroutine gamdevint_s(ia,betar,gamdev,err,message)
! Purpose: returns a deviate drawn from the gamma distribution
! of integer order ia and real (optional) mean betar. ia,betar>0.
! Employs algoritm of Press et al. (1996)
use utilities_dmsl_kit,only:zero,one,two
implicit none
! Dummies
integer(mik), intent(in) :: ia
real(mrk), optional, intent(in) :: betar
real(mrk), intent(out) :: gamdev
integer(mik), intent(out) :: err
character(len=*), intent(out):: message
! Locals
integer(mik),parameter :: armet=10 ! use direct method if ia<=dirmet
real(mrk) :: am, e, h, s, x, y, v(2), arr(armet-1)
! Start procedure
if(ia<1)then  ! error. set output to garbage
  gamdev=-hugeRE;err=-10; message="gamdevint_s/fatal/ia<1"
  return
else  ! looks ok
  err=0; message="gamdevint_s: ok"
endif
if(ia==1)then ! get exponential deviate
  call expdev(x)
elseif(ia<armet)then  ! use direct method,adding waiting times
  call uniran(arr(1:ia))
  x=-log(product(arr(1:ia)))
else  ! use acceptance-rejection method
  do
    call uniran(v)
    v(2)=two*v(2)-one ! generate tangent of a random angle
    if(dot_product(v,v)>one) cycle
    y=v(2)/v(1)
    am=ia-1
    s=sqrt(two*am+one)
    x=s*y + am  ! decide whether to accept or reject x:
    if(x<=zero) cycle ! reject in region of zero probability
    e=(one+y**2)*exp(am*log(x/am)-s*y) ! ratio of target density
    call uniran(h)  ! to comparison function
    if(h<=e) exit  ! reject on basis of second uniform deviate
  enddo
endif
if(present(betar))then  ! mean beta provided as well
  if(betar>zero)then  ! shift to correct mean
    gamdev=betar*x
  else                    ! illegal input
    gamdev=-hugeRE;err=-10; message="gamdevint_s/fatal/betar<0"
    return
  endif
else
  gamdev=x
endif
! End procedure
endsubroutine gamdevint_s
!----------------------------------------------------
subroutine gamdevr_s(alpha,beta,gamdev,err,message)
! Purpose: jacket program: returns a deviate drawn from the
! gamma distribution of real order alpha and scale beta. alpha,beta>0.
! if beta omitted beta=1.
! NB: beta is the scale factor, NOT the inverse scale factor
!     Gelman et al (1995) use the inverse scale factor in Gamma
! Algorithms GS*  (for 0<alpha<1), bounded runtime;
!            exponential sampler (for alpha=1); bounded runtime
!            GKM1 (for 1<alpha<=alpha0), O(alpha^1/2) runtime;
!            GKM2 (for alpha>alpha0), O(1) runtime.
! Ref: Fishman, 1996, pp. 194-200
!----
! For applications where speed is critical, it may be worthwhile
! to check if alpha ~1-10. then, for these cases, it may be worhtwhile
! to exploits a direct method accumulating waiting times.
! Timing analysis shows that, on a serial PII, 350MHz, the a-r and direct
! algorithms are roughly equivalent when alpha<10. then a-r becomes
! substantially faster. also note a-r caters for 0<alpha<infinity
!----
use utilities_dmsl_kit,only:zero,one
implicit none
! Dummies
real(mrk), intent(in) :: alpha
real(mrk), optional, intent(in) :: beta
real(mrk), intent(out) :: gamdev
integer(mik), intent(out) :: err
character(len=*), intent(out):: message
! Locals
real(mrk),parameter::alpha0=2.5_mrk ! change from GKM1 to GKM2 for speed.
! Start procedure
err=0; message="gamdevr_s/ok"
if(alpha<=zero)then     ! illegal argument. alpha must be +ve
  err=-10; message="gamdevr_s/alpha<0"
  gamdev=-hugeRE; return
elseif(alpha<one)then  ! use GS*
  gamdev=gamdev_gs(alpha)
elseif(alpha==one)then ! use exponential sampler
  call expdev_s(gamdev)
elseif(alpha<=alpha0)then  ! use GKM1
  gamdev=gamdev_gkm1(alpha)
else                        ! use GKM2
  gamdev=gamdev_gkm2(alpha)
endif
if(present(beta))then   ! scale factor beta provided as well:
  if(beta>zero)then   ! ok, do scaling
    gamdev=beta*gamdev
  else                    ! illegal input; beta must be +ve
    err=-10; message="gamdevr_s/betar<0"
    gamdev=-hugeRE; return
  endif
endif
! End procedure
endsubroutine gamdevr_s
!----------------------------------------------------
function gamdev_gs(alpha)
! Purpose: implements a fast bounded time algorithm
! for sampling the Gamma distribution Gamma(alpha) for 0<alpha<1.
! ---
! Programmer: Dmitri Kavetski
! History:
!   circa 2000, Station St Waratah
!   20 Feb 2012, Guillaume's suggestion to use log-space
! ---
! Method: acceptance/rejection sampling, algorithm GS*,
! based on Algorithm GS from Ahrens and Dieter (1974).
! Ref: Algorithm GS*, Fishman (1996), pp. 194.
! ---
! Notes:
! 1. The constraint on alpha must be checked by the calling program.
! 2. In Fishman, second acceptance statement has a typographic error:
!    should accept when w>=z, not w<=z)
use utilities_dmsl_kit,only:zero,one,eBase
implicit none
! Dummies
real(mrk),intent(in)::alpha
real(mrk)::gamdev_gs
! Locals
real(mrk)::b,y,z,w,logw
real(mrk),parameter::alphaCrit4Log=0.95_mrk
logical(mlk)::compareInLogSpace
! Start procedure here
if(alpha<zero.or.alpha>one)then ! illegal argument
  gamdev_gs=-hugeRe; return
endif
b=(alpha+eBase)/eBase   ! precompute constant
do ! until accept
  call uniran(y); y=b*y ! generate u ~ U(0,b)
  if(y<=one)then
    z=y**(one/alpha)
    call expdev_s(w)  ! generate exponential deviate E(1)
    if(w>=z)exit      ! accept
  else
    z=-log((b-y)/alpha)
    call uniran(w);   ! uniform call
    compareInLogSpace=(alpha>alphaCrit4Log)
    if(compareInLogSpace)then ! *** Guillaume's enhancement
      logw=(one/(alpha-one))*log(w)
      if(logw>=log(z))exit
    else                      ! *** standard branch from Fishman
      w=w**(one/(alpha-one))
      if(w>=z)exit    ! accept (NB: in Fishman, w<=z, typo?)
    endif
  endif
enddo
gamdev_gs=z
! End procedure here
endfunction gamdev_gs
!----------------------------------------------------
function gamdev_gkm1(alpha)
! Purpose: implements a fast bounded time algorithm
! for sampling the Gamma distribution Gamma(alpha) for
! alpha>1. For efficiency, alpha<=2.5; runtime O(alpha^1/2);
! Method: acceptance/rejection sampling,
!         based on algorithm of Cheng and Feast (1979).
! Ref: Algorithm GKM1, Fishman (1996), pp. 197.
! NB: the constrain on alpha must be checked by the calling program
use utilities_dmsl_kit,only:zero,one,two
implicit none
! Dummies
real(mrk), intent(in) :: alpha
real(mrk) :: gamdev_gkm1
! Locals
real(mrk),parameter::onesixth=one/6._mrk
real(mrk) :: a, b, m, d, x, v
! Start procedure here
if(alpha<=one)then ! illegal argument
  gamdev_gkm1=-hugeRe; return
endif
a=alpha-one; b=(alpha-onesixth/alpha)/a; m=two/a; d=m+two ! setup
do  ! until accept
  call uniran(x); call uniran(v)
  v=b * v / x
  if(m*x-d+v+one/v<=zero)then ! fast acceptance
    exit
  elseif(m*log(x)-log(v)+v-one<=zero)then ! accept
    exit
  endif
enddo
gamdev_gkm1=a * v
! End procedure here
endfunction gamdev_gkm1
!----------------------------------------------------
function gamdev_gkm2(alpha)
! Purpose: implements a fast bounded time algorithm
! for sampling the Gamma distribution Gamma(alpha) for
! alpha>1. For efficiency, alpha>2.5; runtime O(1);
! Method: acceptance/rejection sampling, Algorithm GKM2;
! Ref: Fishman (1996), pp. 200.
! NB: the constrain on alpha must be checked by the calling program
use utilities_dmsl_kit,only:zero,one,two,sqrt2,sqrtEbase
implicit none
! Dummies
real(mrk), intent(in) :: alpha
real(mrk) :: gamdev_gkm2
! Locals
real(mrk),parameter::onesixth=one/6._mrk,c1=one+sqrt2/sqrtEbase
! exact: c1=1+sqrt(2/e)=1.857763884960706796480190_mrk
real(mrk) :: a, b, m, d, f, x, v
! Start procedure here
if(alpha<=one)then ! illegal argument
  gamdev_gkm2=-hugeRe; return
endif
a=alpha-one; b=(alpha-onesixth/alpha)/a ! setup
m=two/a; d=m+two; f=sqrt(alpha)
do    ! until accept
  do  ! until accept
    call uniran(x); call uniran(v)
    x=v+(one-c1*x)/f
    if(x>zero.and.x<one) exit
  enddo
  v=b*v/x
  if(m*x-d+v+one/v<=zero)then ! fast acceptance
    exit
  elseif(m*log(x)-log(v)+v-one<=zero)then ! accept
    exit
  endif
enddo
gamdev_gkm2=a*v
! End procedure here
endfunction gamdev_gkm2
!----------------------------------------------------
subroutine invGamdev_s(alpha,beta,invgamdev,err,message)
! Purpose: Generates numbers from the Inverse Gamma Distribution,
!          with shape alpha and scale beta.
! Method:  Employs the fast acceptance-rejection methods for gamma
!          deviate generation (with real shape parameter);
! NB: alpha and beta are not optional
use utilities_dmsl_kit,only:zero
implicit none
! Dummies
real(mrk),intent(in)::alpha,beta
real(mrk),intent(out)::invgamdev
integer(mik),intent(out)::err
character(*),intent(out)::message
! Locales
real(mrk)::x
! Start procedure
err=0; message="invGamdev_s"
if(alpha<=zero.or.beta<=zero)then ! Check alpha,beta>0
  err=-10; message="invGamdev_s/fatal/alpha,beta<0"
  invgamdev=-hugeRE; return
endif
call gamdev(alpha,gamdev=x,err=err,message=message) ! simplicity itself
if(err/=0.or.x<=zero)then ! some type of error
  err=-10; message="invGamdev_s/fatal/&"//message
  invgamdev=-hugeRE; return
endif
invgamdev=beta/x
! End procedure
endsubroutine invGamdev_s
!----------------------------------------------------
subroutine invChi2dev_s(ndeg,invchi2,err)
! Purpose: Chi^2 sample with integer ndeg degress of freedom
! NB: rely on integer truncation for odd ndeg, ndeg/2=(ndeg-1)/2
use utilities_dmsl_kit,only:zero,half
implicit none
! Dummies
real(mrk),intent(in)::ndeg
real(mrk),intent(out)::invchi2
integer(mik),intent(out)::err
! Locals
real(mrk)::hn
character(1)::jmessage
! Start procedure here
if(ndeg<=zero)then
  err=-100; invchi2=-hugeRE
else
  hn=half*ndeg
  call invGamdev(hn,half,invchi2,err,jmessage)
endif
! End procedure here
endsubroutine invChi2dev_s
!----------------------------------------------------
subroutine invscChi2dev_s(ndeg,s2,invchi2,err)
! Purpose: Chi^2 sample with integer ndeg degress of freedom
! and scale factor s2>0 (Note that s^2 is provided, not s).
! Comments:
! 1. Old comment circa 2002-2003: rely on integer truncation for odd ndeg, ndeg/2=(ndeg-1)/2
!    [DK not sure if this is still applicable]
use utilities_dmsl_kit,only:zero,half
implicit none
! Dummies
real(mrk),intent(in)::ndeg,s2
real(mrk),intent(out)::invchi2
integer(mik),intent(out)::err
! Locals
real(mrk)::hn
character(1)::jmessage ! junk variables
! Start procedure here
if(ndeg<=zero)then  ! negative degrees of freedom, bail out
  invchi2=-hugeRE; err=-100
else                ! use inverse Gamma sampler
  if(s2>zero)then
    hn=half*ndeg
    call invGamdev(hn,hn*s2,invchi2,err,jmessage)
  else
    invchi2=-hugeRE; err=-200
  endif
endif
! End procedure here
endsubroutine invscChi2dev_s
!----------------------------------------------------
subroutine expdev_s(harvest)
! Purpose: returns in harvest an exponentially distributed +ve
! random deviate using the inverse transform method
! NB: assumes uniran is never 0 or 1 exactly
! Comments: for machines with fast arithmetic, log(.) appears
! competitive with the acceptance-rejection EA algorithm coded in
! Fortran-95. When considering much easier parallelisation, the
! inverse transform method looks even more attractive.
implicit none
! Dummies
real(mrk), intent(out) :: harvest
! Start procedure
call uniran(harvest); harvest=-log(harvest)
! End procedure
endsubroutine expdev_s
!----------------------------------------------------
subroutine expdevb_s(beta,harvest,err)
! Purpose: returns in harvest an exponentially distributed +ve
! random deviate with scale beta.
! NB: assumes uniran is never 0 or 1 exactly
use utilities_dmsl_kit,only:zero
implicit none
! Dummies
real(mrk), intent(in) :: beta
real(mrk), intent(out) :: harvest
integer(mik), intent(out) :: err
! Start procedure
if(beta>zero)then
  call uniran(harvest); harvest=-log(harvest)*beta
else
  err=-100; harvest=-hugeRE
endif
! End procedure
endsubroutine expdevb_s
!----------------------------------------------------
subroutine expdev_v(harvest)
! Purpose: returns in harvest an exponentially distributed +ve
! random deviate.
implicit none
! Dummies
real(mrk), dimension(:), intent(out) :: harvest
! Start procedure
call uniran(harvest); harvest=-log(harvest)
! End procedure
endsubroutine expdev_v
!----------------------------------------------------
subroutine expdevb_v(beta,harvest,err)
! Purpose: returns in harvest an exponentially distributed +ve
! random deviate with scale beta.
use utilities_dmsl_kit,only:zero
implicit none
! Dummies
real(mrk), dimension(:), intent(in) :: beta
real(mrk), dimension(:), intent(out) :: harvest
integer(mik), intent(out) :: err
! Start procedure
if(any(beta<=zero))then
  err=-100; harvest=-hugeRE
else
  call uniran(harvest); harvest=-log(harvest)*beta
endif
! End procedure
endsubroutine expdevb_v
!----------------------------------------------------
function poidev2_s(lambda)
! Purpose: returns as a floating point real a sample from
! the Poisson distribution of mean lambda.
! Algorithms: ITR (inverse transform) for lambda<12
!             PRUA (Fishman,1996) for lambda>=12
! Comments: PRUA is somewhat slower than PD (the fastest Poisson
! generator, Fishman,1996), it is considerable simpler, hence
! preferable.
use utilities_dmsl_kit,only:zero,half,one,two,three,four,&
  sqrt2,sqrt3,sqrtEbase,factln
implicit none
! Dummies
real(mrk), intent(in) :: lambda
real(mrk) :: poidev2_s
! Locals
real(mrk),parameter :: lambdacrit=12
real(mrk),parameter::d0=7._mrk,d1=two*sqrt2/sqrtEbase,&
  d2=three-two*sqrt3/sqrtEbase    ! parameters
! exact constants: d0=7.0(9 decimals),d1=2*sqrt(2/e),d2=3-2*sqrt(3/e)
real(mrk),SAVE::lamold=-hugeRE, d3,d4,d5,d6,d7,d8,d9  ! setup constants
real(mrk) :: x, y, z, t, k
! Start procedure here
if(lambda<=zero)then ! illegal input: return junk
  poidev2_s=-hugeRE
elseif(lambda<lambdacrit)then  ! use ITR algorithm
  if(lambda/=lamold)then  ! mean has changed, redo setuo
    d3=exp(-lambda);  lamold=lambda
  endif
  x=d3; y=x; z=zero
  call uniran(t)
  do; if(t<=x) exit               ! do while accumulating discrete cdf
    z=z+one; y=y*lambda/z; x=x+y  ! "integrating/inverting discrete pdf"
  enddo
  poidev2_s=z
else  ! use fast PRUA algorithm
  if(lambda/=lamold)then  ! mean has changed, redo setup
    d3=lambda+half; d4=sqrt(d3); d5=d1*d4+d2; d6=log(lambda)
    d7=aint(lambda,mrk); d8=aint(d3+d0*(d4+one),mrk)
    d9=d6*d7-factln(int(d7,mik)); lamold=lambda
  endif
  do  ! fast acceptance/rejection sampling algorithm
    call uniran(x); call uniran(y)
    z=d3+d5*(y-half)/x
    if(z<zero.or.z>=d8) cycle  ! prelim. test 1: fast rejection
    k=aint(z,mrk); t=k*d6-factln(int(k,mik))-d9
    if(t>=x*(four-x)-three)then ! prelim. test 2: fast acceptance
      exit
    elseif(x*(x-t)>=one)then ! prelim. test 3, fast rejection
      cycle
    elseif(two*log(x)<=t)then ! main test: accept
      exit
    endif
  enddo
  poidev2_s=k
endif
! End procedure
endfunction poidev2_s
!----------------------------------------------------
function bnldev2_s(n,pp)
! Purpose: samples from the Binomial distribution using:
! * ITR algorithm (discrete inverse transform) for n*min(p,1-p)<=12
! * fast BRUA algorithm (Fishman,1996) otherwise (bounded runtime).
! Comments: ln[1!]=0, not 1, silly D.
! NB: Fishman, 1996, appears to be missing undoing the symmetry map.
use utilities_dmsl_kit,only:zero,half,one,two,three,four,sqrt2,sqrt3,sqrtEbase,factln
implicit none
! Dummies
integer(mik), intent(in) :: n
real(mrk), intent(in) :: pp
real(mrk) :: bnldev2_s
! Locals
real(mrk),parameter:: eps=1.e-14_mrk
real(mrk),parameter::d0=7._mrk,d1=two*sqrt2/sqrtEbase,&
  d2=three-two*sqrt3/sqrtEbase   ! parameters
! exact constants: d0=7.0(9 decimals),d1=2*sqrt(2/e),d2=3-2*sqrt(3/e)
real(mrk)::d3,d4,d5,d6,d7,d8,d9,d10,d11,d12  ! setup constants
real(mrk) :: nr, x, y, z, t, w
! Start procedure here
nr=real(n,mrk)
if(n<0.or.pp<zero.or.pp>one+eps)then ! GIGO algorithm!
  bnldev2_s=-hugeRE
elseif(n==0.or.pp<eps.or.pp>one)then  ! Bin(n,pp)==0 by definition
  bnldev2_s=zero            ! if n=0 or pp=0,1
elseif(nr*min(pp,(one-pp))<=12._mrk)then ! use ITR scheme
  t=min(pp,one-pp)              ! symmetry transformation
  x=(one-t)**n; y=x; z=zero
  call uniran(d3)
  do; if(d3<=x) exit  ! do while accumulating discrete cdf
      z=z + one  ! "integrating/inverting discrete pdf"
      y=y * (nr-z+one)*t/(z-z*t)
      x=x + y
  enddo
  bnldev2_s=merge(z,nr-z,pp<=half) ! undo transformation
elseif(nr*min(pp,(one-pp))>=one)then ! use fast BRUA algorithm
  d3=min(pp,one-pp); d4=nr*d3+half
  d5=one-d3; d6=sqrt(nr*d3*d5+half); d7=d1*d6+d2
  d8=d3/d5; d9=log(d8); d10=aint(real(n+1,mrk)*d3,mrk)
  d11=real(min(n+1,floor(d4+d0*d6,mik)),mrk)
  d12=factln(int(d10,mik))+factln(int(n-d10,mik))
  do      ! acceptance/rejection
    call uniran(x); call uniran(y)
    w=d4+d7*(y-half)/x
    if(w<zero.or.w>=d11) cycle      ! fast rejection
    z=aint(w,mrk)
    t=(z-d10)*d9+d12-factln(int(z,mik))-factln(int(nr-z,mik))
    if(x*(four-x)-three<=t)then   ! fast acceptance
      exit
    elseif(x*(x-t)>=one)then      ! fast rejection
      cycle
    elseif(two*log(x)<=t)then     ! accept
      exit
    endif
  enddo
  bnldev2_s=anint(merge(z,nr-z,pp<=0.5),mrk) ! patched by Dmitri
endif
! End procedure here
endfunction bnldev2_s
!----------------------------------------------------
subroutine negbindev_s(alpha,beta,nbdev,err,message)
! Purpose: Samples from the negative binomial distribution using the
! gamma->poisson combination of Gelman et al (1995).
! IN: shape alpha, scale (NOT inverse) beta
! NB: Gelman et al. use bet=inverse scale. Here it is scale.
! Also note NEGBIN scheme of Fishman,1996; NEGBIN uses different arguments,
! NegBin(r,p), with p=1-pr(success);
use utilities_dmsl_kit,only:zero
implicit none
! Dummies
real(mrk), intent(in) :: alpha, beta
real(mrk), intent(out) :: nbdev
integer(mik), intent(out) :: err
character(len=*), intent(out) :: message
! Locals
real(mrk) :: x
! Start procedure
if(alpha<=zero.or.beta<=zero)then
  err=-10; message="negbindev_s/fatal:alpha,beta<0"
  nbdev=-hugeRE; return
else
  call gamdev(alpha,beta,x,err,message);nbdev=poidev(x) ! Gelman's scheme
  err=0; message="negbindev_s:ok"
endif
! End procedure
endsubroutine negbindev_s
!----------------------------------------------------
subroutine znewton(evalFunc,dataIN,dataOUT,x1,x2,f1,f2,tolX,tolF,xscale,fscale,itmax,&
  xroot,froot,fcalls,err,message)
! Purpose: Safeguarged Newton-Raphson-Bisection for solution of scalar nonlinear
! equation f(x)=0.
! Input
!   evalFunc = subroutine for f=f(x) and df=df/dx(x)
!   x1,x2   = solution brackets (inclusive)
!   [f1,f2] = (optional) endpoint values
!   tolX    = scaled tolerance on x
!   tolF    = scaled tolerance of f(x)
!   xscale  = typical scale of x
!   fscale  = typical scale of f
!   itmax   = maximum number of iterations
! Output
!   xroot   = solution f(xroot)~0
!   froot   = f(xroot)
!   fcalls  = number of evalFunc calls
!   err     = error diagnostic
!   message = algorithm performance
! Comments
! - Modified from Press et al. 1996, Numerical Recipes in F-90.
! - See also powerhouse routine zSuperSolver() below if exact derivatives are
!   simple to evaluate or if order of convergence important.
use types_dmsl_kit,only:data_ricz_type
use utilities_dmsl_kit,only:zero,half,two
implicit none
! Dummies
type(data_ricz_type),intent(in),optional::dataIN
type(data_ricz_type),intent(inout),optional::dataOUT
real(mrk),intent(in)::x1,x2,tolX,tolF,xscale,fscale
real(mrk),intent(in),optional::f1,f2
integer(mik),intent(in)::itmax
real(mrk),intent(out)::xroot,froot
integer(mik),intent(out)::fcalls
integer(mik),intent(out)::err
character(*),intent(out)::message
interface
  subroutine evalFunc(dataIN,dataOUT,x,feas,fx,dfdxV,err,message)
  use kinds_dmsl_kit
  use types_dmsl_kit,only:data_ricz_type
  implicit none
  type(data_ricz_type),intent(in),optional::dataIN
  type(data_ricz_type),intent(inout),optional::dataOUT
  real(mrk),intent(in)::x
  logical(mlk),intent(out)::feas
  real(mrk),intent(out),optional::fx,dfdxV(:)
  integer(mik),intent(out)::err
  character(*),intent(out)::message
  endsubroutine evalFunc
endinterface
! Locals
integer(mik)::it
real(mrk)::df,dfV(1),dx,dxold,f,fHi,fLo,xHi,xLo,temp
logical(mlk)::feas
character(100)::msg
! Start procedure here
err=0; fcalls=0; message="znewton/ok"
if(present(f1))then
  fLo=f1
else
  call evalFunc(dataIN,dataOUT,x1,feas,fLo,err=err,message=msg)
  fcalls=fcalls+1
  if(err/=0)then
    err=100;message="f-znewton/userError[x1]/"//msg;return
  elseif(.not.feas)then
    err=150;message="f-znewton/userInfeas[x1]/"//msg;return
  endif
endif
if(abs(fLo)<=tolF*fscale)then
  xroot=x1; froot=fLo; message="znewton/ok/x1"; RETURN
endif
if(present(f2))then
  fHi=f2
else
  call evalFunc(dataIN,dataOUT,x2,feas,fHi,err=err,message=msg)
  fcalls=fcalls+1
  if(err/=0)then
    err=200;message="f-znewton/userError[x2]/"//msg;return
  elseif(.not.feas)then
    err=250;message="f-znewton/userInfeas[x2]/"//msg;return
  endif
endif
if(abs(fHi)<=tolF*fscale)then
  xroot=x2; froot=fHi; message="znewton/ok/x2"; RETURN
endif
if(fLo*fHi>zero)then
  err=10;message="f-znewton/rootNotBracketed"
  xroot=-hugeRE; froot=-hugeRe; RETURN
endif
if(fLo<zero)then ! Orient the search so that f(x1)<0
  xLo=x1; xHi=x2
else
  xHi=x1; xLo=x2
endif ! Initialize the root, the "stepsize before last," and the last step
xroot=xLo+half*(xHi-xLo); dxold=abs(xHi-xLo); dx=dxold
call evalFunc(dataIN,dataOUT,xroot,feas,f,dfV,err,msg)
fcalls=fcalls+1 ; df=dfV(1)
if(err/=0)then
  err=300;message="f-znewton/userError[it=0]/"//msg;return
elseif(.not.feas)then
  err=350;message="f-znewton/userInfeas[it=0]/"//msg;return
endif
do it=1,itmax       ! iteration loop
  if(((xroot-xHi)*df-f)*((xroot-xLo)*df-f)>zero.or.&  ! Bisect if (i) Newton jumps out; or
      abs(two*f)>abs(dxold*df))then                   !           (ii) not collapsing fast enough
    dxold=dx; dx=half*(xHi-xLo); xroot=xLo+dx
    if(xLo==xroot)then
      froot=f; EXIT ! Change in root is negligible.
    endif
  else              ! Newton-Raphson step acceptable. Take it.
    dxold=dx; dx=-f/df; temp=xroot; xroot=xroot+dx
    if(temp==xroot)then
      froot=f; EXIT
    endif
  endif
  if(abs(dx)<tolX*max(abs(xroot),xscale))then ! * Convergence check on dx
    froot=f; message="znewton/ok/tolX";  EXIT
  endif
  call evalFunc(dataIN,dataOUT,xroot,feas,f,dfV,err,msg) ! one call per iteration
  fcalls=fcalls+1; df=dfV(1)
  if(err/=0)then
    err=400;write(message,'(a,i0,a)')"f-znewton/userError[it=",it,"]/"//trim(msg);return
  elseif(.not.feas)then
    err=450;write(message,'(a,i0,a)')"f-znewton/userInfeas[it=",it,"]/"//trim(msg);return
  endif
  if(abs(f)<tolF*fscale)then                  ! * Convergence check on f(x)
    froot=f; message="znewton/ok/tolF";  EXIT
  elseif(f<zero)then                         ! Maintain the bracket on the root.
    xLo=xroot
  else
    xHi=xroot
  endif
enddo
if(it>=itmax+1)then ! exceeded iteration limit: warn user, but keep root estimate
  err=20; message="f-znewton/itmax exceeded"; xroot=xLo+half*(xHi-xLo)
  call evalFunc(dataIN,dataOUT,xroot,feas,froot,dfV,err,msg)
  fcalls=fcalls+1
  if(err/=0)then
    err=500;write(message,'(a,i0,a)')"f-znewton/userError[itmax+1=",it,"]/"//trim(msg);return
  elseif(.not.feas)then
    err=550;write(message,'(a,i0,a)')"f-znewton/userInfeas[itmax+1=",it,"]/"//trim(msg);return
  endif
endif
! End procedure here
endsubroutine znewton
!---------------------------------------
end module numerix_dmsl_kit

