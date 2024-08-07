module utilities_dmsl_kit

use kinds_dmsl_kit    ! numeric kind definitions
implicit none
private

! Publicly-available stuff
public::&
    ! constants
    zero,one,two,three,four,half,ln2,sqrt2,sqrt3,eBase,sqrtEbase,&
    pi,twopi,lnTwoPi,lnPi,lnSqrtPi,lnSqrtTwoPi,sqrtPi,sqrt2Pi,&
    ! strings
    string_number,number_string,changeCase,replacedString,countSubstringInString,&
    ! find location
    ifirstTrueLoc,ifirstFalseLoc,huntLoc,imaxloc,trueLocIndx_f,&
    ! basic math
    assertEq,assertEqLog,quickLinInterp,linearInterp,arthsi,arthsr,bico,factln,&
    getDiffFromMean,fmatmul_mv,&
    gammaln,gammq,gammp,gamf_ser,gamf_cf,betai,betacf,&
    outerprod,getKnorm,upper_rsolv,lower_rsolv,&
    ! I/O
    getSpareUnit,getNumItemsInFile,skipLinesInFile,&
    ! misc
    quickIf,setVars,getXnonMiss,swap,cleanPointers,getDiag,unitmatrix
    

! Mathematical constants
real(mrk),parameter::zero=0._mrk,one=1._mrk,two=2._mrk,three=3._mrk,four=4._mrk,half=0.5_mrk
real(mrk),parameter::ln2=0.69314718055994530941723212145817656807550013436026_mrk
real(mrk),parameter::sqrt2=1.4142135623730950488016887242096980785696718753769_mrk
real(mrk),parameter::sqrt3=1.7320508075688772935274463415058723669428052538104_mrk
real(mrk),parameter::eBase=2.7182818284590452353602874713526624977572470937000_mrk  ! exp(one)
real(mrk),parameter::sqrtEbase=1.6487212707001281468486507878141635716537761007101_mrk
real(mrk),parameter::pi=3.1415926535897932384626433832795028841971693993751_mrk
real(mrk),parameter::twoPi=two*pi
real(mrk),parameter::lnPi=1.1447298858494001741434273513530587116472948129153_mrk
real(mrk),parameter::sqrtPi=1.7724538509055160272981674833411451827975494561224_mrk
real(mrk),parameter::sqrt2pi=sqrt2*sqrtPi
real(mrk),parameter::lnTwoPi=ln2+lnPi
real(mrk),parameter::lnSqrtPi=half*lnPi
real(mrk),parameter::lnSqrtTwoPi=half*lnTwoPi

! Overloads
!--
interface assertEq
  module procedure assertEq_2,assertEq_3,assertEq_4,assertEq_5,assertEq_v
endinterface assertEq
!--
interface assertEqLog
  module procedure lg_assertEq_2,lg_assertEq_3,lg_assertEq_4,lg_assertEq_5,lg_assertEq_v
endinterface assertEqLog
!--
interface huntLoc
  module procedure huntLoc_rsf,huntLoc_rvf
endinterface huntLoc
!--
interface imaxloc
  module procedure imaxloc_r,imaxloc_r_m,imaxloc_i,imaxloc_i_m
endinterface imaxloc
!--
interface hunt
  module procedure hunt_rs,hunt_rv
endinterface hunt
!--
interface number_string
  module procedure int_string1,real_string1,real_string_fmt,int_string_fmt,&
    real_string_fmt_al,int_string_fmt_al,&
    int_string2_v
endinterface number_string
!--
interface rbrakStr
  module procedure rbrakStr_s,rbrakStr_a1
endinterface rbrakStr
!--
interface getSpareUnit
  module procedure getSpareUnit_v,getSpareUnit_1
endinterface getSpareUnit
!--
interface quickLinInterp
  module procedure quickLinInterp_pt1
endinterface quickLinInterp
!--
interface linearInterp
  module procedure linearInterp_eng
endinterface linearInterp
!--
interface gammaln
  module procedure gammaln_s1,gammaln_s2
endinterface gammaln
!--
interface quickIf
  module procedure quickIf_r1,quickIf_i1,quickIf_log1
endinterface quickIf
!--
interface factln
  module procedure factln_2
endinterface factln
!--
interface setVars
  module procedure setVars_r,setVars_i,setVars_l,setVars_ch
endinterface setVars
!--
interface swap
  module procedure swap_r,swap_i,mask_swap_r,mask_swap_i
endinterface swap
!--
interface changeCase
  module procedure changeCase_e
endinterface changeCase
!--
interface replacedString
  module procedure replacedString_trimN,&
    fastReplacedString_trim1_s,fastReplacedString_trim1_v
endinterface replacedString
!--
interface fmatmul_mv
  module procedure fmatmul_m2v1
endinterface fmatmul_mv
!--
interface cleanPointers
  module procedure cleanPointers_x1y2,&
    cleanPointers_funcXY_P1_r,cleanPointers_funcXY_P1_r1,cleanPointers_funcXY_P1_r2,&
    cleanPointers_funcXY_P2_r,cleanPointers_funcXY_P2_r1,cleanPointers_funcXY_P2_r2,&
    cleanPointers_data_ricz
endinterface cleanPointers
!--
interface getDiag
  module procedure getDiag_r2,getDiag_r3
endinterface getDiag
!--
interface unitmatrix
  module procedure unitmatrix_rr,unitmatrix_ii
endinterface unitmatrix
!--
interface getKnorm
  module procedure getKnorm_rv,getKnorm_rm
endinterface getKnorm
!--
interface outerprod
  module procedure outerprod_r,outerprod_r1
endinterface outerprod
!--
interface upper_rsolv
  module procedure upper_rsolv1
endinterface upper_rsolv
!--
interface lower_rsolv
  module procedure lower_rsolv1,lower_rsolv2
endinterface lower_rsolv
!--

! Other
integer(mik),parameter::len_Number_String=50  ! length of number_string character function
character(*),parameter::blankCH=" ",enullCH=""
integer(mik),parameter::lowerToUpperASCII=ichar("A")-ichar("a")
integer(mik),parameter::upperToLowerASCII=ichar("a")-ichar("A")

contains

!----------------------------------------------------
elemental subroutine string_number(string,num,err)
! Purpose: Overloaded for default format (will also attempt
! real read and round to produce integer result
implicit none
! dummies
character(*),intent(in)::string
integer(mik),intent(out)::num
integer(mik),intent(out)::err
! locals
real(mrk)::rtemp
! Start procedure here
read(string,'(bn,i80)',iostat=err)num  ! read as integer
if(err/=EXIT_SUCCESS)then  ! try reading as real and truncating
  read(string,'(bn,f80.0)',iostat=err)rtemp
  if(err==0)num=nint(rtemp,mik)
endif
! End procedure here
endsubroutine string_number
!----------------------------------------------------
elemental function int_string1(num,crit,nFig)
! Purpose: Performs internal write statement for integer values using smart method
!          where if num< crit then the string uses integer number format
!          where if num=>crit then the string uses scientific real number format
!          and nFig is the number of significant figures if using scientific format
implicit none
! dummies
integer(mik),intent(in)::num
integer(mik),intent(in),optional::crit
integer(mik),intent(in),optional::nFig
character(len=len_Number_String)::int_string1
! Start procedure here
integer(mik),parameter::critDef=10000000,nFigDef=9
integer(mik)::crit0,nFig0
! Start procedure here
if(present(crit))then;crit0=crit
else;                 crit0=critDef; endif
if(abs(num)<crit0) then   ! use integer notation
  write(int_string1,'(i0)')num
else                  ! use scientific notation
  if(present(nFig))then;nFig0=nFig
  else;                 nFig0=nFigDef
  endif
  int_string1=real_string_g(real(num,mrk),zero,nFig0)
endif
! End procedure here
endfunction int_string1
!---------------------------------------
elemental function real_string1(num,crit,nFig)
! Purpose: Flexible interface for real_string_g with optional arguments
implicit none
! dummies
real(mrk),intent(in)::num
real(mrk),intent(in),optional::crit
integer(mik),intent(in),optional::nFig
character(len=len_Number_String)::real_string1
! locals
real(mrk),parameter::defCrit=zero   ! default critical
integer(mik),parameter::defFig=3    ! default significant figures
! Start procedure here
if(present(crit).and.present(nFig)) then
  real_string1=real_string_g(num,crit,nFig)
elseif(present(crit).and..not.present(nFig)) then
  real_string1=real_string_g(num,crit,defFig)
elseif(.not.present(crit).and.present(nFig)) then
  real_string1=real_string_g(num,defCrit,nFig)
else  !if(.not.present(crit).and..not.present(nFig)) then
  real_string1=real_string_g(num,defCrit,defFig)
endif
! End procedure here
endfunction real_string1
!---------------------------------------
elemental function real_string_fmt(num,fmt)
! Purpose: Formatted conversion REAL->CHAR.
! Programmer: Dmitri Kavetski
! Kreated: circa 2000, Newcastle-on-Hunter
! History: 10 Mai 2015 AD, Foxhole d'Lyon - support for "*" and "#"
! ---
! Usage
! 1. "fmt" can be either a standard Fortran format string,
!    or take special DMSL-defined values:
!       "*" = full precision exponential untrimmed
!       "#" = full precision adaptive trimmed
! ---
implicit none
! dummies
real(mrk),intent(in)::num
character(*),intent(in)::fmt
character(len_Number_String)::real_string_fmt
! Start procedure here
selectcase(fmt)
case("*")
  write(real_string_fmt,rbrakStr(fmtFullPrecEs(num)))num
case("#")   ! adaptive format
  real_string_fmt=real_string_adapt(num)
case default
  write(real_string_fmt,fmt)num
endselect
! End procedure here
endfunction real_string_fmt
!---------------------------------------
elemental function int_string_fmt(num,fmt)
! Purpose: Overloaded for default format
implicit none
! dummies
integer(mik),intent(in)::num
character(*),intent(in)::fmt
character(len_Number_String)::int_string_fmt
! Start procedure here
write(int_string_fmt,fmt)num
! End procedure here
endfunction int_string_fmt
!---------------------------------------
elemental function real_string_fmt_al(num,fmt,adjl)
! Purpose: Overloaded for default format with left-shift
implicit none
! dummies
real(mrk),intent(in)::num
character(*),intent(in)::fmt
logical(mlk),intent(in)::adjl
character(len_Number_String)::real_string_fmt_al
! Start procedure here
real_string_fmt_al=real_string_fmt(num,fmt)
if(adjl)real_string_fmt_al=adjustl(real_string_fmt_al)
! End procedure here
endfunction real_string_fmt_al
!---------------------------------------
elemental function int_string_fmt_al(num,fmt,adjl)
! Purpose: Overloaded for default format with left-shift
implicit none
! dummies
integer(mik),intent(in)::num
character(*),intent(in)::fmt
logical(mlk),intent(in)::adjl
character(len_Number_String)::int_string_fmt_al
! Start procedure here
write(int_string_fmt_al,fmt)num
if(adjl)int_string_fmt_al=adjustl(int_string_fmt_al)
! End procedure here
endfunction int_string_fmt_al
!---------------------------------------
pure function int_string2_v(num,crit,nFig,lenS)
! Purpose: overloaded with ability specify length of result
implicit none
! dummies
integer(mik),intent(in)::num(:)
integer(mik),intent(in),optional::crit
integer(mik),intent(in),optional::nFig
integer(mik),intent(in)::lenS
character(len=lenS)::int_string2_v(size(num))
! locals
integer(mik),parameter::critDef=100000
integer(mik)::crit0,i,n
! Start procedure here
if(present(crit))then;crit0=crit
else;                 crit0=critDef; endif
n=size(num)
do i=1,n
  if(abs(num(i))<crit0)then ! use integer notation
    write(int_string2_v(i),'(i0)')num(i)
  else                      ! use scientific notation
    int_string2_v(i)=real_string_g(real(num(i),mrk),zero,nFig)
  endif
enddo
! End procedure here
endfunction int_string2_v
!---------------------------------------
pure function ifirstTrueLoc(mask)
! Purpose: returns first .true. ocurrence in a logical array
! Comment: not quite a very elegant implementation.
!          Q: does parallelism justify the means?
implicit none
logical(mlk),intent(in)::mask(:)
integer(mik)::ifirstTrueLoc
integer(mik)::i
!integer(mik),dimension(1)::firstTrueloc
! Start procedure here
!firstTrueloc=maxloc(merge(1,0,mask)); ifirstloc=firstTrueloc(1)
!if(.not.mask(ifirstTrueLoc))ifirstTrueLoc=size(mask)+1
do i=1,size(mask) ! serial implementation
  if(mask(i))exit
enddo
ifirstTrueLoc=i
! End procedure here
endfunction ifirstTrueLoc
!----------------------------------------------------
pure function ifirstFalseLoc(mask)
! Purpose: returns first .false. ocurrence in a logical array
implicit none
logical(mlk),intent(in)::mask(:)
integer(mik)::ifirstFalseLoc
integer(mik)::i
!integer(mik),dimension(1)::firstFalseLoc
! Start procedure here
!firstFalseLoc=minloc(merge(1,0,mask)); ifirstloc=firstFalseLoc(1)
!if(.not.mask(ifirstFalseLoc)) ifirstFalseLoc=size(mask)+1
do i=1,size(mask) ! serial implementation
  if(.not.mask(i))exit
enddo
ifirstFalseLoc=i
! End procedure here
endfunction ifirstFalseLoc
!----------------------------------------------------
pure subroutine assertEq_2(n1,n2,equal,n)
! Purpose: checks that n1-n2 are equal
implicit none
integer(mik),intent(in)::n1,n2
logical(mlk),intent(out)::equal
integer(mik),intent(out)::n
! Start procedure here
if(n1==n2)then
  n=n1;  equal=.true.
else
  equal=.false.
endif
! End procedure here
endsubroutine assertEq_2
!----------------------------------------------------
pure subroutine assertEq_3(n1,n2,n3,equal,n)
! Purpose: checks that n1-n3 are equal
implicit none
integer(mik),intent(in)::n1,n2,n3
logical(mlk),intent(out)::equal
integer(mik),intent(out)::n
! Start procedure here
if(n1==n2.and.n1==n3)then
  n=n1;  equal=.true.
else
  equal=.false.
endif
! End procedure here
endsubroutine assertEq_3
!----------------------------------------------------
pure subroutine assertEq_4(n1,n2,n3,n4,equal,n)
! Purpose: checks that n1-n4 are equal
implicit none
integer(mik),intent(in)::n1,n2,n3,n4
logical(mlk),intent(out)::equal
integer(mik),intent(out)::n
! Start procedure here
if(n1==n2.and.n1==n3.and.n1==n4)then
  n=n1;  equal=.true.
else
  equal=.false.
endif
! End procedure here
endsubroutine assertEq_4
!----------------------------------------------------
pure subroutine assertEq_5(n1,n2,n3,n4,n5,equal,n)
! Purpose: checks that n1-n5 are equal
implicit none
integer(mik),intent(in)::n1,n2,n3,n4,n5
logical(mlk),intent(out)::equal
integer(mik),intent(out)::n
! Start procedure here
if(n1==n2.and.n1==n3.and.n1==n4.and.n1==n5)then
  n=n1;  equal=.true.
else
  equal=.false.
endif
! End procedure here
endsubroutine assertEq_5
!----------------------------------------------------
pure subroutine assertEq_v(iarr,equal,n)
! Purpose: checks that all elements of iarr are equal
implicit none
integer(mik),dimension(:),intent(in)::iarr
logical(mlk),intent(out)::equal
integer(mik),intent(out)::n
! locals
integer(mik)::sizearr
! Start procedure here
sizearr=size(iarr)
selectcase(sizearr)
case(1)
  equal=.true.;   n=iarr(1)
case(2:)
  if(all(iarr(2:)==iarr(1)))then
    equal=.true.; n=iarr(1)
  else
    equal=.false.;n=undefIN
  endif
case default
  equal=.false.;  n=undefIN
endselect
! End procedure here
endsubroutine assertEq_v
!----------------------------------------------------
pure function lg_assertEq_2(n1,n2)
! Purpose: function,true if n1-n2 are equal
implicit none
integer(mik),intent(in)::n1,n2
logical(mlk)::lg_assertEq_2
! Start procedure here
lg_assertEq_2=n1==n2
! End procedure here
endfunction lg_assertEq_2
!----------------------------------------------------
pure function lg_assertEq_3(n1,n2,n3)
! Purpose: function,true if n1-n3 are equal
implicit none
integer(mik),intent(in)::n1,n2,n3
logical(mlk)::lg_assertEq_3
! Start procedure here
lg_assertEq_3=n1==n2.and.n1==n3
! End procedure here
endfunction lg_assertEq_3
!----------------------------------------------------
pure function lg_assertEq_4(n1,n2,n3,n4)
! Purpose: function,true if n1-n4 are equal
implicit none
integer(mik),intent(in)::n1,n2,n3,n4
logical(mlk)::lg_assertEq_4
! Start procedure here
lg_assertEq_4=n1==n2.and.n1==n3.and.n1==n4
! End procedure here
endfunction lg_assertEq_4
!----------------------------------------------------
pure function lg_assertEq_5(n1,n2,n3,n4,n5)
! Purpose: function,true if n1-n5 are equal
implicit none
integer(mik),intent(in)::n1,n2,n3,n4,n5
logical(mlk)::lg_assertEq_5
! Start procedure here
lg_assertEq_5=n1==n2.and.n1==n3.and.n1==n4.and.n1==n5
! End procedure here
endfunction lg_assertEq_5
!----------------------------------------------------
pure function lg_assertEq_v(iarr)
! Purpose: function,true if all elements of iarr are equal
implicit none
integer(mik),dimension(:),intent(in)::iarr
logical(mlk)::lg_assertEq_v
! Start procedure here
lg_assertEq_v=all(iarr(2:)==iarr(1))
! End procedure here
endfunction lg_assertEq_v
!----------------------------------------------------
elemental function quickLinInterp_pt1(x0,f0,x1,f1,xp)
! Purpose: Robust linear interpolation of (x0,f0) and (x1,f1) at xp.
! If x0~x1 then returns (f0+f1)/2
! Programmer: Dmitri Kavetski
! Comments:
! 1. Preserves the monotonicity of the data in floating point arithmetic,
!    eg, if x0<xp<x1 and f0<f1, then f0<fp<f1.
implicit none
! dummies
real(mrk),intent(in)::x0,f0,x1,f1,xp
real(mrk)::quickLinInterp_pt1
! locals
real(mrk)::p
! real(mrk),parameter::tol=2._mrk*epsRe
integer(mik),parameter::nulps=2
! Start procedure here
if(fpCompareSame(a=x0,b=x1,nulps=nulps))then
! if(same(x0,x1,tol))then
  p=half*(f0+f1)
else
  p=(xp-x0)/(x1-x0)
endif
quickLinInterp_pt1=(one-p)*f0+p*f1
! End procedure here
endfunction quickLinInterp_pt1
!----------------------------------------------------
pure subroutine linearInterp_eng(x,xx,yy,guessI,extrapA,extrapB,y,err,message)
! Purpose: Evaluates a tabulated function (arrays xx and yy) at point x.
! Comments:
! 1. use extrap(A,B) to extrapolate if x outside list of abscissas,
!    else will assume constant function outside the known range
implicit none
real(mrk),intent(in)::xx(:),yy(:) ! tabulated function:
real(mrk),intent(in)::x           ! abscissa at which function needs evaluated
integer(mik),intent(in),optional::guessI  ! optional guess for location in table
logical(mlk),intent(in),optional::extrapA,extrapB ! true to extrapolate off the table
real(mrk),intent(out)::y          ! function value at x
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
integer(mik)::i,np
logical(mlk),parameter::extrapDef=.true.
! Start procedure here
np=size(xx)
if(np==1)then
  err=100; message="f-linearInterp/np=1"; return
elseif(size(yy)/=np)then
  err=200; message="f-linearInterp/dimError[xx,yy]"; return
else
  err=0; message="linearInterp/ok"
endif
if(present(guessI))then
  i=huntloc(xx=xx,x=x,jlo=guessI)
else
  i=huntloc(xx=xx,x=x)
endif
if(i==0)then          ! * off the start of the abscissas
  if(quickif(extrapA,extrapDef))then
    i=1               !   do extrapolation
  else
    y=yy(1);  return  !   assume const off the start
  endif
elseif(i==np)then     ! * off the end of the abscissas
  if(quickif(extrapB,extrapDef))then
    i=i-1             !   do extrapolation
  else
    y=yy(np); return  !   assume const off the end
  endif
endif
y=quickLinInterp_pt1(x0=xx(i),f0=yy(i),x1=xx(i+1),f1=yy(i+1),xp=x)
! End procedure here
endsubroutine linearInterp_eng
!----------------------------------------------------
pure function huntLoc_rsf(xx,x,jlo)
! Purpose: Function version of hunt for in-line use
implicit none
! dummies
real(mrk),intent(in)::xx(:),x
integer(mik),intent(in),optional::jlo
integer(mik)::huntLoc_rsf
! locals
integer(mik)::jj
! Start procedure here
if(present(jlo))then
  jj=jlo; call hunt_rs(xx,x,jj)
else
  jj=locate_rs(xx,x)
endif
huntLoc_rsf=jj
! End procedure here
endfunction huntLoc_rsf
!----------------------------------------------------
pure function huntLoc_rvf(xx,x,jlo)
! Purpose: Function version of hunt for in-line use with vector x/jlo
implicit none
! dummies
real(mrk),intent(in)::xx(:),x(:)
integer(mik),intent(in),optional::jlo(:)
integer(mik)::huntLoc_rvf(size(x))
! locals
integer(mik)::n1,n2,jj(size(x))
! Start procedure here
if(present(jlo))then
  n1=size(x); n2=size(jlo)
  if(n1/=n2)then
    huntLoc_rvf=-hugeInt
  else
    jj=jlo; call hunt_rv(xx,x,jj)
  endif
else
  jj=locate_rv(xx,x)
endif
huntLoc_rvf=jj
! End procedure here
endfunction huntLoc_rvf
!----------------------------------------------------
pure function arthsr(n)
! Purpose: simplest real arithmetic progression: 1.0,2.0,...n
implicit none
integer(mik),intent(in)::n
real(mrk),dimension(n)::arthsr
integer(mik)::i
! Start procedure here
forall(i=1:n)arthsr(i)=real(i,mrk)
! End procedure here
endfunction arthsr
!----------------------------------------------------
pure function arthsi(n)
! Purpose: simplest integer arithmetic progression:1,2,3...n
implicit none
integer(mik),intent(in)::n
integer(mik),dimension(n)::arthsi
integer(mik)::i
! Start procedure here
forall(i=1:n)arthsi(i)=i
! End procedure here
endfunction arthsi
!----------------------------------------------------
subroutine getSpareUnit_v(unt,err,message)
! Purpose: returns vector of unused units for file data transfer.
! Unit will be in the range 7->2.2billion (see comment below).
! Comment:
! * Can not be pure as it contains the inquire function,tough life...
! * In Fortran-95,available units range from 0 to 2,147,483,640.
!   Preconnected units 5 (keyboard) and 6 (screen) can be re-connected
!   to a different logical device but to avoid silly errors this is avoided
!   in this procedure.
implicit none
! dummies
integer(mik),dimension(:),intent(out)::unt
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
integer(mik)::i,j,n
logical(mlk)::opened,xist
integer(mik),parameter::minUnits=7,maxUnits=2147483639
! Start procedure here
n=size(unt); j=1
do i=minUnits,maxUnits
  inquire(unit=i,opened=opened,exist=xist) ! check unit status
  if(.not.opened.and.xist)then ! un-opened existing unit found
    unt(j)=i; err=0; message="getSpareUnit/ok"
    j=j+1
    if(j>n)exit
  endif
enddo
if(i>maxUnits)then  ! all units in use
  unt=-1; err=-10; message="getSpareUnit/allUnitsInUse&"//&
      "&(all 2.2billion-u've goda b jokin')"
endif
! End procedure here
endsubroutine getSpareUnit_v
!---------------------------------------
subroutine getSpareUnit_1(unt,err,message)
! Purpose: overloaded for 1 unit
implicit none
! dummies
integer(mik),intent(out)::unt
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
integer(mik)::untv(1)
! Start procedure here
call getSpareUnit_v(untv,err,message)
unt=untv(1)
! End procedure here
endsubroutine getSpareUnit_1
!---------------------------------------
elemental function quickIf_r1(tsource,fsource,mask)
! Purpose: Strict version of merge using IF statements (real version)
! Reasons for existence: 1. allows r1=quickIf(optArg,default,present(optArg).and.mask)
implicit none
! dummies
real(mrk),intent(in),optional::tsource
real(mrk),intent(in)::fsource
logical(mlk),intent(in),optional::mask
real(mrk)::quickIf_r1
! locals
logical(mlk)::mask0,presT
logical(mlk),parameter::maskDef=.true.
! Start procedure here
presT=present(tsource)
if(present(mask))then;mask0=mask
else;                 mask0=maskDef;endif
if(mask0.and.presT)then
  quickIf_r1=tsource
else
  quickIf_r1=fsource
endif
! End procedure here
endfunction quickIf_r1
!----------------------------------------------------
elemental function quickIf_i1(tsource,fsource,mask)
! Purpose: integer version
implicit none
integer(mik),intent(in),optional::tsource
integer(mik),intent(in)::fsource
logical(mlk),intent(in),optional::mask
integer(mik)::quickIf_i1
! locals
logical(mlk)::mask0,presT
logical(mlk),parameter::maskDef=.true.
! Start procedure here
presT=present(tsource)
if(present(mask))then;mask0=mask
else;                 mask0=maskDef;endif
if(mask0.and.presT)then
  quickIf_i1=tsource
else
  quickIf_i1=fsource
endif
! End procedure here
endfunction quickIf_i1
!----------------------------------------------------
elemental function quickIf_log1(tsource,fsource,mask)
! Purpose: logical version.
implicit none
logical(mlk),intent(in),optional::tsource
logical(mlk),intent(in)::fsource
logical(mlk),intent(in),optional::mask
logical(mlk)::quickIf_log1
! locals
logical(mlk)::mask0,presT
logical(mlk),parameter::maskDef=.true.
! Start procedure here
presT=present(tsource)
if(present(mask))then;mask0=mask
else;                 mask0=maskDef;endif
if(mask0.and.presT)then
  quickIf_log1=tsource
else
  quickIf_log1=fsource
endif
! End procedure here
endfunction quickIf_log1
!----------------------------------------------------
elemental function bico(n,k)
! Purpose: Returns the binomial coefficient (n,k)
implicit none
! dummies
integer(mik),intent(in)::n,k
real(mrk)::bico
! Start procedure here
bico=anint(exp(factln(n)-factln(k)-factln(n-k)),mrk)
! End procedure here
endfunction bico
!----------------------------------------------------
elemental function factln_2(n)
! Purpose: fast evaluation of ln[n!],using:
! a. pretabulated values (computed using Maple 5 evalf to 25 decimal places);
! b. fast log(n!) evaluation (Hastings,1955) usign series expansion;
!       Fishman(1996) quotes Stadlober that abs. truncation error<8.e-9;
! [backup: c. use identity n!=Gamma(n+1) (used in NR,Press et al. 1996)];
! NB: factln_2 returns large negative number for n<0
implicit none
! dummies
integer(mik),intent(in)::n
real(mrk)::factln_2
! locals
real(mrk)::nr
real(mrk),parameter::c1=0.9189385332046727417803296_mrk,& ! ln[sqrt(2*pi)]
  c2=one/12._mrk,c3=one/360._mrk
! c4=one/1260._mrk,c5=one/1680._mrk ! for higher order series
!--Table of log-factorials
integer(mik),parameter::ntab=50,ntab1=ntab+1 ! size of table
real(mrk),dimension(0:ntab),parameter::faclogs=(/ &     ! ln(k!) with k:
  zero,zero,0.6931471805599453094172321_mrk, &              ! 0,1,2
  1.791759469228055000812477_mrk,3.178053830347945619646942_mrk,& ! 3,4
  4.787491742782045994247701_mrk,6.579251212010100995060178_mrk,& ! 5,6
  8.525161361065414300165531_mrk,10.60460290274525022841723_mrk,& ! 7,8
  12.80182748008146961120772_mrk,15.10441257307551529522571_mrk,& ! 9,10
  17.50230784587388583928765_mrk,19.98721449566188614951736_mrk,& ! 11,12
  22.55216385312342288557085_mrk,25.19122118273868150009343_mrk,& ! 13,14
  27.89927138384089156608944_mrk,30.67186010608067280375837_mrk,& ! 15,16
  33.50507345013688888400790_mrk,36.39544520803305357621562_mrk,& ! 17,18
  39.33988418719949403622465_mrk,42.33561646075348502965988_mrk,& ! 19,20
  45.38013889847690802616047_mrk,48.47118135183522387963965_mrk,& ! 21,22
  51.60667556776437357044640_mrk,54.78472939811231919009334_mrk,& ! 23,24
  58.00360522298051993929486_mrk,61.26170176100200198476558_mrk,& ! 25,26
  64.55753862700633105895132_mrk,67.88974313718153498289113_mrk,& ! 27,28
  71.25703896716800901007441_mrk,74.65823634883016438548764_mrk,& ! 29,30
  78.09222355331531063141681_mrk,81.55795945611503717850297_mrk,& ! 31,32
  85.05446701758151741396016_mrk,88.58082754219767880362692_mrk,& ! 33,34
  92.13617560368709248333304_mrk,95.71969454214320248495799_mrk,& ! 35,36
  99.33061245478742692932609_mrk,102.9681986145138126987523_mrk,& ! 37,38
  106.6317602606434591262011_mrk,110.3206397147573954290535_mrk,& ! 39,40
  114.0342117814617032329203_mrk,117.7718813997450715388381_mrk,& ! 41,42
  121.5330815154386339623110_mrk,125.3172711493568951252074_mrk,& ! 43,44
  129.1239336391272148825986_mrk,132.9525750356163098828226_mrk,& ! 45,46
  136.8027226373263684696436_mrk,140.6739236482342593987077_mrk,& ! 47,48
  144.5657439463448860089184_mrk,148.4777669517730320675372_mrk/) ! 49,50
!--End Table of log-factorials
! Start procedure here
selectcase(n)
case(0:ntab)     ! use tabulated values
  factln_2=faclogs(n)
case(ntab1:)     ! compute inline,2 options
  nr=real(n,mrk)     ! use fast series expansion algorithm scheme
  factln_2=c1+(nr+half)*log(nr)-nr+c2/nr-c3/nr**3 !+c4/nr**5+c5/nr**7
!  factln2=gammaln(real(n+1,mrk)) ! or use identity: n!=Gamma(n+1)
case default    ! return junk for n<0.
  factln_2=-hugeRE
end select
! End procedure here
endfunction factln_2
!----------------------------------------------------
pure subroutine setVars_r(val,&
   a1, a2, a3, a4, a5, a6, a7, a8, a9,a10,&
  a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,&
  a21,a22,a23,a24,a25,a26,a27,a28,a29,a30)
! Purpose: Sets the supplied variables to 'val' values.
! Comments:
! 1. Useful in trapping un-initialised variables
implicit none
! dummies
real(mrk),intent(in)::val
real(mrk),intent(out),optional:: a1, a2, a3, a4, a5, a6, a7, a8, a9,a10
real(mrk),intent(out),optional::a11,a12,a13,a14,a15,a16,a17,a18,a19,a20
real(mrk),intent(out),optional::a21,a22,a23,a24,a25,a26,a27,a28,a29,a30
! Start procedure here
if(present(a1 ))a1 =val
if(present(a2 ))a2 =val
if(present(a3 ))a3 =val
if(present(a4 ))a4 =val
if(present(a5 ))a5 =val
if(present(a6 ))a6 =val
if(present(a7 ))a7 =val
if(present(a8 ))a8 =val
if(present(a9 ))a9 =val
if(present(a10))a10=val
if(present(a11))a11=val
if(present(a12))a12=val
if(present(a13))a13=val
if(present(a14))a14=val
if(present(a15))a15=val
if(present(a16))a16=val
if(present(a17))a17=val
if(present(a18))a18=val
if(present(a19))a19=val
if(present(a20))a20=val
if(present(a21))a21=val
if(present(a22))a22=val
if(present(a23))a23=val
if(present(a24))a24=val
if(present(a25))a25=val
if(present(a26))a26=val
if(present(a27))a27=val
if(present(a28))a28=val
if(present(a29))a29=val
if(present(a30))a30=val
! End procedure here
endsubroutine setVars_r
!----------------------------------------------------
pure subroutine setVars_i(val,&
   a1, a2, a3, a4, a5, a6, a7, a8, a9,a10,&
  a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,&
  a21,a22,a23,a24,a25,a26,a27,a28,a29,a30)
! Purpose: Overloaded for integer variables.
implicit none
! dummies
integer(mik),intent(in)::val
integer(mik),intent(out),optional:: a1, a2, a3, a4, a5, a6, a7, a8, a9,a10
integer(mik),intent(out),optional::a11,a12,a13,a14,a15,a16,a17,a18,a19,a20
integer(mik),intent(out),optional::a21,a22,a23,a24,a25,a26,a27,a28,a29,a30
! Start procedure here
if(present(a1 ))a1 =val
if(present(a2 ))a2 =val
if(present(a3 ))a3 =val
if(present(a4 ))a4 =val
if(present(a5 ))a5 =val
if(present(a6 ))a6 =val
if(present(a7 ))a7 =val
if(present(a8 ))a8 =val
if(present(a9 ))a9 =val
if(present(a10))a10=val
if(present(a11))a11=val
if(present(a12))a12=val
if(present(a13))a13=val
if(present(a14))a14=val
if(present(a15))a15=val
if(present(a16))a16=val
if(present(a17))a17=val
if(present(a18))a18=val
if(present(a19))a19=val
if(present(a20))a20=val
if(present(a21))a21=val
if(present(a22))a22=val
if(present(a23))a23=val
if(present(a24))a24=val
if(present(a25))a25=val
if(present(a26))a26=val
if(present(a27))a27=val
if(present(a28))a28=val
if(present(a29))a29=val
if(present(a30))a30=val
! End procedure here
endsubroutine setVars_i
!----------------------------------------------------
pure subroutine setVars_l(val,&
   a1, a2, a3, a4, a5, a6, a7, a8, a9,a10,&
  a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,&
  a21,a22,a23,a24,a25,a26,a27,a28,a29,a30)
! Purpose: Overloaded for logical variables.
implicit none
! dummies
logical(mlk),intent(in)::val
logical(mlk),intent(out),optional:: a1, a2, a3, a4, a5, a6, a7, a8, a9,a10
logical(mlk),intent(out),optional::a11,a12,a13,a14,a15,a16,a17,a18,a19,a20
logical(mlk),intent(out),optional::a21,a22,a23,a24,a25,a26,a27,a28,a29,a30
! Start procedure here
if(present(a1 ))a1 =val
if(present(a2 ))a2 =val
if(present(a3 ))a3 =val
if(present(a4 ))a4 =val
if(present(a5 ))a5 =val
if(present(a6 ))a6 =val
if(present(a7 ))a7 =val
if(present(a8 ))a8 =val
if(present(a9 ))a9 =val
if(present(a10))a10=val
if(present(a11))a11=val
if(present(a12))a12=val
if(present(a13))a13=val
if(present(a14))a14=val
if(present(a15))a15=val
if(present(a16))a16=val
if(present(a17))a17=val
if(present(a18))a18=val
if(present(a19))a19=val
if(present(a20))a20=val
if(present(a21))a21=val
if(present(a22))a22=val
if(present(a23))a23=val
if(present(a24))a24=val
if(present(a25))a25=val
if(present(a26))a26=val
if(present(a27))a27=val
if(present(a28))a28=val
if(present(a29))a29=val
if(present(a30))a30=val
! End procedure here
endsubroutine setVars_l
!----------------------------------------------------
pure subroutine setVars_ch(val,&
   a1, a2, a3, a4, a5, a6, a7, a8, a9,a10,&
  a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,&
  a21,a22,a23,a24,a25,a26,a27,a28,a29,a30)
! Purpose: Overloaded for character variables.
implicit none
! dummies
character(*),intent(in)::val
character(*),intent(out),optional:: a1, a2, a3, a4, a5, a6, a7, a8, a9,a10
character(*),intent(out),optional::a11,a12,a13,a14,a15,a16,a17,a18,a19,a20
character(*),intent(out),optional::a21,a22,a23,a24,a25,a26,a27,a28,a29,a30
! Start procedure here
if(present(a1 ))a1 =val
if(present(a2 ))a2 =val
if(present(a3 ))a3 =val
if(present(a4 ))a4 =val
if(present(a5 ))a5 =val
if(present(a6 ))a6 =val
if(present(a7 ))a7 =val
if(present(a8 ))a8 =val
if(present(a9 ))a9 =val
if(present(a10))a10=val
if(present(a11))a11=val
if(present(a12))a12=val
if(present(a13))a13=val
if(present(a14))a14=val
if(present(a15))a15=val
if(present(a16))a16=val
if(present(a17))a17=val
if(present(a18))a18=val
if(present(a19))a19=val
if(present(a20))a20=val
if(present(a21))a21=val
if(present(a22))a22=val
if(present(a23))a23=val
if(present(a24))a24=val
if(present(a25))a25=val
if(present(a26))a26=val
if(present(a27))a27=val
if(present(a28))a28=val
if(present(a29))a29=val
if(present(a30))a30=val
! End procedure here
endsubroutine setVars_ch
!----------------------------------------------------
pure subroutine getDiffFromMean(x,nonmiss,mean,xdiff)
! Purpose: computes packed difference from mean.
! Programmer: Dmitri Kavetski
! Kreated: 16 June 2012 AD, Sonnenthal
implicit none
! dummies
real(mrk),intent(in)::x(:),mean
logical(mlk),intent(in),optional::nonmiss(:)
real(mrk),intent(out)::xdiff(:)
! locals
integer(mik)::i,k,ntot
logical(mlk)::haveNmiss
! Start procedure here
haveNmiss=presentNonZero_log1(nonmiss); ntot=size(x)
if(haveNmiss)then
  k=0; do i=1,ntot
    if(nonmiss(i))then
      k=k+1; xdiff(k)=x(i)-mean
    endif
  enddo
else
  xdiff=x-mean
endif
! End procedure here
endsubroutine getDiffFromMean
!----------------------------------------------------
pure subroutine getXnonMiss(x,nonmiss,haveNmiss,nx,nOk,nMiss,mean,err,message)
! Purpose: Processes missing data situation.
! Programmer: Dmitri Kavetski
! Kreated: 16 June 2012 AD, Sonnenthal
implicit none
! dummies
real(mrk),intent(in)::x(:)
logical(mlk),intent(in),optional::nonmiss(:)
logical(mlk),intent(out)::haveNmiss
integer(mik),intent(out),optional::nx,nOk,nMiss
real(mrk),intent(out),optional::mean
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
character(*),parameter::procnam="getXnonMiss"
integer(mik)::sizeNM,nxL,nOkL,k
! Start procedure here
if(present(mean))mean=undefRN
nxL=size(x); haveNmiss=presentNonZero_log1(nonmiss)
if(haveNmiss)then
  sizeNM=size(nonmiss)
  if(nxL/=sizeNM)then
    write(message,'(a,i0,a,i0,a)')&
      "f-"//procnam//"dimError[nxL=",nxL,"/=sizeNM=",sizeNM,"]"
    err=100; return
  endif
  nOkL=count(nonmiss)
else
  nOkL=nxL
endif
if(present(nx))nx=nxL; if(present(nOk))nOk=nOkL; if(present(nMiss))nMiss=nxL-nOkL
if    (nxL==0)then     ! empty array supplied
  err=100; message="f-"//procnam//"/empty"
elseif(nOkL==0)then        ! all data is missing
  err=200; message="f-"//procnam//"/allMissing"
elseif(nxL==1)then     ! single element
  if(present(mean))mean=x(1)
  err=-100; write(message,'(a,i0,a)')"f-"//procnam//"/dimIssue[size(x)=",nxL,"]"
elseif(nOkL==1)then        ! single non-missing value
  if(present(mean))then
    k=iFirstFalseLoc(nonmiss); mean=x(k)
  endif
  err=-200; write(message,'(a,i0,a)')"f-"//procnam//"/just1nonmiss[k=",k,"]"
else                    ! more usual scenario
  err=0 !; message=procnam//"/ok"
endif
! End procedure here
endsubroutine getXnonMiss
!----------------------------------------------------
elemental subroutine swap_r(a,b)
! Purpose: swaps 2 real numbers
implicit none
real(mrk),intent(inout)::a,b
real(mrk)::temp
! Start procedure here
temp=b; b=a; a=temp
! End procedure here
endsubroutine swap_r
!----------------------------------------------------
elemental subroutine swap_i(a,b)
! Purpose: swaps 2 integer numbers
implicit none
integer(mik),intent(inout)::a,b
integer(mik)::temp
! Start procedure here
temp=b; b=a; a=temp
! End procedure here
endsubroutine swap_i
!----------------------------------------------------
elemental subroutine mask_swap_r(a,b,mask)
! Purpose: masked swap 2 real numbers
implicit none
real(mrk),intent(inout)::a,b
logical(mlk),intent(in)::mask
real(mrk)::temp
! Start procedure here
if(mask)then
  temp=b; b=a; a=temp
endif
! End procedure here
endsubroutine mask_swap_r
!----------------------------------------------------
elemental subroutine mask_swap_i(a,b,mask)
! Purpose: masked swap 2 integer numbers
implicit none
integer(mik),intent(inout)::a,b
logical(mlk),intent(in)::mask
integer(mik)::temp
! Start procedure here
if(mask)then
  temp=b; b=a; a=temp
endif
! End procedure here
endsubroutine mask_swap_i
!----------------------------------------------------
elemental function changeCase_e(string,what)
! Purpose: Changes the case of a character string from
! lower <-> upper or toggle (each letter individually).
implicit none
! Dummy arguments
character(len=*),intent(in)::string,what
character(len=len(string))::changeCase_e
! Local variables
integer(mik)::i,m,lenTrim
! Start procedure here
changeCase_e=string
selectcase(what)
case("u","U","upper","UPPER")
  forall(i=1:len_trim(string),"a"<=string(i:i).and.string(i:i)<="z") &
    changeCase_e(i:i)=char(ichar(string(i:i))+lowerToUpperASCII)
case("l","L","lower","LOWER")
  forall(i=1:len_trim(string),"A"<=string(i:i).and.string(i:i)<="Z") &
    changeCase_e(i:i)=char(ichar(string(i:i))+upperToLowerASCII)
case("t","T","toggle","TOGGLE")
  lenTrim=len_trim(string)
  do i=1,lenTrim
    m=ichar(string(i:i))
    selectcase(m)
    case(ichar("A"):ichar("Z")) ! make upper-case lower-case
      changeCase_e(i:i)=char(m+upperToLowerASCII)
    case(ichar("a"):ichar("z")) ! make lower-case upper-case
      changeCase_e(i:i)=char(m+lowerToUpperASCII)
    endselect
  enddo
endselect
! End procedure here
endfunction changeCase_e
!----------------------------------------------------
subroutine getNumItemsInFile(unt,preRwnd,nskip,nitems,postPos,jchar,err,message)
! Purpose: Loops through a datafile to determine number of items in the file.
! Programmer: Dmitri Kavetski
! ---
! INPUT
!   unt     = file unit to be used
!   preRwnd = .true. to request the file to be pre-rewinded
!   nskip   = number of lines to skip before counting items
!   postPos = after counting items, position file at line 'postPos', usually=nskip
! OUTPUT
!   nitems  = number of items in file
!   jchar   = last line read
! ---
! Comments
! 1. File assumed to be "just opened"
implicit none
! dummies
integer(mik),intent(in)::unt,nskip
logical(mlk),intent(in)::preRwnd
integer(mik),intent(out)::nitems
integer(mik),intent(in)::postPos
character(*),intent(out)::jchar
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
character(*),parameter::procnam="getNumItemsInFile"
integer(mik)::i,jerr
! Start procedure here
jchar=procnam//"/fresh"
call skipLinesInFile(unt=unt,preRwnd=preRwnd,nskip=nskip,&
  lastSkippedLine=jchar,err=err,message=message)
if(err/=EXIT_SUCCESS)then
  err=50;message="f-"//procnam//"/&"//message;return
endif
i=1;do      ! loop through data of interest (after the skip)
  read(unt,'(a)',iostat=jerr)jchar
  selectcase(jerr)
  case(:-1) ! end-of-file
    nitems=i-1; err=0; exit
  case(+1:) ! other error
    nitems=-i; err=50; exit
  case default
    i=i+1
  endselect
enddo
if(err/=EXIT_SUCCESS)then
  write(message,'(a,i0,a)')"f-"//procnam//"/readError[ibad=",nitems,"]"
else
  write(message,'(a,i0,a)')procnam//"/ok[nitem=",nitems,"]"
endif
selectcase(postPos)
case(1:)
  call skipLinesInFile(unt=unt,preRwnd=.true.,nskip=postPos,&
                       lastSkippedLine=jchar,err=err,message=message)
  if(err/=EXIT_SUCCESS)then
    err=50;message="f-"//procnam//"/bug?/strange?/&"//message;return
  endif
case(0)
  rewind(unt)
endselect
! End procedure here
endsubroutine getNumItemsInFile
!----------------------------------------------------
subroutine skipLinesInFile(unt,preRwnd,nskip,lastSkippedLine,err,message)
! Purpose: Skips nskip lines in a file. if(preRwnd) will rewind file before skipping lines.
! ---
! Programmer: Dmitri Kavetski
! Kreated: circa 2009
! History: 7 Sept 2011 AD (added support for nskip<0)
! ---
! Comments
! 1. if pos<0 then will skip backwards
! ---
implicit none
! dummies
integer(mik),intent(in)::unt,nskip
logical(mlk),intent(in),optional::preRwnd
character(*),intent(out),optional::lastSkippedLine
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
character(*),parameter::procnam="skipLinesInFile"
logical(mlk),parameter::preRwndDef=.false.
logical(mlk)::preRwnd0
character(len_stdStrD)::jstr
integer(mik)::i
! Start procedure here
err=0; message=procnam//"/ok"; preRwnd0=quickif(preRwnd,preRwndDef)
if(preRwnd0)then
  rewind(unt,iostat=err)
  if(err/=EXIT_SUCCESS)then
    err=10;message="f-"//procnam//"/rewindError";return
  endif
endif
selectcase(nskip)
case(1:)        ! ** skip forward
  do i=1,nskip
    if(present(lastSkippedLine))then
      read(unt,'(a)',iostat=err)lastSkippedLine
    else
      read(unt,'(a)',iostat=err)jstr
    endif
    selectcase(err)
    case(:-1)   !    end of file
      write(message,'(a,i0,a,i0,a)')"f-"//procnam//"/nskip[",nskip,"]>eof[",i,"]"
      err=20; exit
    case(+1:)   !    read error
      write(message,'(a,i0,a)')"f-"//procnam//"/strange/readError[iskip=",i,"]"
      err=30; exit
    endselect
  enddo
case(:-1)       ! ** skip backwards
  if(present(lastSkippedLine))lastSkippedLine=undefCH
  do i=1,abs(nskip)
    backspace(unt,iostat=err)
    if(err/=EXIT_SUCCESS)then
      write(message,'(a,i0,a)')"f-"//procnam//"/strange/backSpaceError[i=",i,"]"
      err=30; exit
    endif
  enddo
endselect
! End procedure here
endsubroutine skipLinesInFile
!----------------------------------------------------
elemental subroutine replaceSubstrings(oldstring,olds,news,newstring,nrep,err,message)
! Purpose: Replaces each occurence of olds in string with news
! String is then padded with blanks to the original length.
! Comments:
! * note that the entire olds,news (not trim(olds,news)) are replaced.
! * if olds="" then inserts news after each character of oldstring
! * if news="" then removes every ocurrence of olds in oldstring
! * if olds=news="" then does nothing
! * if oldString="" blanks newString
implicit none
! Dummy Arguments
character(*),intent(in)::oldstring,olds,news
character(*),intent(out)::newstring
integer(mik),intent(out),optional::nrep
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
character(*),parameter::procnam="replaceSubstrings"
integer(mik)::lenOldStr,lenTrimOldStr,lenNewStr,lenOlds,lenNews,&
  startScan,startPaste,offset,endPaste,pos,s0
! Start procedure here
lenOldStr=len(oldstring); lenTrimOldStr=len_trim(oldstring); lenNewStr=len(newstring)
lenolds=len(olds); lennews=len(news); offset=0
err=0; message=procnam//"/ok"; if(present(nrep))nrep=0
if(lenTrimOldStr==0)then  ! nothing to modify - return blank
  newString=blankCH
  return
elseif(lenNews==0.and.lenOlds==0)then ! do nothing
  newString=oldString
  return
elseif(lenOlds==0)then  ! insert news between each character of oldString
  offset=1
elseif(lenNewStr==0)then  ! error
  err=10; message="f-"//procnam//"/len(newString)=0"
  return
endif
newString=blankCH ! pre-blank string
! standard replacement of nonzero strings
startScan=1; startPaste=startScan
do  ! scan for next occurence of olds.
  pos=index(oldstring(startScan:lenTrimOldStr),olds)
  if(pos==0)exit  ! no other occurences
  s0=startScan+pos-2+offset ! replace olds with news
  if(s0>lenTrimOldStr)return
  endPaste=startPaste+s0-startScan+lenNews
  if(endPaste>lenNewStr)then
    newString(startPaste:)=oldString(startScan:s0)//news
    err=20; message="f-"//procnam//"/len(newString)TooSmall"
    return
  endif
  newString(startPaste:endPaste)=oldString(startScan:s0)//news
  if(present(nrep))nrep=nrep+1
  startScan=s0+lenOlds+1
  startPaste=endPaste+1
enddo
!DK: simpler just to pre-blank the string
if(startScan<=lenTrimOldStr)then  ! finish up
  endPaste=startPaste+lenTrimOldStr-startScan
  if(endPaste>lenNewStr)then
    newString(startPaste:)=oldString(startScan:endPaste)//news
    err=20; message="f-"//procnam//"/len(newString)TooSmall"
    return
  endif
!   newString(endPaste:)=blankCH  !DK: this is buggy, so just pre-blank it earlier
  newString(startPaste:endPaste)=oldString(startScan:lenTrimOldStr)
endif
! End procedure here
endsubroutine replaceSubstrings
!----------------------------------------------------
pure function replacedString_trimN(string,olds,news,trimS)result(res)
! Purpose: Returns string with every occurence of olds replaced with news.
implicit none
! dummies
character(*),intent(in)::string,olds,news
logical(mlk),intent(in)::trimS
character(merge(len_trim(string),len(string),trimS)+&
          countSubstringInString(string,olds)*(len(news)-len(olds)))::res
! locals
character(*),parameter::procnam="replacedString_trimN"
integer(mik)::nrep,err
character(2*len_trim(string))::message
! Start procedure here
call replaceSubstrings(oldstring=string,olds=olds,news=news,&
  newstring=res,nrep=nrep,err=err,message=message)
if(err/=EXIT_SUCCESS)res="f-"//procnam//"/&"//message
! End procedure here
endfunction replacedString_trimN
!----------------------------------------------------
pure function fastReplacedString_trim1_s(string,olds,news)result(res)
! Purpose: Overloaded for speed: trimmed string replacing single character 'old'
implicit none
! dummies
character(*),intent(in)::string,olds,news
character(len_trim(string))::res
! locals
character(*),parameter::procnam="fastReplacedString_trim1_s"
integer(mik)::nrep,err
character(1)::message ! for speed
! Start procedure here
if(len(olds)>1)then
  err=100; message="f-"//procnam//"/len[olds]>1"
else
  call replaceSubstrings(oldstring=string,olds=olds,news=news,&
    newstring=res,nrep=nrep,err=err,message=message)
endif
if(err/=EXIT_SUCCESS)res=message
! End procedure here
endfunction fastReplacedString_trim1_s
!----------------------------------------------------
pure function fastReplacedString_trim1_v(string,olds,news)result(res)
! Purpose: Overloaded for array olds(:)
implicit none
! dummies
character(*),intent(in)::string,olds(:),news
character(len_trim(string))::res
! locals
character(*),parameter::procnam="fastReplacedString_trim1_v"
integer(mik)::i,n
! Start procedure here
n=size(olds); res=string
do i=1,n
  res=fastReplacedString_trim1_s(res,olds(i),news)
enddo
! End procedure here
endfunction fastReplacedString_trim1_v
!---------------------------------------
elemental function countSubstringInString(string,substring)
! Purpose: counts the number of occurences of substring in string
! Programmer: Dmitri Kavetski
! Kreated: circa 2008
! History:
!   28 aug 2011 AD, Callaghan (debugging argument parsing w Tom the Hammer)
! Comments:
! 1. Checks for full substring, rather than trim(substring), to
!    allow searching for blanks, etc.
implicit none
! dummies
character(*),intent(in)::string
character(*),intent(in)::substring
integer(mik)::countSubstringInString
! locals
integer(mik)::i,lenSub,lenStr,nSub
integer(mik)::err
character(1)::message
! Start procedure here
lenSub=len(substring); lenStr=len_trim(string)
if(lenSub>lenStr)then                             ! substring longer than string
  countSubstringInString=0; return
elseif(lenSub==lenStr.and.string==substring)then  ! exact match
  countSubstringInString=1; return
endif
selectcase(lenSub)
case(1)       ! * trivial case searching for single character
  nSub=0
  do i=1,lenStr
    if(string(i:i)==substring)then
      nSub=nSub+1
    endif
  enddo
case default  ! * slightly trickier case of searching for substring
! The next subroutine provides the full functionality while searching
! for delimiters, but returning the actual vectors is optional 
  call parseStringIntoVector_delS(string=string,delim=substring,&
    narr=nSub,err=err,message=message)
  selectcase(nSub)
  case(1:)      ! adjust in case string did not end with the delimiter
    if(string(lenStr-lenSub+1:lenStr)/=substring)nSub=nSub-1
  case(0)       ! blank string?
  case default  ! this seems unusual
  endselect
endselect
countSubstringInString=nSub
! End procedure here
endfunction countSubstringInString
!----------------------------------------------------
elemental subroutine getGamma_s(x,lnG,DlnG)
! Purpose: Evaluates the natural logarithm of the Gamma function
! ln[Gamma(x)] and its logarithmic derivative (di-gamma function) for x>0.
! Precision: not quite double,|e|<2.e-10,[Nr-77,p.206].
implicit none
! dummies
real(mrk),intent(in)::x
real(mrk),intent(out),optional::lnG,DlnG
! locals
integer(mik),parameter::nCoef=6
real(mrk)::z1,z2,lnZ1,tmp,tmp2,dtmp,dtmp2,sumC,ar(nCoef)
real(mrk),parameter::c1=5.5_mrk,c2=half,&
  c3=2.5066282746310005_mrk,c4=1.000000000190015_mrk
real(mrk),dimension(nCoef),parameter::coef=         &
  (/76.18009172947146_mrk,    -86.50532032941677_mrk,&  
    24.01409824083091_mrk,     -1.231739572450155_mrk,&
     0.1208650973866179e-2_mrk,-0.5395239384953e-5_mrk/)
! Start procedure
if(x<=zero)then         ! illegal arguments
  lnG=-hugeRE; if(present(DlnG))DlnG=-hugeRE; return
endif                   ! common constants
z1=x+c1; z2=x+c2; lnZ1=log(z1); ar=x+arthsr(nCoef); sumC=c4+sum(coef/ar)
if(present(lnG))then    ! ln-gamma function
  tmp=z2*lnZ1-z1; tmp2=log(c3*sumC/x); lnG=tmp+tmp2
endif
if(present(DlnG))then   ! digamma function (logarithmic derivative)
  dtmp=lnZ1+z2/z1-one; dtmp2=-sum(coef/ar**2)/sumC-one/x; DlnG=dtmp+dtmp2
endif
! End procedure
endsubroutine getGamma_s
!----------------------------------------------------
elemental function gammaln_s1(x)
! Purpose: function interface for the lnGamma function
implicit none
! dummies
real(mrk),intent(in)::x
real(mrk)::gammaln_s1
! locals
real(mrk)::lnG
! Start procedure
if(x<=zero)then                 ! illegal arguments
  gammaln_s1=-hugeRE
else                            ! general case
  call getGamma_s(x=x,lnG=lnG); gammaln_s1=lnG
endif
! End procedure
endfunction gammaln_s1
!----------------------------------------------------
elemental function gammaln_s2(x,useBr)
! Purpose: function interface for the lnGamma function.
! if(useBr) allows use of fast branches (good for speed,
! but yields slightly nonuniform accuracy and small
! discontinuities.
implicit none
! dummies
real(mrk),intent(in)::x
logical(mlk),intent(in)::useBr
real(mrk)::gammaln_s2
! locals
real(mrk)::lnG
! Start procedure
if(x<=zero)then                 ! illegal arguments
  gammaln_s2=-hugeRE
elseif(x==half.and.useBr)then   ! special case: G(0.5)
  gammaln_s2=lnSqrtPi
elseif(x==one.and.useBr)then    ! special case: G(1.0)
  gammaln_s2=zero
else                            ! general case
  call getGamma_s(x=x,lnG=lnG); gammaln_s2=lnG
endif
! End procedure
endfunction gammaln_s2
!----------------------------------------------------
elemental function gammq(a,x)
! Purpose: jacket for incomplete gamma function Q(a,x) evaluation
! Selects from two algorithms:
! a. series representation,scheme gamf_ser;
! b. continued fraction representation,scheme gamf_cf;
! Returns large negative number for both illegal input and
! non-converged approximations
implicit none
! dummies
real(mrk),intent(in)::a,x
real(mrk)::gammq
! Start procedure
if(a<=zero.or.x<zero)then   ! early return
  gammq=-hugeRE; return
elseif(x<a+one)then         ! use the series representations
  gammq=one-gamf_ser(a,x)   ! and take its complement
else
  gammq=gamf_cf(a,x)        ! use continued fraction
endif
! End procedure
endfunction gammq
!----------------------------------------------------
elemental function gammp(a,x)
! Purpose: jacket for incomplete gamma function P(a,x) evaluation
! Selects from two algorithms:
! a. series representation,scheme gamf_ser;
! b. continued fraction representation,scheme gamf_cf;
! Returns large negative number for both illegal input and
! non-converged approximations
implicit none
! dummies
real(mrk),intent(in)::a,x
real(mrk)::gammp
! Start procedure
if(a<=zero.or.x<zero)then   ! early return
  gammp=-hugeRE; return
elseif(x<a+one)then         ! use the series representations
  gammp=gamf_ser(a,x)
else                        ! use continued fraction
  gammp=one-gamf_cf(a,x)    ! and take its complement
endif
! End procedure
endfunction gammp
!----------------------------------------------------
elemental function gamf_ser(a,x)
! Purpose: returns the incomplete Gamma function P(a,x)
! evaluated by its series representation.
! NB: currently,no safeguard incorporated for non-convergent cases
! the parameter itmax is changed,from 100 in NR to 1000,000
! to ensure non-convergence does not occur
implicit none
! dummies
real(mrk),intent(in)::a,x
real(mrk)::gamf_ser
! locals
integer(mik),parameter::itmax=1000000
real(mrk),parameter::eps=epsilon(x)
integer(mik)::n
real(mrk)::ap,del,summ
! Start procedure
if(x==zero)then ! early return
  gamf_ser=zero; return
endif
ap=a; summ=one/a; del=summ
do n=1,itmax
  ap=ap+one; del=del*x/ap; summ=summ+del
  if(abs(del)<abs(summ)*eps)exit  ! converged
enddo
if(n<=itmax)then
  gamf_ser=summ*exp(-x+a*log(x)-gammaln(a))
else
  gamf_ser=-hugeRE
endif
! End procedure
endfunction gamf_ser
!----------------------------------------------------
elemental function gamf_cf(a,x)
! Purpose: returns the incomplete Gamma function P(a,x)
! evaluated by its continuous fraction representation.
! Algorithm: modified Lentz' scheme,$5.2 with b0=0
! NB: currently,no safeguard incorporated for non-convergent cases
! (occurs if a too large and itmax too small).
! Currently,itmax is set to 1,000,000 (100 in NR),which
! should avoid most difficulties
implicit none
! dummies
real(mrk),intent(in)::a,x
real(mrk)::gamf_cf
! locals
integer(mik),parameter::itmax=1000000
real(mrk),parameter::eps=epsilon(x),fpmin=tiny(x)/eps
! parameters: itmax is the maximum allowed number of iterations,
! eps is the relative accuracy,fpmin is near the smallest representable
! floating point number
integer(mik)::i
real(mrk)::an,b,c,d,del,h
! Start procedure
if(x==zero)then ! early return
  gamf_cf=one; return
endif
b=x+one-a; c=one/fpmin; d=one/b; h=d ! setup modified...
do i=1,itmax     ! ...Lentz' scheme and iterate to convergence
  an=-i*(i-a); b=b+two; d=an*d+b
  if(abs(d)<fpmin)d=fpmin
  c=b+an/c
  if(abs(c)<fpmin)c=fpmin
  d=one/d; del=d*c; h=h*del
  if(abs(del-one)<=eps)exit
enddo
if(i<=itmax)then
  gamf_cf=exp(-x+a*log(x)-gammaln(a))*h
else
  gamf_cf=-hugeRE
endif
! End procedure
endfunction gamf_cf
!----------------------------------------------------
pure function fmatmul_m2v1(m,v,typeMV)
! Purpose: fast matrix multiplication for matrix * vector
implicit none
! dummies
real(mrk),intent(in)::m(:,:),v(:)
character(*),intent(in)::typeMV
real(mrk),dimension(size(v))::fmatmul_m2v1
! locals
integer(mik)::n,i
! Start procedure here
n=size(v)
selectcase(typeMV)
case default !("GG","gg")         ! 0.  General case
  fmatmul_m2v1=matmul(m,v)
case("LV","lv")                   ! 1.  L*V
  forall(i=1:n)fmatmul_m2v1(i)=dot_product(m(i,1:i),v(1:i))
case("UV","uv")                   ! 2.  U*V
  forall(i=1:n)fmatmul_m2v1(i)=dot_product(m(i,i:n),v(i:n))
case("VL","VTL","vl","vtl","VtL") ! 3.  V(t)*L
  forall(i=1:n)fmatmul_m2v1(i)=dot_product(v(i:n),m(i:n,i))
case("VU","VTU","vu","vtu","VtU") ! 4.  V(t)*U
  forall(i=1:n)fmatmul_m2v1(i)=dot_product(v(1:i),m(1:i,i))
case("LTV","ltv","LtV")           ! 5.  L(t)*V
  forall(i=1:n)fmatmul_m2v1(i)=dot_product(m(i:n,i),v(i:n))
case("UTV","utv","UtV")           ! 6.  U(t)*V
  forall(i=1:n)fmatmul_m2v1(i)=dot_product(m(1:i,i),v(1:i))
case("VLT","VTLT","vlt","vtlt","VtLt","VLt")  ! 7. V(t)*L(t)
  forall(i=1:n)fmatmul_m2v1(i)=dot_product(v(1:i),m(i,1:i))
case("VUT","VTUT","vut","vtut","VUt","VtUt")  ! 8. V(t)*U(t)
  forall(i=1:n)fmatmul_m2v1(i)=dot_product(v(i:n),m(i,i:n))
case("DV","VD","dv","vd")         ! 9.  D*V,or V*D (D=diagonal matrix)
  forall(i=1:n)fmatmul_m2v1(i)=m(i,i)*v(i)
case("SUV","suv")                 ! 10. Symmetric(in upper triangle)*V
  forall(i=1:n)fmatmul_m2v1(i)=&
    dot_product(m(i,i:n),v(i:n))+dot_product(m(1:i-1,i),v(1:i-1))
case("SLV","slv")                 ! 11. Symmetric(in lower triangle)*V
  forall(i=1:n)fmatmul_m2v1(i)=&
    dot_product(m(i,1:i),v(1:i))+dot_product(m(i+1:n,i),v(i+1:n))
case("VSU","vsu")                 ! 12. V(t)*Symmetric(in upper triangle)
  forall(i=1:n)fmatmul_m2v1(i)=&
    dot_product(v(1:i),m(1:i,i))+dot_product(v(i+1:n),m(i,i+1:n))
case("VSL","vsl")                 ! 13. V(t)*Symmetric(in lower triangle)
  forall(i=1:n)fmatmul_m2v1(i)=&
    dot_product(v(i:n),m(i:n,i))+dot_product(v(1:i-1),m(i,1:i-1))
endselect
! End procedure here
endfunction fmatmul_m2v1
!----------------------------------------------------
subroutine cleanPointers_x1y2(x,y,what,err,message)
! Purpose: Garbage collection for pointer data structures.
implicit none
! dummies
real(mrk),pointer::x(:),y(:,:)
integer(mik),intent(in)::what(:)
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
integer(mik)::nw,what0
! Start procedure here
err=0; message="cleanPointers_x1y2/ok"; nw=size(what)
selectcase(nw)
case(1,2)
case default
      err=50;message="f-cleanPointers_x1y2/dimError(what)";return
endselect
what0=what(1)
selectcase(what0)
case(0)   ! do nothing
case(+1)  ! deallocate if associated
  if(associated(x))then
    deallocate(x,stat=err)
    if(err/=EXIT_SUCCESS)then
      err=100;message="f-cleanPointers_x1y2/deallocError(x)";return
    endif
  endif
case(-1)  ! nullify if associated
  if(associated(x))nullify(x)
case default
      err=150;message="f-cleanPointers_x1y2/illegalVal(what)";return
endselect
if(nw==2)what0=what(2)
selectcase(what0)
case(0)   ! do nothing
case(+1)  ! deallocate if associated
  if(associated(y))then
    deallocate(y,stat=err)
    if(err/=EXIT_SUCCESS)then
      err=200;message="f-cleanPointers_x1y2/deallocError(y)";return
    endif
  endif
case(-1)  ! nullify if associated
  if(associated(y))nullify(y)
case default
      err=250;message="f-cleanPointers_x1y2/illegalVal(what)";return
endselect
! End procedure here
endsubroutine cleanPointers_x1y2
!----------------------------------------------------
! pure subroutine cleanPointers_funcXY_P1_r(fvar,what,err,message)
subroutine cleanPointers_funcXY_P1_r(fvar,what,err,message)
! Purpose: Garbage collection for a 'funcXY_P1_r_type' variable.
use types_dmsl_kit,only:funcXY_P1_r_type,funcXY_P1_r_def
implicit none
! dummies
type(funcXY_P1_r_type),intent(inout)::fvar
integer(mik),intent(in)::what
integer(mik),intent(out)::err
character(*),intent(out)::message
! Start procedure here
err=0; message="cleanPointers_funcXY_P1_r/ok"
selectcase(what)
case(0)   ! do nothing
case(+1)  ! deallocate if associated
  if(associated(fvar%x))then
    deallocate(fvar%x,stat=err)
    if(err/=EXIT_SUCCESS)then
      err=100;message="f-cleanPointers_funcXY_P1_r/deallocError(x)";return
    endif
  endif
  if(associated(fvar%y))then
    deallocate(fvar%y,stat=err)
    if(err/=EXIT_SUCCESS)then
      err=200;message="f-cleanPointers_funcXY_P1_r/deallocError(y)";return
    endif
  endif
  fvar=funcXY_P1_r_def
case(-1)  ! nullify if associated
  fvar=funcXY_P1_r_def
!   if(associated(fvar%x))nullify(fvar%x)
!   if(associated(fvar%y))nullify(fvar%y)
case default
      err=300;message="f-cleanPointers_funcXY_P1_r/illegalVal(what)";return
endselect
! End procedure here
endsubroutine cleanPointers_funcXY_P1_r
!----------------------------------------------------
! pure subroutine cleanPointers_funcXY_P1_r1(fvar,what,err,message)
subroutine cleanPointers_funcXY_P1_r1(fvar,what,err,message)
! Purpose: Garbage collection for a 1D array of 'funcXY_P1_r_type' variable.
use types_dmsl_kit,only:funcXY_P1_r_type
implicit none
! dummies
type(funcXY_P1_r_type),intent(inout)::fvar(:)
integer(mik),intent(in)::what
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
integer(mik)::i,n
character(100)::msg
! Start procedure here
err=0; message="cleanPointers_funcXY_P1_r1/ok"; n=size(fvar)
do i=1,n
  call cleanPointers_funcXY_P1_r(fvar(i),what,err,msg)
  if(err/=EXIT_SUCCESS)then
    write(message,'(a,i0,a)')"f-cleanPointers_funcXY_P1_r1/[i=",i,"]/&"//trim(msg)
    err=100; return
  endif
enddo
! End procedure here
endsubroutine cleanPointers_funcXY_P1_r1
!----------------------------------------------------
! pure subroutine cleanPointers_funcXY_P1_r2(fvar,what,err,message)
subroutine cleanPointers_funcXY_P1_r2(fvar,what,err,message)
! Purpose: Garbage collection for a 2D array of 'funcXY_P1_r_type' variable.
use types_dmsl_kit,only:funcXY_P1_r_type
implicit none
! dummies
type(funcXY_P1_r_type),intent(inout)::fvar(:,:)
integer(mik),intent(in)::what
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
integer(mik)::i,n,j,m
character(100)::msg
! Start procedure here
err=0; message="cleanPointers_funcXY_P1_r2/ok"; n=size(fvar,1); m=size(fvar,2)
do j=1,m;do i=1,n
  call cleanPointers_funcXY_P1_r(fvar(i,j),what,err,msg)
  if(err/=EXIT_SUCCESS)then
    write(message,'(a,i0,a,i0,a)')"f-cleanPointers_funcXY_P1_r2/[i=",i,",j=",j,"]/&"//trim(msg)
    err=100; return
  endif
enddo;enddo
! End procedure here
endsubroutine cleanPointers_funcXY_P1_r2
!----------------------------------------------------
! pure subroutine cleanPointers_funcXY_P2_r(fvar,what,err,message)
subroutine cleanPointers_funcXY_P2_r(fvar,what,err,message)
! Purpose: Garbage collection for a 'funcXY_P1_r_type' variable.
use types_dmsl_kit,only:funcXY_P2_r_type,funcXY_P2_r_def
implicit none
! dummies
type(funcXY_P2_r_type),intent(inout)::fvar
integer(mik),intent(in)::what
integer(mik),intent(out)::err
character(*),intent(out)::message
! Start procedure here
err=0; message="cleanPointers_funcXY_P2_r/ok"
selectcase(what)
case(0)   ! do nothing
case(+1)  ! deallocate if associated
  if(associated(fvar%x))then
    deallocate(fvar%x,stat=err)
    if(err/=EXIT_SUCCESS)then
      err=100;message="f-cleanPointers_funcXY_P2_r/deallocError(x)";return
    endif
  endif
  if(associated(fvar%y))then
    deallocate(fvar%y,stat=err)
    if(err/=EXIT_SUCCESS)then
      err=200;message="f-cleanPointers_funcXY_P2_r/deallocError(y)";return
    endif
  endif
  fvar=funcXY_P2_r_def
case(-1)  ! nullify if associated
  fvar=funcXY_P2_r_def
!   if(associated(fvar%x))nullify(fvar%x)
!   if(associated(fvar%y))nullify(fvar%y)
case default
      err=300;message="f-cleanPointers_funcXY_P2_r/illegalVal(what)";return
endselect
! End procedure here
endsubroutine cleanPointers_funcXY_P2_r
!----------------------------------------------------
! pure subroutine cleanPointers_funcXY_P2_r1(fvar,what,err,message)
subroutine cleanPointers_funcXY_P2_r1(fvar,what,err,message)
! Purpose: Garbage collection for a 1D array of 'funcXY_P1_r_type' variable.
use types_dmsl_kit,only:funcXY_P2_r_type
implicit none
! dummies
type(funcXY_P2_r_type),intent(inout)::fvar(:)
integer(mik),intent(in)::what
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
integer(mik)::i,n
character(100)::msg
! Start procedure here
err=0; message="cleanPointers_funcXY_P2_r1/ok"; n=size(fvar)
do i=1,n
  call cleanPointers_funcXY_P2_r(fvar(i),what,err,msg)
  if(err/=EXIT_SUCCESS)then
    write(message,'(a,i0,a)')"f-cleanPointers_funcXY_P2_r1/[i=",i,"]/&"//trim(msg)
    err=100; return
  endif
enddo
! End procedure here
endsubroutine cleanPointers_funcXY_P2_r1
!----------------------------------------------------
! pure subroutine cleanPointers_funcXY_P2_r2(fvar,what,err,message)
subroutine cleanPointers_funcXY_P2_r2(fvar,what,err,message)
! Purpose: Garbage collection for a 2D array of 'funcXY_P2_r_type' variable.
use types_dmsl_kit,only:funcXY_P2_r_type
implicit none
! dummies
type(funcXY_P2_r_type),intent(inout)::fvar(:,:)
integer(mik),intent(in)::what
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
integer(mik)::i,n,j,m
character(100)::msg
! Start procedure here
err=0; message="cleanPointers_funcXY_P2_r2/ok"; n=size(fvar,1); m=size(fvar,2)
do j=1,m;do i=1,n
  call cleanPointers_funcXY_P2_r(fvar(i,j),what,err,msg)
  if(err/=EXIT_SUCCESS)then
    write(message,'(a,i0,a,i0,a)')"f-cleanPointers_funcXY_P2_r2/[i=",i,",j=",j,"]/&"//trim(msg)
    err=100; return
  endif
enddo;enddo
! End procedure here
endsubroutine cleanPointers_funcXY_P2_r2
!----------------------------------------------------
! pure subroutine cleanPointers_data_ricz(d,what,err,message)
subroutine cleanPointers_data_ricz(d,what,err,message)
! Purpose: Garbage collection for a 'data_ricz_type' variable.
use types_dmsl_kit,only:data_ricz_type,data_ricz_def
implicit none
! dummies
type(data_ricz_type),intent(inout)::d
integer(mik),intent(in)::what
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
integer(mik)::i,n
character(100)::msg
! Start procedure here
err=0; message="cleanPointers_data_ricz/ok"
selectcase(what)
case(0)   ! do nothing
case(+1)  ! deallocate if associated
  if(associated(d%rp0))then
    deallocate(d%rp0,stat=err)
    if(err/=EXIT_SUCCESS)then
      err=110;message="f-cleanPointers_data_ricz/deallocError(rp0)";return
    endif
  endif
  if(associated(d%rp1))then
    deallocate(d%rp1,stat=err)
    if(err/=EXIT_SUCCESS)then
      err=110;message="f-cleanPointers_data_ricz/deallocError(rp1)";return
    endif
  endif
  if(associated(d%rp2))then
    deallocate(d%rp2,stat=err)
    if(err/=EXIT_SUCCESS)then
      err=120;message="f-cleanPointers_data_ricz/deallocError(rp2)";return
    endif
  endif
  if(associated(d%rp3))then
    deallocate(d%rp3,stat=err)
    if(err/=EXIT_SUCCESS)then
      err=130;message="f-cleanPointers_data_ricz/deallocError(rp3)";return
    endif
  endif
  if(associated(d%rp4))then
    deallocate(d%rp4,stat=err)
    if(err/=EXIT_SUCCESS)then
      err=140;message="f-cleanPointers_data_ricz/deallocError(rp4)";return
    endif
  endif
  if(associated(d%ip1))then
    deallocate(d%ip1,stat=err)
    if(err/=EXIT_SUCCESS)then
      err=150;message="f-cleanPointers_data_ricz/deallocError(ip1)";return
    endif
  endif
  if(associated(d%zp1))then
    deallocate(d%zp1,stat=err)
    if(err/=EXIT_SUCCESS)then
      err=160;message="f-cleanPointers_data_ricz/deallocError(zp1)";return
    endif
  endif
  if(associated(d%cp1))then
    deallocate(d%cp1,stat=err)
    if(err/=EXIT_SUCCESS)then
      err=170;message="f-cleanPointers_data_ricz/deallocError(cp1)";return
    endif
  endif
  call cleanPointers_funcXY_P2_r(d%fs1,+1,err,msg)
  if(err/=EXIT_SUCCESS)then
    err=180;message="f-cleanPointers_data_ricz/&"//msg;return
  endif
  if(associated(d%rpv1))then
    n=size(d%rpv1)
    do i=1,n
      if(associated(d%rpv1(i)%v))then
        deallocate(d%rpv1(i)%v,stat=err)
        if(err/=EXIT_SUCCESS)then
          write(message,'(a,i0,a)')"f-cleanPointers_data_ricz/deallocError(rpv1(i=",i,")%v"
          err=190; return
        endif
      endif
    enddo
    deallocate(d%rpv1,stat=err)
    if(err/=EXIT_SUCCESS)then
      err=200;message="f-cleanPointers_data_ricz/deallocError(rpv1)";return
    endif
  endif
  d=data_ricz_def
case(-1)  ! nullify if associated
  d=data_ricz_def
!   if(associated(d%rsp1))nullify(d%rsp1)
!   if(associated(d%rsp2))nullify(d%rsp2)
!   if(associated(d%rsp3))nullify(d%rsp3)
!   if(associated(d%rp0)) nullify(d%rp0)
!   if(associated(d%rp1)) nullify(d%rp1)
!   if(associated(d%rp2)) nullify(d%rp2)
!   if(associated(d%rp3)) nullify(d%rp3)
!   if(associated(d%rp4)) nullify(d%rp4)
!   if(associated(d%ip1)) nullify(d%ip1)
!   if(associated(d%zp1)) nullify(d%zp1)
!   if(associated(d%cp1)) nullify(d%cp1)
!   call cleanPointers_funcXY_P2_r(d%fs1,-1,err,msg)
!   if(err/=EXIT_SUCCESS)then
!     err=280;message="f-cleanPointers_data_ricz/&"//msg;return
!   endif
!   [placeholder for cleanPointers for TYPE 'vectorP1_r_type']
case default
    err=300;message="f-cleanPointers_data_ricz/illegalVal(what)";return
endselect
! End procedure here
endsubroutine cleanPointers_data_ricz
!----------------------------------------------------
pure function getDiag_r2(a)
! Purpose: returns the diagonal of real 2D array a
implicit none
! dummies
real(mrk),intent(in)::a(:,:)
real(mrk)::getDiag_r2(1:min(size(a,1),size(a,2)))
integer(mik)::i,n,n1,n2
! Start procedure here
n1=size(a,1); n2=size(a,2); n=min(n1,n2)
forall(i=1:n)getDiag_r2(i)=a(i,i)
! End procedure here
endfunction getDiag_r2
!----------------------------------------------------
pure function getDiag_r3(a)
! Purpose: returns the diagonal of real 3D array a
implicit none
! dummies
real(mrk),intent(in)::a(:,:,:)
real(mrk)::getDiag_r3(1:min(size(a,1),size(a,2),size(a,3)))
integer(mik)::i,n,n1,n2,n3
! Start procedure here
n1=size(a,1); n2=size(a,2); n3=size(a,3); n=min(n1,n2,n3)
forall(i=1:n)getDiag_r3(i)=a(i,i,i)
! End procedure here
endfunction getDiag_r3
!----------------------------------------------------
pure subroutine unitmatrix_rr(a)
! Purpose: overwrites a with a real unit matrix
real(mrk),intent(out)::a(:,:)
integer(mik)::i
! Start procedure here
a(:,:)=zero
forall(i=1:min(size(a,dim=1),size(a,dim=2)))a(i,i)=one
! End procedure here
endsubroutine unitmatrix_rr
!----------------------------------------------------
pure subroutine unitmatrix_ii(a)
! Purpose: overwrites a with an integer unit matrix
integer(mik),intent(out)::a(:,:)
integer(mik)::i
! Start procedure here
a(:,:)=0
forall(i=1:min(size(a,dim=1),size(a,dim=2)))a(i,i)=1
! End procedure here
endsubroutine unitmatrix_ii
!----------------------------------------------------
pure function getKnorm_rv(m,k,avg,sumOnly)
! Purpose: get the Holder "k"-norm of an "n"-vector "m",defined as
! ||v||(k)=[sum(i=1:n)(m(i)**k)]**(1/k),which,in english,
! is the k-th root of the summed k-powered element.
! ---
! Programmer: Dmitri Kavetski
! Kreated: circa 2002      - El Jurasico de Nuevocastillo
! History: 25 June 2014 AD - Eawag Dubendorf FC D36
! ---
! Ref: Golub and van Loan (1983) "Matrix Computations",p.12.
! ---
! Some common norms:
! k=1: Manhattan distance
! k=2: Euclidean norm
! k=inf: max (infinity) norm (here requested as k=0)
! ---
! NB: 1. Note optional averaging. (not classic norm)
!     2. Note 0-length vectors have all norms equal to 0.
! ---
implicit none
! dummies
real(mrk),intent(in)::m(:)
integer(mik),intent(in)::k
logical(mlk),intent(in),optional::avg,sumOnly
real(mrk)::getKnorm_rv
! locals
integer(mik)::n
real(mrk)::nr
logical(mlk),parameter::avgDef=.false.,sumOnlyDef=.false.
logical(mlk)::avg0,sumOnly0
! Start procedure here
n=size(m)
if(n<1)then   ! 0-length vectors have 0 norm
  getKnorm_rv=zero; return
endif
avg0=quickif(avg,avgDef); sumOnly0=quickif(sumOnly,sumOnlyDef)
selectcase(k)
case(0)       ! - infinity norm
  getKnorm_rv=maxval(abs(m))
case(1)       ! - Manhattan distance
  getKnorm_rv=sum(abs(m))
case(2)       ! - Euclidean (uses potentially faster dot_product operator)
  getKnorm_rv=dot_product(m,m)
case(3:)      ! - general case
  getKnorm_rv=sum(abs(m)**k)
case default  ! - illegal value (negative norm index)
  getKnorm_rv=hugeRE; return
endselect
if(avg0)then
  nr=real(n,mrk); getKnorm_rv=getKnorm_rv/nr
endif
if(.not.sumOnly0)then
  selectcase(k)
  case(2)
    getKnorm_rv=sqrt(getKnorm_rv)
  case(3:)
    getKnorm_rv=getKnorm_rv**(one/real(k,mrk))
  endselect
endif
! End procedure here
endfunction getKnorm_rv
!----------------------------------------------------
pure function getKnorm_rm(m,k,avg,sumOnly)
! Purpose: get "k"-norm of an "n"-matrix "m",defined as
! ||m||(k)=sup[x/=0]||m.x||/||x||.
! ---
! Programmer: Dmitri Kavetski
! Kreated: circa 2002      - El Jurasico de Nuevocastillo
! History: 25 June 2014 AD - Eawag Dubendorf FC D36
! ---
! References:
! * Golub and van Loan (1983) "Matrix Computations",p.14.
! * Definitions from Matlab 6.5.
! * Gill,Murray and Wright (1981) Practical Optimization,p.28.
! * Dennis and Schnabel (1996) Numerical methods for optimization,etc.,p.44.
! ---
! Norms supported by this code:
! k=1:    1-norm: max abs column sum
! k=2:    2-norm: generalised Euclidean (Frobenius) norm
! k=inf:  max (infinity) norm (here requested as k=0): max abs row sum
! k>2:    k-norm,here computed as a pseudo-norm paralleling vectors. This is
!         not the classical definition.
! ---
! NB: 1. Note optional averaging. (not classic norm)
!     2. Note 0-length vectors have all norms equal to 0.
!     3. There are different conventions for norms.
!     4. The correct 2-norm is the max eigenvalue of [m(t).m] but is
!        harder to compute,certainly not a "utility" code.
! ---
implicit none
! dummies
real(mrk),intent(in)::m(:,:)
integer(mik),intent(in)::k
logical(mlk),intent(in),optional::avg,sumOnly
real(mrk)::getKnorm_rm
! locals
integer(mik)::n
real(mrk)::nr
logical(mlk),parameter::avgDef=.false.,sumOnlyDef=.false.
logical(mlk)::avg0,sumOnly0
! Start procedure here
n=size(m)
if(n<1)then ! 0-length vectors have 0 norm
  getKnorm_rm=zero; return
endif
nr=real(n,mrk)
avg0=quickif(avg,avgDef); sumOnly0=quickif(sumOnly,sumOnlyDef)
selectcase(k)
case(0)       ! - infinity norm (maximum absolute row sum)
  getKnorm_rm=maxval(sum(abs(m),dim=2))
case(1)       ! - 1-norm        (maximum absolute column sum)
  getKnorm_rm=maxval(sum(abs(m),dim=1))
case(2)       ! - pseudo 2-norm (Frobenius norm)
  getKnorm_rm=sum(abs(m)**2)
case(3:)      ! - pseudo k-norm
  getKnorm_rm=sum(abs(m)**k)
case default  ! - illegal value (negative norm index)
  getKnorm_rm=hugeRE; return
endselect
if(avg0)then
  nr=real(n,mrk); getKnorm_rm=getKnorm_rm/nr
endif
if(.not.sumOnly0)then
  selectcase(k)
  case(2)
    getKnorm_rm=sqrt(getKnorm_rm)
  case(3:)
    getKnorm_rm=getKnorm_rm**(one/real(k,mrk))
  endselect
endif
! End procedure here
endfunction getKnorm_rm
!----------------------------------------------------
pure function outerprod_r(vc,vr)
! Purpose: evaluates the outer product of two vectors
implicit none
real(mrk),dimension(:),intent(in)::vc,vr
real(mrk),dimension(1:size(vc),1:size(vr))::outerprod_r
integer(mik)::i,j
integer(mik)::nc,nr
! Start procedure here
nc=size(vc); nr=size(vr)
forall(i=1:nc,j=1:nr)outerprod_r(i,j)=vc(i)*vr(j)
! End procedure here
endfunction outerprod_r
!----------------------------------------------------
pure function outerprod_r1(vc)
! Purpose: overloaded for self outer-product
implicit none
real(mrk),dimension(:),intent(in)::vc
real(mrk),dimension(1:size(vc),1:size(vc))::outerprod_r1
integer(mik)::i,j,n
! Start procedure here
n=size(vc)
forall(i=1:n,j=1:n,i>j)outerprod_r1(i,j)=vc(i)*vc(j)
forall(i=1:n,j=1:n,i>j)outerprod_r1(j,i)=outerprod_r1(i,j)
forall(i=1:n)outerprod_r1(i,i)=vc(i)**2
! End procedure here
endfunction outerprod_r1
!----------------------------------------------------
pure function imaxloc_r(arr)
! Purpose: returns max location in a real array,much
! more convenient than the intrinsic fortran maxloc/minloc,
! which return location as an array even for a 1D vector
implicit none
real(mrk),dimension(:),intent(in)::arr
integer(mik)::imaxloc_r
integer(mik),dimension(1)::imax
! Start procedure here
imax=maxloc(arr); imaxloc_r=imax(1)
! End procedure here
endfunction imaxloc_r
!----------------------------------------------------
pure function imaxloc_r_m(arr,mask)
! Purpose: returns masked max location in a real array
implicit none
real(mrk),dimension(:),intent(in)::arr
logical(mlk),dimension(:),intent(in)::mask
integer(mik)::imaxloc_r_m
integer(mik),dimension(1)::imax
! Start procedure here
imax=maxloc(arr,mask); imaxloc_r_m=imax(1)
! End procedure here
endfunction imaxloc_r_m
!----------------------------------------------------
pure function imaxloc_i(iarr)
! Purpose: returns max location in an integer array
implicit none
integer(mik),dimension(:),intent(in)::iarr
integer(mik)::imaxloc_i
integer(mik),dimension(1)::imax
! Start procedure here
imax=maxloc(iarr); imaxloc_i=imax(1)
! End procedure here
endfunction imaxloc_i
!----------------------------------------------------
pure function imaxloc_i_m(arr,mask)
! Purpose: returns masked max location in an integer array
implicit none
integer(mik),dimension(:),intent(in)::arr
logical(mlk),dimension(:),intent(in)::mask
integer(mik)::imaxloc_i_m
integer(mik),dimension(1)::imax
! Start procedure here
imax=maxloc(arr,mask); imaxloc_i_m=imax(1)
! End procedure here
endfunction imaxloc_i_m
!----------------------------------------------------
pure subroutine upper_rsolv1(a,d,unitD,b,x,transp,err)
! Solves the set of n linear equations R*x=b,where R is an upper triangular
! matrix stored in upper triangle of a and optionally d (diagonal).
! x is input as the right-hand-side vector of length n
! and is overwritten with the solution vector on output.
! if b provided then assumes RHS in b and puts solution into x.
! Optional transpose requests that U(t)x=b be solved.
! Optional unitD specifies that the diagonal is unity.
implicit none
! Dummies
real(mrk),intent(in)::a(:,:)
real(mrk),intent(in),optional::d(:)
logical(mlk),intent(in),optional::unitD
real(mrk),intent(inout)::x(:)
real(mrk),intent(in),optional::b(:)
logical(mlk),intent(in),optional::transp
integer(mik),intent(out)::err
! Locals
real(mrk)::dd
integer(mik)::i,n
logical(mlk)::ok,transpL,useD,doD
! Start procedure here
call assertEq(size(a,1),size(a,2),size(x),ok,n)
if(ok)then
  err=0
else
  err=100;return
endif
if(n==0)return
useD=present(d)
if(useD)then;if(size(d)/=n)then
  err=200;return
endif;endif
transpL=quickif(transp,.false.); doD=.not.quickif(unitD,.false.)
if(present(b))then    ! *** grab RHS from b
  if(size(b)/=n)then
    err=300;return
  endif
  if(transpL)then
    if(doD)then       ! * general-value-diagonal
      i=1; dd=getDD(); x(1)=b(1)/dd
      do i=2,n       ! ** forward substitution (solve transposed system)
        dd=getDD()
        x(i)=(b(i)-dot_product(a(1:i-1,i),x(1:i-1)))/dd
      enddo
    else              ! * unit-diagonal
      x(1)=b(1)
      do i=2,n
        x(i)=(b(i)-dot_product(a(1:i-1,i),x(1:i-1)))
      enddo
    endif
  else
    if(doD)then
      i=n; dd=getDD(); x(n)=b(n)/d(n)
      do i=n-1,1,-1   ! ** backward substitution
        dd=getDD()
        x(i)=(b(i)-dot_product(a(i,i+1:n),x(i+1:n)))/dd
      enddo
    else
      x(n)=b(n)
      do i=n-1,1,-1
        x(i)=(b(i)-dot_product(a(i,i+1:n),x(i+1:n)))
      enddo
    endif
  endif
else                  ! *** overwrite RHS in x
  if(transpL)then
    if(doD)then
      i=1; dd=getDD(); x(1)=x(1)/dd
      do i=2,n
        dd=getDD()
        x(i)=(x(i)-dot_product(a(1:i-1,i),x(1:i-1)))/dd
      enddo
    else
!      x(1)=x(1)
      do i=2,n
        x(i)=(x(i)-dot_product(a(1:i-1,i),x(1:i-1)))
      enddo
    endif
  else
    if(doD)then
      i=n; dd=getDD(); x(n)=x(n)/d(n)
      do i=n-1,1,-1
        dd=getDD()
        x(i)=(x(i)-dot_product(a(i,i+1:n),x(i+1:n)))/dd
      enddo
    else
!      x(n)=x(n)
      do i=n-1,1,-1
        x(i)=(x(i)-dot_product(a(i,i+1:n),x(i+1:n)))
      enddo
    endif
  endif
endif
! End procedure here
contains
!--
pure function getDD() ! macro to get diagonal pivot
implicit none
real(mrk)::getDD
  if(useD)then
    getDD=d(i)
  else
    getDD=a(i,i)
  endif
endfunction
!--
endsubroutine upper_rsolv1
!----------------------------------------------------
pure subroutine lower_rsolv1(a,d,unitD,b,x,transp,err)
! Solves the set of n linear equations L*x=b,where L is a lower triangular
! matrix stored in lower triangle of a and optionally d (diagonal).
! x is input as the right-hand-side vector of length n
! and is overwritten with the solution vector on output.
! if b provided then assumes RHS in b and puts solution into x
! Optional transpose requests that L(t)x=b be solved
! Optional unitD specifies that the diagonal is unity.
implicit none
! Dummies
real(mrk),intent(in)::a(:,:)
real(mrk),intent(in),optional::d(:)
logical(mlk),intent(in),optional::unitD
real(mrk),intent(inout)::x(:)
real(mrk),intent(in),optional::b(:)
logical(mlk),intent(in),optional::transp
integer(mik),intent(out)::err
! Locals
real(mrk)::dd
integer(mik)::i,n
logical(mlk)::ok,transpL,useD,doD
! Start procedure here
call assertEq(size(a,1),size(a,2),size(x),ok,n)
if(.not.ok)then
  err=100;return
else
  err=0
endif
if(n==0)return
useD=present(d)
if(useD)then;if(size(d)/=n)then
  err=200;return
endif;endif
transpL=quickif(transp,.false.); doD=.not.quickif(unitD,.false.)
if(present(b))then  ! allows keeping RHS
  if(size(b)/=n)then
    err=300;return
  endif
  if(transpL)then   ! solve L'x=b: essentially Ux=b,transposing references to U
    if(doD)then
      i=n; dd=getDD(); x(n)=b(n)/dd
      do i=n-1,1,-1
        dd=getDD()
        x(i)=(b(i)-dot_product(a(i+1:n,i),x(i+1:n)))/dd
      enddo
    else
      x(n)=b(n)
      do i=n-1,1,-1
        x(i)=(b(i)-dot_product(a(i+1:n,i),x(i+1:n)))
      enddo
    endif
  else
    if(doD)then
      i=1; dd=getDD(); x(1)=b(1)/dd
      do i=2,n
        dd=getDD()
        x(i)=(b(i)-dot_product(a(i,1:i-1),x(1:i-1)))/dd
      enddo
    else
      x(1)=b(1)
      do i=2,n
        x(i)=(b(i)-dot_product(a(i,1:i-1),x(1:i-1)))
      enddo
    endif
  endif
else
  if(transpL)then
    if(doD)then
      i=n; dd=getDD(); x(n)=x(n)/dd
      do i=n-1,1,-1
        dd=getDD()
        x(i)=(x(i)-dot_product(a(i+1:n,i),x(i+1:n)))/dd
      enddo
    else
!      x(n)=x(n)
      do i=n-1,1,-1
        x(i)=(x(i)-dot_product(a(i+1:n,i),x(i+1:n)))
      enddo
    endif
  else
    if(doD)then
      i=1; dd=getDD(); x(1)=x(1)/dd
      do i=2,n
        dd=getDD()
        x(i)=(x(i)-dot_product(a(i,1:i-1),x(1:i-1)))/dd
      enddo
    else
!      x(1)=x(1)
      do i=2,n
        x(i)=(x(i)-dot_product(a(i,1:i-1),x(1:i-1)))
      enddo
    endif
  endif
endif
! End procedure here
contains
!--
pure function getDD()
implicit none
real(mrk)::getDD
  if(useD)then
    getDD=d(i)
  else
    getDD=a(i,i)
  endif
endfunction
!--
endsubroutine lower_rsolv1
!----------------------------------------------------
pure SUBROUTINE lower_rsolv2(aBand,D,mBand,b,x,transp,err)
! Purpose: Lower triangular solve for a banded system.
implicit none
! dummies
REAL(mrk),intent(in)::aBand(:,:),D(:),b(:)
integer(mik),intent(in)::mBand
real(mrk),intent(inout)::x(:)
logical(mlk),intent(in)::transp
integer(mik),intent(out)::err
! locals
INTEGER(mik)::n,i,k,iD,kk,kkMax
REAL(mrk)::temp !,summ
! Start procedure here
err=0; n=size(x); iD=mBand+1
if(transp)then
  err=100; x=hugeRe ! not implemented yet
!   do i=n,1,-1
!     sum=x(i)
!     do k=i+1,n
!       summ=summ-A(k,i)*x(k)
!       temp=temp-aBand(k,i)*x(k)
!     enddo
!     x(i)=sum/D(i)
!   enddo
else
  do i=1,n  ! loop over rows in lower-triangular matrix
!     summ=b(i)
    temp=b(i)
    kkMax=(i-1)-max(1,i-mBand)+1
    k=i; do kk=1,kkMax
      k=k-1
!     kk=0
!     do k=i-1,max(1,i-mBand),-1
!       kk=kk+1
!       summ=summ-A(i,k)*x(k)
      temp=temp-aBand(i,iD-kk)*x(k)
    enddo
    if(D(i)==zero)then
      err=10; return
    endif
    x(i)=temp/D(i)
  enddo
endif
! End procedure here
ENDsubroutine lower_rsolv2
!----------------------------------------------------
pure function trueLocIndx_f(mask)
! Purpose: Returns the indices of true elements of mask
implicit none
! dummies
logical(mlk),dimension(:),intent(in)::mask
integer(mik),dimension(count(mask))::trueLocIndx_f
! Start procedure here
trueLocIndx_f=pack(arthsi(size(mask)),mask)
! End procedure here
endfunction trueLocIndx_f
!----------------------------------------------------

! --------
! private
! --------

!----------------------------------------------------
pure function fpCompareSame(a,b,nulps)
! Purpose: True if a and b agree to within nulps trailing units.
implicit none
! dummies
real(mrk),intent(in)::a,b
integer(mik),intent(in)::nulps
logical(mlk)::fpCompareSame
! locals
logical(mlk),parameter::safe=.true.
integer(mik)::fpChk
! Start procedure
call fpCompareSub(a,b,nulps,safe,icomp=fpChk)
fpCompareSame=(fpChk==0)
! End procedure here
endfunction fpCompareSame
!----------------------------------------------------
pure subroutine fpCompareSub(a,b,nulps,safe,icomp,aMinusB)
! Purpose: Compare floating points numbers accounting for FP arithmetic limitations.
! nulps is the number of trailing units that may differ even if a==b.
! Returns: +1 --> a>b, -1 --> a<b, 0 --> a==b
! Optionally returns the difference between the numbers aMinusB=a-b
implicit none
! dummies
real(mrk),intent(in)::a,b
integer(mik),intent(in)::nulps
logical(mlk),intent(in)::safe
real(mrk),intent(out),optional::aMinusB
integer(mik),intent(out)::icomp
! local
real(mrk)::difffp,deltafp,ab,spaceAB
! Start procedure
if(safe.and.nulps<0)then
  if(present(aMinusB))aMinusB=undefRN; icomp=undefIN; return
endif
difffp=a-b;ab=max(abs(a),abs(b));spaceAB=spacing(ab);deltafp=nulps*2*spaceAB
if    (difffp>+deltafp)then
  icomp=+1    ! a>b
elseif(difffp<-deltafp)then
  icomp=-1    ! a<b
else
  icomp=0     ! a==b
endif
if(present(aMinusB))aMinusB=difffp
! End procedure here
endsubroutine fpCompareSub
!----------------------------------------------------
pure subroutine hunt_rs(xx,x,jlo)
! Purpose: given array xx(1:n) and a value "x",returns "jlo" an 
! integer "j" such that "x" lies between xx(jlo) and xx(jlo+1).
! IN: jlo taken as initial guess
! Precondition: xx must be monotonic (either incr or decr).
! Error handling: jlo=0 or jlo=N indicates "x" out of range of "xx"
! Performance: uses expansive-contracting O(logN) bisection algorithm
implicit none
! dummies
real(mrk),intent(in)::xx(:)
real(mrk),intent(in)::x
integer(mik),intent(inout)::jlo
! locals
integer(mik)::n,inc,jhi,jm
logical(mlk)::ascnd     ! true if ascending table
integer(mik),parameter::nulps=2
! Start procedure here
n=size(xx); ascnd=(xx(n)>=xx(1))
if(jlo<=0.or.jlo>n)then ! user's jlo is crap. do fresh bisection
  jlo=0; jhi=n+1
else  ! try to use "jlo" to limit search
  inc=1
  if(x>=xx(jlo).eqv.ascnd)then  ! * hunt up
    do
      jhi=jlo+inc
      if(jhi>n)then ! off the upper edge of table-hunting ended
        jhi=n+1; exit
      else
        if((x<xx(jhi)).eqv.ascnd)exit
        jlo=jhi; inc=inc+inc ! not done hunting: double search area...
      endif
    enddo
  else                          ! * hunt down
    jhi=jlo
    do
      jlo=jhi-inc
      if(jlo<1)then
        jlo=0; exit ! off the lower edge of table-hunting ended
      else
        if(x>=xx(jlo).eqv.ascnd)exit
        jhi=jlo; inc=inc+inc    ! not done hunting yet double the search area
      endif
    enddo
  endif
endif  ! done hunting,value bracketed
do  ! begin final bisection phase
  if(jhi-jlo<=1)then
    if    (fpCompareSame(a=x,b=xx(1),nulps=nulps))then; jlo=1
    elseif(fpCompareSame(a=x,b=xx(n),nulps=nulps))then; jlo=n-1; endif
    exit
  else
    jm=(jhi+jlo)/2
    if((x>=xx(jm)).eqv.ascnd )then; jlo=jm
    else;                           jhi=jm; endif
  endif
enddo
! End procedure here
endsubroutine hunt_rs
!----------------------------------------------------
pure subroutine hunt_rv(xx,x,jlo)
! Purpose: Vector version for multiple x/jlo
implicit none
! dummies
real(mrk),intent(in)::xx(:),x(:)
integer(mik),intent(inout)::jlo(:)
! locals
integer(mik)::i,n1,n2
! Start procedure here
n1=size(x); n2=size(jlo)
if(n1/=n2)then
  jlo=-hugeInt
else
  forall(i=1:n1)jlo(i)=huntLoc_rsf(xx,x(i),jlo=jlo(i))
endif
! End procedure here
endsubroutine hunt_rv
!----------------------------------------------------
pure function locate_rs(xx,x)
! Purpose: Given array xx(1:n) and a value "x",locate returns an 
! integer "j" such that "x" lies between xx(j) and xx(j+1).
! Precondition: xx must be monotonic(either incr or decr).
! Error handling: j=0 or j=N indicates "x" out of range of "xx".
! Comments:
! 1. Uses contracting O(logN) bisection algorithm.
implicit none
! dummies
real(mrk),intent(in)::xx(:),x
integer(mik)::locate_rs
! locals
integer(mik)::n,jl,jm,ju
logical(mlk)::ascnd  ! true is ascending monotonic table
integer(mik),parameter::nulps=2
! Start procedure here
n=size(xx); ascnd=(xx(n)>=xx(1))
jl=0; ju=n+1  ! initialise lower and upper limits
do
  if(ju-jl<=1)exit  ! got it!
  jm=(ju+jl)/2      ! compute a midpoint
  if(ascnd.eqv.(x>=xx(jm)))then
    jl=jm ! replace either the lower limit
  else
    ju=jm ! .. or the upper limit
  endif
enddo ! set the output, being careful with endpoints
if    (fpCompareSame(a=x,b=xx(1),nulps=nulps))then
  locate_rs=1
elseif(fpCompareSame(a=x,b=xx(n),nulps=nulps))then  
  locate_rs=n-1
else
  locate_rs=jl
endif
! End procedure here
endfunction locate_rs
!----------------------------------------------------
pure function locate_rv(xx,x)
! Purpose: vector version of locate
implicit none
! dummies
real(mrk),intent(in)::xx(:),x(:)
integer(mik)::locate_rv(size(x))
! locals
integer(mik)::i
! Start procedure here
forall(i=1:size(x))locate_rv(i)=locate_rs(xx,x(i))
! End procedure here
endfunction locate_rv
!----------------------------------------------------
elemental function real_string_g(num,crit,nFig)
! Purpose: Performs internal write statement for real values using smart method
!          if num<crit then using decimal notation with nFig significant figures
!          otherwise scientific with nFig significant figures.
implicit none
! dummies
real(mrk),intent(in)::num
real(mrk),intent(in)::crit
integer(mik),intent(in)::nFig
character(len=len_Number_String)::real_string_g
! locals
integer(mik),parameter::lenFform=7,nFforms=16,nDecDef=3
character(len=lenFform),dimension(0:nFforms),parameter::formF=&
(/'(i0)   ','(f0.01)','(f0.02)','(f0.03)','(f0.04)','(f0.05)','(f0.06)',&
  '(f0.07)','(f0.08)','(f0.09)','(f0.10)','(f0.11)','(f0.12)','(f0.13)',&
  '(f0.14)','(f0.15)','(f0.16)'/)
integer(mik),parameter::lenEform=11,nEforms=16,nDigDef=3
character(len=lenEform),dimension(1:nEforms),parameter::formE=&
(/'(es08.00e3)','(es09.01e3)','(es10.02e3)','(es11.03e3)','(es12.04e3)',&
  '(es13.05e3)','(es14.06e3)','(es15.07e3)','(es16.08e3)','(es17.09e3)',&
  '(es18.10e3)','(es19.11e3)','(es20.12e3)','(es21.13e3)','(es22.14e3)',&
  '(es23.15e3)'/)
! Start procedure here
if(abs(num)<crit) then  ! ** use fX.Y notation
  selectcase(nFig) ! use specified number of significant figures
  case(0)
    write(real_string_g,formF(0)) nint(num,mik)
  case(1:nFforms)
    write(real_string_g,formF(nFig))   num
  case default
    write(real_string_g,formF(nDecDef))num
  endselect
else                    ! ** use scientific notation
  selectcase(nFig)  ! note trick to avoid the (annoying) fortran habit
  case(1:nEforms)           ! of counting "-" as a significant digit.
    write(real_string_g,formE(nFig))   abs(num)
  case default
    write(real_string_g,formE(nDigDef))abs(num)
  endselect
  real_string_g=adjustl(real_string_g)
  if(num<zero)real_string_g="-"//real_string_g
endif
! End procedure here
endfunction real_string_g
!----------------------------------------------------
elemental function real_string_adapt(num)
! Purpose: Adaptive mag-based REAL->CHAR converter
! Programmer: Dmitri Kavetski
! Kreated: 10 Mai 2015 AD, Foxhole d'Lyon
implicit none
! dummies
real(mrk),intent(in)::num
character(len_Number_String)::real_string_adapt
! locals
real(mrk),parameter::sMin=1.e-4_mrk,sMax=1.e4_mrk ! #: switch f<->es
! Start procedure here
if(sMin<abs(num).and.abs(num)<sMax)then
  write(real_string_adapt,rbrakStr(fmtFullPrecF(num)))num
else
  write(real_string_adapt,rbrakStr(fmtFullPrecEs(num)))num
endif
real_string_adapt=shortestReCH(real_string_adapt)
! End procedure here
endfunction real_string_adapt
!----------------------------------------------------
pure function rbrakStr_s(str)result(res)
! Purpose: Returns str enclosed in round brackets
! Programmer: Dmitri Kavetski, 7 Mai 2015 AD - Foxhole d'Lyon
implicit none
! dummies
character(*),intent(in)::str
character(len_trim(str)+2)::res
! Start procedure here
res="("//trim(str)//")"
! end procedure here
endfunction rbrakStr_s
!----------------------------------------------------
pure function rbrakStr_a1(str)result(res)
! Purpose: overloaded for 1D arrays.
implicit none
! dummies
character(*),intent(in)::str(:)
character(len(str)+2)::res(size(str))
! locals
integer(mik)::i,n
! Start procedure here
n=size(str); do i=1,n; res(i)='('//trim(str(i))//')'; enddo
! end procedure here
endfunction rbrakStr_a1
!----------------------------------------------------
elemental function fmtFullPrecES(a)result(res)
! Purpose: Format for full 'es'-style precision output of real of type 'a'.
! Programmer: Dmitri Kavetski, 7 Mai 2015AD, Foxhole d'Lyon
implicit none
! dummies
integer(mik),parameter::lenFmt=64
real(mrk),intent(in)::a
character(lenFmt)::res
! locals
integer(mik),parameter::xtraEs=8
integer(mik)::iprec,iexpMax,iexpMax10,w,d,e
! Start procedure here
iprec=precision(a); iexpMax=range(a)
iexpMax10=ceiling(log10(real(iexpMax,mrk)))
w=iprec+xtraEs; d=iprec; e=iexpMax10
write(res,'(3(a,i0))')"es",w,".",d,"e",e
! End procedure here
endfunction fmtFullPrecES
!----------------------------------------------------
elemental function fmtFullPrecF(a)result(res)
! Purpose: Format for full 'f'-style precision output of real of type 'a'.
! Programmer: Dmitri Kavetski, 10 Mai 2015AD, Foxhole d'Lyon
implicit none
! dummies
integer(mik),parameter::lenFmt=64
real(mrk),intent(in)::a
character(lenFmt)::res
! locals
integer(mik),parameter::xtraF=3
integer(mik)::iprec,iexpMax10,w,d
! Start procedure here
iprec=precision(a); iexpMax10=ceiling(abs(log10(abs(a))))
if(abs(a)>=one)then
  d=iprec-iexpMax10; w=iprec+xtraF
else
  d=iprec+iexpMax10; w=d+xtraF
endif
write(res,'(2(a,i0))')"f",w,".",d
! End procedure here
endfunction fmtFullPrecF
!----------------------------------------------------
elemental function shortestIntCH(a)result(res)
! Purpose: Shortest complete string representation of an integer:
!          removes leading 0's keeping the rest intact
! Programmer: Dmitri Kavetski, 7 Mai 2015 AD - Foxhole d'Lyon
! Comments:
! 1. Slow thorough check - faster DO-based implementations exist for sure
implicit none
! dummies
character(*),intent(in)::a
character(len(a))::res
! locals
character(*),parameter::procnam="shortestIntCH"
character(len(a))::str
character(1)::signCH
integer(mik)::iFirst,nWhole
! Start procedure here
str=adjustl(a)
selectcase(str(1:1))
case("-","+"); signCH=str(1:1); str=str(2:)
case default;  signCH=blankCH; endselect
nWhole=len_trim(str)
if(nWhole>0)then
  iFirst=verify(trim(str)//"x","0") ! trick to detect multiple 0's
  if(iFirst==0)then;          res=str
  elseif(iFirst>nWhole)then;  res="0"
  else;                       res=str(iFirst:);  endif
else
  res="0"
endif
if(signCH=="-".and.res/="0")res=signCH//res
! End procedure here
endfunction shortestIntCH
!----------------------------------------------------
elemental function shortestReCH(a)result(res)
! Purpose: Shortest complete string representation of a real number:
!          removes trailing decimal 0's keeping the rest intact
! Programmer: Dmitri Kavetski, 7 Mai 2015 AD - Foxhole d'Lyon
! Comments:
! 1. Slow thorough check - faster DO-based implementations exist for sure
implicit none
! dummies
integer(mik),parameter::nXtra=2 ! "." and "0" may be inserted
character(*),intent(in)::a
character(len(a)+nXtra)::res
! locals
character(*),parameter::procnam="shortestReCH"
logical(mlk),parameter::expSP=.true.,chkDcmp=.false.
integer(mik)::n,nDec,iD,iE,iZ
character(8)::expCH
character(1)::signCH
character(len(a)+nXtra)::arec,wholeCH,decCH
! Start procedure here
n=len_trim(a); iD=scan(a,"."); iE=scan(a,"Ee")
if(iE>0)then; expCH=a(iE:);  decCH=a(:iE-1)
else;         expCH=blankCH; decCH=a; endif
if(iD>0)then; wholeCH=a(:iD); decCH=decCH(iD+1:)
else;         wholeCH=decCH;  decCH=blankCH; endif
wholeCH=adjustl(wholeCH)
selectcase(wholeCH(1:1))
case("-","+"); signCH=wholeCH(1:1); wholeCH=wholeCH(2:)
case default;  signCH=blankCH; endselect
if(chkDcmp)then   ! check decomposition
  arec=trim(signCH)//trim(wholeCH)//trim(decCH)//expCH
  if(adjustl(a)/=arec)then
    res=undefCH; return
  endif
endif
nDec=len_trim(decCH)
if(nDec>0)then    ! some decimals present
  iZ=verify(trim(decCH),"0",back=.true.)
  selectcase(iZ)
  case(0);  decCH="0"
  case(1:); decCH=decCH(:iZ); endselect ! trim trailing zeroes
else              ! no decimals at all - too extreme for DK
  decCH="0"
endif
!if(lastChars(wholeCH)==".")wholeCH=trimLastChars(wholeCH)  ! expressive code
iD=len_trim(wholeCH);if(wholeCH(iD:iD)==".")wholeCH(iD:iD)=blankCH  ! fast code
wholeCH=shortestIntCH(wholeCH)
if(wholeCH=="0".and.decCH=="0")then  ! plain zero
  res="0.0"; return
endif
res=trim(wholeCH)//"."//trim(decCH)
if(signCH=="-")res=signCH//res
expCH=expCH(2:); expCH=shortestIntCH(expCH)
if(expCH/="0")then ! only add non-zero exponents
  if(expSP)then
    selectcase(expCH(1:1))
    case("-","+")
    case default; expCH="+"//expCH; endselect
  endif
  res=trim(res)//"E"//expCH ! add non-zero exponent
endif
! End procedure here
endfunction shortestReCH
!----------------------------------------------------
pure function presentNonZero_log1(a)result(res)
! Purpose: overloaded for 1D logical arrays
implicit none
! dummies
logical(mlk),intent(in),optional::a(:)
logical(mlk)::res
! locals
integer(mik)::n
! Start procedure here
res=present(a)
if(res)then
  n=size(a); if(n==0)res=.false.
endif
! End procedure here
endfunction presentNonZero_log1
!----------------------------------------------------
pure subroutine parseStringIntoVector_delS(string,delim,array,narr,narrSeek,&
  cdasIn,blankStdIn,err,message)
! Purpose: Parses trim(string) into array, according to delimitor delim
! Comments:
! 1. The entire delim is used [not trim(delim)], to allow parsing by blanks.
! 2. However, it is trim(string) [not the entire string] that is parsed.
! 3. Uses serial do loop
! 4. If len(delim)==0 then parses each character into elements of array.
! 5. If len_trim(string)=0 returns
! 6. If string starts with delim, inserts blank as first array element
! 7. If string does not end with delim, the last chars form the last element
! 8. cdasIn = consequtive delimiters as single
! 9. blankStdIn = treats consecutive blanks as delimiters (std rules)
implicit none
! dummies
character(*),intent(in)::string
character(*),intent(in)::delim
character(*),intent(out),optional::array(:)
integer(mik),intent(out)::narr              ! number of used elements in arr
integer(mik),intent(in),optional::narrSeek  ! stop after narrSeek delims are found
logical(mlk),intent(in),optional::cdasIn
logical(mlk),intent(in),optional::blankStdIn
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
character(*),parameter::procnam="parseStringIntoVector_delS"
integer(mik)::lenDelim,lenString,sizeArr,startScan,pos,s0,nBlanks,posBlank,adv
logical(mlk)::presArr,pressNarrseek,foundEnuf
logical(mlk),parameter::cdasDef=.false.,blankStdDef=.false.
logical(mlk)::cdas,blankStd,advBlank
! Start procedure here
lenDelim=len(delim);lenString=len_trim(string); presArr=present(array)
pressNarrseek=present(narrSeek); foundEnuf=.false.
cdas=quickif(cdasIn,cdasDef); blankStd=quickif(blankStdIn,blankStdDef)
if(presArr)then;sizeArr=size(array)
else;           sizeArr=hugeInt; endif
narr=0; err=0; message=trim(procnam)//"/ok"
if(lenString==0)then    ! nothing to parse
  err=-10; message="f-"//trim(procnam)//"/len_trim(string)==0"
  return
elseif(lenDelim==0)then ! special case: parse everything
  narr=min(lenString,sizeArr)
  if(pressNarrseek)narr=min(narr,abs(narrSeek))
  if(presArr)forall(pos=1:narr)array(pos)=string(pos:pos)
  if(sizeArr<lenString)then
    err=-20; message="f-"//trim(procnam)//"/size(array)tooSmall"
  endif
  return
endif
startScan=1; posBlank=0
do  ! scan for the provided delimiter
  if(blankStd)then
    nBlanks=verify(string(startScan:lenString),blankCH)-1
    startScan=startScan+nBlanks
  endif
  pos=index(string(startScan:lenString),delim); advBlank=.false.
  if(blankStd)posBlank=scan(string(startScan:lenString),blankCH)
  if(posBlank>0.and.pos>0.and.posBlank<pos-1)then ! blank comes first ...
    if(string(posBlank:pos-1)/=blankCH)advBlank=.true. ! ... and its not followed by delim
  elseif(posBlank>0.and.pos==0)then               ! just a blank (no delim at all)
    advBlank=.true.
  endif
  if(advBlank)then        ! advancing to a blank either way
    pos=posBlank; adv=1   ! the len(blankCH)
  else
    adv=lenDelim
  endif
  selectcase(pos)
  case(0)         ! no other occurences
    exit
  case(1)         ! consecutive delimiter
    if(cdas)then
      startScan=startScan+adv
      cycle
    endif
  endselect
  narr=narr+1     ! found delimiter
  if(narr>sizeArr)then
    err=-30; message="f-"//trim(procnam)//"/size(array)tooSmall"
    narr=sizeArr; goto 6666
  endif
  s0=startScan+pos-2
  if(presArr)array(narr)=string(startScan:s0) ! add to array
  startScan=s0+adv+1
  if(pressNarrseek)then   ! check if found enuf
    if(narr==abs(narrSeek))then
      foundEnuf=.true.; goto 6666
    endif
  endif
enddo
if(.not.foundEnuf.and.narr>=0.and.startScan<=lenString)then  ! finish up array if ...
  narr=narr+1   ! ... the master string did not end with a delimiter
  if(narr>sizeArr)then
    err=-40; message="f-"//trim(procnam)//"/size(array)tooSmall"
    narr=sizeArr; goto 6666
  elseif(presArr)then
    array(narr)=string(startScan:lenString)
  endif
endif
6666 continue
! End procedure here
endsubroutine parseStringIntoVector_delS
!---------------------------------------
elemental function betai(a,b,x)
! Purpose: Returns the incomplete beta function I(x)(a,b),0<=x<=1
! NB: error condition denoted via betai=-hugeRE.
implicit none
! dummies
real(mrk),intent(in)::a,b,x
real(mrk)::betai
! locals
real(mrk)::bt
! Start procedure here
if(x<zero.or.x>one)then
  betai=-hugeRE; return
elseif(x==zero.or.x==one)then
  bt=zero
else                          ! Factors in front of the continued fraction
  bt=exp(gammaln(a+b)-gammaln(a)-gammaln(b)+a*log(x)+b*log(one-x))
endif
if(x <(a+one)/(a+b+two))then   ! Use continued fraction directly
  betai=bt*betacf(a,b,x)/a
else    ! Use continued fraction after making the symmetry transformation
  betai=one-bt*betacf(b,a,one-x)/b
endif
! End procedure here
endfunction betai
!----------------------------------------------------
elemental function betacf(a,b,x)
! Purpose: Engine for betai: Evaluates continued fraction for incomplete
! beta functions using the modified Lentz's method(paragraph 5.2)
! NB: error condition denoted via betacf=-hugeRE.
implicit none
! dummies
real(mrk),intent(in)::a,b,x
real(mrk)::betacf
! locals
integer(mik),parameter::maxit=100000 ! 100 in original NR
real(mrk),parameter::eps=epsRE,fpmin=tinyRE/eps
real(mrk)::aa,c,d,del,h,qab,qam,qap
integer(mik)::m,m2
! Start procedure here
qab=a+b; qap=a+one; qam=a-one ! Factors in the coefficients (6.4.6)
c=one; d=one-qab*x/qap        ! First step of Lentz's method. 
if(abs(d)< fpmin)d=fpmin
d=one/d; h=d
do m=1,maxit
  m2=2*m; aa=m*(b-m)*x/((qam+m2)*(a+m2)); d=one+aa*d ! even step of recurrence
  if(abs(d)<fpmin)d=fpmin
  c=one+aa/c
  if(abs(c)<fpmin)c=fpmin
  d=one/d; h=h*d*c; aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2)); d=one+aa*d
  if(abs(d)<fpmin)d=fpmin                               ! odd step of recurrence
  c=one+aa/c
  if(abs(c)<fpmin)c=fpmin
  d=one/d; del=d*c; h=h*del
  if(abs(del-one)<eps)exit  ! are we done yet?!
enddo
if(m<=maxit)then
  betacf=h
else
  betacf=-hugeRE
endif
! End procedure here
endfunction betacf
!----------------------------------------------------

end module utilities_dmsl_kit