module types_dmsl_kit

use kinds_dmsl_kit    ! numeric kind definitions
implicit none
private

! Publicly-available stuff
public::data_ricz_type,data_ricz_def,&
  funcXY_P1_r_type,funcXY_P1_r_def,&
  funcXY_P2_r_type,funcXY_P2_r_def,&
  vectorP1_r_type,&
  vectorA1_r_type,vectorA2_r_type,vectorA1_i_type


! needed for ricz
type funcXY_P1_r_type
  integer(mik)::n=undefIN
  real(mrk),pointer::x(:)=>null()
  real(mrk),pointer::y(:)=>null()
endtype funcXY_P1_r_type
!--
type(funcXY_P1_r_type),parameter::funcXY_P1_r_def=funcXY_P1_r_type(undefIN,null(),null())
!--
type funcXY_P2_r_type
  integer(mik)::n=undefIN
  integer(mik)::m=undefIN
  real(mrk),pointer::x(:)=>null()
  real(mrk),pointer::y(:,:)=>null()
endtype funcXY_P2_r_type
!--
type(funcXY_P2_r_type),parameter::funcXY_P2_r_def=funcXY_P2_r_type(undefIN,undefIN,null(),null())
!--
type vectorP1_r_type ! type for vector data, with pointer capacity
  real(mrk),pointer::v(:)=>null()
endtype vectorP1_r_type
!-----
type vectorA1_r_type ! type for vector data, Fortran 2000 extension
  real(mrk),allocatable::v(:)
endtype vectorA1_r_type
!--
type vectorA2_r_type ! type for rank-2 array data, Fortran 2000 extension
  real(mrk),allocatable::v(:,:)
endtype vectorA2_r_type
!-----
type vectorA1_i_type
  integer(mik),allocatable::v(:)
endtype vectorA1_i_type
!-----
! ricz all-in-one structure
integer(mik),parameter::chlen_data_ricz=len_vLongStr !100
type data_ricz_type       ! the famous ricz structure that saved dmsl ...
! data values             !DK: make sure to update SUB cleanPointers_data_ricz
  real(mrk)::rs1=undefRN  ! after adding new components. otherwise get a nice
  real(mrk)::rs2=undefRN  ! memory leak
  real(mrk)::rs3=undefRN
  real(mrk)::rs4=undefRN
  integer(mik)::is1=undefIN
  integer(mik)::is2=undefIN
  integer(mik)::is3=undefIN
  integer(mik)::is4=undefIN
  logical(mlk)::ls1=undefLG
  logical(mlk)::ls2=undefLG
  logical(mlk)::ls3=undefLG
  logical(mlk)::ls4=undefLG
  complex(mck)::zs1=undefCZ
  complex(mck)::zs2=undefCZ
  complex(mck)::zs3=undefCZ
  complex(mck)::zs4=undefCZ
  character(chlen_data_ricz)::cs1=undefCH
  character(chlen_data_ricz)::cs2=undefCH
  character(chlen_data_ricz)::cs3=undefCH
  character(chlen_data_ricz)::cs4=undefCH
  character(chlen_data_ricz)::cs5=undefCH
  character(chlen_data_ricz)::cs6=undefCH
  character(chlen_data_ricz)::cs7=undefCH
! 1D pointer components
  real(mrk),pointer::rsp1=>null()
  real(mrk),pointer::rsp2=>null()
  real(mrk),pointer::rsp3=>null()
  real(mrk),pointer::rsp4=>null()
  real(mrk),pointer::rp0(:)=>null()
  real(mrk),pointer::rp1(:)=>null()
  real(mrk),pointer::rp2(:)=>null()
  real(mrk),pointer::rp3(:)=>null()
  real(mrk),pointer::rp4(:)=>null()
  integer(mik),pointer::ip1(:)=>null()
  complex(mck),pointer::zp1(:)=>null()
  character(chlen_data_ricz),pointer::cp1(:)=>null()
! 2D pointer components
  real(mrk),pointer::rpm1(:,:)=>null()
  real(mrk),pointer::rpm2(:,:)=>null()
  real(mrk),pointer::rpm3(:,:)=>null()
  real(mrk),pointer::rpm4(:,:)=>null()
! function table components
  type(funcXY_P2_r_type)::fs1=funcXY_P2_r_def
! more scalable array components
  type(vectorP1_r_type),pointer::rpv1(:)=>null()
endtype data_ricz_type
!---
type(data_ricz_type),parameter::data_ricz_def=data_ricz_type(&
  undefRN,undefRN,undefRN,undefRN,&       ! rs
  undefIN,undefIN,undefIN,undefIN,&       ! is
  undefLG,undefLG,undefLG,undefLG,&       ! ls
  undefCZ,undefCZ,undefCZ,undefCZ,&       ! zs
  undefCH,undefCH,undefCH,undefCH,undefCH,undefCH,undefCH,&       ! cs
  null(),null(),null(),null(),&           ! rsp
  null(),null(),null(),null(),null(),&    ! rp
  null(),&                                ! ip
  null(),&                                ! zp
  null(),&                                ! cp
  null(),null(),null(),null(),&           ! rpm
  funcXY_P2_r_def,&                       ! func
  null())                                 ! rpv
!-----

end module types_dmsl_kit

