module dynamicModelLibrary

!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! WARNING !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!
! this is just an empty shell for the moment, to allow compiling BaM

use kinds_dmsl_kit    ! numeric kind definitions
implicit none
private

! Publicly-available stuff
public::GR4J_DMDL,&
  DMDL_setModel,DMDL_getModelInfo,DMDL_controlModel,DMDL_runModel

! * supported process models
integer(mik),parameter::&   ! IDs: <=0=EXTERNAL (EXE, R, MATLAB), 1-999=common, 1000-1999=CPSM's
  EXTGATE_DMDL  =   -1,&
  COMFUNK_DMDL  =  999,&    ! Common functions
! CPSM's
  EZB_DMDL      = 1000,&    ! Ben's EZ model (simplified locoSPAM)
  HBV_DMDL      = 1001,&    ! HBV model in D's implementation [see tech-note, comm B.Schaefli'07]
  FUSE_DMDL     = 1002,&    ! FUSE model of Martyn Clark, NIWA
  FLEX_DMDL     = 1003,&    ! SUPERFLEX model of Fenicia & Kavetski, CRP-GL & UoN
  GR4J_DMDL     = 1004,&    ! GR4J model of Perrin et al, CEMAGREF
  GR3P_DMDL     = 1005,&    ! GR3P model as described in Berthet' thesis, CEMAGREF
  UBA_DMDL      = 1007      ! UBA model, developed by D
!--

contains

!----------------------------------------------------
subroutine DMDL_setModel(modelID,setupCmd,chvarLibDef,err,message)
! Purpose: This routine sets the model basics.
! At this stage, model parameters are not known.
! Comments:
! 1. Common uses:
!   (a) Determine the actual dimensionality of a model, which,
!       in case of distributed models, may depend on topography data,
!       manual user-specifications, etc.
! 2. For simple, generic lumped models, this routine is seldom used.
implicit none
! dummies
integer(mik),intent(in)::modelID(:)
character(*),intent(in)::setupCmd
character(*),intent(in),optional::chvarLibDef(:,:)
integer(mik),intent(out)::err
character(*),intent(out)::message

! empty
err=666;message='empty subroutine'

endsubroutine DMDL_setModel
!----------------------------------------------------
subroutine DMDL_getModelInfo(modelID,infoCmd,&
  modelName,ninput,nstate,npar,indxName,inputName,stateName,parName,&
  stateLo,stateHi,parLo,parHi,inScal,stateScal,parScal,&
  stateDef,parDef,parSD,parTranDef,parFitDef,&
  err,message)
! Purpose: Returns basic properties of modelID.
! Usage: First call asks for dimensions, second call for supplying specific info.
implicit none
! dummies
integer(mik),intent(in)::modelID(:)
character(*),intent(in)::infoCmd
character(*),intent(out),optional::modelName
integer(mik),intent(out),optional::ninput,nstate,npar
character(*),intent(out),optional::indxName
character(*),intent(out),dimension(:),optional::inputName,stateName,parName
real(mrk),intent(out),dimension(:),optional::stateLo,stateHi,parLo,parHi,&
  inScal,stateScal,parScal,stateDef,parDef,parSD
integer(mik),intent(out),optional::parTranDef(:)
logical(mlk),intent(out),optional::parFitDef(:)
integer(mik),intent(out)::err
character(*),intent(out)::message

! empty
err=666;message='empty subroutine'

endsubroutine DMDL_getModelInfo
!----------------------------------------------------
subroutine DMDL_controlModel(modelID,inittCmd,dataXY,dataProps,parIn,dquanIn,&
  parOut,flexSin,setS0in,stateIn,stateOut,feas,err,message)
! Purpose: Sets/gets model states and parameters. Checks for feasibility of parameters and states.
! Usage:
!  - if(flexSin) then will adjust states to be compatible with parameter values,
!                     eg, if state S exceeds its maximum value Smax, will reset S to Smax.
!  - if(setS0in) then will set all states to initial conditions (either to default or to supplied vals)
implicit none
! dummies
integer(mik),intent(in)::modelID(:)
character(*),intent(in),optional::inittCmd
real(mrk),intent(in),optional::dataXY(:,:),dataProps(:)
real(mrk),intent(in),optional::parIn(:),dquanIn(:),stateIn(:)
logical(mlk),intent(in),optional::flexSin,setS0in
real(mrk),intent(out),optional::parOut(:),stateOut(:)
logical(mlk),intent(out),optional::feas
integer(mik),intent(out)::err
character(*),intent(out)::message

! empty
err=666;message='empty subroutine'

endsubroutine DMDL_controlModel
!----------------------------------------------------
subroutine DMDL_runModel(modelID,runitCmd,iT,dataProps,input,state,feas,err,message)
! Purpose: Performs single step of a system model.
implicit none
! dummies
integer(mik),intent(in)::modelID(:)
character(*),intent(in)::runitCmd
integer(mik),intent(in)::iT
real(mrk),intent(in)::dataProps(:)
real(mrk),intent(in)::input(:)
real(mrk),intent(out)::state(:)
logical(mlk),intent(out)::feas
integer(mik),intent(out)::err
character(*),intent(out)::message

! empty
err=666;message='empty subroutine'

endsubroutine DMDL_runModel
!----------------------------------------------------

end module dynamicModelLibrary
