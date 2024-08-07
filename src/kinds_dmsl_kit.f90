module kinds_dmsl_kit

! global precision
integer,parameter::mrk=selected_real_kind(p=8)   ! global real kind
integer,parameter::mik=selected_int_kind(r=8)    ! global integer kind
integer,     parameter::mlk=kind(.true.)         ! global logical kind
integer,     parameter::mck=kind((1._mrk,1._mrk))     ! global complex kind
real(mrk),   parameter::protoRe=1._mrk                ! prototype of real(mrk) number
integer(mik),parameter::protoInt=1_mik                ! prototype of integer(mik) number
! Large and small numbers
real(mrk),   parameter::tinyRe=tiny(protoRe)          ! smallest real on machine
real(mrk),   parameter::epsRe= epsilon(protoRe)       ! normalised machine accuracy
real(mrk),   parameter::hugeRe=huge(protoRe)          ! largest real on machine
integer(mik),parameter::hugeInt= huge(protoInt)       ! largest integer on machine
! special flag values
real(mrk),   parameter::undefRN=-999999999._mrk       ! flag for undefined real numbers
integer(mik),parameter::undefIN=-999999999            ! flag for undefined integer numbers
character(*),parameter::undefCH="undefined"           ! flag for undefined character strings
logical(mlk),parameter::undefLG=.false.               ! flag for undefined logicals
complex(mck),parameter::undefCZ=(-999999999._mrk,-999999999._mrk) ! flag for undefined complex numbers
integer(mik),parameter::QC_PRESENT=0        ! quality code for present data
integer(mik),parameter::QC_MISSVAL=-6666    ! quality code for missing data
! software aspects
integer(mik),parameter::EXIT_SUCCESS=0,EXIT_FAILURE=+1,EXIT_WARNING=-1
! Platform-specific aspects
character(*),parameter::dirSepCH_win="\"        ! directory name separator
character(*),parameter::dirSepCH_lix="/"
character(*),parameter::dirSepCH=dirSepCH_lix   ! default
! strings
integer(mik),parameter::len_iShortStr=1     ! infrashort string
integer(mik),parameter::len_vShortStr=4     ! 2^2 : usually for quick settings via single args
integer(mik),parameter::len_shortStr=8      ! 2^3 : shortish names
integer(mik),parameter::len_stdStrB=16      ! 2^4 : brief descriptive names
integer(mik),parameter::len_stdStrD=32      ! 2^5 : detailed descriptive name
integer(mik),parameter::len_stdStrL=128     ! 2^7 : typical line-long string
integer(mik),parameter::len_longStr=256     ! 2^8 : long string
integer(mik),parameter::len_vLongStr=1024   ! 2^10: verylong string, usually for file paths or very long messages
integer(mik),parameter::len_uLongStr=8192   ! 2^13: ultralong string
integer(mik),parameter::len_hLongStr=65536  ! 2^15: hyperlong string (use with care)

end module kinds_dmsl_kit
