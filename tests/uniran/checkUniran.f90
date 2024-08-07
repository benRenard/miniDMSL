program main_tests

use kinds_dmsl_kit    ! numeric kind definitions
use uniran1_dmsl_mod
use numerix_dmsl_kit,only:uniran
implicit none
integer(mik),parameter::nsim=10,nsim2= 1048576 ! 1024*1024
integer(mik)::i
real(mrk)::u

write (*,*) '------------'
write (*,*) 'Default seed'
do i=1,nsim
    call uniran(u)
    write(*,*) u
enddo

write (*,*) '------------'
write (*,*) 'seed reseted every iteration'
do i=1,nsim
    call seed_uniran(1)
    call uniran(u)
    write(*,*) u
enddo

write (*,*) '------------'
write (*,*) 'seed reseted every 2 iterations'
do i=1,nsim
    call seed_uniran(mod(i,2))
    call uniran(u)
    write(*,*) u
enddo

open(unit=1,file='uniranDeviates.txt')
do i=1,nsim2
    call uniran(u)
    write(1,'(e14.8)') u
enddo
close(1)
end program main_tests
