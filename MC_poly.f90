program MC_chain
#include "prepro.h"
use com_vars
use ziggurat

implicit none

#ifdef SLIP_LINK
print*, "Slip Link definded"
#endif

! Read input parameters
call read_input()

! Initialize system 
inquire( file='init_positions.dat', exist=init_pos ) 
if ( init_pos )  then
    init_mode = 3
else
    init_mode = 1 !Random Walk
end if    
call init_system(init_mode) ! 1 is random walk;
                    ! 2 is uniformly random (no correlation between beads)
                    ! 3 is read old file with positions 

!DEBUG
print*, "r0"
print*, r0
print*,""
#ifdef SLIP_LINK
print*, "anch_r0"
print*, anch_r0
print*, ""
print*, "anchor"
print*, anchor
print*, ""
print*, "attach"
print*, attach
#endif

print*, "System initialized correctly"

! Meassure Energy
call get_energy()

#ifdef SLIP_LINK
call get_anch_ener()
#endif

print*,"Energy calculated correctly"
print*, ""

!begin time loop
do i_time = 1, n_time
    
!   Try displacement
    call move_attempt()
    
! Accept move with proba = min(1,exp( -delta_E/kT ) )
    call MC_acceptance()
    
#ifdef SLIP_LINK

    call hopp_attempt() ! Try hopping attempt
    
    if (.not. rej_mv ) call MC_hopp_acc() !Accept with probability = min(1,exp( -delta_energy/kT ) )
    
#endif


#if FREE_ENERGY == 1
    call rend_histo()
#endif

! Save data every n_mon trials
    if ( mod(i_time,n_save) .eq. 0 ) then
        call save_data()
    end if

end do ! End time loop

call save_positions()

#if FREE_ENERGY == 1

    call save_rend_histo()
#endif


! Meassure Energy
call get_energy()
print*, "Final bondend energy is: ", energy

#ifdef SLIP_LINK
call get_anch_ener()
print*, "Final Slip_Link energy is: ", sl_sp_ener

#endif

print*, "Finishing run ok!"

#ifdef g3
    call get_g3()
#endif

end program


subroutine rend_histo()
#if FREE_ENERGY == 1

use com_vars
implicit none
integer :: i_x, i_y, i_z

i_x = int( ( r0(1,n_mon) - r0(1,1) ) / x_max  * float(n_bins) )
i_y = int( ( r0(2,n_mon) - r0(2,1) ) / x_max  * float(n_bins) )
i_z = int( ( r0(3,n_mon) - r0(3,1) ) / x_max  * float(n_bins) )

!DEBUG
!print*, i_x, r0(1,n_mon) - r0(1,1),  
 
if( i_x .ge. 0 ) then
    !if( (i_y*i_y .eq. 0) .and. (i_z*i_z .eq.0) ) then !Rend_y = Rend_z = 0
        histo_count = histo_count + 1
        histo_rend(i_x) = histo_rend(i_x) + 1. 
    !end if
end if

#endif
end subroutine



subroutine read_init_pos()
use com_vars
implicit none

open (unit=16, file='init_positions.dat', status='old')

if(n_dim.ne.3) then 
    print*, "If n_dim is not 3, the program does not work"
    call exit()
end if

do i_mon = 1, n_mon
    read(16,"(3f15.8)") r0(1,i_mon), r0(2,i_mon), r0(3,i_mon)
end do

!Read anchor points positions also
#ifdef SLIP_LINK 

do i_anchor = 1, n_anchor
    read(16,"(3f15.8)" ) anch_r0(1,i_anchor), anch_r0(2,i_anchor), anch_r0(3,i_anchor) 
end do

!Read also link information attach, and anchor
do i_anchor = 1, n_anchor
    read(16,*) attach(i_anchor)
end do

do i_mon = 1, n_mon
    read(16,*) anchor(i_mon)
end do

#endif


close(16)

end subroutine

subroutine save_rend_histo()
#include "prepro.h"
#if FREE_ENERGY == 1

use com_vars
implicit none
integer :: i_bin
real(kind=8) :: x_coor

!Take out half the events for Rx=0, because that bin has double size
!This happens because of int(-a) = 0 for |a|<1
histo_rend(0) = histo_rend(0) / 2.

!normalize
!histo_rend(:) = histo_rend(:) / float(histo_count)

histo_rend(:) = histo_rend(:) / ( dx_bin * sum( histo_rend(:) ) )

open (unit=19, file='rend_histo.dat', status='unknown')
open (unit=31, file='Free_Energy_vs_Rx.dat', status='unknown')

do i_bin = 0, n_bins-1
    x_coor = i_bin*dx_bin
    write(19,*) x_coor, histo_rend(i_bin)
    if(histo_rend(i_bin).gt.0) then
        write(31,*) x_coor, -log( histo_rend(i_bin) / histo_rend(0) )
    end if
end do

close(31)
close(19)

print*,"sum(histo_rend(:))",sum(histo_rend(:))
#endif
end subroutine

subroutine save_positions()
#include "prepro.h"
use com_vars
implicit none

open (unit=14, file='last_config.dat', status='unknown')

do i_mon = 1, n_mon
    write(14,"(3f15.8)" ) r0(:,i_mon)
end do

!Save anchor points positions also
#ifdef SLIP_LINK 

do i_anchor = 1, n_anchor
    write(14,"(3f15.8)" ) anch_r0(:,i_anchor)
end do

!Save also link information attach, and anchor
do i_anchor = 1, n_anchor
    write(14,*) attach(i_anchor)
end do

do i_mon = 1, n_mon
    write(14,*) anchor(i_mon)
end do

#endif

close(14)

end subroutine

subroutine get_rcm()
use com_vars

implicit none

r_cm(:) = inv_nmon * sum( r0(:,:), 2)

end subroutine 

subroutine get_rend()
use com_vars

implicit none

r_end(:) = r0(:,n_mon) - r0(:,1) 

end subroutine 



subroutine save_data()
use com_vars
use ziggurat

implicit none
call get_rcm()
call get_rend()

!Save Energy
#ifdef SLIP_LINK
write(60,*) i_time, energy + sl_sp_ener !total energy
write(61,*) i_time, energy ! bonded energy 
write(62,*) i_time, sl_sp_ener !Slip-anchor energy

!DEBUGGING, sl_sp_ener variation  not calculated correctly 
!call get_anch_ener()
!write(62,*) i_time, sl_sp_ener

#else 
write(60,*) i_time, energy

#endif

!Save Chain Center of mass to calculate g3
write(73,"(I15.1,1f15.8)" ) i_time, sum ( ( r_cm(:) - r_cm_init(:) )**2 ) !r_cm(:) - r_cm_init(:) !

!Save mean bead displacement to calculate g1
write(71,"(I15.1,1f15.8)" ) i_time, sum ( ( r0(:,1) - r0_init(:,1) )**2 )

#ifdef g3 
    r_cm_t(:,j_save) = r_cm(:)
    j_save = j_save + 1
#endif


!Save Rend vector
write(81,"(I15.1,3f15.8)" ) i_time, r_end(:)

end subroutine

subroutine MC_acceptance()
#include "prepro.h"
use com_vars
use ziggurat

implicit none
real(kind=8) :: prob_rej

prob_rej = 1. - exp(-delta_energy / Temp)

if( .not. (uni() .le. prob_rej) ) then !If change is not rejected (so, accepted)
    energy = energy + delta_energy !change energy
    r0(:,mv_mon) = r0(:,mv_mon) + dr_mv(:) !perform movement

#ifdef SLIP_LINK
    if ( anchor(mv_mon) .ne. 0 ) then
        energy = energy - delta_sl_sp_ener
        sl_sp_ener = sl_sp_ener + delta_sl_sp_ener
        !print*,"move attempt"
        !call slip_energy_check()
    end if
#endif   
  end if    

end subroutine


subroutine move_attempt()
#include "prepro.h"
use com_vars
use ziggurat

implicit none

!Pick random monomer
mv_mon = int( uni() * float(n_mon) ) + 1 

!Set random move
do i_dim = 1, n_dim
    dr_mv(i_dim) = a_box * ( uni() - .5 )
end do

!get delta_energy
call get_delta_energy()

end subroutine


subroutine hopp_attempt()
#include "prepro.h"
#ifdef SLIP_LINK
use com_vars
use ziggurat

implicit none

!Pick random anchor
mv_anchor = int( uni() * float(n_anchor) ) + 1 

!DEBUG
!print*,"mv_anchor",mv_anchor

!Set random hop
if( uni() .le. 0.5) then
    hop_mv = -1
else
    hop_mv = 1
end if

j_mon = attach(mv_anchor) + hop_mv ! New Monomer to anchor 

!DEBUG
!print*,"j_mon",j_mon

!CHECK IF j_mon = attach(mv_anchor) + hop_mv is 0 or n_mon + 1
if ( (j_mon.eq.0) .or. (j_mon.eq.n_mon+1)) then
    !DEBUG
    !print*, "calling tube_renewal and constraint_release"

    call tube_renewal() !New anchor monomer is out of the chain, perform tube_renewal
    call constraint_release() ! and constraint_release with probability 1

    rej_mv = .True. !Reject Hopping

!CHECK IF j_mon = attach(mv_anchor) + hop_mv is occupied.
else if( anchor(j_mon) .ne. 0 ) then
    
    !DEBUG
    !print*, "anchor(j_mon) .ne. 0", anchor(j_mon)
    
    rej_mv = .True. 
else
   
   !DEBUG
    !print*,"performing hopp"
    
    rej_mv = .False.
end if

!get delta_energy
if (.not. rej_mv) call get_anch_delta_energy()

#endif

end subroutine


subroutine tube_renewal()
#include "prepro.h"
#ifdef SLIP_LINK
use com_vars
use ziggurat
implicit none
integer :: end_mon

!Reclaculate slink spring energy: Substract spring energy
sl_sp_ener = sl_sp_ener - .5 * k_sl_sp * sum( ( anch_r0(:,mv_anchor) - r0(:,attach(mv_anchor)) ) ** 2 )

!Release end bead attachment
anchor( attach(mv_anchor) ) = 0 

!Choose a free chain end monomer (end_mon) 
if( uni() .le. 0.5) then
    end_mon = 1
    if ( anchor(end_mon) .ne. 0) end_mon = n_mon    
else
    end_mon = n_mon
    if ( anchor(end_mon) .ne. 0) end_mon = 1 
end if

!Redo lists
attach(mv_anchor) = end_mon
anchor( end_mon ) = mv_anchor

! create a new anchor point according to
! anch_r0(i_dim,mv_anch) = std_dev * rnor() + r0(i_dim,end_mon)

do i_dim = 1, n_dim
    anch_r0(i_dim,mv_anchor) = std_dev * rnor() + r0(i_dim,end_mon)
end do

!Add new spring energy
sl_sp_ener = sl_sp_ener + .5 * k_sl_sp * sum( ( anch_r0(:,mv_anchor) - r0(:,attach(mv_anchor)) ) ** 2 )

!print*,"tube_renewal"
!call slip_energy_check()
#endif /*anchor*/
end subroutine tube_renewal

subroutine constraint_release()
#include "prepro.h"
#ifdef SLIP_LINK
use com_vars
use ziggurat
implicit none
integer :: new_mon, nei_mon, nei_anch
 
!Release neighbour monomer attachment
nei_anch = anch_neigh( mv_anchor )
nei_mon = attach( nei_anch )

!First take out attach energy
sl_sp_ener = sl_sp_ener - .5 * k_sl_sp * sum( ( anch_r0(:,nei_anch) - r0(:,nei_mon) ) ** 2 )

anchor( nei_mon ) = 0

!Choose a new monomer bead
new_mon = int( float(n_mon) * uni() ) + 1

do while( anchor(new_mon) .ne. 0 ) 
    new_mon = int( float(n_mon) * uni() ) + 1
end do

!Attach to new monomer bead
attach( nei_anch ) = new_mon
anchor( new_mon ) = nei_anch !new_mon is anchored to neighbour of mv_anchor

! Set ach_r0
do i_dim = 1, n_dim
    anch_r0(i_dim,nei_anch) = std_dev * rnor() + r0(i_dim,new_mon)
end do

!Add new Slip_Link spring Energy
sl_sp_ener = sl_sp_ener + .5 * k_sl_sp * sum( ( anch_r0(:,nei_anch) - r0(:,new_mon) ) ** 2 )

!print*,"constraint_release"
!call slip_energy_check()
#endif /*anchor*/
end subroutine constraint_release


subroutine get_anch_delta_energy()
#include "prepro.h"
#ifdef SLIP_LINK
use com_vars
implicit none

delta_energy = 0.

delta_energy = delta_energy - .5 * k_sl_sp * sum( ( anch_r0(:,mv_anchor) - r0(:,attach(mv_anchor)) ) ** 2 )

delta_energy = delta_energy + .5 * k_sl_sp * sum( ( anch_r0(:,mv_anchor) - r0(:,attach(mv_anchor) + hop_mv) ) ** 2 )

#endif
end subroutine

subroutine MC_hopp_acc()
#include "prepro.h"
#ifdef SLIP_LINK
use com_vars
use ziggurat

implicit none
real(kind=8) :: prob_rej

prob_rej = 1. - exp(-delta_energy / Temp)

!DEBUG
!print*,"delta_energy",delta_energy
!print*,"prob_rej",prob_rej

if( .not. (uni() .le. prob_rej) ) then !If change is not rejected (so, accepted)

    !DEBUG
    !print*,"Hopp accepted"
    !print*,"mv_anchor",mv_anchor
    !print*, "release monomer attach(mv_anchor)", attach(mv_anchor)
    !print*, "Attaching mv_anchor to monomer attach(mv_anchor) + hop_mv", attach(mv_anchor) + hop_mv
    
    sl_sp_ener = sl_sp_ener + delta_energy !change energy
    anchor(attach(mv_anchor)) = 0 ! free bead number = attach(mv_anchor)     
    attach(mv_anchor) = attach(mv_anchor) + hop_mv !perform movement
    anchor(attach(mv_anchor)) = mv_anchor
    !print*,"hopp"
    !call slip_energy_check()
end if    

#endif
end subroutine




subroutine get_delta_energy()
#include "prepro.h"
use com_vars
implicit none

delta_energy = 0.
if (mv_mon.ne.1)     delta_energy = delta_energy + .5 * k_spr * sum( ( r0(:,mv_mon) + dr_mv(:) - r0(:,mv_mon-1) ) ** 2 )
if (mv_mon.ne.n_mon) delta_energy = delta_energy + .5 * k_spr * sum( ( r0(:,mv_mon) + dr_mv(:) - r0(:,mv_mon+1) ) ** 2 )  
if (mv_mon.ne.1)     delta_energy = delta_energy - .5 * k_spr * sum( ( r0(:,mv_mon)            - r0(:,mv_mon-1) ) ** 2 )
if (mv_mon.ne.n_mon) delta_energy = delta_energy - .5 * k_spr * sum( ( r0(:,mv_mon)            - r0(:,mv_mon+1) ) ** 2 ) 

#ifdef SLIP_LINK
!Get delta energy from anchor point to monomer mn_mon
if ( anchor(mv_mon) .ne. 0 ) then ! if monomer is attached to an anchor point
    delta_sl_sp_ener = 0.
    delta_sl_sp_ener = delta_sl_sp_ener + .5 * k_sl_sp * sum( ( r0(:,mv_mon) + dr_mv(:) - anch_r0(:,anchor(mv_mon)) ) **2 )
    delta_sl_sp_ener = delta_sl_sp_ener - .5 * k_sl_sp * sum( ( r0(:,mv_mon)            - anch_r0(:,anchor(mv_mon)) ) **2 )   
    delta_energy = delta_energy + delta_sl_sp_ener 
end if

#endif

#if FREE_ENERGY == 2
if (mv_mon .eq. 1) then
    delta_energy = delta_energy + sum( f_ext(:) * ( r0(:,n_mon) - r0(:,1) -  dr_mv(:))) 
    delta_energy = delta_energy - sum( f_ext(:) * ( r0(:,n_mon) - r0(:,1) ))
else if (mv_mon .eq. n_mon) then
    delta_energy = delta_energy + sum( f_ext(:) * ( r0(:,n_mon) +  dr_mv(:) - r0(:,1) )) 
    delta_energy = delta_energy - sum( f_ext(:) * ( r0(:,n_mon) - r0(:,1) )) 
end if
#endif
!print*,delta_energy
end subroutine

subroutine get_energy()
use com_vars
implicit none

energy = 0

!loop over monomers
do i_mon = 2, n_mon 
    energy = energy + k_spr * sum( ( r0(:,i_mon) - r0(:,i_mon-1) ) ** 2 )     
end do

energy = energy * 0.5

#if FREE_ENERGY == 2

energy = energy + sum( f_ext(:) * ( r0(:,n_mon) - r0(:,1)))

#endif
!write(61,*) energy
end subroutine

subroutine get_anch_ener()
#include "prepro.h"
#ifdef SLIP_LINK
use com_vars
implicit none

sl_sp_ener = 0.

!loop over anchor points
do i_anchor = 1, n_anchor
    sl_sp_ener = sl_sp_ener + k_sl_sp * sum( ( anch_r0(:,i_anchor) - r0(:,attach(i_anchor)) ) ** 2 )     
end do

sl_sp_ener = sl_sp_ener * 0.5
#endif
end subroutine


subroutine init_system()
#include "prepro.h"
use com_vars
use ziggurat
implicit none
real(kind=8) signo, dr
real(kind=8),dimension(n_dim) :: r_center

allocate( r0(n_dim,n_mon),r0_init(n_dim,n_mon), boundary(n_dim), dr_mv(n_dim) ) 
allocate( r_end(n_dim), r_cm(n_dim), r_cm_init(n_dim) )
#ifdef g3
    allocate (r_cm_t(n_dim,int(n_time/n_save)+1), g3(int(n_time/n_save))
#endif

k_spr = 3. * ( float( n_mon ) - 1. ) / Rend2

#ifdef SLIP_LINK
    k_sl_sp = 0.5 * k_spr
    std_dev = sqrt( 1. / k_sl_sp )
    n_anchor = int( float (n_mon) / 4. )
    if(mod(n_anchor,2).ne.0) n_anchor = n_anchor + 1 ! Make shure n_anchor is
                                                     !even
    allocate( attach(n_anchor), anchor(n_mon), anch_neigh(n_anchor) )
    allocate( anch_r0(n_dim,n_anchor) )  
#endif anchor

#if FREE_ENERGY == 1
x_max = 4.
n_bins = 60
histo_count = 0
dx_bin = x_max / float(n_bins)
allocate( histo_rend(0:n_bins) )
histo_rend(0:n_bins) = 0.
#endif

!Initiate random seed
call init_rand_seed()

!Set simulation parameters 
boundary(:) = Lx
inv_nmon = 1 / float( n_mon )

select case(init_mode)

case(1) ! Random Walk
    !Set position of initial monomer
    r_center(:) =  boundary(:) / 2.
    r0(:,1) = r_center(:)
    
    !Set position of remaining monomers
    do i_mon = 2, n_mon
        do i_dim = 1, n_dim
            !Set distance to next monomer
            dr = rnor() + a_box
            !Set sign to move
            signo = uni()
            if (signo .lt. 0.5 ) then
                signo = -1.
            else
                signo = 1.
            end if
            r0( i_dim, i_mon ) = r0( i_dim, i_mon - 1 ) + dr * signo  
        end do
    end do

case(2) 
    do i_mon = 1, n_mon
        do i_dim = 1, n_dim
            r0(i_dim, i_mon) = boundary(i_dim) * uni()
        end do
    end do

case(3)
call read_init_pos()

end select

#ifdef SLIP_LINK
!Initiate anchor points

!!Randomly in simulation box:
!do i_anchor = 1, n_anchor
!    do i_dim = 1, n_dim
!        anch_r0(i_dim, i_anchor) = boundary(i_dim) * uni()
!    end do
!end do

if (init_mode .ne. 3) then ! If config is not read from init_positions.dat
    anchor(:) = 0
    !Set which monomers are anchored
    do i_anchor = 1, n_anchor !loop anchor points
    
        j_mon = int( uni() * float(n_mon) ) + 1 !select monomer to be attached to anchor point i_anchor
        
        do while( anchor( j_mon ) .ne. 0 )  !If monomer is not free, try again
            j_mon = int( uni() * float(n_mon) ) + 1
        end do
    
        anchor( j_mon ) = i_anchor ! j_mon is free, anchor j_mon to i_anchor
        attach(i_anchor) = j_mon
       
    end do
    
    !Initiate anchor positions anch_r0, near attached monomers, with gaussian
    !probability distribution
    
    do i_anchor = 1, n_anchor
        do i_dim = 1, n_dim
            anch_r0(i_dim,i_anchor) = std_dev * rnor() + r0(i_dim,attach(i_anchor))
        end do
    end do

end if    !If anchor positions are not read from file init_positions.dat

    !Set neighbors of anchors to perform tube_renewal and coinstraint_release    
    do i_anchor = 1, n_anchor !loop anchor points
        if ( mod(i_anchor,2) .eq. 0) then !if anchor point  is pair
            anch_neigh( i_anchor ) = i_anchor - 1 !neighbour is previous anch
        else ! anchor point is odd
            anch_neigh( i_anchor ) = i_anchor + 1 ! neighbour is next anch
        end if
    end do

#endif


!Center initial position
r_cm_init(:) =  inv_nmon * sum( r0(:,:), 2)

do i_dim = 1, n_dim
    r0(i_dim,:) = r0(i_dim,:) - r_cm_init(i_dim) + boundary(i_dim) / 2.

#ifdef SLIP_LINK
     anch_r0(i_dim,:) = anch_r0(i_dim,:) - r_cm_init(i_dim) + boundary(i_dim) / 2.
#endif

end do

!Save initial position

r_cm_init(:) =  inv_nmon * sum( r0(:,:), 2)
r0_init = r0
end subroutine


subroutine read_input()
use com_vars
implicit none

open (unit = 53, file = "input.dat", status= "old")

read(53,*) n_dim

#if FREE_ENERGY == 2
allocate(f_ext(n_dim))
#endif

read(53,*) Lx !, Ly, Lz Square box
read(53,*) n_time
read(53,*) n_mon
read(53,*) a_box
read(53,*) Rend2
read(53,*) Temp
read(53,*) n_save

#if FREE_ENERGY == 2
read(53,*) f_ext(:)
print*, "Simulating with external force f_ext =", f_ext
#endif

close(53)

end subroutine

      

subroutine init_rand_seed()
use ziggurat
logical :: file_ex
integer :: rand_seed
inquire(file='random_seed.dat',exist=file_ex)
if (file_ex) then
    open(unit = 68, file = "random_seed.dat", status= "old")
    read(68,*) rand_seed
    close(68)
else
    rand_seed = 465178
end if

call zigset( rand_seed )

open(unit = 68, file = "random_seed.dat", status="unknown")
write(68,*) shr3()
close(68)

end subroutine


subroutine get_g3()
use com_vars
#include "prepro.h"
#ifdef g3

implicit none

integer :: n
real (kind=8), dimension(n), intent(out)  :: c
integer ,dimension(:), allocatable :: n_corr
integer :: i,j
real (kind=8) :: x_mean = 0, x_var = 0
n = int(n_time/n_save)
allocate(n_corr(n))
!Initialize counters to  0
g3(:) = 0. 
n_corr = 0

!Now calculate the variance and autocorrelation
do i=1,n-1 
    do j=i+1,n
        g3(j-i) =  sum( (r_cm_t(:,i) - r_cm_t(:,j) )**2 ) !c(j-i+1) + y(i)*y(j) ! Accumulate covariance for each lag
        n_corr(j-i) = n_corr(j-i) + 1 !Count cases for each lag
    end do
end do

do i=1,n-1
    g3(i) = g3(i) / float(n_corr(i)) ! Get mean Covariance
end do
#endif
end subroutine 

subroutine slip_energy_check()
use com_vars
#include "prepro.h"
#ifdef SLIP_LINK
implicit none
real(kind=8) :: old_sl_sp_en

old_sl_sp_en = sl_sp_ener

call get_anch_ener()
if(old_sl_sp_en .ne. sl_sp_ener ) print*, "ERROR: sl_sp_ener badly calculated", sl_sp_ener - old_sl_sp_en
#endif
end subroutine
