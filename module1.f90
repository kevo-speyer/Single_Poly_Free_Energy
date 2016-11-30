module com_vars
!Module with all global
#include "prepro.h"
implicit none
integer :: n_save, i_mon, init_mode, j_mon, n_mon, i_time, n_time, i_dim, n_dim, mv_mon
real(kind=8) , dimension(:,:), allocatable :: r0, r0_init
real(kind=8) , dimension(:), allocatable :: boundary, dr_mv, r_cm,r_cm_init, r_end
real(kind=8) :: Rend2, Temp, a_box, Lx, Ly, Lz, energy, delta_energy, old_energy, k_spr, inv_nmon
logical :: init_pos

#ifdef g3
real(kind=8) , dimension(:,:), allocatable :: r_cm_t
real(kind=8) , dimension(:), allocatable :: g3
integer :: j_save=1
#endif

#ifdef SLIP_LINK
integer , dimension(:), allocatable :: attach, anchor, anch_neigh
!attach(i_anchor) = return the number of the bead attached to anchor i_anchor
!anchor(i_mon) = returns the number of anchor attached to bead i_mon. Returns 0
!if i_mon is not attached to any anchor
real(kind=8) , dimension(:,:), allocatable ::  anch_r0 !anchoring points
!positions

integer :: hop_mv, i_anchor, n_anchor, mv_anchor ! Number of anchor points
real(kind=8) :: sl_sp_ener, delta_sl_sp_ener,  k_sl_sp, std_dev  !Total slip-link energy, slip_spring constant
logical :: rej_mv
#endif

#ifdef  FREE_ENERGY_A
integer :: n_bins, histo_count
real(kind=8) :: x_max, dx_bin
real(kind=8), dimension(:), allocatable :: histo_rend
#endif
end module com_vars
 
