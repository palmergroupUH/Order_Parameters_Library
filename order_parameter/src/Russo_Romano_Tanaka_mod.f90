module Russo_Romano_Tanaka
    use system, only: c_int, c_double, dp, sp
    use order_parameter, only: build_homo_neighbor_list,&
                             & check_nb_list, & 
                             & apply_nearest_neighbor_crit, &  
                             & initialize_Ylm, & 
                             & compute_optimized_q4, & 
                             & compute_optimized_q12, & 
                             & neighbor_averaged_qlm, & 
                             & nnb_dot_product_crystal, & 
                             & select_crystalline_particles  
 
    implicit none 
    ! Default private for all variables and routines in this module
    private 

    ! Choose following subroutines/functions accessible 
    ! by other Fortran modules/routines/functions/
    public :: call_RussoRomanoTanaka


contains

    subroutine call_RussoRomanoTanaka(total_atoms, &
                                    & maxnb,& 
                                    & nnb, &
                                    & cutoff_sqr,&
                                    & connect_cut, &
                                    & crys_cut, & 
                                    & box, &
                                    & xyz, &
                                    & cij) bind(c, name="call_RussoRomanoTanaka")
        implicit none 
    
        ! Passed:
        integer(c_int), intent(in) :: total_atoms
        integer(c_int), intent(in) :: maxnb
        integer(c_int), intent(in) :: nnb
        integer(c_int), intent(in) :: crys_cut
        real(c_double), intent(in) :: cutoff_sqr
        real(c_double), intent(in) :: connect_cut
        real(c_double), intent(in), dimension(1:3) :: box
        real(c_double), intent(in), dimension(1:3, total_atoms) :: xyz

        ! Local:
        integer, dimension(1:total_atoms) :: num_NB_list 
        integer, dimension(1:maxnb, 1:total_atoms) :: NB_list 
        real(dp), dimension(1:4, 1:maxnb, 1:total_atoms) :: Rij
        complex(dp), dimension(-12:12, total_atoms) :: q12 
        !complex(dp), dimension(-4:4, total_atoms) :: q4
        ! coarse-grained:
        complex(dp),dimension(-12:12, total_atoms) :: Q12_cg
        !complex(dp),dimension(-4:4, total_atoms) :: Q4
        logical, dimension(1:total_atoms) :: its_crystal 
        integer :: num_crystal
        integer, dimension(:), allocatable :: crystal_id 
        real(c_double), dimension(1:total_atoms) :: count_bonds 

        ! Return:
        real(c_double), intent(out), dimension(1:nnb, 1:total_atoms) :: cij

        ! get the neighbor list
        call build_homo_neighbor_list(total_atoms, maxnb, cutoff_sqr, xyz, box, num_NB_list, NB_list, Rij) 

        ! check the neighbor list: 
        call check_nb_list(num_NB_list,maxnb,nnb) 

        ! apply the nearest neighbor:
        call apply_nearest_neighbor_crit(nnb, total_atoms, num_NB_list, NB_list, Rij)
         
        ! compute bond orientational order q12
        call compute_optimized_q12(total_atoms, num_NB_list, Rij, q12)

        ! compute the coarse-grained Q12:
        call neighbor_averaged_qlm(total_atoms, 12, num_NB_list, NB_list, q12, Q12_cg) 

        ! compute the dot_product Q12(i) and Q12(j), number bonds, and crystal label
        call nnb_dot_product_crystal(total_atoms, nnb, 12, NB_list, connect_cut, crys_cut, Q12_cg, cij, count_bonds, its_crystal)

        ! assign and select the crystalline phase
        call select_crystalline_particles(total_atoms,its_crystal,num_crystal,crystal_id)

        ! compute q4 and Q4:
        ! call compute_optimized_q4(total_atoms, num_NB_list(crystal_id), Rij(:,:,crystal_id), q4)
        ! call neighbor_averaged_qlm(total_atoms, l, num_NB_list, NB_list, q4, Q4)
        ! call local_bond_order(total_atoms, l, Q4, local_ql) 
        ! compute W4:

        end subroutine


end module  
