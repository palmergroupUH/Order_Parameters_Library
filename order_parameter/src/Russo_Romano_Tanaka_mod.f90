module Russo_Romano_Tanaka
    use system, only: c_int, c_double, dp, sp
    use order_parameter, only: build_homo_neighbor_list,&
                             & check_nb_list, & 
                             & apply_nearest_neighbor_crit, &  
                             & initialize_Ylm, & 
                             & compute_sum_Ylm, & 
                             & neighbor_averaged_qlm, & 
                             & compute_dot_product_nnb 
 
    implicit none 
    ! Default private for all variables and routines in this module
    private 

    ! Choose following subroutines/functions accessible 
    ! by other Fortran modules/routines/functions/
    public :: call_RussoRomanoTanaka


contains

    subroutine initialize_RussoRomanoTanaka(l, sph_const, Plm_const) bind(c,name="initialize_RussoRomanoTanaka")
        implicit none 
        integer(c_int), intent(in) :: l
        real(c_double), intent(out), dimension(-l:l) :: sph_const
        real(c_double), intent(out), dimension(0:l,4) :: Plm_const

        call initialize_Ylm(l, sph_const, Plm_const)

        end subroutine

    subroutine call_RussoRomanoTanaka(sph_const,& 
                                    & Plm_const, &
                                    & total_atoms, &
                                    & maxnb,& 
                                    & nnb, &
                                    & l, &
                                    & cutoff_sqr,&
                                    & box, &
                                    & xyz, &
                                    & cij) bind(c, name="call_RussoRomanoTanaka")
        implicit none 
    
        ! Passed:
        real(c_double), intent(in), dimension(-l:l) :: sph_const 
        real(c_double), intent(in), dimension(0:l,4) :: Plm_const 
        integer(c_int), intent(in) :: total_atoms
        integer(c_int), intent(in) :: maxnb
        integer(c_int), intent(in) :: nnb
        integer(c_int), intent(in) :: l
        real(c_double), intent(in) :: cutoff_sqr
        real(c_double), intent(in), dimension(1:3) :: box
        real(c_double), intent(in), dimension(1:3, total_atoms) :: xyz

        ! Local:
        integer, dimension(1:total_atoms) :: num_NB_list 
        integer, dimension(1:maxnb, 1:total_atoms) :: NB_list 
        real(dp), dimension(1:4, 1:maxnb, 1:total_atoms) :: Rij
        complex(dp),dimension(-l:l, total_atoms) :: Ylm, qlm_nb_ave

        ! Return:
        real(c_double), intent(out), dimension(1:nnb, 1:total_atoms) :: cij
         
        ! get the neighbor list
        call build_homo_neighbor_list(total_atoms, maxnb, cutoff_sqr, xyz, box, num_NB_list, NB_list, Rij) 

        ! check the neighbor list: 
        call check_nb_list(num_NB_list,maxnb,nnb) 

        ! apply the nearest neighbor:
        call apply_nearest_neighbor_crit(nnb, total_atoms, num_NB_list, NB_list, Rij)
         
        ! compute_sum_Ylm  
        call compute_sum_Ylm(sph_const, Plm_const, total_atoms, l, num_NB_list, Rij, Ylm)
         
        ! compute the neighbor_averaged
        call neighbor_averaged_qlm(total_atoms, l, num_NB_list, NB_list, Ylm, qlm_nb_ave) 
    
        ! compute the dot_product:
        call compute_dot_product_nnb(total_atoms, nnb, l, NB_list, qlm_nb_ave, cij) 
        
        end subroutine

    subroutine scalar_product()
        implicit none 
    
        
        end subroutine 
end module  
