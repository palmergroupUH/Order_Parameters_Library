module Russo_Romano_Tanaka
    use system, only: c_int, c_double, dp, sp 
    use order_parameter, only: build_homo_neighbor_list,&
                             & compute_sum_Ylm, & 
                             & neighbor_averaged_qlm, & 
                             & compute_dot_product 
    implicit none 
    ! Default private for all variables and routines in this module
    private 

    ! Choose following subroutines/functions accessible 
    ! by other Fortran modules/routines/functions/
    public :: call_RussoRomanoTanaka


contains

    subroutine call_RussoRomanoTanaka(sph_const,& 
                                    & Plm_const, &
                                    & total_atoms, &
                                    & maxnb,& 
                                    & l, &
                                    & cutoff_sqr,&
                                    & box, &
                                    & xyz) bind(c, name="call_RussoRomanoTanaka")
        implicit none 
    
        ! Passed:
        real(c_double), intent(in), dimension(-l:l) :: sph_const 
        real(c_double), intent(in), dimension(0:l,4) :: Plm_const 
        integer(c_int), intent(in) :: total_atoms
        integer(c_int), intent(in) :: maxnb
        integer(c_int), intent(in) :: l
        real(c_double), intent(in) :: cutoff_sqr
        real(c_double), intent(in), dimension(1:3) :: box
        real(c_double), intent(in), dimension(1:3, total_atoms) :: xyz

        ! Local:
        integer, dimension(1:total_atoms) :: num_NB_list 
        integer, dimension(1:maxnb, 1:total_atoms) :: NB_list 
        real(dp), dimension(1:4, 1:maxnb, 1:total_atoms) :: Rij
        complex(dp),dimension(-l:l, total_atoms) :: Ylm, qlm_nb_ave
        real(dp), dimension(1:maxnb, 1:total_atoms) :: cij 

        ! Return:
         

        ! get the neighbor list
        call build_homo_neighbor_list(total_atoms, maxnb, cutoff_sqr, xyz, box, num_NB_list, NB_list, Rij) 

        ! compute_sum_Ylm  
        call compute_sum_Ylm(sph_const, Plm_const, total_atoms, l, num_NB_list, Rij, Ylm)

        ! compute the neighbor_averaged
        call neighbor_averaged_qlm(total_atoms, l, num_NB_list, NB_list, Ylm, qlm_nb_ave) 
    
        ! compute the dot_product:
        call compute_dot_product(total_atoms, maxnb, l, num_NB_list, NB_list, qlm_nb_ave, cij) 

       
        end subroutine
    

end module  
