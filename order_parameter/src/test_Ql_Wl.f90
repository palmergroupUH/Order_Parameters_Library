module test_Wl_Ql
    use system, only: c_int, c_double, dp, sp
    use order_parameter, only: build_homo_neighbor_list,&
                             & check_nb_list, & 
                             & apply_nearest_neighbor_crit, &  
                             & initialize_Ylm, & 
                             & compute_ave_Ylm, & 
                             & global_bond_order, & 
                             & local_bond_order, & 
                             & calc_Wl  
 
    implicit none 
    ! Default private for all variables and routines in this module
    private 

    ! Choose following subroutines/functions accessible 
    ! by other Fortran modules/routines/functions/
    public :: call_test_Wl_Ql


contains

    subroutine initialize_test_Wl_Ql(l, sph_const, Plm_const) bind(c,name="initialize_test_Wl_Ql")
        implicit none 
        integer(c_int), intent(in) :: l
        real(c_double), intent(out), dimension(-l:l) :: sph_const
        real(c_double), intent(out), dimension(0:l,4) :: Plm_const

        call initialize_Ylm(l, sph_const, Plm_const)

        end subroutine

    subroutine call_test_Wl_Ql(sph_const,& 
                            & Plm_const, &
                            & m1_m2_m3_sum0_mat, & 
                            & wigner_3j_symobol_mat, & 
                            & num_pairs, & 
                            & total_atoms, &
                            & maxnb,& 
                            & nnb, &
                            & l, &
                            & cutoff_sqr,&
                            & box, &
                            & xyz) bind(c, name="call_test_Wl_Ql")
        implicit none 
    
        ! Passed:
        real(c_double), intent(in), dimension(-l:l) :: sph_const 
        real(c_double), intent(in), dimension(0:l,4) :: Plm_const 
        integer(c_int), intent(in) :: num_pairs 
        real(c_double), intent(in), dimension(1:num_pairs) :: wigner_3j_symobol_mat
        integer(c_int), intent(in), dimension(1:3, 1:num_pairs) :: m1_m2_m3_sum0_mat  
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
        complex(dp),dimension(-l:l, total_atoms) :: Ylm 
        real(dp) :: global_op  
        real(dp), dimension(1:total_atoms) :: local_ql, Wl 
        ! Return:
        
        ! get the neighbor list
        
        call build_homo_neighbor_list(total_atoms, maxnb, cutoff_sqr, xyz, box, num_NB_list, NB_list, Rij) 

        ! check the neighbor list: 
        call check_nb_list(num_NB_list,maxnb,nnb) 
        
        ! apply the nearest neighbor:
        call apply_nearest_neighbor_crit(nnb, total_atoms, num_NB_list, NB_list, Rij)
         
        ! compute_sum_Ylm
        call compute_ave_Ylm(sph_const, Plm_const, total_atoms, l, num_NB_list, Rij, Ylm)
        
        ! compute global order parameter 
        call global_bond_order(total_atoms, num_NB_list, l, Ylm, global_op) 

        ! compute local order parameter 
        call local_bond_order(total_atoms, l, Ylm, local_ql)  
        
        ! compute W4, W6, ... 
        call calc_Wl(num_pairs, total_atoms, l, m1_m2_m3_sum0_mat, wigner_3j_symobol_mat, Ylm, Wl)
        print*, global_op        
        print*, Wl(1:5) 
        end subroutine

end module  
