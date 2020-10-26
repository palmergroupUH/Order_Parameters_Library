module bhatia_analysis
    use system, only: c_int, c_double, dp, sp 
    use order_parameter, only: build_homo_neighbor_list,&
                             & check_nb_list, &
                             & apply_nearest_neighbor_crit

    use sorting, only: quick_sort  

    use tetrahedral_order_parameter

    use RDF, only: build_homo_pair_dist_hist, & 
                 & build_hetero_pair_dist_hist, & 
                 & normalize_hist 
       
    use statistics, only: linear_reg 
 
    implicit none 
    private 

    public :: label_high_q_low_q, & 
            & compute_homo_q_pair_corrl, & 
            & compute_hetero_q_pair_corrl  


contains 

    subroutine label_high_q_low_q(total_atoms, & 
                                & maxnb, & 
                                & nnb, & 
                                & xyz, & 
                                & box, & 
                                & cutoff_sqr, & 
                                & sorted_q_indx, & 
                                & q_tetra) bind(c, name="call_label_high_q_low_q")
        implicit none 

        ! Passed 
        integer(c_int), intent(in) :: total_atoms 
        integer(c_int), intent(in) :: maxnb 
        integer(c_int), intent(in) :: nnb 
        real(c_double), intent(in) :: cutoff_sqr 
        real(c_double), intent(in), dimension(1:3) :: box 
        real(c_double), intent(in), dimension(1:3, total_atoms) :: xyz 
    
        ! Local
        integer, dimension(1:total_atoms) :: num_NB_list 
        integer, dimension(1:maxnb, 1:total_atoms) :: NB_list 
        real(dp), dimension(1:4, 1:maxnb, 1:total_atoms) :: Rij 
        integer, dimension(:), allocatable :: tracking_index 

        ! Return
        real(c_double), intent(out), dimension(1:total_atoms) :: q_tetra
        integer, intent(out), dimension(1:total_atoms) :: sorted_q_indx 
        
        call build_homo_neighbor_list(total_atoms, maxnb, cutoff_sqr, xyz, box, num_NB_list, NB_list, Rij) 

        call check_nb_list(num_NB_list,maxnb,nnb) 

        call apply_nearest_neighbor_crit(nnb, total_atoms, num_NB_list, NB_list, Rij) 
 
        call calc_tetrahedral_op(maxnb, nnb, total_atoms, Rij, q_tetra)

        call quick_sort(total_atoms, q_tetra, 1, total_atoms, tracking_index)
   
        sorted_q_indx = tracking_index 

        end subroutine 

    subroutine compute_homo_q_pair_corrl(total_atoms, &
                                       & select_indx, & 
                                       & num_bins, & 
                                       & xyz, & 
                                       & box, & 
                                       & rdf_histogram) bind(c, name="call_homo_q_pair_corrl")
        implicit none 
        
        ! Passed
        integer(c_int), intent(in) :: total_atoms 
        integer(c_int), intent(in) :: num_bins 
        integer(c_int), intent(in), dimension(1:total_atoms) :: select_indx  
        real(c_double), intent(in), dimension(1:3) :: box 
        real(c_double), intent(in), dimension(1:3, total_atoms) :: xyz 
  
        ! Local
        real(dp) :: cutoff
        
        ! Return 
        real(c_double),intent(out),dimension(1:num_bins) :: rdf_histogram  

        ! cutoff is half of the box size and convert the Angstrom to nm 

        cutoff = minval(box)/20.0d0  

        call build_homo_pair_dist_hist(total_atoms,&
                                    & cutoff,&
                                    & num_bins,&
                                    & xyz(:, select_indx)/10.0d0 ,&
                                    & box/10.0d0,&
                                    & rdf_histogram) 

        end subroutine 

    subroutine compute_hetero_q_pair_corrl(total_atoms, & 
                                      & total_atoms_A, & 
                                      & total_atoms_B, &
                                      & select_indx_A, & 
                                      & select_indx_B, &   
                                      & num_bins, &
                                      & xyz, &
                                      & box, &
                                      & rdf_histogram) bind(c, name="call_hetero_q_pair_corrl")
        implicit none 

        ! Passed
        integer(c_int), intent(in) :: total_atoms, total_atoms_A, total_atoms_B 
        integer(c_int), intent(in) :: num_bins 
        integer(c_int), intent(in), dimension(1:total_atoms_A) :: select_indx_A 
        integer(c_int), intent(in), dimension(1:total_atoms_B) :: select_indx_B 
        real(c_double), intent(in), dimension(1:3) :: box 
        real(c_double), intent(in), dimension(1:3, total_atoms) ::  xyz
  
        ! Local
        real(dp) :: cutoff
        
        ! Return 
        real(c_double),intent(out),dimension(1:num_bins) :: rdf_histogram  

        ! cutoff is half of the box size and convert the Angstrom to nm
        cutoff = minval(box)/20.d0

        call build_hetero_pair_dist_hist(total_atoms_A, &
                                       & total_atoms_B, &
                                       & cutoff,&
                                       & num_bins,&
                                       & xyz(:, select_indx_A)/10.0d0, &
                                       & xyz(:, select_indx_B)/10.0d0, &
                                       & box/10.0d0,&
                                       & rdf_histogram) 

        end subroutine

    subroutine compute_SA_SCC_lorentzian(num_SCC, & 
                                       & Q_SCC, & 
                                       & SCC, & 
                                       & num_SA, &
                                       & Q_SA, &
                                       & SA, &
                                       & corrlen_SA, & 
                                       & corrlen_SCC) bind(c, name="call_SA_SCC_lorentzian")
        implicit none 
        ! Passed
        integer(c_int), intent(in) :: num_SCC, num_SA 
        real(c_double), intent(in), dimension(1:num_SCC) :: Q_SCC, SCC
        real(c_double), intent(in), dimension(1:num_SA) :: Q_SA, SA

        ! Local
        real(dp) :: r2_SA, r2_SCC, slope_SA, interc_SA, s_beta_1, s_beta_0, slope_SCC, interc_SCC 
        real(dp), dimension(1:num_SA) :: x_SA, y_SA
        real(dp), dimension(1:num_SCC) :: x_SCC, y_SCC

        ! Return
        
        real(c_double), intent(out) :: corrlen_SA, corrlen_SCC 

        y_SA = 1/SA         

        y_SCC = 1/SCC

        x_SCC = Q_SCC **2

        x_SA = Q_SA **2

        call linear_reg(x_SA, y_SA, num_SA, r2_SA, slope_SA, interc_SA, s_beta_1, s_beta_0)
        
        call linear_reg(x_SCC, y_SCC, num_SCC, r2_SCC, slope_SCC, interc_SCC, s_beta_1, s_beta_0)
        
        corrlen_SA = dsqrt(slope_SA / interc_SA) 

        corrlen_SCC = dsqrt(slope_SCC / interc_SCC) 

        end subroutine 

end module  

