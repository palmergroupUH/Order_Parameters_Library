module structure_factor
    use system, only: c_int, c_double, dp, sp, c_bool  
    use constants, only: pi 
    use histogram, only: init_hist, get_bin_id  
    implicit none 
    private 

    ! Note: 
    ! 1. "sf":  structure factor
    ! 2. "ssf": static structure factor 

    public :: static_sf, & 
            & Ashcroft_Langreth_partial_sf, &   
            & compute_fs_corrected_hr, & 
            & Bhatia_Thornton_ssf 

contains

   

    ! compute N given rho and spherical volume with R radius 
    pure real(dp) function N_R(R, rho)
        implicit none
        real(dp), intent(in) :: R, rho

        N_R = 4.0d0/3.0d0 * pi * rho * R * R * R 

        end function 

    ! u(x) is defined in "An analysis of fluctuations in supercooled TIP4P/2005 water"  
    pure function u_x(num_bins, Q, R)
        implicit none 
        integer, intent(in) :: num_bins
        real(dp), intent(in) :: R 
        real(dp), intent(in), dimension(1:num_bins) :: Q 
        real(dp), dimension(1:num_bins) :: QR 
        real(dp), dimension(1:num_bins) :: u_x 

        QR = Q * R

        u_x = 3.0d0 / (QR * QR * QR) * (dsin(QR) - QR * dcos(QR))

         end function 

    real(dp) function S0(partial, same_species, x_spc_1_2, num_bins_gr, num_bins_q, r_mid, gr, hist_ranges, rho, q_ranges, R_ranges, dr, N) 
        implicit none 
        logical, intent(in) :: partial, same_species  
        integer, intent(in) :: num_bins_gr, num_bins_q 
        real(dp), intent(in) :: rho, dr, N, x_spc_1_2
        real(dp), intent(in), dimension(:) :: r_mid, gr, hist_ranges, R_ranges, q_ranges  
   
        logical :: partial_loc
        logical(c_bool) ::  same_species_loc  
        integer :: i, num_R, bin_pos 
        real(dp), dimension(1:num_bins_q) :: Q_mid, S_Q
        real(dp) :: SN_0 

        ! return
        
        num_R = size(R_ranges)  

        S0 = 0.0d0

        partial_loc = partial 

        same_species_loc = same_species 

        do i = 1, num_R 

            bin_pos = get_bin_id(hist_ranges, dr, R_ranges(i))
    
            if (.not. partial_loc) then  
            
                call static_sf(num_bins_gr, num_bins_q, r_mid(1:bin_pos), gr(1:bin_pos), rho, q_ranges, SN_0, Q_mid, S_Q)

            else
                
                call Ashcroft_Langreth_partial_sf(same_species_loc, &
                                               & num_bins_gr, & 
                                               & num_bins_q, & 
                                               & x_spc_1_2, & 
                                               & r_mid(1:bin_pos), &
                                               & gr(1:bin_pos), & 
                                               & rho, & 
                                               & q_ranges, & 
                                               & SN_0, & 
                                               & Q_mid, & 
                                               & S_Q) 
    
            end if 

            S0 = S0 + SN_0 / ( 1.0d0 - 1.0d0 / N * N_R(R_ranges(i), rho)) 

        end do

        S0 = S0 / num_R

        end function  

    ! For homogeneous and isotropic systems
    ! compute static structure factor from gr
    subroutine static_sf(num_bins_gr, num_bins_q, r_mid, gr, bulk_dens, q_ranges, S0, Q_mid, S_Q) bind(c, name="call_static_sf")
        implicit none 
        ! Passed
        integer(c_int), intent(in) :: num_bins_gr, num_bins_q 
        real(c_double), intent(in) :: bulk_dens 
        real(c_double), intent(in), dimension(1:2)  :: q_ranges 
        real(c_double), intent(in), dimension(1:num_bins_gr) :: r_mid 
        real(c_double), intent(in), dimension(1:num_bins_gr) :: gr  

        ! Local
        integer :: iq
        real(dp) :: dr, Krochecker_delta, Q_interval, integrand
        real(dp), dimension(1:num_bins_q) :: SQ_hist

        ! Return
        real(c_double), intent(out) :: S0 
        real(c_double), intent(out), dimension(1:num_bins_q) :: S_Q
        real(c_double), intent(out), dimension(1:num_bins_q) :: Q_mid

        ! generate Q hsitogram

        call init_hist(q_ranges(1), q_ranges(2), num_bins_q, SQ_hist, Q_interval, Q_mid)

        integrand = 0.0d0 

        dr = r_mid(2) - r_mid(1) 

        S_Q = 0.0d0

        do iq = 1, num_bins_q 

            S_Q(iq) = 1 + 4 * pi * bulk_dens * sum(r_mid * (gr - 1) * dr * dsin(Q_mid(iq) * r_mid)/Q_mid(iq))

        end do

        S0  = 1 + 4 * pi * bulk_dens * sum(r_mid * r_mid * (gr - 1) * dr) 

        end subroutine 

    subroutine Ashcroft_Langreth_partial_sf(same_spec, & 
                                           & num_bins_gr, & 
                                           & num_bins_q, & 
                                           & x_spc_1_2, & 
                                           & r_mid, &
                                           & gr, & 
                                           & bulk_dens, & 
                                           & q_ranges, & 
                                           & S0, & 
                                           & Q_mid, & 
                                           & S_Q) bind(c, name="call_AL_static_sf")
        implicit none

        ! Passed
        logical(c_bool), intent(in) :: same_spec 
        integer(c_int), intent(in) :: num_bins_gr, num_bins_q 
        real(c_double), intent(in) :: bulk_dens 
        real(c_double), intent(in) :: x_spc_1_2 
        real(c_double), intent(in), dimension(1:2)  :: q_ranges 
        real(c_double), intent(in), dimension(1:num_bins_gr) :: r_mid 
        real(c_double), intent(in), dimension(1:num_bins_gr) :: gr  

        ! Local
        integer :: iq 
        real(dp) :: dr, Krochecker_delta, Q_interval, integrand
        real(dp), dimension(1:num_bins_q) :: SQ_hist
       
        ! Return

        real(c_double), intent(out) :: S0
        real(c_double), intent(out), dimension(1:num_bins_q) :: Q_mid
        real(c_double), intent(out), dimension(1:num_bins_q) :: S_Q
    
        call init_hist(q_ranges(1), q_ranges(2), num_bins_q, SQ_hist, Q_interval, Q_mid) 

        if (same_spec) then

            Krochecker_delta = 1.0d0 

        else

            Krochecker_delta = 0.0d0  

        end if 
        
        integrand = 0.0d0 

        S_Q = 0.0d0

        dr = r_mid(2) - r_mid(1) 

        do iq = 1, num_bins_q

            S_Q(iq) = Krochecker_delta  &

                   & + 4 * pi * bulk_dens * x_spc_1_2 &

                   & * sum(r_mid * (gr - 1) * dr * dsin(Q_mid(iq) * r_mid)/Q_mid(iq))
    
        end do

        S0  = Krochecker_delta &
            & + 4 * pi * bulk_dens * x_spc_1_2   &
            & * sum(r_mid * r_mid * (gr - 1) * dr) 


        end subroutine 


    subroutine compute_fs_corrected_hr(corrected, &
                                     & partial, &
                                     & same_species, &
                                     & N1N2, & 
                                     & x_spec_1_2, & 
                                     & N, &
                                     & R_size, & 
                                     & R_ranges, &
                                     & Q_ranges, & 
                                     & num_bins_q, &
                                     & num_bins_gr, &
                                     & rho, & 
                                     & r_mid, &
                                     & gr, &
                                     & r2hr) bind(c, name="call_hr_fs_correct")
        implicit none 

        ! Passed
        logical(c_bool), intent(in) :: corrected, same_species, partial  
        integer(c_int), intent(in) :: num_bins_gr, num_bins_q 
        real(c_double), intent(in) :: N1N2, x_spec_1_2 
        real(c_double), intent(in) :: N 
        integer(c_int), intent(in) :: R_size
        real(c_double), intent(in), dimension(1:2) :: q_ranges 
        real(c_double), intent(in), dimension(1:R_size) :: R_ranges 
        real(c_double), intent(in) :: rho 
        real(c_double), intent(in), dimension(1:num_bins_gr) ::  r_mid, gr

        ! Local
        integer :: i 
        logical :: partial_loc, same_species_loc  
        real(dp) :: denom, corrected_factor, dr
        real(dp), dimension(1:num_bins_q) :: Q_mid, S_Q
        real(dp), dimension(1:2) :: hist_ranges  

        ! Return
        real(c_double), intent(out), dimension(1:num_bins_gr) :: r2hr

        dr = r_mid(2) - r_mid(1) 

        partial_loc = partial

        same_species_loc = same_species 

        hist_ranges = [0.0d0, r_mid(num_bins_gr) + (r_mid(2) - r_mid(1)) / 2 ]  

        if (corrected) then
    
            if (partial_loc) then 

                denom = dsqrt(N1N2)

            else

                denom = N 

            end if 
             
            r2hr = r_mid * r_mid * ((gr - 1) + S0(partial_loc, same_species_loc, x_spec_1_2, num_bins_gr, num_bins_q, r_mid, gr, hist_ranges, rho, q_ranges, R_ranges, dr, N) / denom) 

        else
            
            r2hr = r_mid * r_mid *  (gr - 1 )

        end if 

        end subroutine

    subroutine compute_fs_corrected_ssf(corrected, &
                                     & partial, & 
                                     & same_species, & 
                                     & x_spec_1_2, & 
                                     & N, &
                                     & R_size, & 
                                     & R_ranges, &
                                     & Q_ranges, & 
                                     & num_bins_q, &
                                     & num_bins_gr, &
                                     & rho, & 
                                     & r_mid, &
                                     & gr, &
                                     & S0_fs, &    
                                     & Q_mid, & 
                                     & S_Q_R_corrected) bind(c, name="call_ssf_fs_correct")
        implicit none 

        ! Passed
        logical(c_bool), intent(in) :: partial, same_species, corrected
        integer(c_int), intent(in) :: num_bins_gr, num_bins_q 
        real(c_double), intent(in) :: N, x_spec_1_2 
        integer(c_int), intent(in) :: R_size
        real(c_double), intent(in), dimension(1:2) :: q_ranges 
        real(c_double), intent(in), dimension(1:R_size) :: R_ranges 
        real(c_double), intent(in) :: rho 
        real(c_double), intent(in), dimension(1:num_bins_gr) ::  r_mid, gr

        ! Local
        logical :: partial_loc, same_species_loc 
        real(dp) :: dr, correct_factor, SN_0, end_hist
        integer :: i, num_R 
        real(dp), dimension(1:2) :: hist_ranges  
        real(dp), dimension(1:num_bins_q) :: S_Q

        ! Return
        real(c_double), intent(out) :: S0_fs 
        real(c_double), intent(out), dimension(1:num_bins_q) :: Q_mid 
        real(c_double), intent(out), dimension(1:num_bins_q) :: S_Q_R_corrected
        
        dr = r_mid(2) - r_mid(1)

        end_hist = r_mid(num_bins_gr) + dr / 2
       
        hist_ranges = [0.0d0, end_hist ]  
        
        S_Q_R_corrected = 0

        partial_loc = partial 

        same_species_loc = same_species
    
        if (corrected) then 

            correct_factor = 1.0d0 

        else

            correct_factor = 0.0d0  

        end if 
       
        num_R = size(R_ranges)
 
        S0_fs = S0(partial_loc, same_species_loc, x_spec_1_2, num_bins_gr, num_bins_q, r_mid, gr, hist_ranges, rho, q_ranges, R_ranges, dr, N)
        do i = 1, num_R
    
            if (.not. partial_loc) then

                call static_sf(num_bins_gr, num_bins_q, r_mid, gr, rho, q_ranges, SN_0, Q_mid, S_Q)
    
            else

                call Ashcroft_Langreth_partial_sf(same_species, & 
                                               & num_bins_gr, & 
                                               & num_bins_q, & 
                                               & x_spec_1_2, & 
                                               & r_mid, &
                                               & gr, & 
                                               & rho, & 
                                               & q_ranges, & 
                                               & SN_0, & 
                                               & Q_mid, & 
                                               & S_Q) 

            end if 

            S_Q_R_corrected = S_Q_R_corrected & 
                           & +  S_Q  &
                           & + S0_fs/ N * N_R(R_ranges(i), rho) * u_x(num_bins_q, Q_mid, R_ranges(i)) * correct_factor 
       
        end do  

        S_Q_R_corrected = S_Q_R_corrected / num_R
        
        end subroutine 

    subroutine Bhatia_Thornton_ssf(num_bins_q, &
                                 & x_H_frac, & 
                                 & x_L_frac, &
                                 & SHH, &
                                 & SLL, &
                                 & SHL, &
                                 & SNN, &
                                 & SCC, &
                                 & SNC, & 
                                 & theta_k, &  
                                 & SA) bind(c, name="call_Bhatia_Thornton")
        implicit none 
        ! Passed
        integer(c_int), intent(in) :: num_bins_q 
        real(c_double), intent(in) :: x_H_frac, x_L_frac  

        ! partial structure factor between high q - high q, low q - low q and
        ! high - low q  
        real(c_double), intent(in), dimension(1:num_bins_q) :: SHH, SLL, SHL 
       
        ! Local

        ! Return
        ! structure factor: density-density, concentration - concentration,
        ! density - concentration: 

        real(c_double), intent(out), dimension(1:num_bins_q) :: SNN, SCC, SNC, theta_k, SA 

        ! see Overduin and Patey, "Understanding the Structure Factor and
        ! Isothermal Compressibility of Ambient Water in Terms of Local Structural Environments"
        
        SNN = x_H_frac * SHH + x_L_frac * SLL + 2.0d0 * dsqrt(x_H_frac * x_L_frac) * SHL  

        SCC = x_H_frac * x_L_frac *(x_L_frac * SHH + x_H_frac * SLL - 2.0d0 * dsqrt(x_H_frac * x_L_frac) * SHL)

        SNC = x_H_frac * x_L_frac * (SHH - SLL + (x_L_frac - x_H_frac) / dsqrt(x_H_frac * x_L_frac) * SHL)

        theta_k = x_H_frac * x_L_frac * (SHH * SLL - SHL * SHL) / SCC

        !SA = (-SNC/SCC)**2* SCC 
 
        SA =  SNN - theta_k  

        end subroutine

end module 
