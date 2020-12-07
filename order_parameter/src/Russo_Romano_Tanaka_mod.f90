module Russo_Tanaka
    use system, only: c_int, c_double, dp, sp, convert_c_string_f_string 

    use nb_list, only: build_homo_neighbor_list, &
                     & check_nb_list, &
                     & apply_nearest_neighbor_crit, &
                     & rebuild_neighbor_list_nnb

    use sph_op, only: compute_ave_qlm, & 
                    & coarse_grained_qlm, &
                    & compute_nnb_dij_bond, &
                    & compute_dij_bond, & 
                    & label_crystal_like, & 
                    & local_bond_order, & 
                    & initialize_wigner3j, &
                    & calc_Wl  

    use cluster, only: largest_cluster_high_mem , & 
                     & Allen_Tidesley_cluster 

    implicit none 

    ! Default: all subroutines and functions should be visible only in this module  
    private 

    ! Export following subroutines/functions
    ! "public" make them accessible to other Fortran programs through "use" 
    ! otherwise all subroutines and variables in this module 
    ! can not be used by other programs
    
    public :: Russo_Romano_Tanaka

contains 

    subroutine Russo_Romano_Tanaka(total_atoms, &
                                  & box, &
                                  & xyz, &
                                  & maxnb, & 
                                  & nnb, &
                                  & l, &
                                  & ql_cutoff_sqr, &
                                  & n_bonds, & 
                                  & nxtl, & 
                                  & label, &
                                  & Ql_invar, &
                                  & Wl_invar) bind(c, name="call_Russo_Romano_Tanaka")
        implicit none 

        ! Passed:
        
        ! total number of particles 
        integer(c_int), intent(in) :: total_atoms 
        ! maximum number of neighbors
        integer(c_int), intent(in) :: maxnb 
        ! number of nearest neigbhors
        integer(c_int), intent(in) :: nnb 
        ! the order of spherical harmonics
        integer(c_int), intent(in) :: l 
        ! the number of crystal-like bonds required to be defined as crystal
        integer(c_int), intent(in) :: n_bonds
        ! cutoff for computing spherical harmonics
        real(c_double), intent(in) :: ql_cutoff_sqr
        
        ! system size
        real(c_double), intent(in), dimension(1:3) :: box
        ! coordinates
        real(c_double), intent(in), dimension(1:3, total_atoms) :: xyz

        ! Local:
        character(len=:), allocatable :: bond_crit
        integer, dimension(1:total_atoms) :: num_NB_list 
        integer, dimension(1:maxnb, 1:total_atoms) :: NB_list 
        real(dp), dimension(1:4, 1:maxnb, 1:total_atoms) :: Rij
        integer, dimension(1:total_atoms) :: n_dij_bonds
        integer, dimension(1:total_atoms) :: crystal_id
        real(dp), dimension(1:nnb, 1:total_atoms) :: dij 
        complex(dp), dimension(-l:l, total_atoms) :: qlm, Ql_cg
        real(dp) :: start, finish  
        real(dp) :: interval, upper, lower
        integer :: k, counter, l_2 

        ! Return:
        integer(c_int), intent(out) :: nxtl 
        integer(c_int), intent(inout), dimension(1:total_atoms) :: label 
        real(c_double), intent(inout), dimension(1:total_atoms) :: Ql_invar 
        real(c_double), intent(inout), dimension(1:total_atoms) :: Wl_invar
        
        ! get the neighbor list
        call build_homo_neighbor_list(total_atoms, maxnb, ql_cutoff_sqr, xyz, box, num_NB_list, NB_list, Rij) 
        
        ! check the neighbor list
        call check_nb_list(num_NB_list, maxnb, nnb) 

        ! apply the nearest neighbor
        call apply_nearest_neighbor_crit(nnb, total_atoms, num_NB_list, NB_list, Rij)
        
        ! compute averaged spherical harmonics with its neighbors, qlm
        ! choose the mode "optimized" or "general" 
        ! if optimized version is not implemented, the "general" will be called
        call compute_ave_qlm(total_atoms, l, num_NB_list, Rij, "optimized", qlm)

        ! This style name is used to select the crystal bonds criterion
        ! It can be modified from the source code in the file "sph_mod.f90"
        bond_crit = "RussoTanaka" 

        ! compute the coarse-grained Ql
        call coarse_grained_qlm(total_atoms, l, num_NB_list, NB_list, qlm, Ql_cg)
        
        call compute_nnb_dij_bond(bond_crit, total_atoms, nnb, l, NB_list, n_bonds, Ql_cg, nxtl, dij, n_dij_bonds)

        call label_crystal_like(total_atoms, nxtl, n_bonds, n_dij_bonds, crystal_id) 
       
        ! An order parameter scheme defined by Russo et al.

        l_2 = 4  

        call compute_Ql_Wl(total_atoms,& 
                         & l_2, & 
                         & nxtl, & 
                         & crystal_id, &  
                         & num_NB_list, & 
                         & NB_list, & 
                         & Rij, & 
                         & label, & 
                         & Wl_invar, & 
                         & Ql_invar)
    

        end subroutine 

    subroutine compute_Ql_Wl(total_atoms,& 
                           & l, & 
                           & nxtl, & 
                           & crystal_id, &  
                           & num_NB_list, & 
                           & NB_list, & 
                           & Rij, & 
                           & label, & 
                           & Wl_invar, & 
                           & Ql_invar)
        implicit none 
        ! Passed
        integer, intent(in) :: total_atoms   
        integer, intent(in) :: l   
        integer, intent(in) :: nxtl  
        integer, intent(in), dimension(:) :: num_NB_list 
        integer, intent(in), dimension(:,:) :: NB_list 
        integer, intent(in), dimension(1:total_atoms) :: crystal_id  
        real(dp), intent(in), dimension(:,:,:) :: Rij

        ! Local
        integer  :: num_pairs 
        complex(dp), dimension(-l:l, total_atoms) :: ql 
        complex(dp), dimension(-l:l, total_atoms) :: Ql_cg
        real(dp), dimension(:), allocatable :: wigner_3j_symobol_mat 
        integer, dimension(:, :), allocatable :: m1_m2_m3_sum0_mat 
        integer :: i

        ! Return
        integer, intent(out), dimension(1:total_atoms) :: label
        real(dp), intent(out), dimension(1:total_atoms) :: Wl_invar
        real(dp), intent(out), dimension(1:total_atoms) :: Ql_invar
        
        call initialize_wigner3j(l, num_pairs, m1_m2_m3_sum0_mat, wigner_3j_symobol_mat)  
        
        call compute_ave_qlm(total_atoms, & 
                           & l, &  
                           & num_NB_list, & 
                           & Rij, &
                           & "optimized", &
                           & ql)              

        call coarse_grained_qlm(total_atoms, l, num_NB_list, NB_list, ql, Ql_cg) 

        call local_bond_order(total_atoms, l, Ql_cg, Ql_invar) 
        
        call calc_Wl(num_pairs,total_atoms,l,m1_m2_m3_sum0_mat,wigner_3j_symobol_mat, Ql_cg, Wl_invar)
        
        call label_ice_identity(total_atoms, crystal_id, nxtl, Ql_invar, Wl_invar, label)

        end subroutine 

    subroutine label_ice_identity(total_atoms, crystal_id, num_crystal, Ql, Wl, label)
        implicit none 
        integer, intent(in) :: total_atoms 
        integer, intent(in) :: num_crystal 
        integer, intent(in), dimension(1:total_atoms) :: crystal_id  
        real(dp), intent(in), dimension(1:total_atoms) :: Ql, Wl  
        integer :: i 
        integer, intent(out), dimension(1:total_atoms) :: label

        label = 0
        
        do i = 1, total_atoms 
        
            if (any(i == crystal_id(1:num_crystal))) then
            
                ! clathrate 
                if (Ql(i) < 0.05) then 
        
                    label(i) = 1

                end if

                ! intermediate ice
                if (0.05 < Ql(i) .and. Ql(i) < 0.11) then
        
                    label(i) = 2                 

                end if 

                ! Ice Hexagonal 
                if (Ql(i) > 0.11 .and. Wl(i) > 0) then   
            
                    label(i) = 3

                end if 
                     
                ! Ice 0 
                if (Ql(i) < 0.145 .and. Wl(i) < 0) then 
            
                    label(i) = 4

                end if 

                ! Ice cubic
                if (Ql(i) > 0.145 .and. Wl(i) < 0) then 

                    label(i) = 5

                end if 

            else 

                label(i) = 6 

            end if 

        end do
        
        end subroutine 

end module 
