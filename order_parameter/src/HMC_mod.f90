module HMC
    use system, only: c_int, c_double, dp, sp, convert_c_string_f_string 

    use nb_list, only: build_homo_neighbor_list, &
                     & check_nb_list, &
                     & apply_nearest_neighbor_crit, &
                     & rebuild_neighbor_list_nnb

    use sph_op, only: compute_ave_qlm, & 
                    & coarse_grained_qlm, &
                    & compute_nnb_dij_bond, &
                    & compute_dij_bond, & 
                    & label_crystal_like  

    use cluster, only: largest_cluster_high_mem , & 
                     & Allen_Tidesley_cluster 

    implicit none 

    ! Default: all subroutines and functions should be visible only in this module  
    private 

    ! Export following subroutines/functions
    ! "public" make them accessible to other Fortran programs through "use" 
    ! otherwise all subroutines and variables in this module 
    ! can not be used by other programs
    
    public :: calc_Ql_lcluster

contains 

    subroutine calc_Ql_lcluster(bond_crit_c, &
                              & style_leng, &
                              & total_atoms, & 
                              & box, &
                              & xyz, &
                              & maxnb, & 
                              & nnb, &
                              & l, &
                              & ql_cutoff_sqr, &
                              & cluster_cutoff, &
                              & n_bonds, & 
                              & nxtl, &
                              & largest_cluster) bind(c, name="call_calc_ql_cluster")
        implicit none 

        ! Passed:
        ! the length of keyword
        integer(c_int), intent(in) :: style_leng 
        ! the style of order parameter
        character(len=1), intent(in), dimension(1:style_leng) :: bond_crit_c 
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
        ! cutoff for cluster 
        real(c_double), intent(in) :: cluster_cutoff
        ! system size
        real(c_double), intent(in), dimension(1:3) :: box
        ! coordinates
        real(c_double), intent(in), dimension(1:3, total_atoms) :: xyz

        ! Local:
        character(len=:), allocatable :: bond_crit

        integer, dimension(1:total_atoms) :: num_NB_list 
        integer, dimension(1:maxnb, 1:total_atoms) :: NB_list 
        real(dp), dimension(1:4, 1:maxnb, 1:total_atoms) :: Rij
        real(dp), dimension(1:total_atoms) :: n_dij_bonds
        complex(dp), dimension(-l:l, total_atoms) :: qlm, Ql_cg
        real(dp) :: start, finish  
        real(dp) :: interval, upper, lower
        integer :: k, counter 

        ! Return:
        integer(c_int), intent(out) :: nxtl, largest_cluster 
        
        ! convert C string to Fortran string
        call convert_c_string_f_string(bond_crit_c, style_leng, bond_crit)  

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

        ! choose order parameter style
        
        if (bond_crit == "RussoTanaka") then

            call RussoTanaka(total_atoms, &
                           & l, & 
                           & maxnb, &
                           & nnb, & 
                           & cluster_cutoff, & 
                           & n_bonds, & 
                           & num_NB_list, & 
                           & NB_list, & 
                           & Rij, &
                           & qlm, & 
                           & nxtl, & 
                           & largest_cluster, &
                           & cluster_cutoff, &
                           & xyz, &
                           & box) 

        else if (bond_crit =="ReinhardtDoye") then

            call ReinhardtDoye(total_atoms, &
                             & l, & 
                             & maxnb, &
                             & n_bonds, & 
                             & num_NB_list, & 
                             & NB_list, & 
                             & Rij, & 
                             & qlm, & 
                             & nxtl, & 
                             & largest_cluster, &
                             & cluster_cutoff, & 
                             & xyz, & 
                             & box) 
        
        end if 

        end subroutine 

    subroutine RussoTanaka(total_atoms, &
                         & l, & 
                         & maxnb, &
                         & nnb, &
                         & cluster_cutoff, &
                         & n_bonds, & 
                         & num_NB_list, & 
                         & NB_list, & 
                         & Rij, & 
                         & qlm, & 
                         & nxtl, &
                         & largest_cluster, &
                         & crys_cutoff, &
                         & xyz, &
                         & box)
    
        implicit none 
        ! Passed
        integer, intent(in) :: total_atoms ! total number of particles
        integer, intent(in) :: l ! the order of spherical harmonics
        integer, intent(in) :: maxnb ! maximum number of neighbors
        integer, intent(in) :: nnb ! number of nearest neigbhors
        integer, intent(in) :: n_bonds
        real(dp), intent(in) :: cluster_cutoff ! cutoff
        real(dp), intent(in), dimension(:, :, :) :: Rij
        complex(dp), intent(in), dimension(-l:l, total_atoms) :: qlm

        ! Optional variables below are used in Allen_Tidesley_cluster  
        ! They can be removed from the argument list of called function
        real(dp), intent(in), optional :: crys_cutoff
        real(dp), intent(in), dimension(1:3, 1:total_atoms), optional :: xyz
        real(dp), intent(in), dimension(1:3), optional :: box

        ! Local 
        character(len=:), allocatable :: bond_crit
        complex(dp), dimension(-l:l, total_atoms) :: Ql_cg
        real(dp), dimension(1:nnb, 1:total_atoms) :: dij
        integer, dimension(1:total_atoms) :: n_dij_bonds 
        real(dp) :: crys_cutoff_sqr, start, finish  
        integer, dimension(1:total_atoms) :: crystal_id 

        ! Return:
        integer, intent(inout), dimension(:) :: num_NB_list
        integer, intent(inout), dimension(:, :) :: NB_list
        integer, intent(out) :: nxtl, largest_cluster 

        ! This style name is used to select the crystal bonds criterion
        ! It can be modified from the source code in the file "sph_mod.f90"
        bond_crit = "RussoTanaka" 

        ! compute the coarse-grained Ql
        call coarse_grained_qlm(total_atoms, l, num_NB_list, NB_list, qlm, Ql_cg)

        call compute_nnb_dij_bond(bond_crit, total_atoms, nnb, l, NB_list, n_bonds, Ql_cg, nxtl, dij, n_dij_bonds)
       
        ! To use largest_cluster_high_mem code:
        ! Uncomment the following two routines

        call rebuild_neighbor_list_nnb(total_atoms, cluster_cutoff, maxnb, num_NB_list, NB_list, Rij)
        
        call largest_cluster_high_mem(total_atoms, num_NB_list, NB_list, n_bonds, n_dij_bonds, largest_cluster) 

        ! To use Allen-Tidesely code:
        ! 1. Uncomment the following two routines: "label_crystal_like" and "Allen_Tidesley_cluster"
        ! 2. Also, make sure the xyz, box and crystal_cutoff are passed as an argument 
        
        if (present(crys_cutoff)) then 

            crys_cutoff_sqr = crys_cutoff * crys_cutoff 

        end if 
        
        !call label_crystal_like(total_atoms, nxtl, n_bonds, n_dij_bonds, crystal_id)

        !call Allen_Tidesley_cluster(nxtl, xyz(:, crystal_id(1:nxtl)), box, crys_cutoff_sqr, largest_cluster)

        end subroutine 

                
    subroutine ReinhardtDoye(total_atoms, &
                           & l, & 
                           & maxnb, &
                           & n_bonds, & 
                           & num_NB_list, & 
                           & NB_list, & 
                           & Rij, &
                           & qlm, & 
                           & nxtl, &
                           & largest_cluster,&
                           & crys_cutoff, &
                           & xyz, &
                           & box)
    
        implicit none 
        ! Passed
        integer, intent(in) :: total_atoms ! total number of particles
        integer, intent(in) :: l ! the order of spherical harmonics
        integer, intent(in) :: maxnb ! maximum number of neighbors
        integer, intent(in) :: n_bonds
        real(dp), intent(in), dimension(:, :, :) :: Rij 
        complex(dp), intent(in), dimension(-l:l, total_atoms) :: qlm
    
        ! Passed optional variables below to use Allen_Tidesley_cluster  
        ! They can be removed from the argument of calling function
        real(dp), intent(in), optional :: crys_cutoff
        real(dp), intent(in), dimension(1:3, 1:total_atoms), optional :: xyz
        real(dp), intent(in), dimension(1:3), optional :: box

        ! Local 
        character(len=:), allocatable :: bond_crit
        real(dp), dimension(1:maxnb, 1:total_atoms) :: dij
        integer, dimension(1:total_atoms) :: n_dij_bonds 
        real(dp) :: crys_cutoff_sqr 
        integer, dimension(1:total_atoms) :: crystal_id 
        
        ! Return:
        integer, intent(out) :: nxtl, largest_cluster 
        integer, intent(inout), dimension(:) :: num_NB_list 
        integer, intent(inout), dimension(:, :) :: NB_list 

        ! This style name is used to select the crystal bonds criterion
        ! It can be modified from the source code in the file "sph_mod.f90"

        bond_crit = "ReinhardtDoye" 

        call compute_dij_bond(bond_crit, total_atoms, maxnb, num_NB_list, l, NB_list, n_bonds, qlm, nxtl, dij, n_dij_bonds)
        
        ! To use largest_cluster_high_mem code:
        ! Uncomment the following two routines: "rebuild_neighbor_list_nnb", and
        ! "largest_cluster_high_mem"

        call rebuild_neighbor_list_nnb(total_atoms, crys_cutoff, maxnb, num_NB_list, NB_list, Rij) 

        call largest_cluster_high_mem(total_atoms, num_NB_list, NB_list, n_bonds, n_dij_bonds, largest_cluster)

        ! To use Allen-Tidesely code:
        ! 1. Uncomment the following two routines: "label_crystal_like" and "Allen_Tidesley_cluster"
        ! 2. Also, make sure the "xyz", "box" and "crys_cutoff" are passed as an argument 

        if (present(crys_cutoff)) then 

            crys_cutoff_sqr = crys_cutoff * crys_cutoff 

        end if 
        
        !call label_crystal_like(total_atoms, nxtl, n_bonds, n_dij_bonds, crystal_id)

        !call Allen_Tidesley_cluster(nxtl, xyz(:, crystal_id(1:nxtl)), box, crys_cutoff_sqr, largest_cluster)


        end subroutine

end module 
