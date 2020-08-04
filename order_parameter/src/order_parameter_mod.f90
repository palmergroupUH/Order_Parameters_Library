module order_parameter
    use sorting, only: quick_sort 
    use system, only: c_int, c_double, dp, sp, error_unit 
    use constants, only: pi
    use spherical_harmonics, only: compute_spherical_harmonics, &
                                 & compute_spherical_harmonics_const,&
                                 & compute_associated_legendre_const,& 
                                 & optimized_Y12, & 
                                 & optimized_Y4

    use wigner3j_symbol, only: wigner_3j  
    implicit none
    private
    public :: initialize_Ylm, & 
            ! Compute spherical harmonics
            & compute_ave_Ylm,& 
            & compute_optimized_q4, & 
            & compute_optimized_q12, & 
            ! Compute the second/third order invariants 
            & local_bond_order, & 
            & global_bond_order, & 
            & build_homo_neighbor_list, &
            & rebuild_neighbor_list_nnb, & 
            & apply_nearest_neighbor_crit, &
            & check_nb_list, &
            & neighbor_averaged_qlm, &
            & nnb_dot_product, &
            & nnb_dot_product_crystal, &
            & select_crystalline_particles, & 
            & calc_Wl, & 
            & convert_c_string_f_string
              
contains

! -------------------------------------------------------------------------------------------------
!                                  Neighbor List
! -------------------------------------------------------------------------------------------------

    subroutine build_homo_neighbor_list(total_atoms, maxnb, cutoff_sqr_nb, xyz, box, num_NB_list, NB_list, Rij)
        implicit none 
        ! Passed
        integer,intent(in) :: maxnb
        integer,intent(in) :: total_atoms
        real(dp),intent(in) :: cutoff_sqr_nb
        real(dp),intent(in),dimension(:,:) :: xyz
        real(dp),intent(in),dimension(:) :: box

        ! Local
        integer :: iatom,jatom
        real(dp) :: distance_ij, rij_sqr
        real(dp),dimension(1:size(box)) :: xyz_temp, xyz_separate

        ! Return
        integer,intent(out),dimension(1:total_atoms) :: num_NB_list
        integer,intent(out),dimension(1:maxnb,1:total_atoms) ::  NB_list
        real(dp),intent(out),dimension(1:4,1:maxnb,1:total_atoms) ::Rij

        num_NB_list = 0 
    
        NB_list = 0 

        Rij = 0.0d0

        do iatom = 1, total_atoms - 1

            xyz_temp = xyz(:,iatom)

            do jatom = iatom + 1, total_atoms

                xyz_separate = xyz_temp - xyz(:,jatom )

                ! periodical boundary condition
                xyz_separate = xyz_separate - box*dnint(xyz_separate/box)

                rij_sqr = xyz_separate(1)*xyz_separate(1) & 
                      & + xyz_separate(2)*xyz_separate(2) &
                      & + xyz_separate(3)*xyz_separate(3)

                if ( rij_sqr < cutoff_sqr_nb ) then

                    distance_ij = dsqrt(rij_sqr)

                    num_NB_list(iatom) = num_NB_list(iatom) + 1

                    num_NB_list(jatom) = num_NB_list(jatom) + 1

                    NB_list(num_NB_list(iatom),iatom) = jatom

                    NB_list(num_NB_list(jatom),jatom) = iatom

                    Rij(1:3,num_NB_list(iatom),iatom) = xyz_separate

                    Rij(4,num_NB_list(iatom),iatom) = distance_ij

                    Rij(1:3,num_NB_list(jatom),jatom) = -xyz_separate

                    Rij(4,num_NB_list(jatom),jatom) = distance_ij

                end if

            end do

        end do
 
        end subroutine

    subroutine rebuild_neighbor_list_nnb(total_atoms, &
                                       & cutoff, &
                                       & maxnb, &
                                       & num_NB_list, &
                                       & NB_list, &
                                       & Rij)
        implicit none 

        ! Passed 
        integer, intent(in) :: total_atoms 
        integer, intent(in) :: maxnb 
        real(dp), intent(in) :: cutoff 
        real(dp), intent(in), dimension(:,:,:) :: Rij 

        ! Local
        integer, dimension(1:maxnb) :: temp_nb_lst 
        integer :: count_nb
        integer :: iatom, jnb 
        real(dp) :: dist 

        ! Return
        integer, intent(inout), dimension(1:total_atoms) :: num_NB_list
        integer, intent(inout), dimension(1:maxnb, 1:total_atoms) :: NB_list

        do iatom = 1, total_atoms

            count_nb = 0

            temp_nb_lst = 0 

            do jnb = 1, maxnb

                dist = Rij(4,jnb, iatom)

                if (dist <= cutoff .and. dist /= 0.0d0) then

                    count_nb = count_nb + 1

                    temp_nb_lst(count_nb) = NB_list(jnb, iatom)

                end if

            end do 

            num_NB_list(iatom) = count_nb

            NB_list(:, iatom) = temp_nb_lst
    
        end do 
        
        end subroutine 

    subroutine check_nb_list(num_NB_list,maxnb,nnb)
        implicit none
        integer,intent(in),dimension(:) :: num_NB_list
        integer,intent(in) :: maxnb, nnb

        ! if Nothing wrong, continues
        if (all(num_NB_list <= maxnb) .AND. & 
          & all(num_NB_list >= nnb)   .AND. & 
          & all(num_NB_list > 0)) then

            continue

        else if ( any(num_NB_list > maxnb) ) then

            write(error_unit,*) "Error! Maximum numbers of neighbors is larger than maxnb"
            write(error_unit,*) "Check the followings: "
            write(error_unit,*) "1. Cutoff is too much ? "
            write(error_unit,*) "2. Input maxnb is too small ? "
            write(error_unit,*) "3. The atoms are read from trajectories & 
                                & correctly (different kinds of species & 
                                & used properly) ? "

            STOP "Decrease the cutoff or increase maxnb"

        else if (any(num_NB_list < nnb)) THEN

            write(error_unit,*) "Minimum number of atoms is:", & 
                                & MinVAL(num_NB_list), & 
                                & "The required nearest neighbors is:",nnb

            STOP "Some atoms may not have enough neighbors to meet the nearest neighbors requirement"

        end if

        if (any(num_NB_list ==0)) then

            write(error_unit,*) "Some atoms does not have any neighbors at this cutoff "
            stop

        end if

        IF ( any(maxnb - num_NB_list  < 5) ) THEN

            write(error_unit,*) "WARNinG: Maxnb is only less than 5 & 
                                & bigger than maximum numberof neighbors"

        end IF

        end subroutine

    subroutine use_nearest_neighbors(nnb, total_atoms, num_NB_list, NB_list, Rij)
        implicit none

        ! Passed
        integer, intent(in) :: nnb,total_atoms

        ! Local
        integer :: iatom,jatom,num_NB
        real(dp),dimension(1:3)   :: xyz
        integer,dimension(:),allocatable :: tracking_index

        ! Return
        integer, intent(inout), dimension(:) :: num_NB_list
        integer, intent(inout),dimension(:,:) :: NB_list 
        real(dp), intent(inout),dimension(:,:,:) :: Rij 

        DO iatom = 1,total_atoms

            ! search the nearest neighbors and return index

            num_NB = num_NB_list(iatom)
            
            call quick_sort(num_NB, Rij(4,1:num_NB,iatom),1,num_NB, tracking_index)
            
            NB_list(1:num_NB,iatom) = NB_list(tracking_index,iatom)

            Rij(1:3,1:num_NB,iatom) = Rij(1:3,tracking_index,iatom)
    
        end do

        num_NB_list = nnb

        end subroutine

    subroutine check_enough_neighbors(nnb, num_NB_list)
        implicit none 
        integer, intent(in) :: nnb 
        integer, intent(in), dimension(:) :: num_NB_list
        integer :: iatom

        if (any(num_NB_list - nnb < 0)) then 

            write(error_unit,*) "Some atoms may not have enough neighbors & 
                                & to meet the nearest neighbors requirement"

            stop    

        end IF 

        end subroutine 


    subroutine apply_nearest_neighbor_crit(nnb, total_atoms, num_NB_list, NB_list, Rij)
        implicit none 
        integer,intent(in) :: nnb,total_atoms  
        integer, intent(inout), dimension(:) :: num_NB_list
        integer, intent(inout),dimension(:,:) :: NB_list 
        real(dp), intent(inout),dimension(:,:,:) :: Rij 

        select case(nnb)   

            case(0) 

                ! No need to nearest neighbors scheme and return 
                continue

            case(1:) 

                call check_enough_neighbors(nnb, num_NB_list) 
                
                ! apply nearest neighbors scheme    
                call use_nearest_neighbors(nnb,total_atoms, num_NB_list, NB_list, Rij) 

            case(:-1) 

                stop "number of nearest neighbors can not be negative "
    
        end select  

        end subroutine 

! -------------------------------------------------------------------------------------------------
!                                  Spherical Harmonics
! -------------------------------------------------------------------------------------------------


    subroutine initialize_Ylm(l, sph_const, Plm_const)
        implicit none 
        integer, intent(in) :: l
        real(dp), intent(out), dimension(-l:l) :: sph_const
        real(dp), intent(out), dimension(0:l,4) :: Plm_const

        call compute_spherical_harmonics_const(l, sph_const)

        call compute_associated_legendre_const(l, Plm_const) 

        end subroutine 

    ! General spherical harmonics for any order 
    ! Optimized Y12 and Y4 are automatically used
    ! at l = 12 or l = 4
    subroutine compute_ave_Ylm(sph_const,&
                              & Plm_const,&
                              & total_atoms,& 
                              & l, &
                              & num_NB_list, &
                              & Rij,& 
                              & Ylm_ave)
        implicit none 
        ! Passed:
        integer, intent(in) :: l, total_atoms 
        real(dp), intent(in), dimension(-l:l) :: sph_const  
        real(dp),intent(in),dimension(0:l,4) :: Plm_const 
        integer,intent(in),dimension(:) :: num_NB_list 
        real(dp),intent(in),dimension(:,:,:) :: Rij
        
        ! Local:
        integer :: iatom, jnb, num_NB
        complex(dp), dimension(-l:l) :: Ylm, Ylm_sum_nb
        real(dp) :: rcx,rcy,rcz,rl
        real(dp) :: sin_theta, cos_theta, sin_phi, cos_phi
      
        ! Return:
        complex(dp), intent(out), dimension(-l:l,total_atoms) :: Ylm_ave

        do iatom = 1, total_atoms
    
            Ylm_sum_nb = 0.0d0

            num_NB = num_NB_list(iatom)

            do jnb = 1, num_NB  

                rcx = Rij(1,jnb,iatom) 
    
                rcy = Rij(2,jnb,iatom)

                rcz = Rij(3,jnb,iatom) 

                rl = Rij(4,jnb,iatom)

                cos_theta = rcz/rl 

                sin_theta = dsqrt(1-cos_theta*cos_theta)

                if (sin_theta /= 0) then 

                    cos_phi = rcx/(rl*sin_theta)

                    sin_phi = rcy/(rl*sin_theta)

                else 
    
                    cos_phi = 1.0d0
        
                    sin_phi = 0.0d0

                end if 

                if (l /= 12 .and. l /=4 ) then
                
                    call compute_spherical_harmonics(sph_const, Plm_const, l, cos_theta, cos_phi, sin_phi, Ylm)
    
                else if (l==12) then
                    
                    call optimized_Y12(sin_theta, cos_theta, cos_phi, sin_phi, Ylm)
            
                else if (l==4) then 
                    
                    call optimized_Y4(sin_theta, cos_theta, cos_phi, sin_phi, Ylm)

                end if
    
                Ylm_sum_nb = Ylm_sum_nb + Ylm

            end do

            Ylm_ave(:,iatom) = Ylm_sum_nb/num_NB

        end do 

        end subroutine

    ! Optimized spherical harmonics 
    ! Currently, only l = 12 and l = 4 is supported:
    ! These should give the same results as "compute_ave_Ylm"
    ! at the given order. 
    subroutine compute_optimized_q12(total_atoms, & 
                                   & num_NB_list, &
                                   & Rij,& 
                                   & Ylm_ave)
        implicit none 
        ! Passed:
        integer, intent(in) :: total_atoms 
        integer,intent(in),dimension(:) :: num_NB_list 
        real(dp),intent(in),dimension(:,:,:) :: Rij
        
        ! Local:
        integer :: iatom, jnb, num_NB
        complex(dp), dimension(-12:12) :: Ylm, Ylm_sum_nb
        real(dp) :: rcx,rcy,rcz,rl
        real(dp) :: sin_theta, cos_theta, sin_phi, cos_phi
      
        ! Return:
        complex(dp), intent(out), dimension(-12:12,total_atoms) :: Ylm_ave

        do iatom = 1, total_atoms
    
            Ylm_sum_nb = 0.0d0

            num_NB = num_NB_list(iatom)

            do jnb = 1, num_NB  

                rcx = Rij(1,jnb,iatom) 
    
                rcy = Rij(2,jnb,iatom)

                rcz = Rij(3,jnb,iatom) 

                rl = Rij(4,jnb,iatom)

                cos_theta = rcz/rl 

                sin_theta = dsqrt(1-cos_theta*cos_theta)

                if (sin_theta /= 0) then 

                    cos_phi = rcx/(rl*sin_theta)

                    sin_phi = rcy/(rl*sin_theta)

                else 
    
                    cos_phi = 1.0d0
        
                    sin_phi = 0.0d0

                end if 

                call optimized_Y12(sin_theta, cos_theta, cos_phi, sin_phi, Ylm) 

                Ylm_sum_nb = Ylm_sum_nb + Ylm 

            end do 

            Ylm_ave(:,iatom) = Ylm_sum_nb/num_NB 

        end do

        end subroutine 
                
    subroutine compute_optimized_q4(total_atoms, & 
                                  & num_NB_list, &
                                  & Rij,& 
                                  & Ylm_ave)
        implicit none 
        ! Passed:
        integer, intent(in) :: total_atoms 
        integer,intent(in),dimension(:) :: num_NB_list 
        real(dp),intent(in),dimension(:,:,:) :: Rij
        
        ! Local:
        integer :: iatom, jnb, num_NB
        complex(dp), dimension(-4:4) :: Ylm, Ylm_sum_nb
        real(dp) :: rcx,rcy,rcz,rl
        real(dp) :: sin_theta, cos_theta, sin_phi, cos_phi
      
        ! Return:
        complex(dp), intent(out), dimension(-4:4,total_atoms) :: Ylm_ave

        do iatom = 1, total_atoms
    
            Ylm_sum_nb = 0.0d0

            num_NB = num_NB_list(iatom)

            do jnb = 1, num_NB  

                rcx = Rij(1,jnb,iatom) 
    
                rcy = Rij(2,jnb,iatom)

                rcz = Rij(3,jnb,iatom) 

                rl = Rij(4,jnb,iatom)

                cos_theta = rcz/rl 

                sin_theta = dsqrt(1-cos_theta*cos_theta)

                if (sin_theta /= 0) then 

                    cos_phi = rcx/(rl*sin_theta)

                    sin_phi = rcy/(rl*sin_theta)

                else 
    
                    cos_phi = 1.0d0
        
                    sin_phi = 0.0d0

                end if 

                call optimized_Y4(sin_theta, cos_theta, cos_phi, sin_phi, Ylm) 

                Ylm_sum_nb = Ylm_sum_nb + Ylm 

            end do 

            Ylm_ave(:,iatom) = Ylm_sum_nb/num_NB 

        end do

        end subroutine         

    subroutine local_bond_order(total_atoms, l, Ylm_ave, local_ql) 
        implicit none 

        ! Passed 
        integer,intent(in) :: total_atoms,l  
        complex(dp), intent(in), dimension(-l:l,total_atoms) :: Ylm_ave         

        ! Local:
        integer :: i
        real(dp) ::  norm_const 
        complex(dp), dimension(-l:l) :: local_Ylm  

        ! Return:
        real(dp), intent(out), dimension(1:total_atoms) :: local_ql 
                
        norm_const = (4.0d0*pi)/(2*l+1) 
        
        do i = 1, total_atoms

            local_Ylm = Ylm_ave(:,i)
    
            local_ql(i) = zsqrt(norm_const*sum(local_Ylm*dconjg(local_Ylm))) 
            
        end do 

        end subroutine 

    subroutine global_bond_order(total_atoms, num_NB_list, l, Ylm_ave, global_op)
        implicit none 

        ! Passed 
        integer,intent(in) :: total_atoms,l  
        integer, intent(in), dimension(:) :: num_NB_list
        complex(dp), intent(in), dimension(-l:l,total_atoms) :: Ylm_ave         

        ! Local:
        integer :: i,  num_NB 
        real(dp) :: norm_const 
        complex(dp), dimension(-l:l) :: sum_Ylm

        ! Return:
        real(dp), intent(out) :: global_op

        sum_Ylm = 0.0d0

        do i = 1, total_atoms

            sum_Ylm = sum_Ylm + Ylm_ave(:,i)*num_NB_list(i)  

        end do
    
        norm_const = 4.0d0*pi/(2*l+1)

        sum_Ylm = sum_Ylm/sum(num_NB_list)

        global_op = zsqrt(norm_const*sum(sum_Ylm*dconjg(sum_Ylm)))

        end subroutine 

    subroutine neighbor_averaged_qlm(total_atoms,l,num_NB_list,NB_list,Ylm,Qlm)
        implicit none

        ! Passed:

        integer,intent(in) :: total_atoms,l
        integer,intent(in),dimension(:) :: num_NB_list
        integer,intent(in),dimension(:,:) :: NB_list
        complex(dp), intent(in), dimension(:,:) :: Ylm 

        ! Local:
        complex(dp),dimension(-l:l) :: sum_neighbor_qlm
        integer :: i,j,j_nb,num_NB

        ! Return
        complex(dp),intent(out),dimension(-l:l,1:total_atoms) :: Qlm

        do i = 1, total_atoms

            num_NB = num_NB_list(i)

            sum_neighbor_qlm = Ylm(:,i)

            do j = 1, num_NB

                ! NB_list(j,i) return the j neighbor (index) of i particle
                sum_neighbor_qlm = sum_neighbor_qlm + Ylm(:,NB_list(j,i))

            end do

            Qlm(:,i) = sum_neighbor_qlm/num_NB

        end do

        end subroutine 

    subroutine non_nnb_product(total_atoms, num_NB_list, maxnb, l, NB_list, qlm, cij)
        implicit none

        ! Passed
        integer, intent(in) :: total_atoms
        integer, intent(in), dimension(:) :: num_NB_list 
        integer, intent(in) :: l
        integer, intent(in) :: maxnb
        complex(dp),intent(in),dimension(:,:) :: qlm
        integer,intent(in),dimension(:,:) :: NB_list

        ! Local
        real(dp) :: qlm_norm_i, qlm_norm_j
        complex(dp),dimension(-l:l) :: Ylm_i, Ylm_j
        integer :: i,j

        ! Return
        real(dp),intent(out),dimension(1:maxnb,1:total_atoms) :: cij

        cij = 0.0d0

        do i = 1, total_atoms

            Ylm_i = qlm(:,i)

            qlm_norm_i = zsqrt(sum(Ylm_i*dconjg(Ylm_i)))
            
            do j = 1, num_NB_list(i)
    
                ! Find neighbor j, and extract qlm 
                Ylm_j = qlm(:,NB_list(j,i))

                qlm_norm_j = zsqrt(sum(Ylm_j*dconjg(Ylm_j)))
                
                cij(j,i) = dot_product((Ylm_i)/qlm_norm_i, (Ylm_j)/qlm_norm_j)  !SUM(Ylm_i/qlm_norm_i*DCONJG(Ylm_j/qlm_norm_j)) 
                
            end do

        end do

        end subroutine 
    

    subroutine nnb_dot_product(total_atoms, nnb, l, NB_list, qlm, cij)
        implicit none

        ! Passed
        integer, intent(in) :: total_atoms
        integer, intent(in) :: l
        integer, intent(in) :: nnb
        complex(dp),intent(in),dimension(:,:) :: qlm
        integer,intent(in),dimension(:,:) :: NB_list

        ! Local
        real(dp) :: qlm_norm_i, qlm_norm_j, sum_ylm
        complex(dp),dimension(-l:l) :: Ylm_i, Ylm_j
        integer :: i,j, m 

        ! Return
        real(dp),intent(out),dimension(1:nnb,1:total_atoms) :: cij

        cij = 0.0d0

        do i = 1, total_atoms

            Ylm_i = qlm(:,i)

            qlm_norm_i = zsqrt(sum(Ylm_i*dconjg(Ylm_i)))
            
            do j = 1,nnb
    
                ! Find neighbor j, and extract qlm 
                Ylm_j = qlm(:,NB_list(j,i))

                qlm_norm_j = zsqrt(sum(Ylm_j*dconjg(Ylm_j)))
    
                sum_ylm = 0.0d0

                do m  = -l, l
               
                    sum_ylm = sum_ylm + Ylm_i(m)*dconjg(Ylm_j(m))  

                end do 
                
                !cij(j,i) = dot_product((Ylm_i)/qlm_norm_i, (Ylm_j)/qlm_norm_j)  !SUM(Ylm_i/qlm_norm_i*DCONJG(Ylm_j/qlm_norm_j)) 
                cij(j,i) = sum_ylm/(qlm_norm_i*qlm_norm_j)  

            end do

        end do

        end subroutine 

    subroutine nnb_dot_product_crystal(total_atoms,& 
                                     & nnb,&
                                     & l,&
                                     & NB_list,&
                                     & conntect_cut,&
                                     & crys_bonds_cut,&
                                     & qlm, & 
                                     & cij, &
                                     & count_bonds, & 
                                     & its_crystal)
        implicit none 

        ! Passed
        integer, intent(in) :: total_atoms
        integer, intent(in) :: l
        integer, intent(in) :: nnb 
        real(dp), intent(in) :: conntect_cut
        integer, intent(in) :: crys_bonds_cut 
        complex(dp), intent(in),dimension(:,:) :: qlm
        integer, intent(in),dimension(:,:) :: NB_list

        ! Local
        real(dp) :: qlm_norm_i, qlm_norm_j, scalar_prod, sum_ylm 
        complex(dp),dimension(-l:l) :: Ylm_i, Ylm_j
        integer :: i,j, bonds_per_mol, m  

        ! Return
        real(dp),intent(out),dimension(1:nnb,1:total_atoms) :: cij
        real(dp), intent(out), dimension(1:total_atoms) :: count_bonds 
        logical, intent(out), dimension(1:total_atoms) :: its_crystal

        cij = 0.0d0

        count_bonds = 0.0d0 

        its_crystal = .False. 

        do i = 1, total_atoms

            Ylm_i = qlm(:,i)

            qlm_norm_i = zsqrt(sum(Ylm_i*dconjg(Ylm_i)))

            bonds_per_mol = 0
            
            do j = 1,nnb
    
                ! Find neighbor j, and extract qlm 
                Ylm_j = qlm(:,NB_list(j,i))

                qlm_norm_j = zsqrt(sum(Ylm_j*dconjg(Ylm_j)))

                ! compute the dot product of two complex vectors

                sum_ylm = 0.0d0

                do m = -l, l
    
                    sum_ylm = sum_ylm + Ylm_i(m)*dconjg(Ylm_j(m))

                end do
                
                !scalar_prod = dot_product((Ylm_i)/qlm_norm_i, (Ylm_j)/qlm_norm_j)  !SUM(Ylm_i/qlm_norm_i*DCONJG(Ylm_j/qlm_norm_j)) 
                scalar_prod = sum_ylm/(qlm_norm_i*qlm_norm_j) 

                cij(j,i) = scalar_prod

                if (scalar_prod > conntect_cut) then

                    bonds_per_mol = bonds_per_mol + 1 
    
                end if 

            end do

            ! if number of crystal bonds of a particle exceeds certain threshold
            if (bonds_per_mol >= crys_bonds_cut) then 
   
                its_crystal(i) = .True.  
                
            end if 

            count_bonds(i) = bonds_per_mol 

        end do
                
        end subroutine 

! -------------------------------------------------------------------------------------------------
!                                    Third-Order Invariant Wl 
! -------------------------------------------------------------------------------------------------

    pure subroutine initialize_wigner3j(l, num_pairs, m1_m2_m3_sum0_mat, wigner_3j_symobol_mat)
        implicit none 
        ! Passed
        integer, intent(in) :: l 
    
        ! Local 
        
        ! Retrun
        integer, intent(out) :: num_pairs  
        integer, intent(out), dimension(:,:), allocatable :: m1_m2_m3_sum0_mat
        real(dp), intent(out), dimension(:), allocatable :: wigner_3j_symobol_mat
         
        num_pairs = count_pairs_m1_m2_m3(l)

        allocate(m1_m2_m3_sum0_mat(1:3, 1:num_pairs))
        allocate(wigner_3j_symobol_mat(1:num_pairs))

        call wigner_3j_symbol_coefficients(l, num_pairs, m1_m2_m3_sum0_mat,wigner_3j_symobol_mat)         

        end subroutine 

    pure integer function count_pairs_m1_m2_m3(l)
        implicit none
        integer,intent(in) :: l
        integer :: m1,m2,m3,counter

        count_pairs_m1_m2_m3 = 0

        do m1 = -l,l

            do m2 = -l,l

                do m3 = -l,l

                    if ((m1+m2+m3) == 0) then

                        count_pairs_m1_m2_m3 = count_pairs_m1_m2_m3 + 1

                    end if

                end do

            end do

        end do

        end function

    pure subroutine wigner_3j_symbol_coefficients(l, num_pairs, m1_m2_m3_sum0_mat,wigner_3j_symobol_mat)
        implicit none

        ! Passed 
        integer,intent(in) :: l
        integer, intent(in) :: num_pairs

        ! Local
        integer :: m1, m2, m3, counter

        ! Return
        integer, intent(inout), dimension(1:3,1:num_pairs) :: m1_m2_m3_sum0_mat
        real(dp),intent(inout), dimension(1:num_pairs) :: wigner_3j_symobol_mat

        counter = 0

        do m1 = -l,l

            do m2 = -l,l

                do m3 = -l,l

                    ! only m1 + m2 + m3 = 0 will be summed together 

                    if (( m1+m2+m3 ) == 0)  Then

                        counter = counter + 1

                        wigner_3j_symobol_mat(counter) = wigner_3j(l,l,l,m1,m2,m3)

                        m1_m2_m3_sum0_mat(:,counter) = [m1,m2,m3]

                    end if

                end do

            end do

        end do

        end subroutine

    subroutine calc_Wl(num_pairs,total_atoms, l,m1_m2_m3_sum0_mat,wigner_3j_symobol_mat,Qlm,calc_W)
        implicit none
        ! Passed
        integer,intent(in) :: total_atoms,num_pairs,l
        integer,intent(in),dimension(:,:) :: m1_m2_m3_sum0_mat
        real(dp),intent(in),dimension(:) :: wigner_3j_symobol_mat
        complex(dp),intent(in),dimension(-l:l,1:total_atoms) :: Qlm

        ! Local
        integer :: i,j
        real(dp) :: calc_w_inter 
    
        ! Return
        real(dp),intent(inout),dimension(:) :: calc_W

        calc_w = 0.0d0

        do i = 1,total_atoms

            calc_w_inter = 0.0d0

            do j = 1, num_pairs
            
                
                calc_w_inter = calc_w_inter + product(Qlm(m1_m2_m3_sum0_mat(:,j),i)) * wigner_3j_symobol_mat(j)

            end do
    
            calc_w(i) = calc_w_inter/((sum(Qlm(:, i)*dconjg(Qlm(:, i))))**(1.5d0))

        end do

        end subroutine

    subroutine select_crystalline_particles(total_atoms,its_crystal,num_crystal, crystal_id)
        implicit none 
        ! Passed
        integer, intent(in) :: total_atoms 
        logical, intent(in), dimension(1:total_atoms) :: its_crystal 

        ! Local
        integer :: i, counter 
        
        ! Return
        integer, intent(out) :: num_crystal 
        integer, intent(out), dimension(:), allocatable :: crystal_id 

        counter = 0 

        num_crystal = count(its_crystal) 

        allocate(crystal_id(1:num_crystal))

        do i = 1, total_atoms

            counter = counter + 1  

            if (its_crystal(i)) then  
       
                crystal_id(counter) = i  

            end if 
    
        end do 
        
        end subroutine 

    subroutine convert_c_string_f_string(str,strlen,f_string)
        implicit none
        integer ,intent(in):: strlen
        character(len=1),intent(in),dimension(1:strlen) :: str
        character(len=:),allocatable,intent(out) :: f_string
        integer :: i

        f_string = ""

        do i = 1,strlen

            f_string = f_string//str(i)

        end do

        end subroutine

end module  
