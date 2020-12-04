module sph_op
    use system, only: c_int, c_double, dp, sp, error_unit 
    use constants, only: pi
    use spherical_harmonics, only: compute_spherical_harmonics, &
                                 & compute_spherical_harmonics_const,&
                                 & compute_associated_legendre_const,& 
                                 & optimized_Y12, & 
                                 & optimized_Y8, & 
                                 & optimized_Y6, & 
                                 & optimized_Y4

    use wigner3j_symbol, only: wigner_3j  

    implicit none
    ! Default: all subroutines and functions should be visible only 
    private

    ! Export following subroutines/functions
    ! "public" make them accessible to other Fortran programs through "use" 
    ! otherwise other subroutines and variables can not be used by other programs

    public :: compute_ave_qlm, & 
            & local_bond_order, & 
            & global_bond_order, & 
            & coarse_grained_qlm, &
            & compute_nnb_dij, &
            & compute_dij, &
            & compute_nnb_dij_bond, &
            & compute_dij_bond, & 
            & lablel_crystal_like, & 
            & calc_Wl 

contains 

! -------------------------------------------------------------------------------------------------
!                                  Modifiable
! -------------------------------------------------------------------------------------------------

    logical function its_crys_bonds(keyword, dij)
        ! crystal bonds criterion 
        ! modify if needed
        implicit none
        character(len=*), intent(in) :: keyword
        real(dp), intent(in) :: dij 

        its_crys_bonds = .False.

        if (keyword == "RussoTanaka") then 

            if (dij > 0.75) then 

                its_crys_bonds = .True. 
    
            end if 

        else if (keyword == "ReinhardtDoye") then
    
            if (dij < -0.82d0 .or. ((dij > -0.145d0) .and. (dij < -0.065d0))) then  

                its_crys_bonds = .True.
       
            end if 

        end if

        end function

! -------------------------------------------------------------------------------------------------
!                                  Spherical Harmonics
! -------------------------------------------------------------------------------------------------


    subroutine initialize_Ylm(l, sph_const, Plm_const)
        implicit none 
        integer, intent(in) :: l
        real(dp), intent(out), dimension(-l:l) :: sph_const
        real(dp), intent(out), dimension(0:l,4) :: Plm_const

        call compute_associated_legendre_const(l, Plm_const) 

        call compute_spherical_harmonics_const(l, sph_const)

        end subroutine 

    ! General spherical harmonics for any order 
    ! Optimized Y12, Y8, Y6 and Y4 are automatically used
    ! at l = 12, l=8, l=6 or l = 4
    subroutine compute_ave_qlm(total_atoms,&
                             & l, &
                             & num_NB_list, &
                             & Rij,& 
                             & mode, &
                             & Ylm_ave)
        implicit none 
        ! Passed:
        integer, intent(in) :: l, total_atoms 
        integer,intent(in),dimension(:) :: num_NB_list 
        real(dp),intent(in),dimension(:,:,:) :: Rij
        character(len=*), intent(in) :: mode 
        
        ! Local:
        real(dp), dimension(-l:l) :: sph_const  
        real(dp), dimension(0:l,4) :: Plm_const 
        integer :: iatom, jnb, num_NB
        complex(dp), dimension(-l:l) :: Ylm, Ylm_sum_nb
        real(dp) :: rcx,rcy,rcz,rl
        real(dp) :: sin_theta, cos_theta, sin_phi, cos_phi
      
        ! Return:
        complex(dp), intent(out), dimension(-l:l,total_atoms) :: Ylm_ave
        
        call initialize_Ylm(l, sph_const, Plm_const)

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
                
                if (mode=="optimized") then  

                    if (l == 12) then

                        call optimized_Y12(sin_theta, cos_theta, cos_phi, sin_phi, Ylm)

                    else if (l == 8) then

                        call optimized_Y8(sin_theta, cos_theta, cos_phi, sin_phi, Ylm)

                    else if (l == 6) then 

                        call optimized_Y6(sin_theta, cos_theta, cos_phi, sin_phi, Ylm)

                    else if (l == 4) then 

                        call optimized_Y4(sin_theta, cos_theta, cos_phi, sin_phi, Ylm)
    
                    else

                        call compute_spherical_harmonics(sph_const, Plm_const, l, cos_theta, cos_phi, sin_phi, Ylm)
                        
                    end if

                else if (mode =="general") then 

                    call compute_spherical_harmonics(sph_const, Plm_const, l, cos_theta, cos_phi, sin_phi, Ylm)
                    
                else

                    stop "Please check the spelling and choose the keyword: 'optimized' & 
                      &  or 'general' for the mode variable in subroutine 'compute_ave_qlm'"
    
                end if
                
                Ylm_sum_nb = Ylm_sum_nb + Ylm

            end do

            Ylm_ave(:,iatom) = Ylm_sum_nb / num_NB

        end do
        
        
        end subroutine

    subroutine norm_sph_cmplx_v(l, Ylm, norm_2)
        implicit none 
        ! Passed
        integer, intent(in) :: l
        complex(dp), intent(in), dimension(-l:l) :: Ylm

        ! Local
        integer :: m 

        ! Return
        real(dp), intent(out) :: norm_2

        norm_2 = 0.0d0 

        do m = -l, l 
           
            norm_2 = norm_2 + Ylm(m) * dconjg(Ylm(m)) 
    
        end do 

        norm_2 = dsqrt(norm_2) 

        end subroutine

    subroutine local_bond_order(total_atoms, l, Ylm_ave, local_ql)
        implicit none

        ! Passed
        integer,intent(in) :: total_atoms, l
        complex(dp), intent(in), dimension(-l:l,total_atoms) :: Ylm_ave

        ! Local:
        integer :: m, i
        real(dp) ::  norm_const, sum_y
        complex(dp), dimension(-l:l) :: local_Ylm

        ! Return:
        real(dp), intent(out), dimension(1:total_atoms) :: local_ql 
                
        norm_const = (4.0d0*pi)/(2*l+1)
        
        do i = 1, total_atoms

            local_Ylm = Ylm_ave(:,i)
        
            sum_y = 0.0d0

            do m = -l, l 

                sum_y = sum_y + local_Ylm(m) * dconjg(local_Ylm(m)) 
    
            end do 
    
            local_ql(i) = norm_const * sum_y  !zsqrt(norm_const*sum(local_Ylm*dconjg(local_Ylm))) 
            
        end do 

        local_ql = dsqrt(local_ql) 
    
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


    subroutine coarse_grained_qlm(total_atoms, l, num_NB_list, NB_list, Ylm, Qlm)

        !Lechner, W., & Dellago, C. (2008). Accurate determination of crystal
        !structures based on averaged local bond order parameters. The Journal of
        !Chemical Physics, 129(11), 114707. https://doi.org/10.1063/1.2977970 
    
        implicit none

        ! Passed:

        integer,intent(in) :: total_atoms,l
        integer,intent(in),dimension(:) :: num_NB_list
        integer,intent(in),dimension(:,:) :: NB_list
        complex(dp), intent(in), dimension(:,:) :: Ylm 

        ! Local:
        complex(dp),dimension(-l:l) :: sum_neighbor_qlm
        integer :: i,j,j_nb,num_NB, j_atom 

        ! Return
        complex(dp),intent(out),dimension(-l:l,1:total_atoms) :: Qlm

        Qlm = 0.0d0 
        
        do i = 1, total_atoms

            num_NB = num_NB_list(i)

            sum_neighbor_qlm = Ylm(:,i)

            do j = 1, num_NB

                ! NB_list(j,i) return the j neighbor (index) of i particle
                j_atom = NB_list(j,i)

                sum_neighbor_qlm = sum_neighbor_qlm + Ylm(:,j_atom)
                
            end do
             
            Qlm(:,i) = sum_neighbor_qlm/num_NB

        end do

        end subroutine 


    subroutine compute_dij(total_atoms, num_NB_list, maxnb, l, NB_list, qlm, dij)
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
        real(dp),intent(out),dimension(1:maxnb,1:total_atoms) :: dij

        dij = 0.0d0

        do i = 1, total_atoms

            Ylm_i = qlm(:,i)

            qlm_norm_i = zsqrt(sum(Ylm_i*dconjg(Ylm_i)))
            
            do j = 1, num_NB_list(i)
    
                ! Find neighbor j, and extract qlm 
                Ylm_j = qlm(:,NB_list(j,i))

                qlm_norm_j = zsqrt(sum(Ylm_j*dconjg(Ylm_j)))
                
                dij(j,i) = dot_product((Ylm_i)/qlm_norm_i, (Ylm_j)/qlm_norm_j)  !SUM(Ylm_i/qlm_norm_i*DCONJG(Ylm_j/qlm_norm_j)) 
                
            end do

        end do

        end subroutine 
    

    subroutine compute_nnb_dij(total_atoms, nnb, l, NB_list, qlm, dij)
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
        real(dp),intent(out),dimension(1:nnb,1:total_atoms) :: dij

        dij = 0.0d0

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
                dij(j,i) = sum_ylm/(qlm_norm_i*qlm_norm_j)  

            end do

        end do

        end subroutine 

    subroutine compute_dij_bond(style, &
                              & total_atoms, &
                              & maxnb, &
                              & num_NB_list, &
                              & l, &
                              & NB_list, &
                              & n_bonds, &
                              & qlm, & 
                              & nxtl, & 
                              & dij, &
                              & count_bonds)
        implicit none 

        ! Passed
        character(len=*), intent(in) :: style
        integer, intent(in) :: total_atoms
        integer, intent(in) :: maxnb
        integer, intent(in) :: l
        integer, intent(in) :: n_bonds
        integer, intent(in), dimension(:) :: num_NB_list 
        complex(dp), intent(in), dimension(:,:) :: qlm
        integer, intent(in),dimension(:,:) :: NB_list

        ! Local
        real(dp) :: qlm_norm_i, qlm_norm_j, scalar_prod, sum_ylm 
        complex(dp),dimension(-l:l) :: Ylm_i, Ylm_j
        integer :: i,j, count_crys_bonds, m  

        ! Return
        integer, intent(out) :: nxtl
        real(dp), intent(out), dimension(1:maxnb,1:total_atoms) :: dij
        integer, intent(out), dimension(1:total_atoms) :: count_bonds 

        dij = 0.0d0

        count_bonds = 0.0d0

        nxtl = 0

        do i = 1, total_atoms

            Ylm_i = qlm(:,i)

            qlm_norm_i = zsqrt(sum(Ylm_i*dconjg(Ylm_i)))
            
            count_crys_bonds = 0
            
            do j = 1, num_NB_list(i)
    
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
               
                dij(j,i) = scalar_prod
                
                if (its_crys_bonds(style, scalar_prod)) then

                    count_crys_bonds = count_crys_bonds + 1 
    
                end if 

            end do

            count_bonds(i) = count_crys_bonds 
            
            if (count_bonds(i) >= n_bonds) then 
    
                nxtl = nxtl + 1 
    
            end if 

        end do
        
        end subroutine 

    subroutine compute_nnb_dij_bond(style, &
                                  & total_atoms, &
                                  & nnb, &
                                  & l, &
                                  & NB_list, &
                                  & n_bonds, &
                                  & qlm, & 
                                  & nxtl, & 
                                  & dij, &
                                  & count_bonds)
        implicit none 

        ! Passed
        character(len=*), intent(in) :: style
        integer, intent(in) :: total_atoms
        integer, intent(in) :: l
        integer, intent(in) :: nnb 
        integer, intent(in) :: n_bonds
        complex(dp), intent(in),dimension(:,:) :: qlm
        integer, intent(in),dimension(:,:) :: NB_list

        ! Local
        real(dp) :: qlm_norm_i, qlm_norm_j, scalar_prod, sum_ylm 
        complex(dp),dimension(-l:l) :: Ylm_i, Ylm_j
        integer :: i,j, count_crys_bonds, m  

        ! Return
        integer, intent(out) :: nxtl
        real(dp), intent(out), dimension(1:nnb,1:total_atoms) :: dij
        integer, intent(out), dimension(1:total_atoms) :: count_bonds 

        dij = 0.0d0

        count_bonds = 0.0d0

        nxtl = 0
    
        do i = 1, total_atoms

            Ylm_i = qlm(:,i)
            
            !qlm_norm_i =  zsqrt(sum(Ylm_i * dconjg(Ylm_i)))

            call norm_sph_cmplx_v(l, Ylm_i, qlm_norm_i)

            count_crys_bonds = 0

            do j = 1,nnb
    
                ! Find neighbor j, and extract qlm 
                Ylm_j = qlm(:,NB_list(j,i))

                !qlm_norm_j = zsqrt(sum(Ylm_j*dconjg(Ylm_j)))
                call norm_sph_cmplx_v(l, Ylm_j, qlm_norm_j) 

                ! Three different ways to compute the dot product 
                ! method 1: scalar_prod = SUM(Ylm_i/qlm_norm_i*DCONJG(Ylm_j/qlm_norm_j))  
                ! method 2: scalar_prod = dot_product((Ylm_i)/qlm_norm_i, (Ylm_j)/qlm_norm_j) 
                
                ! method 3: compute the dot product of two complex vectors
                sum_ylm = 0.0d0

                do m = -l, l
    
                    sum_ylm = sum_ylm + Ylm_i(m)*dconjg(Ylm_j(m))

                end do

                scalar_prod = sum_ylm/(qlm_norm_i * qlm_norm_j) 
                
                dij(j,i) = scalar_prod
                
                if (its_crys_bonds(style, scalar_prod)) then

                    count_crys_bonds = count_crys_bonds + 1 
    
                end if 

            end do

            count_bonds(i) = count_crys_bonds 
        
            if (count_bonds(i) >= n_bonds) then 
    
                nxtl = nxtl + 1 
    
            end if 


        end do
                
        end subroutine 

    subroutine lablel_crystal_like(total_atoms, nxtl, n_bonds, n_dij_bonds, crystal_id)
        implicit none 

        ! Passed
        integer, intent(in) :: total_atoms
        integer, intent(in) :: nxtl
        integer, intent(in) :: n_bonds
        integer, intent(in), dimension(:) :: n_dij_bonds 

        ! Local
        integer :: counter, i

        ! Return
        integer, intent(out), dimension(:), allocatable  :: crystal_id
    
        counter = 0

        allocate(crystal_id(1:nxtl))

        crystal_id = 0 

        do i = 1, total_atoms 

            if (n_dij_bonds(i) >= n_bonds) then

                counter = counter + 1 

                crystal_id(counter) = i

            end if 

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

        do i = 1, total_atoms

            calc_w_inter = 0.0d0

            do j = 1, num_pairs
            
                
                calc_w_inter = calc_w_inter + product(Qlm(m1_m2_m3_sum0_mat(:,j),i)) * wigner_3j_symobol_mat(j)

            end do
    
            calc_w(i) = calc_w_inter/((sum(Qlm(:, i)*dconjg(Qlm(:, i))))**(1.5d0))

        end do

        end subroutine

end module 
