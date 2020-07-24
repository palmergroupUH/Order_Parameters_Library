module order_parameter
    use sorting, only: selection_sort 
    use system, only: c_int, c_double, dp, sp, error_unit 
    use constants, only: pi
    use spherical_harmonics, only: compute_spherical_harmonics, &
                                 & compute_spherical_harmonics_const,&
                                 & compute_associated_legendre_const,& 
                                 & optimized_Q12
    implicit none
    private
    public :: initialize_Ylm, & 
            & compute_sum_Ylm,& 
            & build_homo_neighbor_list, &
            & apply_nearest_neighbor_crit, & 
            & check_nb_list, &             
            & neighbor_averaged_qlm, &
            & compute_dot_product_nnb, & 
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
            
            call selection_sort(num_NB,Rij(4,1:num_NB,iatom),tracking_index, nnb)
            
            NB_list(1:nnb,iatom) = NB_list(tracking_index(1:nnb),iatom)

            Rij(1:3,1:nnb,iatom) = Rij(1:3,tracking_index(1:nnb),iatom)
    
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

    subroutine compute_sum_Ylm(sph_const,&
                              & Plm_const,&
                              & total_atoms,& 
                              & l, &
                              & num_NB_list, &
                              & Rij,& 
                              & Ylm_sum)
        implicit none 
        ! Passed:
        integer, intent(in) :: l, total_atoms 
        real(dp), intent(in), dimension(-l:l) :: sph_const  
        real(dp),intent(in),dimension(0:l,4) :: Plm_const 
        integer,intent(in),dimension(:) :: num_NB_list 
        real(dp),intent(in),dimension(:,:,:) :: Rij
        
        ! Local:
        integer :: iatom, jnb
        complex(dp), dimension(-l:l) :: Ylm, Ylm_sum_nb
        real(dp) :: rcx,rcy,rcz,rl
        real(dp) :: sin_theta, cos_theta, sin_phi, cos_phi
      
        ! Return:
        complex(dp), intent(out), dimension(-l:l,total_atoms) :: Ylm_sum

        do iatom = 1, total_atoms
    
            Ylm_sum_nb = 0.0d0

            do jnb = 1, num_NB_list(iatom)

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

                if (l /= 12) then
                
                    call compute_spherical_harmonics(sph_const, Plm_const, l, cos_theta, cos_phi, sin_phi, Ylm)
    
                else if (l==12) then

                    call optimized_Q12(sin_theta, cos_theta, cos_phi, sin_phi, Ylm)

                end if
    
                Ylm_sum_nb = Ylm_sum_nb + Ylm

            end do

            Ylm_sum(:,iatom) = Ylm_sum_nb/num_NB_list(iatom)  

        end do 

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
        complex(dp),intent(out),dimension(:,:) :: Qlm

        do i = 1, total_atoms

            num_NB = num_NB_list(i)

            sum_neighbor_qlm = Ylm(:,i)

            do j = 1, num_NB

                ! NB_list(j,i) return the j neighbor (index) of i particle
                sum_neighbor_qlm = sum_neighbor_qlm + Ylm(:,NB_list(j,i))

            end do

            Qlm(:,i) = sum_neighbor_qlm/(num_NB)

        end do

        end subroutine 

    subroutine compute_dot_product_nnb(total_atoms, nnb, l, NB_list, qlm, cij)
        implicit none

        ! Passed
        integer, intent(in) :: total_atoms
        integer, intent(in) :: l
        integer, intent(in) :: nnb 
        complex(dp),intent(in),dimension(:,:) :: qlm
        integer,intent(in),dimension(:,:) :: NB_list

        ! Local
        real(dp) :: qlm_norm_i, qlm_norm_j
        complex(dp),dimension(-l:l) :: Ylm_i, Ylm_j
        integer :: i,j,num_comp

        ! Return
        real(dp),intent(out),dimension(1:nnb,1:total_atoms) :: cij

        cij = 0.0d0

        do i = 1,total_atoms

            Ylm_i = qlm(:,i)

            qlm_norm_i = zsqrt(sum(Ylm_i*dconjg(Ylm_i)))
            
            do j = 1,nnb
    
                ! Find neighbor j, and extract qlm 
                Ylm_j = qlm(:,NB_list(j,i))

                qlm_norm_j = zsqrt(sum(Ylm_j*dconjg(Ylm_j)))
                
                cij(j,i) = dot_product((Ylm_i)/qlm_norm_i, (Ylm_j)/qlm_norm_j)  !SUM(Ylm_i/qlm_norm_i*DCONJG(Ylm_j/qlm_norm_j)) 
                
            end do

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
