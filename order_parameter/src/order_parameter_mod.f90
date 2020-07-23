module order_parameter
    use system, only: c_int, c_double, dp, sp 
    use spherical_harmonics, only: compute_spherical_harmonics, &
                                 & optimized_Q12   
    implicit none
    private
    public :: compute_sum_Ylm,& 
            & build_homo_neighbor_list, &
            & neighbor_averaged_qlm, &
            & compute_dot_product 
              
contains

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

            Ylm_sum(:,iatom) = Ylm_sum_nb 

        end do 

        end subroutine

    subroutine neighbor_averaged_qlm(total_atoms,l,num_NB_list,NB_list,Ylm,qlm_nb_ave)
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
        complex(dp),intent(out),dimension(:,:) :: qlm_nb_ave

        do i = 1, total_atoms

            num_NB = num_NB_list(i)

            sum_neighbor_qlm = Ylm(:,i)

            do j = 1, num_NB 

                ! NB_list(j,i) return the j neighbor (index) of i particle
                sum_neighbor_qlm = sum_neighbor_qlm + Ylm(:,NB_list(j,i)) 

            end do 

            qlm_nb_ave(:,i) = sum_neighbor_qlm/(num_NB)

        end do
        
        end subroutine 

    subroutine compute_dot_product(total_atoms, maxnb, l, num_NB_list, NB_list, qlm, cij)
        implicit none

        ! Passed
        integer, intent(in) :: total_atoms
        integer, intent(in) :: l
        integer, intent(in) :: maxnb 
        complex(dp),intent(in),dimension(:,:) :: qlm
        integer,intent(in),dimension(:) :: num_NB_list
        integer,intent(in),dimension(:,:) :: NB_list

        ! Local
        real(dp) :: qlm_norm_i, qlm_norm_j
        complex(dp),dimension(-l:l) :: Ylm_i, Ylm_j
        integer :: i,j,num_NB,nb,num_comp

        ! Return
        real(dp),intent(out),dimension(1:maxnb,1:total_atoms) :: cij

        cij = 0.0d0

        do i = 1,total_atoms

            Ylm_i = qlm(:,i)

            num_NB = num_NB_list(i)

            qlm_norm_i = ZSQRT(SUM(Ylm_i*DCONJG(Ylm_i)))

            do j = 1,num_NB

                nb = NB_list(j,i)

                Ylm_j = qlm(:,nb)

                qlm_norm_j =ZSQRT(SUM(Ylm_j*DCONJG(Ylm_j)))

                cij(j,i) = SUM(Ylm_i/qlm_norm_i*DCONJG(Ylm_j)/qlm_norm_j)

            end do

        end do

        end subroutine 


end module  
