module nb_list
    use system, only: c_int, c_double, dp, sp, error_unit 
    use constants, only: pi
    use sorting, only: quick_sort 

    implicit none
    private
    public :: build_homo_neighbor_list, & 
            & rebuild_neighbor_list_nnb, & 
            & check_nb_list, &
            & apply_nearest_neighbor_crit

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
        integer, intent(inout), dimension(:) :: num_NB_list
        integer, intent(inout), dimension(:, :) :: NB_list

        do iatom = 1, total_atoms

            count_nb = 0

            temp_nb_lst = 0 

            do jnb = 1, num_NB_list(iatom)

                dist = Rij(4,jnb, iatom)

                if (dist <= cutoff) then

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

                ! No need to use nearest neighbors scheme and return
                continue

            case(1:) 

                call check_enough_neighbors(nnb, num_NB_list) 
                
                ! apply nearest neighbors scheme    
                call use_nearest_neighbors(nnb,total_atoms, num_NB_list, NB_list, Rij) 

            case(:-1) 

                stop "number of nearest neighbors can not be negative "
    
        end select  

        end subroutine 


end module 
