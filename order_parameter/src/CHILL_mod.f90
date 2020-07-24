module CHILL
    use system, only: c_int, c_double, dp, sp
    use order_parameter, only: build_homo_neighbor_list,&
                             & check_nb_list, & 
                             & apply_nearest_neighbor_crit, &  
                             & initialize_Ylm, & 
                             & compute_sum_Ylm, & 
                             & neighbor_averaged_qlm, & 
                             & compute_dot_product_nnb, & 
                             & convert_c_string_f_string
    implicit none 
    private 


contains 

    subroutine initialize_CHILL(l, sph_const, Plm_const) bind(c, name="initialize_CHILL")
        implicit none
        ! Passed
        integer(c_int), intent(in) :: l
        
        ! Return
        real(c_double), intent(out), dimension(-l:l) :: sph_const
        real(c_double), intent(out), dimension(0:l,4) :: Plm_const

        call initialize_Ylm(l, sph_const, Plm_const)

        end subroutine 

    subroutine call_CHILL(sph_const,& 
                        & Plm_const, &
                        & CHILL_keyword, &
                        & strlength, &
                        & total_atoms, &
                        & maxnb,& 
                        & nnb, &
                        & l, &
                        & cutoff_sqr,&
                        & box, &
                        & xyz, &
                        & cij, &
                        & chill_id_list) bind(c, name="call_CHILL")
        implicit none 

        ! Passed:
        integer(c_int), intent(in) :: strlength
        character(len=1), intent(in), dimension(1:strlength) :: CHILL_keyword
        real(c_double), intent(in), dimension(-l:l) :: sph_const 
        real(c_double), intent(in), dimension(0:l,4) :: Plm_const 
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
        character(len=:), allocatable :: CHILL_key 

        ! Return:
        real(c_double), intent(out), dimension(1:nnb, 1:total_atoms) :: cij
        integer(c_int), intent(out),dimension(1:total_atoms) :: chill_id_list 

        call convert_c_string_f_string(CHILL_keyword, strlength, CHILL_key)

        ! get the neighbor list
        
        call build_homo_neighbor_list(total_atoms, maxnb, cutoff_sqr, xyz, box, num_NB_list, NB_list, Rij) 
        
        ! check the neighbor list: 
        call check_nb_list(num_NB_list,maxnb,nnb) 
        
        ! apply the nearest neighbor:
        call apply_nearest_neighbor_crit(nnb, total_atoms, num_NB_list, NB_list, Rij)
       
        ! compute_sum_Ylm  
        call compute_sum_Ylm(sph_const, Plm_const, total_atoms, l, num_NB_list, Rij, Ylm)
        
        ! compute the dot_product:
        call compute_dot_product_nnb(total_atoms, nnb, l, NB_list, Ylm, cij)
        
        ! apply the CHILL:
        call apply_CHILL_or_CHILL_plus(CHILL_key, total_atoms, nnb, cij, NB_list, chill_id_list)
        
        end subroutine

    subroutine apply_CHILL_or_CHILL_plus(CHILL_keyword, total_atoms, nnb, cij, NB_list, chill_id_list)
        implicit none
        ! Passed 
        character(len=*), intent(in) :: CHILL_keyword
        integer, intent(in) :: total_atoms,nnb
        real(dp), intent(in),dimension(:,:)  :: cij
        integer, intent(in),dimension(:,:) :: NB_list
        
        ! Local
        integer, dimension(1:2,1:total_atoms) :: cij_staggered_eclipsed 
        real(dp), dimension(1:2) :: eclipsed_bonds

        ! Return
        integer, intent(out),dimension(1:total_atoms) :: chill_id_list
         
        if (CHILL_keyword == "CHILL") then

            eclipsed_bonds = [-0.2, -0.05]

            call staggered_or_eclipsed(total_atoms, nnb, cij, eclipsed_bonds, cij_staggered_eclipsed)
            
            call apply_CHILL(total_atoms, cij_staggered_eclipsed, nnb, NB_list, chill_id_list)
            
        else if (CHILL_keyword == "CHILL+") then 

            eclipsed_bonds = [-0.35, 0.25]

            call staggered_or_eclipsed(total_atoms, nnb, cij, eclipsed_bonds, cij_staggered_eclipsed)

            call apply_CHILL_plus(total_atoms, cij_staggered_eclipsed, nnb, NB_list, chill_id_list)
       
        else
        
            stop "CHILL keyword not recgonized: Please choose 'CHILL' or 'CHILL+' "

        end if

        end subroutine

    pure subroutine staggered_or_eclipsed( total_atoms,nnb,cij, eclispsed_bounds, cij_staggered_eclipsed)
        implicit none
        ! Passed
        real(dp),intent(in), dimension(:) :: eclispsed_bounds 
        integer, intent(in) :: total_atoms,nnb
        real(dp), intent(in), dimension(:,:)  :: cij

        ! Local
        integer :: i,j,stagger,eclipse
        real(dp) :: bond, lower_bound, upper_bound 
    
        ! Return
        integer, intent(inout), dimension(:,:) :: cij_staggered_eclipsed

        cij_staggered_eclipsed = 0

        lower_bound = eclispsed_bounds(1) 

        upper_bound = eclispsed_bounds(2)

        do i = 1, total_atoms

            stagger = 0; eclipse = 0

            do j = 1, nnb

                bond = cij(j,i)

                if ( bond < -0.8d0 ) THEN

                    stagger = stagger + 1

                else if ( bond < upper_bound .and. bond > lower_bound ) then

                    eclipse = eclipse + 1

                end if

            end do

            cij_staggered_eclipsed(1,i) = stagger

            cij_staggered_eclipsed(2,i) = eclipse

        end do

        end subroutine

    pure subroutine apply_CHILL(total_atoms,cij_staggered_eclipsed,nnb,NB_list,chill_id_list)
        implicit none
        integer,intent(in) :: total_atoms,nnb
        integer,intent(in),dimension(:,:) :: cij_staggered_eclipsed
        integer,intent(in),dimension(:,:) :: NB_LIST
        integer :: i,j,nb,num_Ic,num_ifC,num_Ih,num_LDA,counter,num_liquid,num_ice
        integer,intent(out),dimension(1:total_atoms) :: chill_id_list

        counter = 0

        chill_id_list = 0

        ! In cij_staggered_eclipsed, 1st row is staggered, 2nd row is eclipsed
        do i = 1, total_atoms

            ! Define Ice cubic:
            if (cij_staggered_eclipsed(1,i) == 4 &
               .and. cij_staggered_eclipsed(2,i) == 0) then

                chill_id_list(i) = 1

            end if

            ! Define the Ice Hexagonal:
            if (cij_staggered_eclipsed(1,i) == 3 &
               & .and. cij_staggered_eclipsed(2,i) ==1) then

                chill_id_list(i) = 2

            end if
    
            ! Define the Interfacial Ice:
            ! For a particle with exactly 2 staggered bond 
            if (cij_staggered_eclipsed(1,i) == 2) then

                do j = 1,nnb

                    nb = NB_list(j,i)

                    ! At least one neighbor with more than two staggered bonds
                    if ( cij_staggered_eclipsed(1,nb) > 2) then

                        chill_id_list(i) = 3

                        exit

                    end if

                end do

            ! For a particle with 3 staggered bond and 0 eclipsed bond 
            else if ( cij_staggered_eclipsed(1,i) == 3  &
                    & .and. cij_staggered_eclipsed(2,i) ==0 ) then

                do j = 1,nnb

                    nb = NB_list(j,i)

                    ! at least one neighbor with two staggered bonds
                    if (cij_staggered_eclipsed(1,nb) == 2) then

                        chill_id_list(i) = 3

                        exit

                    end if

                end do

            end if

            ! Define the liquid or amorphous phase
            if (chill_id_list(i) == 0) then

                chill_id_list(i) = 4

            end if

        end do

    end subroutine

    pure subroutine apply_CHILL_plus(total_atoms,cij_staggered_eclipsed,nnb,NB_list,chill_id_list)
        implicit none
        integer,intent(in) :: total_atoms,nnb
        integer,intent(in),dimension(:,:) :: cij_staggered_eclipsed
        integer,intent(in),dimension(:,:) :: NB_LIST
        integer :: i,j,nb,num_Ic,num_ifC,num_Ih,num_LDA,counter,num_liquid,num_ice
        integer,intent(out),dimension(1:total_atoms) :: chill_id_list

        counter = 0

        chill_id_list = 0

        ! In cij_staggered_eclipsed, 1st row is staggered, 2nd row is eclipsed
        do i = 1, total_atoms

            ! Define Ice cubic:
            if (cij_staggered_eclipsed(1,i) == 4 &
               .and. cij_staggered_eclipsed(2,i) == 0) then

                chill_id_list(i) = 1

            end if

            ! Define the Ice Hexagonal 
            if (cij_staggered_eclipsed(1,i) == 3 &
               & .and. cij_staggered_eclipsed(2,i) ==1) then

                chill_id_list(i) = 2

            end if

            ! Define the Interfacial Ice
    
            ! For a particle with exactly 2 staggered bond 
            if (cij_staggered_eclipsed(1,i) == 2) then

                do j = 1,nnb

                    nb = NB_list(j,i)

                    ! At least one neighbor with more than two staggered bonds 
                    if ( cij_staggered_eclipsed(1,nb) > 2) then

                        chill_id_list(i) = 3

                        exit

                    end if

                end do

            ! For a particle with 3 staggered bond and 0 eclipised bond
            else if ( cij_staggered_eclipsed(1,i) == 3  &
                    & .and. cij_staggered_eclipsed(2,i) ==0 ) then

                do j = 1,nnb

                    nb = NB_list(j,i)

                    ! At least one neighbor with more than one staggered bond
                    if (cij_staggered_eclipsed(1,nb) > 1) then

                        chill_id_list(i) = 3

                        exit

                    end if

                end do

            end if

            ! Define the clathhydrate 
            ! All 4 bonds are eclipsed 
            if (cij_staggered_eclipsed(2,i) == 4) then

                chill_id_list(i) = 5
                
            end if 

            ! Define the interfacial clathhydrate 

            ! 3 eclipsed bond and any staggered bond
            if (cij_staggered_eclipsed(2,i) == 3) then

                chill_id_list(i) = 6 

            end if 

            ! Define the liquid
            if (chill_id_list(i) == 0) then

                chill_id_list(i) = 4

            end if

        end do

    end subroutine

end module 
