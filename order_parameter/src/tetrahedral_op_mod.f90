module tetrahedral_order_parameter
    use sorting, only: selection_sort 
    use system, only: c_int, c_double, dp, sp, error_unit 

    implicit none

    private
    public :: calc_tetrahedral_op

contains 

    subroutine calc_tetrahedral_op(maxnb, nnb, total_atoms, Rij, q_tetrahedral)
        implicit none

        ! Passed
        integer(c_int), intent(in) :: nnb,total_atoms,maxnb
        real(c_double), intent(in), dimension(1:4,1:maxnb,1:total_atoms) :: Rij
    
        ! Local
        real(dp) :: rik_val, ril_val,dotproduct,sumt,costheta,ql, one_third, three_over_eight  
        real(dp),dimension(1:3) :: rik, ril 
        integer :: i, k, l

        ! Return 
        real(c_double), intent(inout), dimension(1:total_atoms) :: q_tetrahedral

        q_tetrahedral = 0.0d0

        one_third = 1.0d0/3.0d0

        three_over_eight = 3.0d0/8.0d0

        do i = 1,total_atoms

            sumt = 0.0

            do k = 1,nnb-1

                rik = Rij(1:3,k,i)

                rik_val =  Rij(4,k,i)

                do l = k+1,nnb

                    ril = Rij(1:3,l,i)

                    ril_val = Rij(4,l,i)

                    costheta = dot_product(rik, ril) / (rik_val * ril_val)

                    sumt = sumt + (costheta + one_third) * (costheta + one_third)

                end do

            end do

            q_tetrahedral(i) = 1 - three_over_eight * sumt

        end do

        end subroutine

end module 
