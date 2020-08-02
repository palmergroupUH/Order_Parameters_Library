!----------------------------------------------- Program Descriptions ----------------------------------------------- 

! This module includes subroutines and functions to compute Wigner3j symobl: 

! The implementation can be found in [1]: 

! Reference:

! [1] S.-T. Lai. Computation of algebraic formulas for wigner 3-j, 6-j, and 9-j
! symbols by maple.International Journal of Quantum Chemistry, 52(3):593–607,1994

! [2] https://mathworld.wolfram.com/Wigner3j-Symbol.html 

! Note:
! Implementation of wigner3j symbol:
! Reference Wigner3j symbol values for comparison
! j1 j2 j3 m1 m2 m3 
! 4  4  4  0 0 0  
! 4  4  4  0 -1 1   

! 4  4  4  -1 1 0   
! 4  4  4  1 -1 0   
! 4  4  4  0 2 -2   
! 4  4  4  0 -3 3   
! 4  4  4  -4 4 0   

! Date composed by Jingxiang Guo : 8/1/2020 

module wigner3j_symbol
	use system, only: c_int, c_double, dp, error_unit ! use c double precision and units for error stream
	use factorial, only: compute_factorial 
	implicit none 
	! all variables and subroutines/functions are private  
	private 


	! Export these subroutines/functions ( other subroutines/functions are still private) : 	
	public :: compute_wigner3j_symbol 

	! Global variables

contains

! -------------------------------------------------------------------------------------------------
!                                        Rules Satisfied
! -------------------------------------------------------------------------------------------------

    logical function satisfy_rule1(j1,j2,j3,m1,m2,m3)
        implicit none
        integer,intent(in) :: j1,j2,j3,m1,m2,m3

        if ( abs(m1) <= j1 .and. abs(m2) <= j2 .and. abs(m3) <= j3 ) then

            satisfy_rule1 = .TRUE.

        else

            write(error_unit,*) "rule 1 does not satisfy for Wigner 3j calculations: "

            write(error_unit,*) "abs(m1) < j1, abs(m2) < j2, abs(m3) < j3 "

            satisfy_rule1 = .FALSE.

        end if

        end function

    logical function satisfy_rule2(m1,m2,m3)
        implicit none
        integer,intent(in) :: m1,m2,m3

        if ( m1+m2+m3 == 0 ) then

            satisfy_rule2 = .TRUE.

        else

            write(error_unit,*) "rule 2 does not satisfy for Wigner 3j calculations: "

            write(error_unit,*) "m1+m2+m3 == 0"

            satisfy_rule2 = .FALSE.

        end if

        end function

    logical function satisfy_rule3(j1,j2,j3)
        implicit none
        integer,intent(in) :: j1,j2,j3

        ! triangular inequalities:

        if ( abs(j1-j2) <= j3 .and. j3 <= j1+j2) then

            satisfy_rule3 = .TRUE.

        else

            write(error_unit,*) "rule 3 does not satisfy for Wigner 3j calculations"

            write(error_unit,*) "Does not satisfy the triangluar inequalities "

            satisfy_rule3 = .FALSE.

        end if

        end function

    logical function satisfy_rule4(jarray)
        implicit none
        real(dp),intent(in),dimension(:) :: jarray
        real(dp) :: remainder

        remainder = dmod(sum(jarray),1.0d0)

        if ( remainder  == 0.0d0 ) THEN

            satisfy_rule4 =.TRUE.

        ELSE

            write(error_unit,*) "rule 4 does not satisfy for Wigner 3j calculations"

            write(error_unit,*) "j1+j2+j3 must be integer"

            satisfy_rule4 = .FALSE.

        end if

        end function

    pure real(dp) function wigner_3j(j1,j2,j3,m1,m2,m3)
        implicit none 
        ! Passed
        integer, intent(in) :: j1,j2,j3,m1,m2,m3
       
        ! Local
        integer ::  z, firstterm,u,lower_z,upper_z 
        real(dp) :: delta, term2,sumterm,secondterm, second_const 
        ! Return
        
    
        ! Ref. [1]   
        delta = dsqrt((compute_factorial(j1+j2-j3) * & 
                     & compute_factorial(j1-j2+j3) * & 
                     & compute_factorial(-j1+j2+j3)) &
                     & /compute_factorial(j1+j2+j3+1))

        second_const = dsqrt(compute_factorial(j3-m3)*compute_factorial(j3+m3) & 
                          & /(compute_factorial(j2+m2) * & 
                            & compute_factorial(j2-m2) * & 
                            & compute_factorial(j1-m2-m3) * & 
                            & compute_factorial(j1+m2+m3)))

        u = -j1 + j3 + j2 

        ! z lower and upper only shown in the sample program section of Ref. [1] 
        lower_z = max(-j1-m2-m3,u-j3-m3,0) 

        upper_z = min(j2+j3-m2-m3,j3-m3,u)  

        wigner_3j = 0.0d0 

        sumterm = 0.0d0 

        do z = lower_z,upper_z 
    
            firstterm = (-1)**(2*j2-j1-m1+z) 

            secondterm = compute_factorial(j2+j3-m2-m3-z)*compute_factorial(j1+m2+m3+z) & 
                        & /( compute_factorial(z)* & 
                          & compute_factorial(j3-m3-z) * & 
                          & compute_factorial(u-z) * & 
                          & compute_factorial(j3+m3-u+z))

            sumterm = sumterm + firstterm*secondterm 
        
        end do 

        wigner_3j = delta*second_const*sumterm

        end function

    real(c_double) function compute_wigner3j_symbol(j1,j2,j3,m1,m2,m3) bind(c, name="call_wigner3j_symobol")
        implicit none
        integer(c_int), intent(in) :: j1,j2,j3,m1,m2,m3

        ! Check if all 3 rules satisfied (rule1 - rule 3 are independent)  
        if (    &  
           & satisfy_rule3(j1,j2,j3) .AND. &

           & satisfy_rule2(m1,m2,m3) .AND. &

           & satisfy_rule1(j1,j2,j3,m1,m2,m3)) THEN

            compute_wigner3j_symbol = wigner_3j(j1,j2,j3,m1,m2,m3)

        ! if not, wigner3j symobol = 0.0
        else

            compute_wigner3j_symbol = 0.0d0

        end if

        end function

end module	
