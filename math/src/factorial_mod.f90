!----------------------------------------------- Program Descriptions ----------------------------------------------- 
! This program computes the factorial or double factorial for an integer number: n  
! Input : n  
! Output: n! or n!! 
! By Jingxiang Guo on 11/27/2019 

module factorial 
	use system, only: dp
	implicit none 
	! all variables and subroutines/functions are private  
	private 

	! Global variables: 	

	! Export these subroutines/functions ( other subroutines/functions are still private) : 	
	public :: compute_factorial,&
			& compute_factorial_recur,&
		    & compute_doublefactorial  

	! all subroutines/functions return double precision 

contains

	! elemental: perform factorial calculations for each element of an array. So, the
	! input n, can be both a scalar or an array  
	elemental real(dp) function compute_factorial(n) 
		implicit none 
		! Passed 
		integer,intent(in) :: n 	
		integer :: i 

		compute_factorial = 1.0d0  

		do i =  n,1,-1 

			compute_factorial = compute_factorial*i 	

		end do 	

		end function 

	! Method 2: compute factorial by recursive functions c interoperable 
	recursive function compute_factorial_recur(n) result(fact) 
		implicit none 
    
        ! Passed
		integer,intent(in) :: n 
		integer :: i 
        real(dp) :: fact

		if ( n <= 1) then   

			fact = 1.0d0  

		else 
		
			fact = n*compute_factorial_recur(n-1) 
		
		end if 	

		end function 	
	
	! compute double factorial: n!! = n(n-2)(n-4).......  
	elemental real(dp) function compute_doublefactorial(n) result( double_fact) 
		implicit none 
        ! Passed
		integer,intent(in) :: n 	
		integer :: i 

		double_fact = 1.0d0  

		do i = n,1,-2 

			double_fact = double_fact*i 
		
		end do 

		end function 

end module 
