!----------------------------------------------- Program Descriptions ----------------------------------------------- 

! This library includes subroutines and functions to compute: 
	!  Wigner3j symobl  


! The implementation can be found in [1]: 

! Reference:
! [1] S.-T. Lai. Computation of algebraic formulas for wigner 3-j, 6-j, and 9-j
! symbols by maple.International Journal of Quantum Chemistry, 52(3):593–607,1994

! Instructions: 
! 1. How to compute wigner3j symbol:  
! 	a. 
!	b. 


! Date composed by Jingxiang Guo : 11/27/2019 

module wigner3j_symbol
	use constants,only: pi 
	use system, only: dp ! use double precision kind number 
	use factorial  
    use spherical_harmonics 
	implicit none 
	! all variables and subroutines/functions are private  
	private 


	! Export these subroutines/functions ( other subroutines/functions are still private) : 	
	public :: compute_associated_legendre_const, & 
			& copmute_associated_legendre_poly, & 
			& compute_spherical_harmonics_const, & 			

	! Global variables

contains


end module	
