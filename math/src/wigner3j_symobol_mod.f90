!----------------------------------------------- Program Descriptions ----------------------------------------------- 

! This library includes subroutines and functions to compute: 
	!  Wigner3j symobl  

! The algorithm can be found in:  
! https://earlelab.rit.albany.edu/documents/Lai369j.pdf 
! Eq (7), Eq (10),Eq (11) 
! http://mathworld.wolfram.com/SphericalHarmonic.html 
! Eq (6)  
! Instructions: 
! 1. How to compute wigner3j symbol:  
! 	a. first, call function compute_associated_legendre_const first to get preconst 
!	b. then, call function copmute_associated_legendre_poly by passing preconst as argument  


! Date composed by Jingxiang Guo : 11/27/2019 

module wigner3j_symbol
	use constants,only: pi 
	use system, only: dp ! use double precision kind number 
	use factorial  
	implicit none 
	! all variables and subroutines/functions are private  
	private 


	! Export these subroutines/functions ( other subroutines/functions are still private) : 	
	public :: compute_associated_legendre_const, & 
			& copmute_associated_legendre_poly, & 
			& compute_spherical_harmonics_const, & 			
			& compute_spherical_harmonics  

	! Global variables

contains
	
	pure function compute_associated_legendre_const(l) result( preconst)  
		! compute associated legendre polynomial constant
		implicit none 

		! Passed: 
		integer,intent(in) :: l 

		! Local: 
		integer,dimension(0:l) :: m_list 
		integer :: i 

		! Return: 
		real(dp),dimension(0:l,4)  :: preconst

		m_list = [(i,i=0,l)] 	

		preconst = 0.0d0 

		preconst(:,1) =  (-1)**m_list*compute_doublefactorial(2*m_list -1 ) 

		preconst(:,2) = m_list/2.0d0
		
		preconst(:,3) = (2*m_list(0:l-1)+1) 
	
		preconst(:,4) = (-1)**m_list*compute_factorial(l-m_list)/compute_factorial(l+m_list) 

		end function 

	pure function compute_postive_associated_legendre_poly(preconst,l,x) result(Plm_pos)  
		! Implement Eq (7) , Eq (10) and Eq (11) 
		implicit none 

		! Passed:  
		real(dp),intent(in) :: x 
		integer,intent(in) :: l 
		real(dp),intent(in),dimension(0:l,4) :: preconst

		! Local:  	
		real(dp),dimension(0:l) :: P_ll,P_l_l1
		integer :: i,current_l,current_m 
		real(dp) :: p_lm_2, p_lm,p_lm_1,norm_1 
		
		! Return: 
		real(dp),dimension(0:l)  :: Plm_pos	
	
		! Identity: Eq (10)  
		P_ll = preconst(:,1)*(1-x*x)**preconst(:,2)
		
		! Identity: Eq (11) 	

		P_l_l1(0:l-1) = x*preconst(:,3)*P_ll(0:l-1) 
		
		! Put last P_ll in P_l_l1 since P_l_l1 is one element less than P_ll
		P_l_l1(l) = P_ll(l) 

		! The last two elements of Plm are calculated directly from P_ll and P_l_l1
		Plm_pos(l-1:l) = P_l_l1(l-1:l) 

		! loop backward from l-2 to 0 ( horizontal move )  
		do current_l = 0,l-2  

			p_lm = P_ll(current_l) 

			p_lm_1 = P_l_l1(current_l) 

			! get current l from previous l-1, and l-2 results ( vertical move)  
			do current_m = current_l + 2 , l 					

				norm_1 = 1.0d0/dble((current_m-current_l)) 

				! Eq (7) 
				p_lm_2 = x*(2*current_m-1)*norm_1*p_lm_1 - p_lm*(current_m+current_l-1)*norm_1
		
				p_lm = p_lm_1 

				p_lm_1 = p_lm_2	

			end do 
	
			Plm_pos(current_l) = p_lm_2 
											
		end do 
			
		end function 

	pure function compute_negative_associated_legendre_poly(l,Plm_pos,preconst) result( Plm_neg) 
		! compute negative part of associated legendre polynomial 
		implicit none 
		! Passed:  
		integer,intent(in) :: l 
		real(dp),intent(in),dimension(0:l) :: Plm_pos
		real(dp),intent(in),dimension(0:l,4) :: preconst

		! Local:  	
		integer :: i,indx 	

		! Return:
		real(dp),dimension(-l:-1) :: Plm_neg 			
	
		do i = -l,-1 	

			indx = i + 2*(-i) 
		
			Plm_neg(i) = Plm_pos(indx)*preconst(indx,4) 			

		end do 

		end function 

	pure function copmute_associated_legendre_poly(preconst,x,l) result(Plm) 

		! combine both positive and negative part of associated legendre  polynomial 
		implicit none 
		! Passed 

		real(dp),intent(in) :: x 
		integer,intent(in) :: l 
		real(dp),intent(in),dimension(0:l,4) :: preconst

		! Local 

		real(dp),dimension(-l:-1) :: Plm_neg
		real(dp),dimension(0:l) :: Plm_pos

		! Return 	

		real(dp),dimension(-l:l) :: Plm
		
		Plm_pos =  compute_postive_associated_legendre_poly(preconst,l,x) 

		! put positive part inside Plm_pos 

		Plm(0:l) = Plm_pos 

		! add negative part  
		Plm(-l:-1) = Plm_pos(l:1:-1)*preconst(l:1:-1,4)  
		
		end function 

!---------------------------------------------------------------------------------------
!-------------------------------- Spherical Harmonics ----------------------------------
!---------------------------------------------------------------------------------------

	pure function compute_spherical_harmonics_const(l) result(preconst)  
	 	!compute spherical harmonics normalization constant(l,m) as prefactor  
		implicit none 
		! Passed 
		integer,intent(in) :: l 

		! Local 
		integer :: i 
		integer,dimension(-l:l) :: m_list 

		! Return 	
		real(dp),dimension(-l:l)  :: preconst

		m_list(-l:l) = [(i,i=-l,l,1)] 

		preconst = dsqrt((2*l+1)*compute_factorial(l-m_list)/(4*pi*compute_factorial(l+m_list))) 

		end function 
	
	pure function compute_exp_phi(l,cos_phi,sin_phi) result(Exp_phi)  
		implicit none 

		! Passed 
		integer,intent(in) :: l 	
		real(dp),intent(in) :: sin_phi,cos_phi
		
		! Local:  
		integer :: i
		real(dp),dimension(0:l) :: angle_lst
		real(dp),dimension(-l:l,2) :: sin_cos_phi 

		! Return:  
		complex(dp),dimension(-l:l) :: Exp_phi ! first column is cos and second column is sin

		sin_cos_phi(0,1) = 1  

		sin_cos_phi(0,2) = 0  

		sin_cos_phi(1,1) = cos_phi 

		sin_cos_phi(1,2) = sin_phi 

		do i = 2, l
		
			sin_cos_phi(i,1) = sin_cos_phi(i-1,1)*cos_phi - sin_cos_phi(i-1,2)*sin_phi 			
	
			sin_cos_phi(i,2) = sin_cos_phi(i-1,2)*cos_phi + sin_cos_phi(i-1,1)*sin_phi 	

		end do 		

		sin_cos_phi(-l:-1,1) = sin_cos_phi(l:1:-1,1)  

		sin_cos_phi(-l:-1,2) = sin_cos_phi(l:1:-1,2)  

		Exp_phi = cmplx(sin_cos_phi(:,1), -sin_cos_phi(:,2))  	

		end function 

	pure function compute_spherical_harmonics(sph_const,Plm_const,l,costheta,sin_cos_phi) result(Ylm) 
		! compute spherical harmonics ; constant (l,m)*associated_legendre_polynomial*expotential(mphi) 	
		! Implemet Eq(6) in http://mathworld.wolfram.com/SphericalHarmonic.html 
		implicit none 	
        !Passed 
        integer,intent(in) :: l 
        real(dp),intent(in),dimension(-l:l) :: sph_const
		real(dp),intent(in),dimension(2) :: sin_cos_phi 
        real(dp),intent(in),dimension(0:l,4) :: Plm_const  
        real(dp),intent(in) :: costheta

        ! Local 
		real(dp) :: start,finish
        real(dp),dimension(-l:l)  :: Plm
        complex(dp),dimension(-l:l) ::  Exp_phi

        ! Return 
        complex(dp),dimension(-l:l) :: Ylm 

		Plm = copmute_associated_legendre_poly(Plm_const,costheta,l) 

		Exp_phi = compute_exp_phi(l,sin_cos_phi(1),sin_cos_phi(2)) 
		
		!Exp_phi = compute_exp_phi_backup(l,0.5d0) 
		
		Ylm = sph_const*Plm*Exp_phi 

		end function 

end module 
