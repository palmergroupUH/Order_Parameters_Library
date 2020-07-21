!----------------------------------------------- Program Descriptions ----------------------------------------------- 

! This library includes subroutines and functions to compute: 
	!  Associated legendre polynomial 
	!  Spherical harmonics ( complex number ) for any given order l, m  

! The algorithm implemented for computing associated legendre polynomial and spherical harmonics can be found in:  
! http://mathworld.wolfram.com/AssociatedLegendrePolynomial.html 
! Eq (7), Eq (10),Eq (11) 
! http://mathworld.wolfram.com/SphericalHarmonic.html 
! Eq (6)  
! Instructions: 
! 1. How to compute associated legendre polynomial: 
! 	a. first, call function "compute_associated_legendre_const" first to get preconst 
!	b. then, call function "copmute_associated_legendre_poly" by passing preconst as argument  
! 2. How to compute spherical harmonics
!   a. first call function "compute_spherical_harmonics_const" to get preconst  
!	b. then, call function "compute_spherical_harmonics" by passing preconst as argument 
!   to get complex Ylm with dimension( -l : l ) 


! Date composed by Jingxiang Guo : 11/27/2019 

module spherical_harmonics
	use system, only: dp, c_double, c_int ! use double precision kind number
	use constants,only: pi 
	use factorial  
	implicit none 
	! all variables and subroutines/functions are private  
	private 


	! Export these subroutines/functions ( other subroutines/functions are still private) : 	
	public :: compute_associated_legendre_const, & 
			& copmute_associated_legendre_poly, & 
			& compute_spherical_harmonics_const, & 			
			& compute_spherical_harmonics,& 
            & optimized_Q12 

	! Global variables

contains
	
! ----------------------------------------------------------------------------------------------------
!                      Associated Legendre Polynomial 
! ----------------------------------------------------------------------------------------------------

	pure subroutine compute_associated_legendre_const(l, preconst) bind(c, name="call_associated_legendre_const")
		! compute the constants of associated legendre polynomial 
		implicit none 

		! Passed: 
		integer(c_int),intent(in) :: l 

		! Local: 
		integer,dimension(0:l) :: m_list 
		integer :: i 

		! Return: 
		real(c_double),intent(out), dimension(0:l,4)  :: preconst

		m_list = [(i,i=0,l)] 	

		preconst = 0.0d0 

		preconst(:,1) =  (-1)**m_list*compute_doublefactorial(2*m_list -1 ) 

		preconst(:,2) = m_list/2.0d0
		
		preconst(:,3) = (2*m_list(0:l-1)+1) 
	
		preconst(:,4) = (-1)**m_list*compute_factorial(l-m_list)/compute_factorial(l+m_list) 

		end subroutine  

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

	pure subroutine copmute_associated_legendre_poly(preconst,x,l, Plm) bind(c, name="call_associated_legendre_poly")

		! combine both positive and negative part of associated legendre  polynomial 

		implicit none 
		! Passed 

		real(c_double),intent(in) :: x 
		integer(c_int),intent(in) :: l 
		real(c_double),intent(in),dimension(0:l,4) :: preconst

		! Local 

		real(dp),dimension(-l:-1) :: Plm_neg
		real(dp),dimension(0:l) :: Plm_pos

		! Return 	

		real(c_double),intent(out),dimension(-l:l) :: Plm
		
		Plm_pos =  compute_postive_associated_legendre_poly(preconst,l,x) 

		! put positive part inside Plm_pos 

		Plm(0:l) = Plm_pos 

		! add negative part  
		Plm(-l:-1) = Plm_pos(l:1:-1) * preconst(l:1:-1,4)  
		
		end subroutine 

! ------------------------------------------------------------------------------------------------- 
!                                 Spherical Harmonics
!------------------------------------------------------------------------------------------------- 

	pure subroutine compute_spherical_harmonics_const(l, preconst) bind(c, name="call_spherical_harmonics_const")
	 	!compute spherical harmonics normalization constant(l,m) as prefactor  
		implicit none 
		! Passed 
		integer(c_int),intent(in) :: l 

		! Local 
		integer :: i 
		integer,dimension(-l:l) :: m_list 

		! Return 	
		real(c_double),intent(out), dimension(-l:l)  :: preconst

		m_list(-l:l) = [(i,i=-l,l,1)]

		preconst = dsqrt((2*l+1)*compute_factorial(l-m_list)/(4*pi*compute_factorial(l+m_list)))

		end subroutine 
	
	pure function compute_exp_phi(l,cos_phi,sin_phi) result(Exp_phi)  
		implicit none 

		! Passed 
		integer,intent(in) :: l 	
		real(dp),intent(in) :: cos_phi, sin_phi
		
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

		sin_cos_phi(-l:-1,2) = -sin_cos_phi(l:1:-1,2)  

		Exp_phi = cmplx(sin_cos_phi(:,1), sin_cos_phi(:,2))  	

		end function 

	subroutine compute_spherical_harmonics(sph_const,Plm_const,l,costheta,sin_cos_phi, Ylm)
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
        complex(dp),intent(out), dimension(-l:l) :: Ylm

		call copmute_associated_legendre_poly(Plm_const,costheta,l, Plm) 

		Exp_phi = compute_exp_phi(l,sin_cos_phi(1),sin_cos_phi(2)) 
		
		!Exp_phi = compute_exp_phi_backup(l,0.5d0) 
        
		Ylm = sph_const*Plm*Exp_phi 

		end subroutine 

    ! This subroutine is exposed as c callable 
    subroutine call_spherical_harmonics(sph_const,&
                                           & Plm_const,&
                                           & l, &
                                           & cos_theta, &
                                           & sin_cos_phi,&
                                           & Ylm_complex) bind(c, name="call_spherical_harmonics")

		implicit none 	
        !Passed 
        integer(c_int),intent(in) :: l 
        real(c_double),intent(in),dimension(-l:l) :: sph_const
        real(c_double),intent(in),dimension(0:l,4) :: Plm_const  
        real(c_double),intent(in) :: cos_theta
		real(c_double),intent(in),dimension(2) :: sin_cos_phi

        ! Local: 
        complex(dp), dimension(-l:l) :: Ylm

        ! Return:
        real(c_double),intent(out), dimension(-l:l,2) :: Ylm_complex 
       
        call compute_spherical_harmonics(sph_const,Plm_const,l,cos_theta,sin_cos_phi, Ylm)
        
        Ylm_complex(:,1) = real(Ylm, dp) 

        Ylm_complex(:,2) = dimag(Ylm)
        
        end subroutine 

    pure subroutine optimized_Q12(sin_theta, cos_theta, cos_phi, sin_phi, Ylm)
            implicit none 
            ! Passed
            real(dp), intent(in) :: sin_theta, cos_theta, cos_phi, sin_phi 

            ! Local
            real(dp),dimension(1:12) :: sint, cost, sinp, cosp 

            ! Return
            complex(dp), intent(out), dimension(-12:12) :: Ylm
       
            ! assign the inital parameters 
            sint(1) = sin_theta

            cost(1) = cos_theta

            sinp(1) = sin_phi

            cosp(1) = cos_phi

            ! calculate powers of sin(theta)
            sint(2) = sint(1)*sint(1)
            sint(3) = sint(2)*sint(1)
            sint(4) = sint(3)*sint(1) 
            sint(5) = sint(4)*sint(1)
            sint(6) = sint(5)*sint(1)
            sint(7) = sint(6)*sint(1) 
            sint(8) = sint(7)*sint(1)
            sint(9) = sint(8)*sint(1) 
            sint(10) = sint(9)*sint(1) 
            sint(11) = sint(10)*sint(1) 	
            sint(12) = sint(11)*sint(1)

            ! calculate powers of cos(theta)
            cost(2) = cost(1)*cost(1)
            cost(3) = cost(2)*cost(1)
            cost(4) = cost(3)*cost(1)
            cost(5) = cost(4)*cost(1)
            cost(6) = cost(5)*cost(1)
            cost(7) = cost(6)*cost(1)
            cost(8) = cost(7)*cost(1)
            cost(9) = cost(8)*cost(1)
            cost(10) = cost(9)*cost(1)
            cost(11) = cost(10)*cost(1)
            cost(12) = cost(11)*cost(1)


            ! calculate sin/cos (2phi)
            sinp(2) = 2.0d0*sinp(1)*cosp(1)
            cosp(2) = 2.0d0*cosp(1)**2 - 1.0d0

            ! calculate sin/cos (3phi)
            sinp(3) = sinp(2)*cosp(1) + cosp(2)*sinp(1)
            cosp(3) = cosp(2)*cosp(1) - sinp(2)*sinp(1)


            ! calculate sin/cos (4phi) 

            sinp(4) = sinp(3)*cosp(1) + cosp(3)*sinp(1) 
            cosp(4) = cosp(3)*cosp(1) - sinp(3)*sinp(1) 

            ! calculate sin/cos (5phi)

            sinp(5) = sinp(4)*cosp(1) + cosp(4)*sinp(1) 
            cosp(5) = cosp(4)*cosp(1) - sinp(4)*sinp(1) 

            ! calculate sin/cos (6phi) 
            sinp(6) = sinp(5)*cosp(1) + cosp(5)*sinp(1) 
            cosp(6) = cosp(5)*cosp(1) - sinp(5)*sinp(1) 

            ! calculate sin/cos (7phi) 

            sinp(7) = sinp(6)*cosp(1) + cosp(6)*sinp(1) 
            cosp(7) = cosp(6)*cosp(1) - sinp(6)*sinp(1) 

            ! calculate sin/cos (8phi)

            sinp(8) = sinp(7)*cosp(1) + cosp(7)*sinp(1) 
            cosp(8) = cosp(7)*cosp(1) - sinp(7)*sinp(1) 

            ! calculate sin/cos (9phi)

            sinp(9) = sinp(8)*cosp(1) + cosp(8)*sinp(1) 
            cosp(9) = cosp(8)*cosp(1) - sinp(8)*sinp(1) 

            ! Calculate sin/cos (10phi) 
            sinp(10) = sinp(9)*cosp(1) + cosp(9)*sinp(1) 
            cosp(10) = cosp(9)*cosp(1) - sinp(9)*sinp(1) 

            ! Calculate sin/cos (11phi) 

            sinp(11) = sinp(10)*cosp(1) + cosp(10)*sinp(1) 
            cosp(11) = cosp(10)*cosp(1) - sinp(10)*sinp(1) 

            ! Calculate sin/cos (12phi)

            sinp(12) = sinp(11)*cosp(1) + cosp(11)*sinp(1) 
            cosp(12) = cosp(11)*cosp(1) - sinp(11)*sinp(1) 

		    ! Update the harmonic sum
		
		    Ylm(-12) = 0.566266663742191d0*cmplx(cosp(12),-sinp(12))*sint(12)
		    Ylm(-11) = 2.77412876903310d0*CMPLX(cosp(11),-sinp(11))*sint(11)*cost(1)
		    Ylm(-10) = 0.409022972331817d0*CMPLX(cosp(10),-sinp(10))*sint(10)*(-1+23*cost(2))
		    Ylm(-9) = 1.10763944520068d0*CMPLX(cosp(9),-sinp(9))*sint(9)*cost(1)*(-3+23*cost(2))
		    Ylm(-8) = 0.362560114310785d0*CMPLX(cosp(8),-sinp(8))*sint(8)*(1-42*cost(2)+161*cost(4))
		    Ylm(-7) = 0.725120228621570d0*CMPLX(cosp(7),-sinp(7))*cost(1)*(5-70*cost(2)+161*cost(4))*sint(7)
		    Ylm(-6) = 0.06791373178178367d0*CMPLX(cosp(6),-sinp(6))*(-5+285*cost(2)-1995*cost(4)+3059*cost(6))*sint(6)
		    Ylm(-5) = 0.762329748554085d0*CMPLX(cosp(5),-sinp(5))*(-5+95*cost(2)-399*cost(4)+437*cost(6))*sint(5)*cost(1)
		    Ylm(-4) = 0.06536923664453508d0*CMPLX(cosp(4),-sinp(4))*(5-340*cost(2)+3230*cost(4)-9044*cost(6)+7429*cost(8))*sint(4)
		    Ylm(-3) = 0.08715898219271344d0*CMPLX(cosp(3),-sinp(3))*cost(1)*(45-1020*cost(2)+5814*cost(4)-11628*cost(6)+7429*cost(8))*sint(3)
		    Ylm(-2) = 0.106747516436237d0*CMPLX(cosp(2),-sinp(2))*(-3+225*cost(2)-2550*cost(4)+9690*cost(6)-14535*cost(8)+7429*cost(10))*sint(2)
		    Ylm(-1) = 0.01720392001939924d0*CMPLX(cosp(1),-sinp(1))*cost(1)*(-231+5775*cost(2)-39270*cost(4)+106590*cost(6)-124355*cost(8)+52003*cost(10))*sint(1)
		    Ylm(0) = 0.001377415975458390d0*(231-18018*cost(2)+225225*cost(4)-1021020*cost(6)+2078505*cost(8)-1939938*cost(10)+676039*cost(12))
		    Ylm(1) = -0.01720392001939924d0*CMPLX(cosp(1),sinp(1))*cost(1)*(-231+5775*cost(2)-39270*cost(4)+106590*cost(6)-124355*cost(8)+52003*cost(10))*sint(1)
		    Ylm(2) = 0.106747516436237d0*CMPLX(cosp(2),sinp(2))*(-3+225*cost(2)-2550*cost(4)+9690*cost(6)-14535*cost(8)+7429*cost(10))*sint(2)
		    Ylm(3) = -0.08715898219271344d0*CMPLX(cosp(3),sinp(3))*cost(1)*(45-1020*cost(2)+5814*cost(4)-11628*cost(6)+7429*cost(8))*sint(3)
		    Ylm(4) = 0.06536923664453508d0*CMPLX(cosp(4),sinp(4))*(5-340*cost(2)+3230*cost(4)-9044*cost(6)+7429*cost(8))*sint(4)
		    Ylm(5) = -0.762329748554085d0*CMPLX(cosp(5),sinp(5))*(-5+95*cost(2)-399*cost(4)+437*cost(6))*sint(5)*cost(1)
		    Ylm(6) = 0.06791373178178367d0*CMPLX(cosp(6),sinp(6))*(-5+285*cost(2)-1995*cost(4)+3059*cost(6))*sint(6)
		    Ylm(7) = -0.725120228621570d0*CMPLX(cosp(7),sinp(7))*cost(1)*(5-70*cost(2)+161*cost(4))*sint(7)
		    Ylm(8) = 0.362560114310785d0*CMPLX(cosp(8),sinp(8))*sint(8)*(1-42*cost(2)+161*cost(4))
		    Ylm(9) = -1.10763944520068d0*CMPLX(cosp(9),sinp(9))*sint(9)*cost(1)*(-3+23*cost(2))
		    Ylm(10) = 0.409022972331817d0*CMPLX(cosp(10),sinp(10))*sint(10)*(-1 + 23*cost(2))
		    Ylm(11) = -2.77412876903310d0*CMPLX(cosp(11),sinp(11))*sint(11)*cost(1)
		    Ylm(12) = 0.566266663742191d0*CMPLX(cosp(12),sinp(12))*sint(12)

        end subroutine 

   ! This subroutine is exposed as c callable 
   pure subroutine call_optimized_Q12(sin_theta, cos_theta, cos_phi, sin_phi, Ylm_complex) bind(c, name="call_optimized_Y12")
            implicit none 
            ! Passed
            real(c_double), intent(in) :: sin_theta, cos_theta, cos_phi, sin_phi 

            ! Local
            real(dp),dimension(1:12) :: sint, cost, sinp, cosp 
            complex(dp), dimension(-12:12) :: Ylm 

            ! Return
            real(c_double), intent(out), dimension(-12:12,2) :: Ylm_complex

            call optimized_Q12(sin_theta, cos_theta, cos_phi, sin_phi, Ylm)
                   
            Ylm_complex(:,1) = real(Ylm, dp) 

            Ylm_complex(:,2) = dimag(Ylm)

            end subroutine 

! -------------------------------------------------------------------------------------------------
!                                           Wigner 3j
! -------------------------------------------------------------------------------------------------
	 
end module 
