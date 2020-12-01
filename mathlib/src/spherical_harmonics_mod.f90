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
    ! To other Fortran subroutines/functions
	public :: compute_associated_legendre_const, & 
			& copmute_associated_legendre_poly, & 
			& compute_spherical_harmonics_const, & 			
			& compute_spherical_harmonics,& 
            & optimized_Y12, &
            & optimized_Y8, & 
            & optimized_Y6, & 
            & optimized_Y4 

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

		Exp_phi = dcmplx(sin_cos_phi(:,1), sin_cos_phi(:,2))  	

		end function 

	subroutine compute_spherical_harmonics(sph_const,Plm_const,l,costheta, cos_phi, sin_phi, Ylm)
		! compute spherical harmonics ; constant (l,m)*associated_legendre_polynomial*expotential(mphi) 	
		! Implemet Eq(6) in http://mathworld.wolfram.com/SphericalHarmonic.html 
		implicit none 	
        !Passed 
        integer,intent(in) :: l 
        real(dp),intent(in),dimension(-l:l) :: sph_const
		real(dp),intent(in) :: cos_phi, sin_phi 
        real(dp),intent(in),dimension(0:l,4) :: Plm_const  
        real(dp),intent(in) :: costheta

        ! Local 
		real(dp) :: start,finish
        real(dp),dimension(-l:l)  :: Plm
        complex(dp),dimension(-l:l) ::  Exp_phi

        ! Return 
        complex(dp),intent(out), dimension(-l:l) :: Ylm

		call copmute_associated_legendre_poly(Plm_const,costheta,l, Plm) 

		Exp_phi = compute_exp_phi(l,cos_phi,sin_phi) 
		
		!Exp_phi = compute_exp_phi_backup(l,0.5d0) 
        
		Ylm = sph_const*Plm*Exp_phi 

		end subroutine 

    ! This subroutine is exposed as c callable 
    subroutine call_spherical_harmonics(sph_const,&
                                           & Plm_const,&
                                           & l, &
                                           & cos_theta, &
                                           & cos_phi,&
                                           & sin_phi,&
                                           & Ylm_complex) bind(c, name="call_spherical_harmonics")

		implicit none 	
        !Passed 
        integer(c_int),intent(in) :: l 
        real(c_double),intent(in),dimension(-l:l) :: sph_const
        real(c_double),intent(in),dimension(0:l,4) :: Plm_const  
        real(c_double),intent(in) :: cos_theta
		real(c_double),intent(in) :: cos_phi, sin_phi

        ! Local: 
        complex(dp), dimension(-l:l) :: Ylm

        ! Return:
        real(c_double),intent(out), dimension(-l:l,2) :: Ylm_complex 
       
        call compute_spherical_harmonics(sph_const,Plm_const,l,cos_theta, cos_phi, sin_phi, Ylm)
        
        Ylm_complex(:,1) = real(Ylm, dp) 

        Ylm_complex(:,2) = dimag(Ylm)
        
        end subroutine 

    pure subroutine optimized_Y12(sin_theta, cos_theta, cos_phi, sin_phi, Y)
            implicit none 
            ! Passed
            real(dp), intent(in) :: sin_theta, cos_theta, cos_phi, sin_phi 

            ! Local
            real(dp),dimension(1:12) :: sint, cost, sinp, cosp 
            real(dp), dimension(-12:12) :: pre_and_post 

            ! Return
            complex(dp), intent(out), dimension(-12:12) :: Y
       
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
	        ! From Mathmatica function: "SphericalHarmonicsY[12,m,theta, phi]"	

            pre_and_post(-12) = 0.566266663742191d0 * sint(12)
            pre_and_post(-11) = 2.77412876903310d0 * sint(11)*cost(1)
            pre_and_post(-10) = 0.409022972331817d0 * sint(10)*(-1.0d0+23.0d0*cost(2))
            pre_and_post(-9) = 1.10763944520068d0 * sint(9)*(-3.0d0*cost(1) +23.0d0*cost(3))
            pre_and_post(-8) = 0.362560114310785d0 * sint(8)*(1.0d0-42.0d0*cost(2)+161.0d0*cost(4))
            pre_and_post(-7) = 0.725120228621570d0 * sint(7) * (5.0d0*cost(1)-70.0d0*cost(3)+161.0d0*cost(5))
            pre_and_post(-6) = 0.06791373178178367d0 * sint(6)*(-5.0d0+285.0d0*cost(2)-1995.0d0*cost(4)+3059.0d0*cost(6))
            pre_and_post(-5) = 0.762329748554085d0 * sint(5) * (-5.0d0*cost(1) + 95.0d0 *cost(3)-399.0d0 *cost(5)+437.0d0 *cost(7))
            pre_and_post(-4) = 0.06536923664453508d0 * sint(4) * (5.0d0 - 340.0d0 *cost(2)+3230.0d0 *cost(4)-9044.0d0 *cost(6)+7429.0d0 *cost(8))
            pre_and_post(-3) = 0.08715898219271344d0 * sint(3) * (45.0d0*cost(1)-1020.0d0 *cost(3)+ 5814.0d0 *cost(5)-11628.0d0 *cost(7)+ 7429.0d0 *cost(9))
            pre_and_post(-2) = 0.106747516436237d0 * sint(2) *(-3.0d0 + 225.0d0 *cost(2)-2550.0d0 *cost(4)+ 9690.0d0 *cost(6)-14535.0d0 *cost(8)+7429.0d0 *cost(10))
            pre_and_post(-1) = 0.01720392001939924d0 * sint(1) * (-231.0d0 * cost(1) +5775.0d0 *cost(3)-39270.0d0 *cost(5)+106590.0d0 *cost(7)-124355.0d0 *cost(9)+52003.0d0 *cost(11))
            pre_and_post(0) = 0.001377415975458390d0 * (231.0d0 -18018.0d0 *cost(2)+225225.0d0 *cost(4)-1021020.0d0 *cost(6)+2078505.0d0 *cost(8)-1939938.0d0 *cost(10)+676039.0d0 *cost(12))

            pre_and_post(1) = -pre_and_post(-1)
            pre_and_post(2) = pre_and_post(-2)
            pre_and_post(3) = -pre_and_post(-3)
            pre_and_post(4) = pre_and_post(-4)
            pre_and_post(5) = -pre_and_post(-5)
            pre_and_post(6) = pre_and_post(-6)
            pre_and_post(7) = -pre_and_post(-7)
            pre_and_post(8) = pre_and_post(-8)
            pre_and_post(9) = -pre_and_post(-9)
            pre_and_post(10) = pre_and_post(-10)
            pre_and_post(11) = -pre_and_post(-11)
            pre_and_post(12) = pre_and_post(-12)


            Y(-12) = dCMPLX(cosp(12),-sinp(12))* pre_and_post(-12)
            Y(-11) = dCMPLX(cosp(11),-sinp(11))* pre_and_post(-11)
            Y(-10) = dCMPLX(cosp(10),-sinp(10))* pre_and_post(-10)
            Y(-9) = dCMPLX(cosp(9),-sinp(9))* pre_and_post(-9)
            Y(-8) = dCMPLX(cosp(8),-sinp(8))* pre_and_post(-8)
            Y(-7) = dCMPLX(cosp(7),-sinp(7))* pre_and_post(-7)
            Y(-6) = dCMPLX(cosp(6),-sinp(6))* pre_and_post(-6)
            Y(-5) = dCMPLX(cosp(5),-sinp(5))* pre_and_post(-5)
            Y(-4) = dCMPLX(cosp(4),-sinp(4))* pre_and_post(-4)
            Y(-3) = dCMPLX(cosp(3),-sinp(3))* pre_and_post(-3)
            Y(-2) = dCMPLX(cosp(2),-sinp(2))* pre_and_post(-2)
            Y(-1) = dCMPLX(cosp(1),-sinp(1))* pre_and_post(-1)
            Y(0) = pre_and_post(0)
            Y(1) = dCMPLX(cosp(1),sinp(1))* pre_and_post(1)
            Y(2) = dCMPLX(cosp(2),sinp(2))* pre_and_post(2)
            Y(3) = dCMPLX(cosp(3),sinp(3))* pre_and_post(3)
            Y(4) = dCMPLX(cosp(4),sinp(4))* pre_and_post(4)
            Y(5) = dCMPLX(cosp(5),sinp(5))* pre_and_post(5)
            Y(6) = dCMPLX(cosp(6),sinp(6))* pre_and_post(6)
            Y(7) = dCMPLX(cosp(7),sinp(7))* pre_and_post(7)
            Y(8) = dCMPLX(cosp(8),sinp(8))* pre_and_post(8)
            Y(9) = dCMPLX(cosp(9),sinp(9))* pre_and_post(9)
            Y(10) = dCMPLX(cosp(10),sinp(10))* pre_and_post(10)
            Y(11) = dCMPLX(cosp(11),sinp(11))* pre_and_post(11)
            Y(12) = dCMPLX(cosp(12),sinp(12))* pre_and_post(12)
                            
          end subroutine 

    pure subroutine optimized_Y8(sin_theta, cos_theta, cos_phi, sin_phi, Y)
        implicit none 
        ! Passed
        real(dp), intent(in) :: sin_theta, cos_theta, cos_phi, sin_phi 

        ! Local
        real(dp), dimension(1:8) :: sint, cost, sinp, cosp 
        real(dp), dimension(-8:8) :: pre_and_post

        ! Return
        complex(dp), intent(out), dimension(-8:8) :: Y

        Y = 0.0d0       

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

        ! calculate powers of cos(theta)
        cost(2) = cost(1)*cost(1)
        cost(3) = cost(2)*cost(1)
        cost(4) = cost(3)*cost(1)
        cost(5) = cost(4)*cost(1)
        cost(6) = cost(5)*cost(1)
        cost(7) = cost(6)*cost(1)
        cost(8) = cost(7)*cost(1)

        ! calculate sin/cos (2phi)
        sinp(2) = 2.0d0*sinp(1)*cosp(1)
        cosp(2) = 2.0d0*cosp(1)*cosp(1) - 1.0d0

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

        pre_and_post(-8) = 0.515428984397284d0 * sint(8)
        pre_and_post(-7) = 2.06171593758914d0 * cost(1) * sint(7) 
        pre_and_post(-6) = 0.376416108728495d0 * (-1.0d0 + 15.0d0 * cost(2) ) * sint(6) 
        pre_and_post(-5) = 2.43945519537307d0 * (cost(1) * (-1.0d0 + 5.0d0 * cost(2))*sint(5)) 
        pre_and_post(-4) = 0.338291568889025d0 * (1.0d0 - 26.0d0 * cost(2) + 65.0d0 * cost(4)) * sint(4)  
        pre_and_post(-3) = 0.873465074979714d0 * cost(1) * (3.0d0 - 26.0d0 * cost(2) + 39.0d0 * cost(4)) * sint(3) 
        pre_and_post(-2) = 0.322548355192883d0 * (-1.0d0 + 33.0d0 * cost(2) - 143.0d0 * cost(4) + 143.0d0 * cost(6)) * sint(2)  
        pre_and_post(-1) = 0.07710380440405712d0 * cost(1) * (-35.d0 + 385.0d0 * cost(2) - 1001d0*cost(4) + 715d0*cost(6)) * sint(1) 
        pre_and_post(0) = 0.009086770491564996d0 * (35.0d0 - 1260d0* cost(2) + 6930d0*cost(4) - 12012d0*cost(6) + 6435d0*cost(8))
        pre_and_post(1) = -pre_and_post(-1)
        pre_and_post(2) = pre_and_post(-2)
        pre_and_post(3) = -pre_and_post(-3)
        pre_and_post(4) = pre_and_post(-4)
        pre_and_post(5) = -pre_and_post(-5)
        pre_and_post(6) = pre_and_post(-6)
        pre_and_post(7) = -pre_and_post(-7)
        pre_and_post(8) = pre_and_post(-8)

        Y(-8) = dcmplx(cosp(8),-sinp(8)) * pre_and_post(-8)  
        Y(-7) = dcmplx(cosp(7),-sinp(7)) * pre_and_post(-7)  
        Y(-6) = dcmplx(cosp(6),-sinp(6)) * pre_and_post(-6)  
        Y(-5) = dcmplx(cosp(5),-sinp(5)) * pre_and_post(-5) 
        Y(-4) = dcmplx(cosp(4),-sinp(4)) * pre_and_post(-4)  
        Y(-3) = dcmplx(cosp(3),-sinp(3)) * pre_and_post(-3) 
        Y(-2) = dcmplx(cosp(2),-sinp(2)) * pre_and_post(-2) 
        Y(-1) = dcmplx(cosp(1),-sinp(1)) * pre_and_post(-1)
        Y(0) = pre_and_post(0)  
        Y(1) = dcmplx(cosp(1),sinp(1)) * pre_and_post(1) 
        Y(2) = dcmplx(cosp(2),sinp(2)) * pre_and_post(2)  
        Y(3) = dcmplx(cosp(3),sinp(3)) * pre_and_post(3) 
        Y(4) = dcmplx(cosp(4),sinp(4)) * pre_and_post(4)
        Y(5) = dcmplx(cosp(5),sinp(5)) * pre_and_post(5) 
        Y(6) = dcmplx(cosp(6),sinp(6)) * pre_and_post(6) 
        Y(7) = dcmplx(cosp(7),sinp(7)) * pre_and_post(7) 
        Y(8) = dcmplx(cosp(8),sinp(8)) * pre_and_post(8) 
        
        end subroutine 
 
    pure subroutine optimized_Y6(sin_theta, cos_theta, cos_phi, sin_phi, Y)
        implicit none 
        ! Passed
        real(dp), intent(in) :: sin_theta, cos_theta, cos_phi, sin_phi 

        ! Local
        real(dp), dimension(1:6) :: sint, cost, sinp, cosp 
        real(dp), dimension(-6:6) :: pre_and_post

        ! Return
        complex(dp), intent(out), dimension(-6:6) :: Y

        Y = 0.0d0       

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

        ! calculate powers of cos(theta)
        cost(2) = cost(1)*cost(1)
        cost(3) = cost(2)*cost(1)
        cost(4) = cost(3)*cost(1)
        cost(5) = cost(4)*cost(1)
        cost(6) = cost(5)*cost(1)


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

        ! Update the harmonic sum
        ! From Mathmatica function: "SphericalHarmonicsY[4,m,theta, phi]"	
    
        ! pre and post factor: 

        pre_and_post(-6) = 0.483084113580066d0 * sint(6)
        pre_and_post(-5) = 1.67345245810010d0 * sint(5) * cost(1)
        pre_and_post(-4) = 0.356781262853998d0 * (11.0d0 * cost(2)- 1.0d0) * sint(4)
        pre_and_post(-3) = 0.651390485867716d0 * (-3.0d0*cost(1) + 11.0d0 * cost(3)) * sint(3)
        pre_and_post(-2) = 0.325695242933858d0 * (1.0d0 - 18.0d0 * cost(2) + 33.0d0* cost(4)) * sint(2)
        pre_and_post(-1) = 0.411975516301141d0 * (5.0d0*cost(1) - 30.0d0 * cost(3) + 33.0d0*cost(5)) * sint(1)
        pre_and_post(0) = 0.06356920226762842d0 * (105.0d0 * cost(2) - 315.0d0 * cost(4) + 231.0d0 * cost(6) - 5.0d0)
        pre_and_post(1) = -pre_and_post(-1)
        pre_and_post(2) = pre_and_post(-2)
        pre_and_post(3) = -pre_and_post(-3)
        pre_and_post(4) = pre_and_post(-4)
        pre_and_post(5) = -pre_and_post(-5)
        pre_and_post(6) = pre_and_post(-6)


        Y(-6) = dcmplx(cosp(6),-sinp(6))*pre_and_post(-6)  
        Y(-5) = dcmplx(cosp(5),-sinp(5))*pre_and_post(-5) 
        Y(-4) = dcmplx(cosp(4),-sinp(4))*pre_and_post(-4)  
        Y(-3) = dcmplx(cosp(3),-sinp(3))*pre_and_post(-3) 
        Y(-2) = dcmplx(cosp(2),-sinp(2))*pre_and_post(-2) 
        Y(-1) = dcmplx(cosp(1),-sinp(1))*pre_and_post(-1)
        Y(0) = pre_and_post(0) 
        Y(1) = dcmplx(cosp(1),sinp(1))* pre_and_post(1) 
        Y(2) = dcmplx(cosp(2),sinp(2))* pre_and_post(2)  
        Y(3) = dcmplx(cosp(3),sinp(3))* pre_and_post(3) 
        Y(4) = dcmplx(cosp(4),sinp(4))* pre_and_post(4)
        Y(5) = dcmplx(cosp(5),sinp(5))* pre_and_post(5) 
        Y(6) = dcmplx(cosp(6),sinp(6))* pre_and_post(6) 

        end subroutine 

    pure subroutine optimized_Y4(sin_theta, cos_theta, cos_phi, sin_phi, Ylm)
            implicit none  
            ! Passed
            real(dp), intent(in) :: sin_theta, cos_theta, cos_phi, sin_phi 

            ! Local
            real(dp), dimension(1:4) :: sint, cost, sinp, cosp 
            real(dp), dimension(-4:4) :: pre_and_post

            ! Return
            complex(dp), intent(out), dimension(-4:4) :: Ylm
       
            ! assign the inital parameters 
            sint(1) = sin_theta

            cost(1) = cos_theta

            sinp(1) = sin_phi

            cosp(1) = cos_phi

            ! calculate powers of sin(theta)
            sint(2) = sint(1)*sint(1)
            sint(3) = sint(2)*sint(1)
            sint(4) = sint(3)*sint(1) 

            ! calculate powers of cos(theta)
            cost(2) = cost(1)*cost(1)
            cost(3) = cost(2)*cost(1)
            cost(4) = cost(3)*cost(1)


            ! calculate sin/cos (2phi)
            sinp(2) = 2.0d0*sinp(1)*cosp(1)
            cosp(2) = 2.0d0*cosp(1)**2 - 1.0d0

            ! calculate sin/cos (3phi)
            sinp(3) = sinp(2)*cosp(1) + cosp(2)*sinp(1)
            cosp(3) = cosp(2)*cosp(1) - sinp(2)*sinp(1)


            ! calculate sin/cos (4phi) 

            sinp(4) = sinp(3)*cosp(1) + cosp(3)*sinp(1) 
            cosp(4) = cosp(3)*cosp(1) - sinp(3)*sinp(1) 

		    ! Update the harmonic sum
	        ! From Mathmatica function: "SphericalHarmonicsY[4,m,theta, phi]"	

		    pre_and_post(-4) = 0.442532692444983d0 *sint(4)
		    pre_and_post(-3) = 1.25167147089835d0 * cost(1)*(sint(3)) 
		    pre_and_post(-2) = 0.334523271778645d0 * (-1 + 7.0d0*cost(2))*sint(2) 
		    pre_and_post(-1) = 0.473087347878780d0 * (-3.0d0 * cost(1) + 7.0d0*cost(3))*sint(1) 
		    pre_and_post(0) =  0.105785546915204d0 * (3.0d0 -30.0d0*cost(2) + 35.0d0*cost(4))
		    pre_and_post(1) = -pre_and_post(-1) 
		    pre_and_post(2) = pre_and_post(-2)  
		    pre_and_post(3) = -pre_and_post(-3) 
		    pre_and_post(4) = pre_and_post(-4) 

		    Ylm(-4) = dCMPLX(cosp(4),-sinp(4)) * pre_and_post(-4)  
		    Ylm(-3) = dCMPLX(cosp(3),-sinp(3))* pre_and_post(-3) 
		    Ylm(-2) = dCMPLX(cosp(2),-sinp(2))* pre_and_post(-2)
		    Ylm(-1) = dCMPLX(cosp(1),-sinp(1))* pre_and_post(-1) 
		    Ylm(0) = pre_and_post(0) 
		    Ylm(1) = dCMPLX(cosp(1),sinp(1))* pre_and_post(1)
		    Ylm(2) = dCMPLX(cosp(2),sinp(2)) * pre_and_post(2) 
		    Ylm(3) = dCMPLX(cosp(3),sinp(3)) * pre_and_post(3)  
		    Ylm(4) = dCMPLX(cosp(4),sinp(4)) * pre_and_post(4) 
        
            end subroutine

    ! This subroutine is exposed as c callable 
    pure subroutine call_optimized_Y12(sin_theta, cos_theta, cos_phi, sin_phi, Ylm_complex) bind(c, name="call_optimized_Y12")
            implicit none 
            ! Passed
            real(c_double), intent(in) :: sin_theta, cos_theta, cos_phi, sin_phi 

            ! Local
            real(dp),dimension(1:12) :: sint, cost, sinp, cosp 
            complex(dp), dimension(-12:12) :: Ylm 

            ! Return
            real(c_double), intent(out), dimension(-12:12,2) :: Ylm_complex

            call optimized_Y12(sin_theta, cos_theta, cos_phi, sin_phi, Ylm)
                   
            Ylm_complex(:,1) = real(Ylm, dp) 

            Ylm_complex(:,2) = dimag(Ylm)

            end subroutine 

    ! This subroutine is exposed as c callable 
    pure subroutine call_optimized_Y8(sin_theta, cos_theta, cos_phi, sin_phi, Ylm_complex) bind(c, name="call_optimized_Y8")
            implicit none 
            ! Passed
            real(c_double), intent(in) :: sin_theta, cos_theta, cos_phi, sin_phi 

            ! Local
            real(dp),dimension(1:8) :: sint, cost, sinp, cosp 
            complex(dp), dimension(-8:8) :: Ylm 

            ! Return
            real(c_double), intent(out), dimension(-8:8,2) :: Ylm_complex

            call optimized_Y8(sin_theta, cos_theta, cos_phi, sin_phi, Ylm)
                   
            Ylm_complex(:,1) = real(Ylm, dp) 

            Ylm_complex(:,2) = dimag(Ylm)

            end subroutine 

    ! This subroutine is exposed as c callable 
    pure subroutine call_optimized_Y6(sin_theta, cos_theta, cos_phi, sin_phi, Ylm_complex) bind(c, name="call_optimized_Y6")
            implicit none 
            ! Passed
            real(c_double), intent(in) :: sin_theta, cos_theta, cos_phi, sin_phi 

            ! Local
            real(dp),dimension(1:6) :: sint, cost, sinp, cosp 
            complex(dp), dimension(-6:6) :: Ylm 

            ! Return
            real(c_double), intent(out), dimension(-6:6,2) :: Ylm_complex

            call optimized_Y6(sin_theta, cos_theta, cos_phi, sin_phi, Ylm)
                   
            Ylm_complex(:,1) = real(Ylm, dp) 

            Ylm_complex(:,2) = dimag(Ylm)

            end subroutine 

    ! This subroutine is exposed as c callable 
    pure subroutine call_optimized_Y4(sin_theta, cos_theta, cos_phi, sin_phi, Ylm_complex) bind(c, name="call_optimized_Y4")
            implicit none 
            ! Passed
            real(c_double), intent(in) :: sin_theta, cos_theta, cos_phi, sin_phi 

            ! Local
            real(dp),dimension(1:4) :: sint, cost, sinp, cosp 
            complex(dp), dimension(-4:4) :: Ylm 

            ! Return
            real(c_double), intent(out), dimension(-4:4,2) :: Ylm_complex

            call optimized_Y4(sin_theta, cos_theta, cos_phi, sin_phi, Ylm)
                   
            Ylm_complex(:,1) = real(Ylm, dp) 

            Ylm_complex(:,2) = dimag(Ylm)

            end subroutine 
end module 
