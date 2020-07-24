!----------------------------------------------- Program Descriptions ----------------------------------------------- 
! This modules provide subroutines to perform trial translate 
! Date composed by Jingxiang Guo : 04/15/2020 
! The subroutine/function names with "call_xxx" can be called from Python or C 

module lj_pot
	use system
	implicit none 
	! all variables and subroutines/functions are private  
	private 
    
    real(dp),parameter :: pi = 3.14159265358979d0


	! Export these subroutines/functions as public ( other subroutines/functions are still private) : 	
	public :: compute_isite_LJ_cut,& 
            & compute_lj_cut_properties,&
            & correct_lj_cut_lrc_pe,&
            & correct_lj_cut_press,& 
            & correct_lj_cut_lrc_press 

	! Global varialbes: 
	
contains

!----------------------------------------------------------------------------------------------------------------- 
!                                               Interfacing with C                                                 
!----------------------------------------------------------------------------------------------------------------- 
	
    subroutine compute_isite_LJ_cut(&
                            & sigma,& 
                            & epslion,& 
                            & total_atoms,& 
                            & xyz,& 
                            & box,& 
                            & icoord,& 
                            & iselect,&
                            & cutoff,& 
                            & potential_energy,&
                            & sum_virial) &
                            & bind(c,name="compute_isite_LJ_cut")
        implicit none 

        !Passed: 
        integer(c_int),intent(in) :: iselect,total_atoms  
        real(c_double),intent(in) :: cutoff,sigma,epslion  
        real(c_double),intent(in),dimension(1:3,total_atoms) :: xyz  
        real(c_double),intent(in),dimension(1:3) :: box,icoord 

        !Local: 
        integer :: iatom 
        real(dp) :: cutoff_sqr, sum_diff_sqr, & 
                  & twelves_terms, six_terms, & 
                  & inv_sqr, sigma2 
        real(dp),dimension(1:3) :: r_sep,fij 

        !Output: 
        real(c_double),intent(out) :: potential_energy,sum_virial 
        
        potential_energy = 0.0d0 

        cutoff_sqr = cutoff*cutoff 

        sigma2 = sigma*sigma

        sum_virial = 0 

        do iatom = 1,total_atoms 

            if (iatom == iselect) cycle 

            r_sep =  xyz(:,iatom) - icoord  

            ! periodical boundary condition and minimum image convention 
            r_sep = r_sep - box*dnint(r_sep/box)  

            sum_diff_sqr = r_sep(1)*r_sep(1) + r_sep(2)*r_sep(2) + r_sep(3)*r_sep(3)  
            
            ! sum_diff_sqr is 0 for a self term, the potential energy is 0 (no
            ! impact on the total potential)  
            if (sum_diff_sqr < cutoff_sqr) then            

                inv_sqr = sigma2/sum_diff_sqr 
        
                six_terms = inv_sqr*inv_sqr*inv_sqr  

                twelves_terms = six_terms*six_terms 

                fij = r_sep*inv_sqr*(twelves_terms - 0.5d0*six_terms) 

                sum_virial = sum_virial + dot_product(fij,r_sep) 
    
                potential_energy = potential_energy + (twelves_terms - six_terms)  

            end if 

        end do 
        
        potential_energy = 4 * epslion * potential_energy
  
        end subroutine

    subroutine compute_lj_cut_properties(& 
                                        & sigma,& 
                                        & epslion,& 
                                        & total_atoms,& 
                                        & xyz,& 
                                        & box,& 
                                        & cutoff,& 
                                        & total_pe,&
                                        & sum_virial) &
                                        & bind(c,name="compute_lj_cut_properties")

        implicit none 

        !Passed: 
        integer(c_int),intent(in) :: total_atoms  
        real(c_double),intent(in) :: cutoff,sigma,epslion  
        real(c_double),intent(in),dimension(1:3,total_atoms) :: xyz  
        real(c_double),intent(in),dimension(1:3) :: box

        !Local: 
        integer :: iatom,jatom  
        real(dp),dimension(1:3) :: xyz_i 
        real(dp) :: cutoff_sqr, sum_diff_sqr, & 
                  & twelves_terms, six_terms, & 
                  & inv_sqr, sigma2 
        real(dp),dimension(1:3) :: r_sep,fij

        !Output: 
        real(c_double),intent(out) :: total_pe,sum_virial 

        total_pe = 0.0d0 

        cutoff_sqr = cutoff*cutoff 
    
        sigma2 = sigma*sigma

        sum_virial = 0 

        do iatom = 1,total_atoms-1  

            xyz_i = xyz(:,iatom) 

            do jatom = iatom+1,total_atoms 

                r_sep =  xyz(:,jatom) - xyz_i  

                ! periodical boundary condition and minimum image convention 
                r_sep = r_sep - box*dnint(r_sep/box)  

                sum_diff_sqr = r_sep(1)*r_sep(1) + r_sep(2)*r_sep(2) + r_sep(3)*r_sep(3)  

                ! sum_diff_sqr is 0 for a self term, the potential energy is 0 (no
                ! impact on the total potential)

                if (sum_diff_sqr < cutoff_sqr) then

                    inv_sqr = sigma2/sum_diff_sqr 
            
                    six_terms = inv_sqr*inv_sqr*inv_sqr 

                    twelves_terms = six_terms*six_terms 

                    fij = r_sep*inv_sqr*(twelves_terms - 0.5d0*six_terms)

                    sum_virial = sum_virial + dot_product(fij,r_sep)

                    total_pe = total_pe + (twelves_terms - six_terms)

                end if 

            end do 

        end do 
            
        total_pe = 4 * epslion * total_pe         

        sum_virial = 48 * sum_virial 
      
        end subroutine 

    ! "Frenkel & Smit": eq 3.2.5
    subroutine correct_lj_cut_lrc_pe(sigma, epslion, total_atoms, volume, r_cut, u_tail) bind(c,name="compute_lj_lrc")
        implicit none 
        
        ! Passed: 
        integer(c_int),intent(in) :: total_atoms  
        real(c_double),intent(in) :: r_cut,sigma,epslion,volume  
        
        ! Local:
        real(dp) :: power3_term,power9_term
           
        ! Output:
        real(c_double),intent(out) :: u_tail

        power3_term = (sigma/r_cut)*(sigma/r_cut)*(sigma/r_cut)
        
        power9_term = (1.0d0/3.0d0)*power3_term*power3_term*power3_term

        u_tail =8.0d0/3.0d0*pi*(total_atoms/volume)*epslion*sigma*sigma*sigma*(power9_term-power3_term)

        end subroutine    

    subroutine correct_lj_cut_press(sigma, epslion, total_atoms, volume, r_cut, u_tail) bind(c,name="correct_lj_cut_press")
        implicit none 
        
        ! Passed: 
        integer(c_int),intent(in) :: total_atoms  
        real(c_double),intent(in) :: r_cut,sigma,epslion,volume  
        
        ! Local:
        real(dp) :: power3_term,power9_term
           
        ! Output:
        real(c_double),intent(out) :: u_tail

        power3_term = (sigma/r_cut)*(sigma/r_cut)*(sigma/r_cut)
        
        power9_term = (1.0d0/3.0d0)*power3_term*power3_term*power3_term

        u_tail =8.0d0/3.0d0*pi*(total_atoms/volume)*epslion*sigma*sigma*sigma*(power9_term-power3_term)

        end subroutine
    
    subroutine correct_lj_cut_lrc_press(sigma, epslion, total_atoms, volume, r_cut, u_tail) bind(c,name="correct_lj_cut_press")
        implicit none 
        
        ! Passed: 
        integer(c_int),intent(in) :: total_atoms  
        real(c_double),intent(in) :: r_cut,sigma,epslion,volume  
        
        ! Local:
        real(dp) :: power3_term,power9_term
           
        ! Output:
        real(c_double),intent(out) :: u_tail

        power3_term = (sigma/r_cut)*(sigma/r_cut)*(sigma/r_cut)
        
        power9_term = (1.0d0/3.0d0)*power3_term*power3_term*power3_term

        u_tail =8.0d0/3.0d0*pi*(total_atoms/volume)*epslion*sigma*sigma*sigma*(power9_term-power3_term)

        end subroutine


end module 
