!----------------------------------------------- Program Descriptions ----------------------------------------------- 
! This modules provide subroutines to perform trial translate 
! Date composed by Jingxiang Guo : 04/15/2020 
! The subroutine/function names with "call_xxx" can be called from Python or C 

module lj_potential
	use system
	implicit none 
	! all variables and subroutines/functions are private  
	private 

	! Export these subroutines/functions as public ( other subroutines/functions are still private) : 	
	public :: compute_isite_LJ_cut

	! Global varialbes: 
	
contains

!----------------------------------------------------------------------------------------------------------------- 
!                                               Interfacing with C                                                 
!----------------------------------------------------------------------------------------------------------------- 
	
    subroutine compute_isite_LJ_cut(&
                            & total_atoms,& 
                            & xyz,& 
                            & box,& 
                            & icoord,&
                            & cutoff,& 
                            & potential_energy) &
                            & bind(c,name="compute_isite_LJ_cut")
        implicit none 

        !Passed: 
        integer(c_int),intent(in) :: total_atoms  
        real(c_double),intent(in) :: cutoff 
        real(c_double),intent(in),dimension(1:3,total_atoms) :: xyz  
        real(c_double),intent(in),dimension(1:3) :: box,icoord  

        !Local: 
        integer :: iatom 
        real(dp) :: cutoff_sqr, sum_diff_sqr, twelves_terms, six_terms  
        real(dp),dimension(1:3) :: r_sep  

        !Output: 
        real(c_double),intent(out) :: potential_energy   

        potential_energy = 0.0d0 

        cutoff_sqr = cutoff*cutoff 

        do iatom = 1,total_atoms 

            r_sep =  xyz(:,iatom) - icoord  

            ! periodical boundary condition and minimum image convention 
            r_sep = r_sep - box*dnint(r_sep/box)  

            sum_diff_sqr = r_sep(1)*r_sep(1) + r_sep(2)*r_sep(2) + r_sep(3)*r_sep(3)  

            ! sum_diff_sqr is 0 for a self term, the potential energy is 0 (no
            ! impact on the total potential)  
            if (sum_diff_sqr < cutoff_sqr) then
    
                inv_sqr = 1/sum_diff_sqr
        
                six_terms = inv_sqr*sum_diff_sqr*sum_diff_sqr  

                twelves_terms = six_terms*six_terms 

                potential_energy = potential_energy + (twelves_terms - six_terms)  

            end if 

        end do 
        
        potential_energy = 4 * potential_energy         
  
        end subroutine 


end module 
