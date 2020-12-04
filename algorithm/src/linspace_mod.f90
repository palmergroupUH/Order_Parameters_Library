module matlab_linspace
    use system 
    implicit none 
    private 


    public :: linspace 

! overload the generic function overloading

interface linspace 

	module procedure :: linspace_int, linspace_sp, linspace_dp 

end interface 

contains 

	pure function linspace_int(first,last,num) 
		implicit none 
        ! Passed
		integer,intent(in) :: num 
		integer,intent(in) :: first,last 
        ! Local
		real(dp) :: interval 
		integer :: i 
		real(dp),dimension(:),allocatable :: linspace_int 	

		allocate(linspace_int(1:num)) 

		if ( mod(last-first,num) == 0 ) then 
		
			do i = 1,num-1 

				linspace_int = first + (i-1) 	

			end do 
				
		else 

			interval = (last-first)/dble((num-1))  
		
			DO i = 1,num
				
				linspace_int(i) = first + interval*(i-1) 

			end DO 	

		end if 
		
		end function  

	pure function linspace_sp(first,last,num ) 
		implicit none 
		integer,intent(in) :: num 
		real(sp),intent(in) :: first,last 
		real(dp) :: interval 
		integer :: i
		real(dp),dimension(:),allocatable :: linspace_sp  	

		interval = (last-first)/(num-1) 

		allocate(linspace_sp(1:num)) 
		
		DO i = 1,num
			
			linspace_sp(i) = first + interval*(i-1) 

		end DO 	
		
		end function 

	pure function linspace_dp(first,last,num ) 
		implicit none 
		integer,intent(in) :: num 
		real(dp),intent(in) :: first,last 
		real(dp) :: interval 
		integer :: i
		real(dp),dimension(:),allocatable :: linspace_dp  	

		interval = (last-first)/(num-1) 

		allocate(linspace_dp(1:num)) 
		
		do i = 1,num
			
			linspace_dp(i) = first + interval*(i-1) 

		end do 	
		
		end function 


end module 
