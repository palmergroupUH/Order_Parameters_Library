!----------------------------------------------- Program Descriptions ----------------------------------------------- 

! This program contains programming environment information used by current compiler 

! The variables defined could be used by all modules in the library 

! Date composed by Jingxiang Guo : 11/27/2019 

module system 
	use iso_c_binding ,only: c_int,c_double,c_char,c_float,c_short,c_long,c_ptr 
	use iso_fortran_env,only : input_unit,output_unit,error_unit
	implicit none 
	! all variables and subroutines/functions are private  
	private 
	
	! Placeholders variables (separating input and output in subroutine argument):
	integer :: O__O

	! To make code portable, select appropriate precision:
	integer,parameter :: sp = selected_real_kind(6, 37) ! Use Single Precision
	integer,parameter :: dp = selected_real_kind(15,307) ! Use Double Precision
	integer,parameter :: qp = selected_real_kind(33, 4931) ! Use Quadruple-Precision

	! Current precision 
	integer,parameter :: cp=dp ! ( double precision for calculations)  

	! all command_line_argument ( at least 200 bytes )  
	integer :: num_arg 
	character(len=200),dimension(:),allocatable :: arg 

	! Export these environmental variables: 
	public :: input_unit,output_unit,error_unit, & 
			  ! ctypes variable 
			  & c_int,c_double,c_char,c_float,c_short,c_long,c_ptr,&
			  ! precision kind parameters: 
			  & sp,dp,qp,cp,& 
			  ! command line variables:
			  & num_arg,arg,& 
			  & use_command_line,& 
			  ! placeholders variables 
			  & O__O,& 
              & convert_c_string_f_string 
		
			

contains
		
    subroutine use_command_line()
        implicit none
        integer :: i

        num_arg = command_argument_count()

        allocate(arg(1:num_arg))

        do i = 1,num_arg

            call get_command_argument(i,arg(i))

        end  do

        end  subroutine

    subroutine convert_c_string_f_string(str,strlen,f_string) 
        implicit none 
        integer ,intent(in):: strlen 
        character(len=1),intent(in),dimension(1:strlen) :: str 
        character(len=:),allocatable,intent(out) :: f_string         
        integer :: i 

        f_string = " "

        do i = 1,strlen

            f_string = f_string//str(i) 
           
        end do 

        end subroutine 

end module 
