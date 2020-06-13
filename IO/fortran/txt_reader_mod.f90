module txt_reader
    use system
    implicit none 
    private 

    public :: get_txt_lines_columns,&
              & loadtxt 

contains 

    subroutine get_txt_lines_columns(txtfile,&
                                   & strleng,& 
                                   & num_lines,& 
                                   & num_columns) &
                                   & bind(c,name="get_txt_lines_columns") 
        implicit none 

        ! Passed: 
        integer(c_int),intent(in) :: strleng 
        character(kind=c_char,len=1),intent(in),dimension(1:strleng) :: txtfile 

        ! Local: 
        character(len=:),allocatable :: filename 
        integer :: i,IOstatus, unit_number
        integer,parameter :: max_columns=10000 
        character(len=10000) :: line 
        real(dp),dimension(1:max_columns) :: test_array 

        ! Output:
        integer(c_int),intent(out) ::  num_lines 
        integer(c_int),intent(out) :: num_columns

        call convert_c_string_f_string(txtfile,strleng,filename) 
       
        open(newunit=unit_number,file=filename,action="read",form="formatted") 

            num_lines = 0 

            do 
        
                read(unit_number,'(A)',IOSTAT=IOstatus) line  

                if ( IOstatus /= 0) then 

                    exit

                end if 

                num_lines = num_lines + 1 
                
            end do 

        !call determine_columns(filename,num_columns) 

        do i = 1, max_columns 

            read(line,*,iostat=IOstatus) test_array(1:i)  
            
            if (IOstatus == -1) then 

                exit

            end if 

        end do 

        num_columns = i -1  

        close(unit_number) 

        end subroutine 

    subroutine loadtxt(txtfile,&
                     & strleng,&
                     & num_lines,&
                     & num_col,& 
                     & skiprows,& 
                     & loaded_data)& 
                     & bind(c,name="load_txt")
        implicit none 

        ! Passed 
        integer(c_int),intent(in) :: strleng 
        character(kind=c_char,len=1),dimension(1:strleng),intent(in) :: txtfile
        integer(c_int),intent(in) :: num_lines,num_col,skiprows

        ! Local 
        character(len=:),allocatable :: filename 
        integer :: i,unit_number	

        ! Output: 
        real(c_double),intent(out),dimension(1:num_lines-skiprows,1:num_col) :: loaded_data

        call convert_c_string_f_string(txtfile,strleng,filename) 
        
        open(newunit=unit_number,file=filename,action="read",form="formatted") 

            do i = 1,skiprows
                
                read(unit_number,*)

            end do 

            do i = 1, num_lines - skiprows  

                read(unit_number,*) loaded_data(i,1:num_col) 
                
            end do 
        
        close(unit_number) 

        end subroutine 


end module 

