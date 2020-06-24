module xyz_reader
    use system
    implicit none 
    private 

    public :: call_read_xyz_header,&
            & call_read_xyz_with_box,& 
            & call_read_xyz_no_box,&
            & call_read_xyz_no_box_chunk,& 
            & call_read_xyz_with_box_chunk

contains 

    subroutine call_read_xyz_header(xyzfile,& 
                             & strleng, & 
                             & total_atoms,&  
                             & total_frames) & 
                             & bind(c,name="call_read_xyz_header") 
        implicit none 

        ! Passed: 
        integer(c_int),intent(in) :: strleng 
        character(kind=c_char,len=1),intent(in),dimension(1:strleng) :: xyzfile 

        ! Local: 
        character(len=:),allocatable :: filename 
        integer :: unit_number,num_lines,IOstatus 

        ! Output:
        integer(c_int),intent(out) ::  total_atoms 
        integer(c_int),intent(out) ::  total_frames   

        call convert_c_string_f_string(xyzfile,strleng,filename) 
       
        open(newunit=unit_number,file=filename,action="read",form="formatted") 


            read(unit_number,*) total_atoms  

            num_lines = 1 

            do 

                read(unit_number,'(A)',IOSTAT=IOstatus) 

                if (IOstatus /= 0) then  

                    exit 

                end if 

                num_lines = num_lines + 1 

            end do 
       
            total_frames = num_lines/(total_atoms + 2)  
                
        close(unit_number) 

        end subroutine 

    subroutine call_read_xyz_with_box(&
                                  & xyzfile,& 
                                  & strleng, & 
                                  & current_frame,& 
                                  & total_atoms, & 
                                  & box, & 
                                  & xyz) &
                                  & bind(c,name="call_read_xyz_with_box") 
        implicit none 

        ! Passed: 
        integer(c_int),intent(in) :: strleng,total_atoms, current_frame 
        character(kind=c_char,len=1),intent(in),dimension(1:strleng) :: xyzfile 

        ! Local: 
        character(len=:),allocatable :: filename
        integer :: i,iatom ,IOstatus, unit_number,atomID
        character(len=10) :: atom_type 

        ! Output:
        real(c_double),intent(out),dimension(1:3) :: box   
        real(c_double),intent(out),dimension(1:3,total_atoms) :: xyz 

        call convert_c_string_f_string(xyzfile,strleng,filename) 
       
        open(newunit=unit_number,file=filename,action="read",form="formatted") 

            do i = 1, current_frame - 1

                read(unit_number,*)   
                read(unit_number,*) 

                do iatom = 1,total_atoms 

                    read(unit_number,*) 

                end do 

            end do 
           
            read(unit_number,*)   
            read(unit_number,*) box  
            
            do iatom = 1,total_atoms 

                read(unit_number,*) atom_type, xyz(1:3,iatom) 
                
            end do 

        close(unit_number) 
        
        end subroutine  

    subroutine call_read_xyz_no_box(&
                          & xyzfile,& 
                          & strleng,& 
                          & current_frame,& 
                          & total_atoms,& 
                          & xyz) &
                          & bind(c,name="call_read_xyz_no_box") 
        implicit none 

        ! Passed: 
        integer(c_int),intent(in) :: strleng,total_atoms,current_frame 
        character(kind=c_char,len=1),intent(in),dimension(1:strleng) :: xyzfile 

        ! Local: 
        character(len=:),allocatable :: filename 
        integer :: i,iatom,IOstatus, unit_number,atomID
        character(len=10) :: atom_type 

        ! Output:
        real(c_double),intent(out),dimension(1:3,total_atoms) :: xyz 

        call convert_c_string_f_string(xyzfile,strleng,filename) 
            
        open(newunit=unit_number,file=filename,action="read",form="formatted") 
            
            do i = 1, current_frame - 1

                read(unit_number,*)   
                read(unit_number,*) 

                do iatom = 1,total_atoms 

                    read(unit_number,*) 

                end do 

            end do 
           
            
            read(unit_number,*)   
            read(unit_number,*)   
            
            do iatom = 1,total_atoms 

                read(unit_number,*) atom_type, xyz(1:3,iatom) 
                
            end do 

        close(unit_number) 
        
        end subroutine  

    subroutine call_read_xyz_no_box_chunk(xyzfile,& 
                                    & strleng,& 
                                    & total_atoms,&
                                    & start,&
                                    & num_frames,& 
                                    & xyz)& 
                                    & bind(c,name="call_read_xyz_no_box_chunk") 
        implicit none 

        ! Passed: 
        integer(c_int),intent(in) :: strleng,total_atoms,start,num_frames 
        character(kind=c_char,len=1),intent(in),dimension(1:strleng) :: xyzfile 

        ! Local: 
        character(len=:),allocatable :: filename 
        integer :: i,iatom,IOstatus, unit_number,atomID,iframe 
        character(len=10) :: atom_type 

        ! Output:
        real(c_double),intent(out),dimension(1:3,total_atoms,num_frames) :: xyz 

        call convert_c_string_f_string(xyzfile,strleng,filename) 
        
        open(newunit=unit_number,file=filename,action="read",form="formatted") 

            do i = 1, start - 1

                read(unit_number,*)   

                read(unit_number,*) 

                do iatom = 1,total_atoms 

                    read(unit_number,*) 

                end do 

            end do 
           
            do iframe = start, start + num_frames -1  

                read(unit_number,*)   

                read(unit_number,*)   

                do i = 1,total_atoms 

                    read(unit_number,*) atom_type, xyz(1:3,i,iframe) 
                    
                end do 

            end do 

        close(unit_number) 
        
        end subroutine  

    subroutine call_read_xyz_with_box_chunk(xyzfile,& 
                                    & strleng,& 
                                    & total_atoms,&
                                    & start,&
                                    & num_frames,& 
                                    & box, & 
                                    & xyz)& 
                                    & bind(c,name="call_read_xyz_with_box_chunk") 
        implicit none 

        ! Passed: 
        integer(c_int),intent(in) :: strleng,total_atoms,start,num_frames 
        character(kind=c_char,len=1),intent(in),dimension(1:strleng) :: xyzfile 

        ! Local: 
        character(len=:),allocatable :: filename 
        integer :: i,iatom,IOstatus, unit_number,atomID,iframe 
        character(len=10) :: atom_type 

        ! Output:
        
        real(c_double),intent(out),dimension(1:3,num_frames) :: box 
        real(c_double),intent(out),dimension(1:3,total_atoms,num_frames) :: xyz 

        call convert_c_string_f_string(xyzfile,strleng,filename) 
       
        open(newunit=unit_number,file=filename,action="read",form="formatted") 

            do i = 1, start - 1

                read(unit_number,*)   

                read(unit_number,*) 

                do iatom = 1,total_atoms 

                    read(unit_number,*) 

                end do 

            end do 
           
            do iframe = 1, num_frames 

                read(unit_number,*)   

                read(unit_number,*) box(:,iframe)  

                do iatom = 1,total_atoms 

                    read(unit_number,*) atom_type, xyz(1:3,iatom,iframe) 
                    
                end do 

            end do 

        close(unit_number) 
        
        end subroutine  

end module 


