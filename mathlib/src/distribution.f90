module histogram
    use system, only: c_int, c_double, dp, sp   
    implicit none 

    private 

    public :: init_hist, combine_hist,array_into_hist, normalize_hist, get_bin_id 

contains 

    subroutine init_hist(lower, upper, num_bins, hist, interval, r_mid) bind(c, name="call_init_hist")
        implicit none 
    
        ! Passed
        real(c_double), intent(in) :: lower, upper
        integer(c_int), intent(in) :: num_bins

        ! Local
        integer :: i  

        ! Return
        real(c_double), intent(out), dimension(1:num_bins) :: hist, r_mid  
        real(c_double), intent(out) :: interval  
        
        interval = (upper - lower) / num_bins 

        hist = 0.0d0 

        r_mid = 0.0d0 

        do i = 1, num_bins

            r_mid(i) = interval * (i-1) + 0.5d0 * interval + lower 

        end do 

        end subroutine 

    subroutine combine_hist(num_bins, new_hist, hist) bind(c, name="call_combine_hist")
        implicit none 

        ! Passed
        integer(c_int), intent(in) :: num_bins
        real(c_double), intent(in), dimension(1:num_bins) :: new_hist 

        ! Return 
        real(c_double), intent(inout), dimension(1:num_bins) :: hist     
       
        hist = hist + new_hist  

        end subroutine 

    subroutine array_into_hist(N, array, interval, upper, lower, num_bins, hist) bind(c, name="call_array_into_hist")
        implicit none 

        ! Passed
        integer(c_int), intent(in) :: N 
        integer(c_int), intent(in) :: num_bins 
        real(c_double), intent(in) :: interval, upper, lower 
        real(c_double), intent(in), dimension(1:N) :: array

        ! Local
        integer :: i, bin_index 

        ! Return
        real(c_double), intent(inout), dimension(1:num_bins) :: hist
        
      
        ! check any values are out of bounds 
        if (any(array > upper) .or. any(array < lower)) then  

            stop "Out of bound error: some values are outside of histogram ranges"

        end if 
 
        do i = 1, N
    
            bin_index = int((array(i) - lower) / interval) + 1  

            hist(bin_index) = hist(bin_index) + 1.0d0

        end do 

        end subroutine 

    subroutine normalize_hist(num_bins, interval, hist) bind(c, name="call_normalize_hist") 
        implicit none
            
        ! Passed
        integer(c_int), intent(in) :: num_bins
        real(c_double), intent(in) :: interval 

        ! Return
        real(c_double), intent(inout), dimension(1:num_bins) :: hist         

        ! normalize by its account

        hist = hist / (sum(hist) * interval)

        end subroutine 

    pure integer function get_bin_id(ranges, bin_interval, pos) 
        implicit none 
        real(dp), intent(in) :: pos, bin_interval 
        real(dp), intent(in), dimension(1:2) :: ranges
        real(dp) :: lower, upper     

        lower = ranges(1) 

        upper = ranges(2) 

        if (pos < upper .and. pos > lower) then  

            get_bin_id = int((pos-lower)/bin_interval) + 1

        else if (pos < ranges(1) ) then 

            get_bin_id = 0  

        else if (pos == ranges(2)) then  
        
            get_bin_id = int((pos-lower)/bin_interval)

        end if            
        
        end function 
end module 
