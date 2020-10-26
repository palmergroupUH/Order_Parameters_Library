!----------------------------------------------- Program Descriptions ----------------------------------------------- 

! This module includes subroutines and functions to compute Wigner3j symobl: 

! The implementation can be found in [1]: 

! Note:
! Implementation of basics statisics function 

! Date composed by Jingxiang Guo : 8/28/2020 

module statistics
    use system 
    implicit none 
    private 

    public :: mean, stdv, std, linear_reg  


interface mean  

    module procedure :: mean_int, mean_sp, mean_dp

end interface 

interface stdv

    module procedure :: stdv_int, stdv_sp, stdv_dp 

end interface 

interface std

    module procedure :: std_int, std_sp, std_dp 

end interface 

interface linear_reg

    module procedure :: linear_reg_sp, linear_reg_dp

end interface 


contains

    pure real(dp) function mean_int(array, num_data) 
        implicit none 
        integer, intent(in) :: num_data 
        integer, intent(in), dimension(:) :: array 
        integer :: i

        mean_int = 0.0d0  

        do i = 1, num_data  

            mean_int = mean_int + array(i) 
              
        end do 

        mean_int = mean_int/num_data 
        
        end function 

    pure real(dp) function mean_sp(array, num_data) 
        implicit none 
        integer, intent(in) :: num_data 
        real(sp), intent(in), dimension(:) :: array 
        integer :: i

        mean_sp = 0.0d0  

        do i = 1, num_data  

            mean_sp = mean_sp + array(i)
              
        end do 
    
        mean_sp = mean_sp/num_data 

        end function 

    pure real(dp) function mean_dp(array, num_data) 
        implicit none 
        integer, intent(in) :: num_data 
        real(dp), intent(in), dimension(:) :: array 
        integer :: i

        mean_dp = 0.0d0  

        do i = 1, num_data  

            mean_dp = mean_dp + array(i) 
              
        end do 
    
        mean_dp = mean_dp/num_data 

        end function 

    pure real(dp) function stdv_int(array, num_data, dof) 
        implicit none 
        integer, intent(in) :: num_data, dof  
        integer, intent(in), dimension(:) :: array 
        integer :: element 
        real(dp) :: mean_val 
        integer :: i

        stdv_int = 0.0d0

        mean_val = mean(array, num_data)

        do i = 1, num_data

            element = array(i)

            stdv_int = stdv_int + (element - mean_val)*(element - mean_val)

        end do

        stdv_int = stdv_int/(num_data - dof)
 
        end function  

    pure real(dp) function stdv_sp(array, num_data, dof) 
        implicit none 
        integer, intent(in) :: num_data, dof  
        real(sp), intent(in), dimension(:) :: array 
        real(sp) :: element 

        real(dp) :: mean_val 
        integer :: i

        stdv_sp = 0.0d0

        mean_val = mean(array, num_data)

        do i = 1, num_data

            element = array(i)

            stdv_sp = stdv_sp + (element - mean_val)*(element - mean_val)

        end do

        stdv_sp = stdv_sp/(num_data - dof)
 
        end function  

    pure real(dp) function stdv_dp(array, num_data, dof) 
        implicit none 
        integer, intent(in) :: num_data, dof  
        real(dp), intent(in), dimension(:) :: array 
        real(dp) :: element 
        real(dp) :: mean_val 
        integer :: i

        stdv_dp = 0.0d0

        mean_val = mean(array, num_data)

        do i = 1, num_data

            element = array(i)

            stdv_dp = stdv_dp + (element - mean_val)*(element - mean_val)

        end do

        stdv_dp = stdv_dp/(num_data - dof)
 
        end function  

    pure real(dp) function std_int(array, num_data, dof) 
        implicit none 
        integer, intent(in) :: num_data, dof  
        integer, intent(in), dimension(:) :: array 

        std_int = dsqrt(stdv(array, num_data, dof))  

        end function  

    pure real(dp) function std_sp(array, num_data, dof) 
        implicit none 
        integer, intent(in) :: num_data, dof  
        real(sp), intent(in), dimension(:) :: array 

        std_sp = dsqrt(stdv(array, num_data, dof)) 
 
        end function  

    pure real(dp) function std_dp(array, num_data, dof)
        implicit none 
        integer, intent(in) :: num_data, dof
        real(dp), intent(in), dimension(:) :: array 

        std_dp = dsqrt(stdv(array, num_data, dof))  

        end function  

    pure subroutine linear_reg_dp(X_measure, Y_measure, num_data, r2, beta_1, beta_0, s_beta_1, s_beta_0)
        implicit none 
        ! Passed 
        integer, intent(in) :: num_data 
        real(dp), intent(in), dimension(:) :: X_measure, Y_measure 

        ! Local
        integer :: i 
        real(dp) :: x_bar, y_bar, diff_x, diff_y 
        real(dp) :: sum_x_sqr, sum_y_sqr, sum_corss 
        real(dp) :: s
        real(dp) :: r_corr 

        ! Return
        real(dp), intent(out) :: r2, beta_0, beta_1, s_beta_0, s_beta_1  

        x_bar = mean(X_measure, num_data)

        y_bar = mean(Y_measure, num_data)
     
        sum_x_sqr = 0.0d0 

        sum_corss = 0.0d0 

        sum_y_sqr = 0.0d0 
 
        do i = 1, num_data 

            diff_x = X_measure(i) - x_bar  

            diff_y = Y_measure(i) - y_bar 

            sum_x_sqr = sum_x_sqr + diff_x * diff_x  

            sum_y_sqr = sum_y_sqr + diff_y * diff_y 

            sum_corss = sum_corss + diff_x * diff_y

        end do 

        r_corr = sum_corss/(dsqrt(sum_x_sqr) * dsqrt(sum_y_sqr)) 
        
        ! slope of estimated coeff

        beta_1 = sum_corss / sum_x_sqr  

        ! intercept of the estimated coeff

        beta_0 = y_bar - x_bar * beta_1  

        ! error of standard deviation:

        s = dsqrt(((1 - r_corr*r_corr) * sum_y_sqr) / (num_data - 2) ) 
        
        ! stndard deviation of estimate beta_1 

        s_beta_1 = s / dsqrt(sum_x_sqr) 

        ! stndard deviation of estimate beta_0
        
        s_beta_0 = s * dsqrt(1.0d0/num_data + x_bar*x_bar /(sum_x_sqr))

        r2 = r_corr * r_corr 

        end subroutine  

    pure subroutine linear_reg_sp(X_measure, Y_measure, num_data, r2, beta_1, beta_0, s_beta_1, s_beta_0)
        implicit none 
        ! Passed
        integer, intent(in) :: num_data 
        real(sp), intent(in), dimension(:) :: X_measure, Y_measure 

        ! Local
        integer :: i 
        real(dp) :: x_bar, y_bar, diff_x, diff_y 
        real(dp) :: sum_x_sqr, sum_y_sqr, sum_corss 
        real(dp) :: s
        real(dp) :: r_corr

        ! Return
        real(dp), intent(out) :: r2, beta_0, beta_1, s_beta_0, s_beta_1  

        x_bar = mean(X_measure, num_data)

        y_bar = mean(Y_measure, num_data)
     
        sum_x_sqr = 0.0d0 

        sum_corss = 0.0d0 

        sum_y_sqr = 0.0d0 
 
        do i = 1, num_data 

            diff_x = X_measure(i) - x_bar  

            diff_y = Y_measure(i) - y_bar 

            sum_x_sqr = sum_x_sqr + diff_x * diff_x  

            sum_y_sqr = sum_y_sqr + diff_y * diff_y 

            sum_corss = sum_corss + diff_x * diff_y

        end do 

        r_corr = sum_corss / (dsqrt(sum_x_sqr) * dsqrt(sum_y_sqr)) 

        ! slope of estimated coeff

        beta_1 = sum_corss / sum_x_sqr  

        ! intercept of the estimated coeff

        beta_0 = y_bar - x_bar * beta_1  

        ! error of standard deviation:

        s = dsqrt(((1 - r_corr*r_corr) * sum_y_sqr) / num_data - 2 ) 

        ! stndard deviation of estimate beta_1 

        s_beta_1 = s / dsqrt(sum_x_sqr) 

        ! stndard deviation of estimate beta_0
        
        s_beta_0 = s * dsqrt(1.0d0/num_data + x_bar*x_bar /(sum_x_sqr))

        r2 = r_corr * r_corr 

        end subroutine 

    subroutine linear_regression_for_C(X_measure, & 
                                     & Y_measure, & 
                                     & num_data, & 
                                     & results) bind(c, name="call_linear_regression")

        implicit none 
        ! Passed 
        integer(c_int), intent(in) :: num_data 
        real(c_double), intent(in), dimension(1:num_data) :: X_measure, Y_measure 

        ! Local
        integer :: i 
        real(dp) :: x_bar, y_bar, diff_x, diff_y 
        real(dp) :: sum_x_sqr, sum_y_sqr, sum_corss 
        real(dp) :: s
        real(dp) :: r_corr 
        real(dp) :: r2, beta_0, beta_1, s_beta_0, s_beta_1

        ! Return
        real(c_double), intent(out), dimension(1:5) :: results

        call linear_reg_dp(X_measure, Y_measure, num_data, r2, beta_1, beta_0, s_beta_1, s_beta_0)

        results = [r2, beta_1, beta_0, s_beta_1, s_beta_0]

        end subroutine

end module 
