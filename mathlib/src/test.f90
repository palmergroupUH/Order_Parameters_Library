module function_test
    use iso_c_binding
    implicit none 
    private

    public :: call_fortran
    

contains 

    function call_fortran(a) result(gg) bind(c, name="call_fortran_like_c")
        implicit none
        real(c_double), intent(in) :: a
        real(c_double) :: gg 
        integer :: i

        gg = a * 1.5d0
        print*, gg

        end function 

end module 

