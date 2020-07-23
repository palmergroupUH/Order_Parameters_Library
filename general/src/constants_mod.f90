!----------------------------------------------- Program Descriptions ----------------------------------------------- 

! This program define many physical constants used by other modules 
! Date composed by Jingxiang Guo : 11/27/2019 

module constants
    ! import the precision (machine dependent) defined for each data type.
	use system,only: dp,sp 
	implicit none 
	! Global Variables:

	real(dp),parameter,public :: pi = 3.14159265358979323d0 
	real(dp),parameter,public :: Bohr = 0.529177d0 
	real(dp),parameter,public :: Mole = 6.02214179E+23 
	real(dp),parameter,public :: water_molar_mass = 18.01528 ! g/mol  
	
	! TIP4P:
	real(dp),parameter,public :: OM = 0.1577d0 ! Angstrom 
	real(dp),parameter,public :: boltzman = 1.38064852E-23 !   
	real(dp),parameter,public :: gas_constant = 8.3144598 
	integer,dimension(1),public :: seed=(/197394210/) 

	! Global varialbes:

contains
	


end module 
