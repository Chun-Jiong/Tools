!-- the statistic variables -------------------------------------------------
MODULE stat_vrbls
	IMPLICIT NONE

	!-- constant ----------------------------------------------------------
	real(8), parameter          :: eps    = 1.d-14    ! very small number
	real(8), parameter          :: tol    = 0.2d0     ! tolerance for Cor

	!-- statistic parameters ----------------------------------------
	! THIS IS ALMOST PROJECT-INDEPENDENT 
	integer(4), parameter       :: MxBlck = 2**12     ! maximum number of blocks for statistics
	integer(4), parameter       :: MnBlck = 2**6      ! minimum number of blocks

	integer(4)                  :: NBlck              ! # blocks
	logical                     :: prt
	!-----------------------------------------------------------------


	! GO TO MODIFY PARAMETER LIST in initialize.f90
	integer(4)              :: NObs_b                 ! #basic     observables
	integer(4)              :: NObs_c                 ! #composite observables
	integer(4)              :: NObs                   ! = NObs_b+NObs_c
	!-----------------------------------------------------------------


	!-- Statistics ---------------------------------------------------
	! THIS IS PROJECT-DEPENDENT 
	real(8), allocatable    :: Quan(:)                ! Measured quantities
	real(8), allocatable    :: Obs(:,:)               ! 1st--#quan.  2nd--#block
	real(8), allocatable    :: Ave(:), Dev(:), Cor(:) ! average, error bars, and correlation of observables
	!-----------------------------------------------------------------

	!-- Output file --------------------------------------------------
	character(100), parameter :: output="out.dat"

END MODULE stat_vrbls
