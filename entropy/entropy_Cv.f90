#include "Modules/mdl_vrbls.f90"
#include "Modules/stat_vrbls.f90"
#include "Modules/rng_vrbls.f90"

PROGRAM main
	use mdl_vrbls
	use stat_vrbls
	use rng_vrbls
	implicit none

	integer(4)                  :: narg
	character(200)              :: filename
	character(100)              :: chartemp
	integer(4)                  :: method
	logical                     :: existyn
	integer(4)                  :: iblck
	integer(4)                  :: iobs
	real(8)                     :: beta, c
	real(8)                     :: entropy_infinity

	real(8), allocatable        :: binObs(:,:)
	real(8), allocatable        :: entropy(:,:)
	integer(4)                  :: bins

	integer(4)                  :: i

	narg = command_argument_count()

	select case (narg)
	case (2)
		call getarg(1,filename)
		call getarg(2,chartemp)
		read(chartemp,*) entropy_infinity
		method = 0
	case (4)
		call getarg(1,filename)
		call getarg(2,chartemp)
		read(chartemp,*) entropy_infinity
		call getarg(3,chartemp)
		read(chartemp,*) method
		call getarg(4,chartemp)
		read(chartemp,*) bins
	case default
		write(*,*) "Please enter the file name (and method, bins)."
		write(*,*) "OR if you want help, please use the parameter -h/--help"
		stop
	end select

	entropy_infinity = log(entropy_infinity)

	if( method/=0 .and. method/=1 ) then
		write(*,*) "method must be 0/1"
		stop
	end if

	inquire(file=trim(filename), exist=existyn)
	if( .not. existyn ) then
		write(*,*) "Err: No ", trim(filename)
		stop
	end if

	open(1,file=trim(filename), action="read")
	NBlck = 0
	read(1,*) NObs
	do
		read(1,*,end=10) c
		do iobs=1, NObs
			read(1,*) beta, c
		end do
		NBlck = NBlck+1
	end do
	10 close(1)

	allocate(Obs(NObs,NBlck))
	allocate(temperature(NObs))
	Obs = 0.d0

	open(1,file=trim(filename),action="read")
	read(1,*) NObs
	do iblck=1, NBlck
		read(1,*) c
		do iobs=1, NObs
			read(1,*) beta, c
			Obs(iobs,iblck) = c
			temperature(iobs) = 1.d0/beta
		end do
	end do
	close(1)
	!--- read data done --------------------------------------------------!
	
	if( method==1 ) then
		!write(*,'(A16,I8)') "Initial bins  = ", NBlck
		!write(*,'(A16,I8)') "Random  bins  = ", bins
		call date_and_time(date, time, zone, tval)
		Seed = tval(5)*3600+tval(6)*60+tval(7)
		call set_RNG

		allocate(binObs(NObs,bins))
		do iobs=1, NObs
			do iblck=1, bins
				i = int(rn()*NBlck)+1
				binObs(iobs,iblck) = Obs(iobs,i)
			end do
		end do
		NBlck = bins
		deallocate(Obs)
		allocate(Obs(NObs,NBlck))
		Obs = binObs
		deallocate(binObs)
	end if

	!--- calculate entropy -----------------------------------------------!
	allocate(entropy(NObs,NBlck))
	do iblck = 1, NBlck
		entropy(1,iblck) = entropy_infinity
		do i = 2, NObs
			entropy(i,iblck) = entropy(i-1,iblck) - &
				& (Obs(i,iblck)+Obs(i-1,iblck))/2.d0 /(temperature(i)+temperature(i-1))*2.d0 &
				& *abs(temperature(i)-temperature(i-1))
		end do 
	end do 
	Obs = entropy
	deallocate(entropy)
	!--- calculate entropy done ------------------------------------------!


	allocate(Ave(NObs))
	allocate(Dev(NObs))
	allocate(Cor(NObs))

	Ave = 0.d0
	Dev = 0.d0
	Cor = 0.d0

	call stat_analy(NBlck)
	call write2file(NObs)

	stop

CONTAINS
#include "Statistics/statistics.f90"
#include "WriteRead/wr_cnf_stat.f90"
#include "RNG/rng.f90"

END PROGRAM main
