PROGRAM MAIN
	use, intrinsic :: iso_c_binding
	implicit none
#include <fftw3.f03>
	type(C_PTR)                 :: plan
	real(8), parameter          :: PI=3.14159265358979d0
	real(8), parameter          :: eps = 1.d-14    ! very small number
	real(8), parameter          :: tol = 0.2d0     ! tolerance for Cor

	complex(C_DOUBLE_COMPLEX), allocatable :: fftw_in(:), fftw_out(:)

	integer(4)                  :: narg                 ! the number of input argument
	character(100)              :: filename
	integer(4)                  :: deltat               ! the time unit
	real(8), allocatable        :: a(:,:)               ! sample record
	real(8), allocatable        :: apqamq(:)            ! A(q)A(-q)
	real(8), allocatable        :: autocf(:,:)          ! autocorrelation function
	integer(4)                  :: NBlck                ! the number of blcks
	integer(4)                  :: Nsamp                ! sample length in a block
	real(8)                     :: c
	real(8)                     :: nor
	real(8), allocatable        :: ave(:), dev(:), cor(:)
	real(8)                     :: devp, devn
	integer(4)                  :: iblck
	integer(4)                  :: isamp

	narg = command_argument_count()

	if( narg == 0 ) then
		write(*,*) "Please enter the file name."
		write(*,*) "OR if you want help, please use the parameter -h/--help"
		stop
	end if

	!--- get the file name ---!
	call getarg(1,filename)

	if( trim(filename) == "-h" .or. trim(filename) == "--help" ) then
		call help()
		stop
	end if

	!--- open file ---!
	open(1,file=trim(filename),action="read")
	!--- check data ---!
	read(1,*) deltat, Nsamp
	NBlck = 0
	do
		read(1,*,end=10) iblck
		do isamp = 1, Nsamp
			read(1,*) c
		end do
		NBlck = NBlck+1
	end do
	10 rewind(1)
	allocate(a(NBlck,Nsamp))
	allocate(autocf(NBlck,Nsamp))
	!--- read data ---!
	read(1,*) deltat, Nsamp
	do iblck = 1, NBlck
		read(1,*) c
		do isamp = 1, Nsamp
			read(1,*) c
			a(iblck,isamp) = c
		end do
	end do
	close(1)
	!--- read data done ---!

	!print *, NBlck, Nsamp

	!--- calculate deltaA ---!
	do iblck = 1, NBlck
		c = sum(a(iblck,1:Nsamp))/Nsamp
		a(iblck,:) = a(iblck,:) - c
	end do


	!--- the fftw setting --!
	allocate(fftw_in(2*Nsamp))
	allocate(fftw_out(2*Nsamp))
	allocate(apqamq(2*Nsamp))

	plan = fftw_plan_dft_1d(2*Nsamp, fftw_in, fftw_out, FFTW_FORWARD, FFTW_ESTIMATE)


	!--- calculate autocorrelation function ---!
	do iblck = 1, NBlck
		do isamp = 1, Nsamp
			fftw_in(isamp) = CMPLX(a(iblck,isamp),0.d0)
		end do
		fftw_in(Nsamp+1:2*Nsamp) = CMPLX(0.d0,0.d0)
		call fftw_execute_dft(plan, fftw_in, fftw_out)
		do isamp = 1, 2*Nsamp
			apqamq(isamp) = Real( conjg(fftw_out(isamp)) * fftw_out(isamp) )
		end do

		do isamp = 1, 2*Nsamp
			fftw_in(isamp) = CMPLX(apqamq(isamp),0.d0)
		end do
		call fftw_execute_dft(plan, fftw_in, fftw_out)
		autocf(iblck,1:Nsamp) = Real(fftw_out(1:Nsamp))
		autocf(iblck,1:Nsamp) = autocf(iblck,1:Nsamp) / (2.d0*Nsamp)

		nor = 0.d0
		do isamp = Nsamp, 1, -1
			nor = nor + a(iblck,Nsamp-isamp+1)**2
			autocf(iblck, isamp) = autocf(iblck, isamp) / nor
		end do
	end do


	!--- calculate ave, dev ---!
	allocate(ave(Nsamp),dev(Nsamp),cor(Nsamp))
	if( NBlck>1 ) then
		nor = 1.d0/NBlck
		do isamp = 1, Nsamp
			ave(isamp) = sum(autocf(1:NBlck,isamp)) / NBlck
			devp = 0.d0;  cor(isamp) = 0.d0;  dev(isamp) = 0.d0
			do iblck = 1,  NBlck
				devn   = autocf(iblck,isamp)-ave(isamp)
				dev(isamp) = dev(isamp)+devn*devn
				cor(isamp) = cor(isamp)+devn*devp
				devp   = devn
			enddo 
			dev(isamp) = dev(isamp)*nor
			cor(isamp) = cor(isamp)*nor
			if(dev(isamp)>eps)   cor(isamp) = cor(isamp)/dev(isamp)
			dev(isamp) = dsqrt(dev(isamp)/(NBlck-1.d0))
			!if(dabs(cor(isamp))>tol) prt = .false.
		end do
	else if( NBlck==1 ) then
		ave(1:Nsamp) = autocf(NBlck,1:Nsamp)
		dev = 0.d0
		cor = 0.d0
	else
		write(*,*) "Err: NBlck<1", NBlck
		stop
	end if

	write(*,'(A12,I6)') "# Blocks = ", NBlck
	do isamp = 1, Nsamp
		write(*,'(I6, 2ES16.8)') isamp-1, ave(isamp), dev(isamp)
		!write(*,'(ES24.1, 2ES16.8)') (isamp-1)*deltat, ave(isamp), dev(isamp)
	end do


	!--- destroy allocated arrarys ---!
	deallocate(a,autocf,apqamq)
	deallocate(ave,dev,cor)
	deallocate(fftw_in,fftw_out)
	!--- distroy plane ---!
	call fftw_destroy_plan(plan)

CONTAINS
	SUBROUTINE help
		implicit none

		write(*,*) "The program is used to calculate the autocorrelation function using FFTW library"
		write(*,*) "Usage:  ./AutocorrelationFunction filename"
		write(*,*) "The data format:"
		write(*,*) "\Delta t     Length(Len)"
		write(*,*) "0   -- used to seperate different bins"
		write(*,*) "A1"
		write(*,*) "A2"
		write(*,*) "..."
		write(*,*) "$A_{Len}$"
		write(*,*) "0"
		write(*,*) "A1"
		write(*,*) "A2"
		write(*,*) "..."
		write(*,*) "A_{Len}"
		write(*,*) "0"
		write(*,*) "A1"
		write(*,*) "A2"
		write(*,*) "..."
		write(*,*) "A_{Len}"
		write(*,*) "."
		write(*,*) "."
		write(*,*) "."
	END SUBROUTINE help

END PROGRAM MAIN
