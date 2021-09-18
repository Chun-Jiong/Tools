PROGRAM MAIN
	use, intrinsic :: iso_c_binding
	implicit none
#include <fftw3.f03>
	type(C_PTR)                    :: plan
	real(8), parameter             :: PI=3.14159265358979d0

	complex(C_DOUBLE_COMPLEX), allocatable :: fftw_in(:,:), fftw_out(:,:)


	integer(4)                  :: maxLx
	integer(4)                  :: Lx, Ly, Lz
	integer(4)                  :: ix, iy, iz
	integer(4)                  :: i, j, k

	character(100)            :: arg
	character(100)            :: fname
	real(8)                   :: ave, sigma2
	real(8), allocatable      :: rhot(:)
	real(8)                   :: c1, c2
	real(8)                   :: tauint
	integer(4)                :: maxt
	complex(C_DOUBLE_COMPLEX) :: c

	call getarg(1,fname)
	call getarg(2,arg)
	read(arg,*) maxLx

	!--- open file ---!
	open(1,file=trim(fname),action="read")
	Lx=0
	do
		read(1,*,end=10) c1, c2
		Lx = Lx+1
	end do
	10 rewind(1)
	if( Lx>maxLx ) Lx=maxLx

	!--- important ---!
	Lx = Lx*2
	!-----------------!
	Ly = 1
	allocate(fftw_in(Lx,Ly))
	allocate(fftw_out(Lx,Ly))
	plan = fftw_plan_dft_2d(Ly,Lx, fftw_in, fftw_out, FFTW_FORWARD, FFTW_ESTIMATE)

	allocate(rhot(0:Lx/2-1))
	rhot = 0.d0

	
	!--- setting fftw_in ---!
	ave = 0.d0
	iy=1
	do ix=1, Lx/2
		read(1,*) c1, c2
		fftw_in(ix,iy) = CMPLX(c2,0.d0)
		ave = ave + c2
		!write(*,*) c2
	end do
	close(1)
	ave = ave/(Lx/2)
	fftw_in(1:Lx/2,iy) = fftw_in(1:Lx/2,iy) - CMPLX(ave,0.d0)
	fftw_in(Lx/2+1:Lx,iy) = CMPLX(0.d0,0.d0)

	!===========================================!
	!--- cal rho(t) ---!
	rhot = 0.d0
	do i=0, Lx/2-1
		do j=0, Lx/2-1-i
			rhot(i) = rhot(i) + REAL(fftw_in(j +1, 1))*REAL(fftw_in(j+i +1, 1))
		end do
		!write(*,*) rhot(i)
	end do
	!write(*,*) "----------------------------"

	!--- DO FFT ---!
	call fftw_execute_dft(plan, fftw_in, fftw_out)

	!!--- print the result ---!
	!do iy=1, Ly
	!do ix=1, Lx
	!	write(*,*) fftw_out(ix,iy)
	!end do
	!end do

	!!===========================================!
	!!--- normal FT ---!
	!write(*,*) "----------------------------"
	!do iy=1, Ly
	!do ix=1, Lx
	!	c = cmplx(0.d0,0.d0)
	!	do j=1, Ly
	!	do i=1, Lx
	!		c = c + fftw_in(i,j)*exp(-cmplx(0.d0,1.d0)*2.d0*PI*(1.d0*(ix-1)*(i-1)/Lx+1.d0*(iy-1)*j/Ly))
	!	end do
	!	end do
	!	write(*,*) c
	!end do
	!end do

	!--- DO FFT AGAIN ---!
	iy = 1
	do ix=1, Lx
		c2 = abs(fftw_out(ix,iy))**2
		fftw_in(ix,iy) = CMPLX(c2,0.d0)
	end do
	call fftw_execute_dft(plan, fftw_in, fftw_out)
	!!!!!!!!!
	fftw_out = fftw_out/Lx
	!!!!!!!!!
	!--- print the result ---!
	!write(*,*) REAL(fftw_out(1,1))/(Lx/2)
	!do iy=1, Ly
	!do ix=1, Lx/2
	!	!write(*,*) fftw_out(ix,iy)
	!	write(*,*) rhot(ix-1), REAL(fftw_out(ix,iy))
	!end do
	!end do

	!--- BEGIN TO CALCULATE \tau ---!
	sigma2 = REAL(fftw_out(1,1))/(Lx/2)
	maxt = 10000
	if( Lx/2<maxt ) maxt=Lx/2

	tauint = -0.5d0
	do i=1, maxt
		c2 = REAL(fftw_out(i,1))/(Lx/2 - (i-1))/sigma2
		if( c2>0.d0 ) then
			tauint = tauint + c2
			!write(*,'(2ES16.8)') c2, tauint
		else
			!write(*,'(A2,ES16.8)') "# ", tauint
			write(*,'(ES16.8)') tauint
			exit
		end if
	end do
	!write(*,*) "# ", tauint

	!!--- use rhot to calculate \tau ---!
	!tauint = 0.d0
	!tauint = tauint+0.5d0
	!sigma2 = rhot(0)/(Lx/2)
	!do i=2, Lx/2
	!	tauint = tauint +  rhot(i-1)/(Lx/2 - (i-1))/sigma2
	!end do
	!write(*,*) tauint

	!--- destroy allocated arrarys ---!
	deallocate(fftw_in,fftw_out)
	!--- distroy plane ---!
	call fftw_destroy_plan(plan)
END PROGRAM MAIN
