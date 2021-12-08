! image data
! infile's format
!
! 0  !also could be other integer number
! real(SqSq(0,0,0))  aimag(SqSq(0,0,0))
! real(SqSq(1,0,0))  aimag(SqSq(1,0,0))
! real(SqSq(2,0,0))  aimag(SqSq(2,0,0))
! ...
! 0
! real(SqSq(0,0,0))  aimag(SqSq(0,0,0))
! real(SqSq(1,0,0))  aimag(SqSq(1,0,0))
! real(SqSq(2,0,0))  aimag(SqSq(2,0,0))
! ...

program structure_factor
	implicit none

	integer(4)                  :: narg                 ! the number of input argument

	character(50)               :: filename
	integer(4)                  :: GN
	complex(8), allocatable     :: bin(:,:,:,:,:,:,:)
	complex(8), allocatable     :: SqSq(:,:,:,:,:,:)
	real(8)                     :: q(3), dq(3)
	integer(4)                  :: Lx, Ly, Lz
	integer(4)                  :: subl
	integer(4)                  :: h, l, m
	integer(4)                  :: h1, h2, l1, l2, m1, m2
	integer(4)                  :: fm
	integer(4)                  :: NBlck
	real(8)                     :: reclatvec(3,3)
	real(8)                     :: sublatvec(3,3)
	!integer(4)                  :: nperiod
	real(8)                     :: nperiod
	character(50)               :: charnum, charformat
	integer(4)                  :: iGN
	integer(4)                  :: ix, iy, iz
	integer(4)                  :: sb1, sb2
	real(8)                     :: clx(2)
	complex(8)                  :: c
	real(8)                     :: br
	integer(4)                  :: ndata
	integer(4)                  :: i

	narg = command_argument_count()
	if( narg<7 ) then
		write(*,*) "Please enter: filename, GN, Lx, Ly, Lz, subl, nperiod; fm"
		stop
	end if

	call getarg(1, filename)
	call getarg(2, charnum)
	read(charnum,*) GN
	call getarg(3, charnum)
	read(charnum,*) Lx
	call getarg(4, charnum)
	read(charnum,*) Ly
	call getarg(5, charnum)
	read(charnum,*) Lz
	call getarg(6, charnum)
	read(charnum,*) subl
	call getarg(7, charnum)
	read(charnum,*) nperiod
	if( narg>7 ) then
		call getarg(8, charnum)
		read(charnum,*) fm      ! 1: keep q=0;  0: set S(q=0)=0
	else
		fm=1
	end if

	allocate(SqSq(0:GN,subl,subl,0:Lx-1,0:Ly-1,0:Lz-1))

	reclatvec = 0.d0
	open(1,file="reclatvec.dat",action="read")
	do i=1, 3
		read(1,*) reclatvec(i,1:3)
	end do
	close(1)

	sublatvec = 0.d0
	open(1,file="sublatvec.dat",action="read")
	do i=1, subl
		read(1,*) sublatvec(i,1:3)
	end do
	close(1)

	NBlck = 0
	!SqSq = 0.d0
	SqSq = CMPLX(0.d0,0.d0)
	ndata = 0
	!--- START READ DATA -------------------------------------------------!
	open(1,file=trim(filename)//".txt",action="read")
	!--- output bin ---!
	open(2,file=trim(filename)//"-bin.dat",action="write")
	do
		read(1,*,end=10) i
		do iz=0, Lz-1; do iy=0, Ly-1; do ix=0, Lx-1
			write(2,'(4I4)') ix+iy*Lx+iz*Lx*Ly+1, ix, iy, iz
			q = reclatvec(1,1:3)/Lx*ix + reclatvec(2,1:3)/Ly*iy + reclatvec(3,1:3)/Lz*iz
			do iGN=0, GN
				c = CMPLX(0.d0,0.d0)
				do sb2=1, subl; do sb1=1,subl      !!!!! sb2, sb1
					br = -dot_product(q,sublatvec(sb1,1:3)-sublatvec(sb2,1:3))
					read(1,*) clx(1), clx(2)
					c = c + CMPLX(clx(1),clx(2))*CMPLX(dcos(br),dsin(br))
					ndata = ndata+1
				end do; end do
				write(2,'(ES16.8)') real(c)
			end do
		end do; end do; end do
		NBlck = NBlck+1
	end do
	10 rewind(1)
	close(2)
	if( GN==0 ) call system("rm -f "//trim(filename)//"-bin.dat")
	if( ndata/=(Lx*Ly*Lz*subl*subl)*NBlck*(GN+1) .or. &
		ix/=Lx .or. iy/=Ly .or. iz/=Lz .or. &
		sb1/=subl+1 .or. sb2/=subl+1 .or. &
		iGN/=GN+1 ) then
		write(*,*) "Err: the data is incomplete."
		stop
	end if
	!write(*,*) ndata, ix, iy, iz, sb1, sb2, iGN
	write(*,'(A10,I6)') "NBlck = ", NBlck

	do
		read(1,*,end=20) i
		do iz=0, Lz-1; do iy=0, Ly-1; do ix=0, Lx-1
			do iGN=0, GN
				do sb2=1, subl; do sb1=1,subl      !!!!! sb2, sb1
					read(1,*) clx(1), clx(2)
					SqSq(iGN,sb1,sb2,ix,iy,iz) = SqSq(iGN,sb1,sb2,ix,iy,iz) + CMPLX(clx(1),clx(2))
					!ndata = ndata+1
				end do; end do
			end do
		end do; end do; end do
		!NBlck = NBlck+1
	end do
	20 close(1)
	!--- READ DATA END ---------------------------------------------------!

	!if( fm==0 ) then
	!	do sb1=1, subl
	!		SqSq(:,sb1,sb1,0,0,0)  = CMPLX(0.d0,0.d0)
	!	end do
	!end if

	SqSq = SqSq/(NBlck*1.d0)

	open(1,file=trim(filename)//"-sum.dat",action="write")
	open(2,file=trim(filename)//"-trace.dat",action="write")
	if( GN>0 ) open(3,file=trim(filename)//"-tau.dat",action="write")
	if( Lz==1 ) then
		m1=0; m2=0
	else
		!m1=-nperiod*Lz+1; m2=nperiod*Lz-1
		m1=-floor(nperiod*Lz+1.d-6); m2=floor(nperiod*Lz+1.d-6)
	end if
	l1=-floor(nperiod*Ly+1.d-6); l2=floor(nperiod*Ly+1.d-6)
	h1=-floor(nperiod*Lx+1.d-6); h2=floor(nperiod*Lx+1.d-6)
	do m=m1, m2
	do l=l1, l2
		do h=h1, h2
	!do m=0, Lz-1
	!do l=0, Ly-1
	!	do h=0, Lx-1
			iGN=0
			q = reclatvec(1,1:3)/Lx*h + reclatvec(2,1:3)/Ly*l + reclatvec(3,1:3)/Lz*m
			ix = mod(h+floor(nperiod+1.d0)*Lx,Lx)
			iy = mod(l+floor(nperiod+1.d0)*Ly,Ly)
			iz = mod(m+floor(nperiod+1.d0)*Lz,Lz)
			c = CMPLX(0.d0,0.d0)
			dq = reclatvec(1,1:3)/Lx*(h-ix) + reclatvec(2,1:3)/Ly*(l-iy) + reclatvec(3,1:3)/Lz*(m-iz)
			do sb2=1,subl; do sb1=1,subl
				!br = -dot_product(dq,sublatvec(sb1,1:3)-sublatvec(sb2,1:3))
				br = -dot_product(q,sublatvec(sb1,1:3)-sublatvec(sb2,1:3))
				c = c + SqSq(iGN,sb1,sb2,ix,iy,iz)*CMPLX(dcos(br),dsin(br))
			end do; end do
			!-------------------
			if( fm==0 ) then
				if( mod(h,Lx)==0 .and. mod(l,Ly)==0 .and. mod(m,Lz)==0 ) c= CMPLX(0.d0,0.d0)
			end if
			!-------------------
			write(1,'(5ES12.4)') q(1),q(2),q(3),real(c),aimag(c)

			c = CMPLX(0.d0,0.d0)
			do sb1=1, subl
				c = c + SqSq(iGN,sb1,sb1,ix,iy,iz)
			end do
			!c = c/subl
			!-------------------
			if( fm==0 ) then
				if( mod(h,Lx)==0 .and. mod(l,Ly)==0 .and. mod(m,Lz)==0 ) c= CMPLX(0.d0,0.d0)
			end if
			!-------------------
			write(2,'(5ES12.4)') q(1),q(2),q(3),real(c),aimag(c)

			ix=h; iy=l; iz=m
			!if( GN>0 ) then
			if(   GN>0  .and.             &
				  ix>=0 .and. ix<Lx .and. &
				  iy>=0 .and. iy<Ly .and. &
				  iz>=0 .and. iz<Lz       &
				 ) then
				write(3,'(A2, 4I4)') "#", ix+iy*Lx+iz*Lx*Ly+1, ix, iy, iz
				do iGN=0, GN
					c = CMPLX(0.d0,0.d0)
					do sb2=1,subl; do sb1=1,subl
						!br = -dot_product(dq,sublatvec(sb1,1:3)-sublatvec(sb2,1:3))
						br = -dot_product(q,sublatvec(sb1,1:3)-sublatvec(sb2,1:3))
						c = c + SqSq(iGN,sb1,sb2,ix,iy,iz)*CMPLX(dcos(br),dsin(br))
					end do; end do
					write(3,'(ES12.4)') real(c)
				end do
				write(3,'()'); write(3,'()')
			end if
		end do
		write(1,'()')
		write(2,'()')
	end do
	end do
	close(1)
	close(2)
	if( GN>0 ) close(3)
END program structure_factor
