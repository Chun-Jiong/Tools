program main
	implicit none
	real(8), allocatable      :: initial(:,:)
	real(8)                   :: s1, s2
	integer(8)                :: maxl
	real(8)                   :: c1, c2

	integer(8)                :: istart
	real(8)                   :: unitp, p
	real(8)                   :: is1, is2
	real(8)                   :: ss1, ss2

	character(100)            :: arg
	character(100)            :: fname
	integer(8)                :: i, i1, i2

	call getarg(1,fname)
	call getarg(2,arg)
	read(arg,*) istart
	call getarg(3,arg)
	read(arg,*) unitp
	write(*,'(A12,A50)')  "  name = ", trim(fname)
	write(*,'(A12,I8)')   "istart = ", istart
	write(*,'(A12,F8.3)') " unitp = ", unitp
	!stop

	open(1,file=trim(fname),action="read")
	do
		read(1,*,end=10) maxl, c1, c2
	end do
	10 rewind(1)
	!write(*,*) maxl
	allocate(initial(maxl,2))
	s1=0.d0; s2=0.d0
	do
		read(1,*,end=11) maxl, c1, c2
		initial(maxl,1:2) = (/c1,c2/)
		s1 = s1 + c1
		s2 = s2 + c2
	end do
	11 close(1)

	open(2,file=trim(fname)//".dat",action="write")
	ss1=0.d0; ss2=0.d0
	do i=1, istart
		write(2,'(F16.1,2ES16.8)') i*1.d0, initial(i,1)/s1, initial(i,2)/s2
		write(*,'(F16.1)') i*1.d0
		ss1 = ss1 + initial(i,1)
		ss2 = ss2 + initial(i,2)
	end do
	p=0.d0
	i2=istart

	do
		p=p+unitp
		if( floor(10**p-10**(p-unitp))>0 ) exit
	end do
	p=p-unitp

	is1=0.d0; is2=0.d0
	p=p+unitp
	i1=i2
	i2=i1+floor(10**p-10**(p-unitp))
	do i=istart+1, maxl
		if( i<=i2 ) then
			is1 = is1 + initial(i,1)
			is2 = is2 + initial(i,2)
		else
			if( ss1<s1-0.1d0 .or. ss2<s2-0.1d0 ) then
				ss1 = ss1 + is1
				ss2 = ss2 + is2
				write(2,'(F16.1,2ES16.8)') (i1+i2)/2.d0, is1/s1/(i2-i1), is2/s2/(i2-i1)
				write(*,'(F16.1)') (i1+i2)/2.d0
			else
				exit
			end if
			is1=initial(i,1)
			is2=initial(i,2)
			p=p+unitp
			i1=i2
			i2=i1+floor(10**p-10**(p-unitp))
		end if
	end do
	write(2,'(A2,2F16.1)') "# ", s1,  s2
	write(2,'(A2,2F16.1)') "# ", ss1, ss2
	close(2)
end program main
