program main
    implicit none
    real(8)                     :: ln0o5, p
    character(50)               :: filename, char1, charp
    integer(4)                  :: D
    integer(4)                  :: L, RL, subl, Vol
    real(8), allocatable        :: betaenergy(:,:)
    integer(4), parameter       :: maxb=1000
    real(8)                     :: list(maxb,2)
    integer(4)                  :: nbeta
    real(8)                     :: delta(2)
    real(8)                     :: c(2)
    real(8)                     :: slope
    integer(4)                  :: i, j

    write(*,*) "enter: filename  dimension  L  subl  p"
    call getarg(1,filename)
    call getarg(2, char1)
    read(char1,*) D
    call getarg(3, char1)
    read(char1,*) RL
    call getarg(4, char1)
    read(char1,*) subl
    call getarg(5, charp)
    read(charp,*) p

    Vol = RL**D * subl

    if( p<0.15 ) then
        write(*,*) "p<0.15"
        stop
    end if

    ln0o5 = log(p)

    open(1,file=trim(filename),action="read")
    nbeta = 0
    do
        read(1,*,end=10) L, c(1:2)
        nbeta = nbeta+1
    end do
    10 rewind(1)
    write(*,*) "nbeta = ", nbeta
    allocate(betaenergy(nbeta,2))
    do i=1, nbeta
        read(1,*) L, c(1:2)
        betaenergy(i,1) = c(1)
        betaenergy(i,2) = c(2)/c(1)
    end do
    close(1)

    betaenergy(:,2) = betaenergy(:,2) * Vol
    list(1,1:2) = betaenergy(1,1:2)
    i = 2
    j = 1
    do while ( i<=nbeta .and. j<maxb )
        delta(1) = betaenergy(i,1) - list(j,1)
        delta(2) = betaenergy(i,2) - list(j,2)
        !if( exp(delta(1)*delta(2))>0.15d0 ) then
        if( abs(exp(delta(1)*delta(2))-1.d0)<1.d-5 ) then ! if energy curve is flat
            list(j+1,1) = betaenergy(i,1)
            list(j+1,2) = betaenergy(i,2)
            j = j+1
            i = i+1
        else if( exp(delta(1)*delta(2))>p ) then
            i = i+1
        else
            slope = (betaenergy(i,2)-list(j,2))/(betaenergy(i,1)-list(j,1))
            delta(1) = sqrt( ln0o5/slope )
            list(j+1,1) = list(j,1) + delta(1)
            list(j+1,2) = list(j,2) + slope*delta(1)
            j = j+1
        end if
    end do

    open(1,file="betalist.txt",action="write")
    do i=1, j-1
        !write(1,'(2ES16.8)') list(i,1:2)
        !write(1,'(I4,3ES16.8)') i-1, list(i,1), list(i,2)/Vol*list(i,1), exp( (list(i,1)-list(i+1,1))*(list(i,2)-list(i+1,2)) )
        write(1,'(I4,ES16.8,3ES16.8)') i-1, list(i,1), &
            & 1.d0/list(i,1), list(i,2)/Vol, exp( (list(i,1)-list(i+1,1))*(list(i,2)-list(i+1,2)) )
    end do
    close(1)
    write(*,*) "length = ", j

    deallocate(betaenergy)
end program main
