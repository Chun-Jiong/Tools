!============== Write to files =====================================
!! THIS IS PROJECT-DEPENDENT 
!SUBROUTINE write2file(iblck)
!	IMPLICIT NONE
!	integer(4), intent(in)      :: iblck
!	integer(4)                  :: i, j, k
!
!	!-- open file ----------------------------------------------------
!	!open (id,file=trim(datafile),  access="append")
!	open (id,file=trim(datafile),  action="write")
!	write(id, *) "===================================================="
!
!	!-- write to data file ------------------------------------------------
!	!-- quantities --
!	write(id,"(A16,I6,4F16.8,I9,I3,I9)") trim(ident), Lx, beta, Jcp(1:3), Totsamp, nw, Seed
!
!	do i = 1, Nobs
!		write(id,'(I5,2ES20.8,F12.5)') i, Ave(i), Dev(i), Cor(i)
!		write(* ,'(I5,2ES20.8,F12.5)') i, Ave(i), Dev(i), Cor(i)
!	enddo
!
!END SUBROUTINE write2file
!===================================================================


SUBROUTINE write2file(iobs)
	IMPLICIT NONE
	integer(4), intent(in)      :: iobs
	integer(4)                  :: i, j, k

	!-- open file ----------------------------------------------------
	open (1,file=trim(output), action="write")
	do i = 1, iobs
		write(1,'(I4,3ES20.8)') i, Ave(i), Dev(i), Cor(i)
	enddo
	close(1)
END SUBROUTINE write2file
