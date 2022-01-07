
!==============Statistics ==========================================
!! THIS IS PROJECT-INDEPENDENT 
SUBROUTINE stat_analy(iblck)
	implicit none
	integer          :: iblck
	integer          :: j, k, k0
	double precision :: devn, devp, nor

	NObs_b = NObs
	NObs_c = 0

	! -- calculate average -------------------------------------------
	nor  = 1.d0/(iblck*1.d0)
	do j = 1, NObs_b
		Ave(j) = nor*Sum(Obs(j,1:iblck))
	enddo

	!Coarsen: do
		! -- calculate error and t=1 correlation for basics obs.--------
		prt = .true.
		DO j = 1, NObs_b
			devp = 0.d0;  Cor(j) = 0.d0;  Dev(j) = 0.d0
			do k = 1,  iblck
				devn   = Obs(j,k)-Ave(j)
				Dev(j) = Dev(j)+devn*devn
				Cor(j) = Cor(j)+devn*devp
				devp   = devn
			enddo 
			Dev(j)   = Dev(j)*nor;        Cor(j) = Cor(j)*nor
			if(Dev(j)>eps)                Cor(j) = Cor(j)/Dev(j)
			Dev(j)   = dsqrt(Dev(j)/(iblck-1.d0))
			if(dabs(Cor(j))>tol) prt = .false.
		ENDDO 

		!IF(prt)                         EXIT Coarsen 
		!IF(NBlck<=64)    THEN
		!	prt = .false.;                EXIT Coarsen 
		!ENDIF

	!	! -- coarsen blocking ------------------------------------------
	!	nor   = nor*2.d0
	!	call coarsen_data
	!enddo Coarsen 

	!! -- define auxillary variables and average of composite obs.-----
	if( NObs_c>0 ) then
		call cal_Obs_comp(iblck)
	end if

	!! -- calculate error and t=1 correlation for composite obs.-----
	!do j = 1+NObs_b, NObs
	!	devp = 0.d0;  Cor(j) = 0.d0;  Dev(j) = 0.d0
	!	DO k = 1,  iblck
	!		devn   = Obs(j,k)-Ave(j)
	!		Dev(j) = Dev(j)+devn*devn
	!		Cor(j) = Cor(j)+devn*devp
	!		devp   = devn
	!	ENDDO
	!	Dev(j)   = Dev(j)*nor;        Cor(j) = Cor(j)*nor
	!	IF(Dev(j)>eps)                Cor(j) = Cor(j)/Dev(j)
	!	Dev(j)   = dsqrt(Dev(j)/(iblck-1.d0))
	!enddo

	return
END SUBROUTINE stat_analy
!===================================================================

!==============Calculate fluctuation quantity ==================
!! THIS IS PROJECT-INDEPENDENT 
!! C = V(Ave(b2)-Ave(b1)^2)
SUBROUTINE cal_fluct(ib, jb,b1,b2)
	implicit none
	integer(4), intent(in)      :: ib
	integer(4), intent(in)      :: jb, b1, b2
	integer(4)                  :: k

	!-- Average ----------------------------------------------------
	!Ave(jb) = Vol*(Ave(b2)-Ave(b1)**2.d0)

	!-- Obs(j,k) series --------------------------------------------
	do k = 1, ib
		Obs(jb,k) = Obs(b2,k)-Obs(b1,k)**2.d0
		Obs(jb,k) = Obs(jb,k)
	end do
	Ave(jb) = SUM(Obs(jb,1:ib))/ib

END SUBROUTINE cal_fluct


!==============Calculate composite observables =====================
!! THIS IS PROJECT-INDEPENDENT 
!! call in 'stat_alan'
SUBROUTINE cal_Obs_comp(ib)
	implicit none
	integer(4), intent(in)      :: ib
	integer(4)                  :: jb, b2, b3, b1

	!-- calculate the average ----------------------------------------

	jb = NObs_b+1;   call cal_fluct(ib, jb, 1, 2)  ! Cv
	jb = NObs_b+2;   call cal_fluct(ib, jb, 3, 4)  ! X0
END SUBROUTINE cal_Obs_comp
!===================================================================
