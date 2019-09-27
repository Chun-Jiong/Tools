      program mrq
c     version jun 1995  
c     levenberg-marquardt method for nonlinear fits (num. recipes)
      implicit real*8(a-h,o-z)   
      external qcalc
      parameter (maxdata=400,maxiter=2000,maxp=20) 
      integer idim,ndata,n,i,j,m,ii,qi,qj,iter
      integer l(maxdata),lr(maxdata),lista(maxp)  
      character*8  aa(maxp) 
      character*73 text 
      real*8 cc,wkc,sccr,yi,yc,dchisq,chisq,ochisq,conv,alamda
      real*8 erra(maxp),dyda(maxp),a(maxp),alpha(maxp,maxp) 
      real*8 covar(maxp,maxp) 
      real*8 c1(maxdata)
      real*8 q1(maxdata),dq1(maxdata) 
      real*8 x(maxdata),y(maxdata),sig(maxdata) 
      call init(npar,mfit,mins,maxs,aa,a,cc,yi,yt,sccr,conv,lista)
      open (9,file='ou0.tau',access='append') 
      open (8,file='cmp.tau') 
      open (10,file='lmq.conv') 
      ndata=0 
      n=1 
      call input(n,l,q1,dq1,c1,ndata,x,lr,sig,y,yt,cc,sccr,mins,maxs)
      alamda=-1.d0  
      call mrqmin (lr,x,y,sig,ndata,a,npar,lista,mfit,
     *covar,alpha,maxp,chisq,qcalc,alamda)  
c     fit until convergence has been reached
      ochisq=chisq  
      dchisq=-1.d0 
      write(10,'(''convergence data mrq''/''maxiter='',i5
     */''conv='',d18.10//'' iter    dchisq''/)') maxiter,conv
      iter=0
  122 continue
      if (((dchisq.gt.0).and.(dchisq.le.conv)).or.(iter.ge.maxiter)) 
     *goto 123
      call mrqmin (lr,x,y,sig,ndata,a,npar,lista,mfit,
     *       covar,alpha,maxp,chisq,qcalc,alamda)  
      dchisq=ochisq-chisq
      ochisq=chisq 
      iter = iter+1
      write(10,'(i5,d18.10)') iter,dchisq
      goto 122
  123 continue
c     calculate covariances 
      alamda = 0
      call mrqmin (lr,x,y,sig,ndata,a,npar,lista,mfit, 
     *     covar,alpha,maxp,chisq,qcalc,alamda)  
c     write results to file xxx.outp
      do i=1,npar 
         erra(i)=dsqrt(covar(i,i)) 
      enddo 
      write (9,*) '****************************************'
      write(9,*) 
     *'q(k,l)=l**(2*a(1)-2)*(a(2)+a(3)*l**a(4)+a(5)*l**(2-2*a(1)))'
      write (9,*)
      open (7,file='lmq.par') 
      do i=1,npar 
         write (7,'(i4,f22.10)') i,a(i)
         write (9,'(i3,2x,a8,"=",f16.10," +-",f16.10)')
     *       i,aa(i),a(i),erra(i)
      enddo 
      write (9,'(/a13,d16.10)')'convergence: ',dchisq
      write (9,'(a13,d16.10)') 'chi squared: ',chisq 
      write (9,'(a13,i4)')    '# degr. fr.: ',(ndata-mfit)
      write (9,'(a13,i4)')    'iterations:  ',iter 
      write (9,'(a20,f8.4)')  'scaling-width crit. ',sccr     

      write (9,'(A4,A15,A20,A4,3A20)') 'l','c1', 'input', '    ',
     * 'error', 'computed', 'rel. dev.'
      do i=1,ndata 
         call qcalc (lr(i),x(i),a,yc,dyda,npar) 
c        write (9,'(2i4,f12.6,f10.6," +-",3f10.6)') i,lr(i),x(i)
         write (9,'(i4,f15.10,f20.10," +- ",3f20.10)') lr(i),x(i),
     *    y(i),  sig(i),yc,((y(i)-yc)/sig(i)) 
      enddo
c     if (4*iter.ge.maxiter) call dtest(lr,x,a,dyda,mfit)
      call dtest(lr,x,a,dyda,mfit)
      close (7)
      close (8)
      close (9)
      close (10)
      end 
      subroutine dtest(lr,x,a,dyda1,mfit)
      implicit real*8(a-h,o-z)   
      parameter (maxdata=400,maxp=20)
      real*8 dyda1(maxp),dyda2(maxp),a(maxp)
      real*8 x(maxdata)
      integer lr(maxdata)
      save
      da=1.d-4
      write(10,100) da
  100 format('slow convergence: test derivatives for da=',d12.4//
     *'par    analytic        numerical     difference'/)
      call qcalc (lr(1),x(1),a,y1,dyda1,npar) 
      do 101 i=1,mfit
      a(i)=a(i)+da
      call qcalc (lr(1),x(1),a,y2,dyda2,npar) 
      a(i)=a(i)-da
      dydaa=0.5d0*(dyda1(i)+dyda2(i))
      dydan=(y2-y1)/da
      dif=dydaa-dydan
      write(10,102)i,dydaa,dydan,dif
  102 format(i3,2d16.8,d12.4)
  101 continue 
      return
      end 
c  subroutine marquardt-methode
      subroutine mrqmin(l,x,y,sig,ndata,a,npar,lista,mfit,    
     *    covar,alpha,mxp,chisq,funcs,alamda) 
      implicit real*8(a-h,o-z)   
      external funcs
      parameter (maxp=20) 
      real*8 x(ndata),y(ndata),sig(ndata),a(mxp),atry(maxp)
      real*8 alpha(mxp,mxp),beta(maxp)
      real*8  covar(mxp,mxp)
      real*8 da(maxp)   
      integer lista(mxp),l(ndata),k,j,kk,ihit
      save
      if(alamda.lt.0.)then
        kk=mfit+1 
        do 12 j=1,npar
          ihit=0
          do 11 k=1,mfit
            if(lista(k).eq.j)ihit=ihit+1
11        continue
          if (ihit.eq.0) then 
            lista(kk)=j 
            kk=kk+1 
          else if (ihit.gt.1) then
            stop 'improper permutation in lista' 
          endif 
12      continue
        if (kk.ne.(npar+1)) stop 'improper permutation in lista' 
        alamda=0.001d0
        call mrqcof(l,x,y,sig,ndata,a,npar,lista,mfit,alpha,beta,mxp, 
     *    chisq,funcs)
        ochisq=chisq
        do 13 j=1,npar
          atry(j)=a(j)
13      continue
      endif 
      do 15 j=1,mfit
        do 14 k=1,mfit
          covar(j,k)=alpha(j,k) 
14      continue
        covar(j,j)=alpha(j,j)*(1.d0+alamda) 
        da(j)=beta(j) 
15    continue
      call gaussj(covar,mfit,mxp,da,1,1)
      if(alamda.eq.0.)then
        call covsrt(covar,mxp,npar,lista,mfit)
        return
      endif 
      do 16 j=1,mfit
        atry(lista(j))=a(lista(j))+da(j)
16    continue
      call mrqcof(l,x,y,sig,ndata,atry,npar,lista,mfit,covar,da,mxp,
     *   chisq,funcs) 
      if(chisq.lt.ochisq)then 
        alamda=0.1d0*alamda 
        ochisq=chisq
        do 18 j=1,mfit
          do 17 k=1,mfit
            alpha(j,k)=covar(j,k) 
17        continue
          beta(j)=da(j) 
          a(lista(j))=atry(lista(j))
18      continue
      else
        alamda=10.d0*alamda 
        chisq=ochisq
      endif 
      return
      end 
      subroutine mrqcof(l,x,y,sig,ndata,a,npar,lista,mfit,alpha,beta, 
     *    nalp,chisq,funcs) 
      implicit real*8(a-h,o-z)   
      external funcs
      parameter (maxp=20) 
      real*8 x(ndata),y(ndata),sig(ndata),
     *    dyda(maxp),a(nalp)
      integer lista(mfit),l(ndata)
      real*8 alpha(nalp,nalp),beta(nalp)
      save
      do 12 j=1,mfit
        do 11 k=1,j 
          alpha(j,k)=0. 
11      continue
        beta(j)=0.
12    continue
      chisq=0.
      do 15 i=1,ndata 
        call funcs(l(i),x(i),a,ymod,dyda,npar)  
        sig2i=1./(sig(i)*sig(i))
        dy=y(i)-ymod
        do 14 j=1,mfit
          wt=dyda(lista(j))*sig2i 
          do 13 k=1,j 
            alpha(j,k)=alpha(j,k)+wt*dyda(lista(k)) 
13        continue
          beta(j)=beta(j)+dy*wt 
14      continue
        chisq=chisq+dy*dy*sig2i 
15    continue
      do 17 j=2,mfit
        do 16 k=1,j-1 
          alpha(k,j)=alpha(j,k) 
16      continue
17    continue
      return
      end 
      subroutine gaussj(a,n,np,b,m,mp)
      implicit real*8 (a-h,o-z) 
      parameter (nmax=50) 
      dimension a(np,np),b(np,mp),ipiv(nmax),indxr(nmax),indxc(nmax)
      save
      do 11 j=1,n 
        ipiv(j)=0 
11    continue
      do 22 i=1,n 
        big=0.
        do 13 j=1,n 
          if(ipiv(j).ne.1)then
            do 12 k=1,n 
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then 
                  big=abs(a(j,k)) 
                  irow=j
                  icol=k
                endif 
              else if (ipiv(k).gt.1) then 
                stop 'singular matrix' 
              endif 
12          continue
          endif 
13      continue
        ipiv(icol)=ipiv(icol)+1 
        if (irow.ne.icol) then
          do 14 l=1,n 
            dum=a(irow,l) 
            a(irow,l)=a(icol,l) 
            a(icol,l)=dum 
14        continue
          do 15 l=1,m 
            dum=b(irow,l) 
            b(irow,l)=b(icol,l) 
            b(icol,l)=dum 
15        continue
        endif 
        indxr(i)=irow 
        indxc(i)=icol 
        if (a(icol,icol).eq.0.) stop 'singular matrix.'
        pivinv=1./a(icol,icol)
        a(icol,icol)=1. 
        do 16 l=1,n 
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m 
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0. 
            do 18 l=1,n 
              a(ll,l)=a(ll,l)-a(icol,l)*dum 
18          continue
            do 19 l=1,m 
              b(ll,l)=b(ll,l)-b(icol,l)*dum 
19          continue
          endif 
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n 
            dum=a(k,indxr(l)) 
            a(k,indxr(l))=a(k,indxc(l)) 
            a(k,indxc(l))=dum 
23        continue
        endif 
24    continue
c        write(17,*)'ex',a
      return
      end 
      subroutine covsrt(covar,ncvm,npar,lista,mfit) 
      implicit real*8(a-h,o-z)   
      real*8 covar(ncvm,ncvm) 
      dimension lista(mfit) 
      save
      do 12 j=1,npar-1
        do 11 i=j+1,npar
          covar(i,j)=0. 
11      continue
12    continue
      do 14 i=1,mfit-1
        do 13 j=i+1,mfit
          if(lista(j).gt.lista(i)) then 
            covar(lista(j),lista(i))=covar(i,j) 
          else
            covar(lista(i),lista(j))=covar(i,j) 
          endif 
13      continue
14    continue
      swap=covar(1,1) 
      do 15 j=1,npar
        covar(1,j)=covar(j,j) 
        covar(j,j)=0. 
15    continue
      covar(lista(1),lista(1))=swap 
      do 16 j=2,mfit
        covar(lista(j),lista(j))=covar(1,j) 
16    continue
      do 18 j=2,npar
        do 17 i=1,j-1 
          covar(i,j)=covar(j,i) 
17      continue
18    continue
      return
      end 
      subroutine inpar(a,npar)
      implicit real*8(a-h,o-z)   
      real*8 a(1),b
      open (7,file='lmq.par') 
11    read (7,*,end=10) i,b
      if (i.gt.npar) then
         write(9,'("i=",i4," exceeds npar")') i
         stop
      endif
      a(i)=b
      goto 11
10    continue
      close (7) 
      return
      end 
      subroutine init(npar,mfit,mins,maxs,aa,a,cc,yi,yt,sccr,conv,ip)
      implicit real*8(a-h,o-z)   
      parameter (maxp=20) 
      character*8  aa(maxp) 
      real*8 a(maxp),cc,sccr,conv
	  real*8 yi,yt,yh,y2
      integer ip(maxp)
      save
      sccr=2.00d0
      conv=0.2d-7
c     number of parameters to be fitted
      mfit = 3
c     number of parameters in expression
      npar = 5

c     y=l**(2*a(1)-2)*(a(2)+a(3)*l**a(4)+a(5)*l**(2-2*a(1)))
c     define 6-char identification for each parameter
      aa(1) =' yh ' 
      aa(2) =' a0 ' 
      aa(3) =' b1 ' 
      aa(4) =' yi ' 
      aa(5) =' b2 ' 

c     additional input from terminal
      write (*,'(a)') 'minimum system size?'
      read (*,*) mins
c     write (*,'(a)') 'irrelevant exponent?'
c     read (*,*) yi
c     take interval for relevant data-points 
c            write (*,'(a)') 'value of critical coupling: _'
c            read (*,*) cc
c            write (*,'(a)') 'value scaling-width criterion: _' 
c            read (*,*) sccr

c     define initial values of a(i),yt,cc
c     y=l**(2*a(1)-2)*(a(2)+a(3)*l**a(4)+a(5)*l**(2-2*a(1)))
      yh =  1.875d0
      yi = -2.000d0
      y2 = 2.d0-2*yh
      do ik=1,npar
         a(ik)=0.0d0
      enddo
      a(1) =  yh
      a(2) =  0.653d0
      a(3) = -0.118d0
      a(4) =  yi
      a(5) =  0.d0   

      call inpar(a,npar)
c     define suitable permutation ip of parameters
      do i=1,maxp 
         ip(i)=0
         if (i.le.mfit) ip(i)=i   
      enddo
c     ip(i)=j means that:
c     parameter #j is put on the i-th place of the list of
!      ip(4 )=6
c     parameters to be fitted (length mfit)
      return
      end
      subroutine qcalc (l,x,a,y,dyda,npar)
c     calculate q(k,l) and its first derivatives
      implicit real*8(a-h,o-z)   
      parameter (maxp=20) 
      real*8 a(maxp),dyda(maxp),y,x,wk
      integer l 

c     y=l**(2*a(1)-2)*(a(2)+a(3)*l**a(4)+a(5)*l**(2-2*a(1)))
      w1=l**(2*a(1)-2)
      w2=l**a(4)
      w3=l**(2-2*a(1))
      dl=dlog(1.0d0*l)
      w1dl=w1*dl*2.d0
      w2dl=w2*dl
      w3dl=w3*dl*(-2.d0)

      y=w1*(a(2)+a(3)*w2+a(5)*w3)

      dyda(1) = w1dl * (a(2)+a(3)*w2+a(5)*w3)
     *        + w1   * a(5) * w3dl
      dyda(2) = w1
      dyda(3) = w1 * w2
      dyda(4) = w1 * a(3) * w2dl
      dyda(5) = w1 * w3
      return
      end 
      subroutine input(n,li,q,erq,c1,ndata,x,l,ery,y,yt,cc,sccr,m1,m2)
c     read input and reject data if necessary
c     define:   ndata # accepted data points
c     define:   x(ndata) input for independent variable
c     define:   y(ndata) input for   dependent variable
c     define: ery(ndata) error in y
c     define:   l(ndata) input for finite-size parameter 
      implicit real*8(a-h,o-z)   
      parameter (maxdata=400) 
      integer l(1),li(1)
      real*8 q(1),erq(1),c1(1),x(1),y(1),ery(1) 
      real*8 yt,cc,wkc,sccr
      save
      do 20 i=1,maxdata
         read (8,*,end=10) li(n),c1(n),q(n),erq(n)
c         print *, li(n),q(n),erq(n)
         n=n+1 
         if(n+1.gt.maxdata)then
            write(*,*)' number of datapoints > max=',
     *         maxdata 
         endif 
  20  continue
  10  continue
      do i=1, n-1 
        if(li(i).ge.m1) then
            ndata=ndata+1 
            l(ndata)=li(i)
            x(ndata)=c1(i)
            y(ndata)=q(i)
            ery(ndata)=erq(i) 
        endif
      enddo
      return
      end
