*
*---  sample program illustrating the use of DGEXPV ...
*     Non-symmetric problem (Example 6.3 in the Expokit report) ...
*
      implicit none
      external dgcoov

      double precision tic, tac, clock

*---  matrix data ... 
*---  BEWARE: these values should match those in dgmatv.f
      integer n, nz, nmax, nzmax
      parameter( nmax = 5000, nzmax = 600000 )
      integer ia(nzmax), ja(nzmax)
      double precision a(nzmax)
      common /RMAT/ a, ia, ja, nz, n

*---  arguments variables ...
      integer m, mmax, lwsp, liwsp
      parameter( mmax = 3 )
      parameter( lwsp = nmax*(mmax+2)+5*(mmax+2)**2+7, liwsp = nmax )
      integer iwsp(liwsp)
      double precision t, tol, anorm
      double precision v(nmax), w(nmax), wsp(lwsp)

      integer i, itrace, iflag
      double precision ZERO, ONE
      parameter( ZERO=0.0d0, ONE=1.0d0 )

      intrinsic ABS
*
*---  Executable statements ...

*---  load a Harwell-Boeing matrix ...
      n = 9
      nz = 6
      
      ia(1) = 1
      ia(2) = 1
      ia(3) = 2
      ia(4) = 2
      ia(5) = 3
      ia(6) = 3
      ja(1) = 1
      ja(2) = 3
      ja(3) = 1
      ja(4) = 2
      ja(5) = 1
      ja(6) = 3
      a(1) = 0.6
      a(2) = 0.2
      a(3) = 0.6
      a(4) = 0.6
      a(5) = 0.6
      a(6) = 0.6

      do i = 1,nz
         print*, ia(i)
      enddo
      do i = 1,nz
         print*, ja(i)
      enddo
      do i = 1,nz
         print*, a(i)
      enddo


*---  compute the infinite norm of A ...
      do i = 1,n
         wsp(i) = ZERO
      enddo
      do i = 1,nz
         wsp(ia(i)) = wsp(ia(i)) + ABS( a(i) )
      enddo
      anorm = wsp(1)
      do i = 2,n
         if ( anorm.lt.wsp(i) ) anorm =  wsp(i)
      enddo
      write(UNIT=*,FMT='(A,E8.2)') '||A||_inf= ',anorm

*---  the operand vector v is set to (1, ..., 1)^T ...
      v(1) = ONE
      do i = 2,n
         v(i) = ZERO
      enddo


*---  set other input arguments ...
      t = 1.0d0
      tol = 1000.0d0
      m = 3
      itrace = 0

*---  compute w = exp(t*A)v with COO format ...
      tic = clock()
      call DGEXPV( n, m, t,v,w, tol, anorm,
     .             wsp,lwsp, iwsp,liwsp, dgcoov, itrace, iflag ) 
      tac = clock()

*---
      print 9001,'----------------------------------------------------'
      print 9001,'DGEXPV (COO) has completed:'
      print 9001,'----------------------------------------------------'
      print 9001,'w(1:4) ='
      do i = 1,10
         print*,w(i)
      enddo

*---  display some statistics if desired ...
      print 9001,'final report----------------------------------------'
      print 9002,'runtime   = ',tac-tic
      print 9002,'||A||_inf = ',anorm
      print 9003,'nz        =',nz
      print 9003,'n         =',n
      print 9003,'m         =',m
      print 9003,'itrace    =',itrace
      print 9003,'iflag     =',iflag
      print 9003,'ibrkflag  =',iwsp(6)
      print 9003,'mbrkdwn   =',iwsp(7)
      print 9003,'nstep     =',iwsp(4)
      print 9003,'nreject   =',iwsp(5)
      print 9003,'nmult     =',iwsp(1)
      print 9003,'nexph     =',iwsp(2)
      print 9003,'nscale    =',iwsp(3)

      print 9002,'tol       = ',tol
      print 9002,'t         = ',t
      print 9002,'tbrkdwn   = ',wsp(7)
      print 9002,'step_min  = ',wsp(1)
      print 9002,'step_max  = ',wsp(2)
      print 9002,'max_round = ',wsp(3)
      print 9002,'sum_round = ',wsp(4)
      print 9002,'max_error = ',wsp(5)
      print 9002,'sum_error = ',wsp(6)
      print 9002,'hump      = ',wsp(9)
      print 9002,'scale-norm= ',wsp(10)

*----------------------------------------------------------------------|
*----------------------------------------------------------------------|

 9001 format(A)
 9002 format(A,E8.2)
 9003 format(A,I9)
 9004 format( 4(1X,D11.4) )
      END
*----------------------------------------------------------------------|




