*
*---  sample program illustrating the computation of small matrix
*     exponentials in full with Expokit. Refer to the Expokit
*     documentation for more details about the methods, and 
*     especially the domain of applicability of the Chebyshev scheme.
*
      implicit none
      integer m,i,j,k,ideg,mprint,lda,ldh,lwsp,liwsp,ns,iflag,iexp
      parameter ( ideg=6, lda=50, ldh=lda )
      parameter ( lwsp=4*ldh*ldh+ideg+1, liwsp=ldh )
*
      integer iwsp(liwsp), iseed(4)
      double precision t, A(lda,lda), H(ldh,ldh), y(ldh)
      double precision wsp(lwsp), s1, s2
      complex*16 Hc(ldh,ldh), yc(ldh), wspc(lwsp)
*
      intrinsic CMPLX, CONJG, MIN, DBLE, IMAG

*-------------------------
*     REAL CASE
*-------------------------

*---  set A = random symmetric negative define matrix ...
      t = 1.0d0
      m = 2
      do j = 1,m
         do i = j,m
            A(i,j) = 0.6
            A(j,i) = A(i,j)
         enddo
*         A(j,j) = -2.5d0 + A(j,j)
      enddo

*---  maximum number of rows/columns to be printed
      mprint = MIN(2,m)

      print*,"t ="
      print*,t
      print 9000,"REAL SYMMETRIC CASE","*******************************"
      print*,"A = "
      print 9001,( (A(i,j), j=1,mprint), i=1,mprint )
*
 9000 format( /,A,/,A )
 9001 format( 2(1X,D11.4) )
*---  Some compliers (e.g., g77) generate 'Unsupported FORMAT specifier'
*     with the specification above. In this case, simply use this form:
* 9001 format( 5(1X,D11.4) )


*---  Pade ...
      call DGPADM( ideg, m, t, A,lda, wsp,lwsp, iwsp,iexp,ns, iflag )
      print 9000,"With DGPADM:","exp(t*A) ="
      print 9001,( (wsp(iexp+(j-1)*m+i-1), j=1,mprint), i=1,mprint )

      call DSPADM( ideg, m, t, A,lda, wsp,lwsp, iwsp,iexp,ns, iflag )
      print 9000,"With DSPADM:","exp(t*A) ="
      print 9001,( (wsp(iexp+(j-1)*m+i-1), j=1,mprint), i=1,mprint )

*---  Chebyshev
      do i = 1,m
         y(i) = 0.0d0
      enddo
      y(1) = 1.0d0
      call DGCHBV( m,t , A,lda, y, wsp, iwsp, iflag )
      print 9000,"With DGCHBV:","exp(t*A)e_1 ="
      do i = 1,mprint
         print*,y(i)
      enddo

      do i = 1,m
         y(i) = 0.0d0
      enddo
      y(1) = 1.0d0
      call DSCHBV( m,t , A,lda, y, wsp, iwsp, iflag )
      print 9000,"With DSCHBV:","exp(t*A)e_1 ="
      do i = 1,mprint
         print*,y(i)
      enddo

*---  set H = upper Hessenberg part of A ...
      do j = 1,m
         do i = 1,MIN(j+1,m)
            H(i,j) = A(i,j)
         enddo
         do k = i,m
            H(k,j) = 0.0d0
         enddo
      enddo

      print 9000,"REAL UPPER HESSENBERG CASE","************************"
      print*,"H ="
      print 9001,( (H(i,j), j=1,mprint), i=1,mprint )

*---  Pade ...
      call DGPADM( ideg, m, t, H,ldh, wsp,lwsp, iwsp,iexp,ns, iflag )
      print 9000,"With DGPADM:","exp(t*H) ="
      print 9001,( (wsp(iexp+(j-1)*m+i-1), j=1,mprint), i=1,mprint )

*---  Chebyshev
      do i = 1,m
         y(i) = 0.0d0
      enddo
      y(1) = 1.0d0
      call DNCHBV( m,t, A,lda, y, wsp, iflag )
      print 9000,"With DNCHBV:","exp(t*A)e_1 ="
      do i = 1,mprint
         print*,y(i)
      enddo


      END
