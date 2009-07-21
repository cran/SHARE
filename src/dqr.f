c 
c    dqr is a fortran wrapper to convert the array to matrix and feed into dqrls
c
c
c
c

      subroutine dqr(xvec,n,p,y,ny,tol,b,rsd,qty,k,jpvt,qraux,work)
      integer n,p,ny,k,jpvt(p)
      double precision xvec(n*p),y(n,ny),tol,b(p,ny),rsd(n,ny),
     .                 qty(n,ny),qraux(p),work(p)
           
c
c     internal variables.
c
      integer i,j,count
      double precision x(n,p)

      count=1
      DO 20 i=1,n
          DO 10 j=1,p
             x(i,j)=xvec(count)
             count=count+1
   10     CONTINUE
   20 CONTINUE

      call dqrls(x,n,p,y,ny,tol,b,rsd,qty,k,jpvt,qraux,work)

      return
      END      







c
c     dqrfit is a subroutine to compute least squares solutions
c     to the system
c
c     (1)               x * b = y
c
c     which may be either under-determined or over-determined.
c     the user must supply a tolerance to limit the columns of
c     x used in computing the solution.  in effect, a set of
c     columns with a condition number approximately bounded by
c     1/tol is used, the other components of b being set to zero.
c
c     on entry
c
c        x      double precision(n,p).
c               x contains n-by-p coefficient matrix of
c               the system (1), x is destroyed by dqrfit.
c
c        n      the number of rows of the matrix x.
c
c        p      the number of columns of the matrix x.
c
c        y      double precision(n,ny)
c               y contains the right hand side(s) of the system (1).
c
c        ny     the number of right hand sides of the system (1).
c
c        tol    double precision
c               tol is the nonnegative tolerance used to
c               determine the subset of columns of x included
c               in the solution.  columns are pivoted out of
c               decomposition if
c
c        jpvt   integer(p)
c               the values in jpvt are permuted in the same
c               way as the columns of x.  this can be useful
c               in unscrambling coefficients etc.
c
c        work   double precision(2*p)
c               work is an array used by dqrdc2 and dqrsl.
c
c     on return
c
c        x      contains the output array from dqrdc2.
c               namely the qr decomposition of x stored in
c               compact form.
c
c        b      double precision(p,ny)
c               b contains the solution vectors with rows permuted
c               in the same way as the columns of x.  components
c               corresponding to columns not used are set to zero.
c
c        rsd    double precision(n,ny)
c               rsd contains the residual vectors y-x*b.
c
c        qty    double precision(n,ny)     t
c               qty contains the vectors  q y.   note that
c               the initial p elements of this vector are
c               permuted in the same way as the columns of x.
c
c        k      integer
c               k contains the number of columns used in the
c               solution.
c
c        jpvt   has its contents permuted as described above.
c
c        qraux  double precision(p)
c               qraux contains auxiliary information on the
c               qr decomposition of x.
c
c
c     on return the arrays x, jpvt and qraux contain the
c     usual output from dqrdc, so that the qr decomposition
c     of x with pivoting is fully available to the user.
c     in particular, columns jpvt(1), jpvt(2),...,jpvt(k)
c     were used in the solution, and the condition number
c     associated with those columns is estimated by
c     abs(x(1,1)/x(k,k)).
c
c     dqrfit uses the linpack routines dqrdc and dqrsl.
c
      subroutine dqrls(x,n,p,y,ny,tol,b,rsd,qty,k,jpvt,qraux,work)
      integer n,p,ny,k,jpvt(p)
      double precision x(n,p),y(n,ny),tol,b(p,ny),rsd(n,ny),
     .                 qty(n,ny),qraux(p),work(p)
c
c     internal variables.
c
      integer info,j,jj,kk
c
c     reduce x.
c
      call dqrdc2(x,n,n,p,tol,k,qraux,jpvt,work)
c
c     solve the truncated least squares problem for each rhs.
c
      if(k .gt. 0) then
         do 20 jj=1,ny
   20       call dqrsl(x,n,n,k,qraux,y(1,jj),rsd(1,jj),qty(1,jj),
     1           b(1,jj),rsd(1,jj),rsd(1,jj),1110,info)
      else
         do 30 i=1,n
            do 30 jj=1,ny
   30           rsd(i,jj) = y(i,jj)
      endif
c
c     set the unused components of b to zero.
c
      kk = k + 1
      do 50 j=kk,p
         do 40 jj=1,ny
            b(j,jj) = 0.d0
   40    continue
   50 continue
      return
      end













c
c     dqrdc2 uses householder transformations to compute the qr
c     factorization of an n by p matrix x.  a limited column
c     pivoting strategy based on the 2-norms of the reduced columns
c     moves columns with near-zero norm to the right-hand edge of
c     the x matrix.  this strategy means that sequential one
c     degree-of-freedom effects can be computed in a natural way.
c
c     i am very nervous about modifying linpack code in this way.
c     if you are a computational linear algebra guru and you really
c     understand how to solve this problem please feel free to
c     suggest improvements to this code.
c
c     Another change was to compute the rank.
c
c     on entry
c
c        x       double precision(ldx,p), where ldx .ge. n.
c                x contains the matrix whose decomposition is to be
c                computed.
c
c        ldx     integer.
c                ldx is the leading dimension of the array x.
c
c        n       integer.
c                n is the number of rows of the matrix x.
c
c        p       integer.
c                p is the number of columns of the matrix x.
c
c        tol     double precision
c                tol is the nonnegative tolerance used to
c                determine the subset of the columns of x
c                included in the solution.
c
c        jpvt    integer(p).
c                integers which are swapped in the same way as the
c                the columns of x during pivoting.  on entry these
c                should be set equal to the column indices of the
c                columns of the x matrix (typically 1 to p).
c
c        work    double precision(p,2).
c                work is a work array.
c
c     on return
c
c        x       x contains in its upper triangle the upper
c                triangular matrix r of the qr factorization.
c                below its diagonal x contains information from
c                which the orthogonal part of the decomposition
c                can be recovered.  note that if pivoting has
c                been requested, the decomposition is not that
c                of the original matrix x but that of x
c                with its columns permuted as described by jpvt.
c
c        k       integer.
c                k contains the number of columns of x judged
c                to be linearly independent.
c
c        qraux   double precision(p).
c                qraux contains further information required to recover
c                the orthogonal part of the decomposition.
c
c        jpvt    jpvt(k) contains the index of the column of the
c                original matrix that has been interchanged into
c                the k-th column.
c
c     this version dated 22 august 1995
c     ross ihaka
c
c     bug fixes 29 September 1999 BDR (p > n case, inaccurate ranks)
c
c
c     dqrdc uses the following functions and subprograms.
c
c     blas daxpy,ddot,dscal,dnrm2
c     fortran dabs,dmax1,min0,dsqrt
c
      subroutine dqrdc2(x,ldx,n,p,tol,k,qraux,jpvt,work)
      integer ldx,n,p
      integer jpvt(*)
      double precision x(ldx,*),qraux(*),work(p,2),tol
c
c     internal variables
c
      integer i,j,l,lp1,lup,k
      double precision dnrm2,tt,ttt
      double precision ddot,nrmxl,t
c
c
c     compute the norms of the columns of x.
c
      do 70 j = 1, p
         qraux(j) = dnrm2(n,x(1,j),1)
         work(j,1) = qraux(j)
         work(j,2) = qraux(j)
         if(work(j,2) .eq. 0.0d0) work(j,2) = 1.0d0
   70 continue
c
c     perform the householder reduction of x.
c
      lup = min0(n,p)
      k = p + 1
      do 200 l = 1, lup
c
c     previous version only cycled l to lup
c
c     cycle the columns from l to p left-to-right until one
c     with non-negligible norm is located.  a column is considered
c     to have become negligible if its norm has fallen below
c     tol times its original norm.  the check for l .le. k
c     avoids infinite cycling.
c
   80    continue
         if (l .ge. k .or. qraux(l) .ge. work(l,2)*tol) go to 120
            lp1 = l+1
            do 100 i=1,n
               t = x(i,l)
               do 90 j=lp1,p
                  x(i,j-1) = x(i,j)
   90          continue
               x(i,p) = t
  100       continue
            i = jpvt(l)
            t = qraux(l)
            tt = work(l,1)
            ttt = work(l,2)
            do 110 j=lp1,p
               jpvt(j-1) = jpvt(j)
               qraux(j-1) = qraux(j)
               work(j-1,1) = work(j,1)
               work(j-1,2) = work(j,2)
  110       continue
            jpvt(p) = i
            qraux(p) = t
            work(p,1) = tt
            work(p,2) = ttt
            k = k - 1
            go to 80
  120    continue
         if (l .eq. n) go to 190
c
c           compute the householder transformation for column l.
c
            nrmxl = dnrm2(n-l+1,x(l,l),1)
            if (nrmxl .eq. 0.0d0) go to 180
               if (x(l,l) .ne. 0.0d0) nrmxl = dsign(nrmxl,x(l,l))
               call dscal(n-l+1,1.0d0/nrmxl,x(l,l),1)
               x(l,l) = 1.0d0 + x(l,l)
c
c              apply the transformation to the remaining columns,
c              updating the norms.
c
               lp1 = l + 1
               if (p .lt. lp1) go to 170
               do 160 j = lp1, p
                  t = -ddot(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)
                  call daxpy(n-l+1,t,x(l,l),1,x(l,j),1)
                  if (qraux(j) .eq. 0.0d0) go to 150
                     tt = 1.0d0 - (dabs(x(l,j))/qraux(j))**2
                     tt = dmax1(tt,0.0d0)
                     t = tt
c
c modified 9/99 by BDR. Re-compute norms if there is large reduction
c The tolerance here is on the squared norm
c In this version we need accurate norms, so re-compute often.
c  work(j,1) is only updated in one case: looks like a bug -- no longer used
c
c                     tt = 1.0d0 + 0.05d0*tt*(qraux(j)/work(j,1))**2
c                     if (tt .eq. 1.0d0) go to 130
                     if (dabs(t) .lt. 1e-6) go to 130
                        qraux(j) = qraux(j)*dsqrt(t)
                     go to 140
  130                continue
                        qraux(j) = dnrm2(n-l,x(l+1,j),1)
                        work(j,1) = qraux(j)
  140                continue
  150             continue
  160          continue
  170          continue
c
c              save the transformation.
c
               qraux(l) = x(l,l)
               x(l,l) = -nrmxl
  180       continue
  190    continue
  200 continue
      k = min0(k - 1, n)
      return
      end









c
c     dqrsl applies the output of dqrdc to compute coordinate
c     transformations, projections, and least squares solutions.
c     for k .le. min(n,p), let xk be the matrix
c
c            xk = (x(jpvt(1)),x(jpvt(2)), ... ,x(jpvt(k)))
c
c     formed from columnns jpvt(1), ... ,jpvt(k) of the original
c     n x p matrix x that was input to dqrdc (if no pivoting was
c     done, xk consists of the first k columns of x in their
c     original order).  dqrdc produces a factored orthogonal matrix q
c     and an upper triangular matrix r such that
c
c              xk = q * (r)
c                       (0)
c
c     this information is contained in coded form in the arrays
c     x and qraux.
c
c     on entry
c
c        x      double precision(ldx,p).
c               x contains the output of dqrdc.
c
c        ldx    integer.
c               ldx is the leading dimension of the array x.
c
c        n      integer.
c               n is the number of rows of the matrix xk.  it must
c               have the same value as n in dqrdc.
c
c        k      integer.
c               k is the number of columns of the matrix xk.  k
c               must nnot be greater than min(n,p), where p is the
c               same as in the calling sequence to dqrdc.
c
c        qraux  double precision(p).
c               qraux contains the auxiliary output from dqrdc.
c
c        y      double precision(n)
c               y contains an n-vector that is to be manipulated
c               by dqrsl.
c
c        job    integer.
c               job specifies what is to be computed.  job has
c               the decimal expansion abcde, with the following
c               meaning.
c
c                    if a.ne.0, compute qy.
c                    if b,c,d, or e .ne. 0, compute qty.
c                    if c.ne.0, compute b.
c                    if d.ne.0, compute rsd.
c                    if e.ne.0, compute xb.
c
c               note that a request to compute b, rsd, or xb
c               automatically triggers the computation of qty, for
c               which an array must be provided in the calling
c               sequence.
c
c     on return
c
c        qy     double precision(n).
c               qy conntains q*y, if its computation has been
c               requested.
c
c        qty    double precision(n).
c               qty contains trans(q)*y, if its computation has
c               been requested.  here trans(q) is the
c               transpose of the matrix q.
c
c        b      double precision(k)
c               b contains the solution of the least squares problem
c
c                    minimize norm2(y - xk*b),
c
c               if its computation has been requested.  (note that
c               if pivoting was requested in dqrdc, the j-th
c               component of b will be associated with column jpvt(j)
c               of the original matrix x that was input into dqrdc.)
c
c        rsd    double precision(n).
c               rsd contains the least squares residual y - xk*b,
c               if its computation has been requested.  rsd is
c               also the orthogonal projection of y onto the
c               orthogonal complement of the column space of xk.
c
c        xb     double precision(n).
c               xb contains the least squares approximation xk*b,
c               if its computation has been requested.  xb is also
c               the orthogonal projection of y onto the column space
c               of x.
c
c        info   integer.
c               info is zero unless the computation of b has
c               been requested and r is exactly singular.  in
c               this case, info is the index of the first zero
c               diagonal element of r and b is left unaltered.
c
c     the parameters qy, qty, b, rsd, and xb are not referenced
c     if their computation is not requested and in this case
c     can be replaced by dummy variables in the calling program.
c     to save storage, the user may in some cases use the same
c     array for different parameters in the calling sequence.  a
c     frequently occuring example is when one wishes to compute
c     any of b, rsd, or xb and does not need y or qty.  in this
c     case one may identify y, qty, and one of b, rsd, or xb, while
c     providing separate arrays for anything else that is to be
c     computed.  thus the calling sequence
c
c          call dqrsl(x,ldx,n,k,qraux,y,dum,y,b,y,dum,110,info)
c
c     will result in the computation of b and rsd, with rsd
c     overwriting y.  more generally, each item in the following
c     list contains groups of permissible identifications for
c     a single callinng sequence.
c
c          1. (y,qty,b) (rsd) (xb) (qy)
c
c          2. (y,qty,rsd) (b) (xb) (qy)
c
c          3. (y,qty,xb) (b) (rsd) (qy)
c
c          4. (y,qy) (qty,b) (rsd) (xb)
c
c          5. (y,qy) (qty,rsd) (b) (xb)
c
c          6. (y,qy) (qty,xb) (b) (rsd)
c
c     in any group the value returned in the array allocated to
c     the group corresponds to the last member of the group.
c
c     linpack. this version dated 08/14/78 .
c     g.w. stewart, university of maryland, argonne national lab.
c
c     dqrsl uses the following functions and subprograms.
c
c     BLAS      daxpy,dcopy,ddot
c     Fortran   dabs,min0,mod
c
      subroutine dqrsl(x,ldx,n,k,qraux,y,qy,qty,b,rsd,xb,job,info)
      integer ldx,n,k,job,info
      double precision x(ldx,*),qraux(*),y(*),qy(*),qty(*),b(*),rsd(*),
     *                 xb(*)
c
c     internal variables
c
      integer i,j,jj,ju,kp1
      double precision ddot,t,temp
      logical cb,cqy,cqty,cr,cxb
c
c
c     set info flag.
c
      info = 0
c
c     determine what is to be computed.
c
      cqy = job/10000 .ne. 0
      cqty = mod(job,10000) .ne. 0
      cb = mod(job,1000)/100 .ne. 0
      cr = mod(job,100)/10 .ne. 0
      cxb = mod(job,10) .ne. 0
      ju = min0(k,n-1)
c
c     special action when n=1.
c
      if (ju .ne. 0) go to 40
         if (cqy) qy(1) = y(1)
         if (cqty) qty(1) = y(1)
         if (cxb) xb(1) = y(1)
         if (.not.cb) go to 30
            if (x(1,1) .ne. 0.0d0) go to 10
               info = 1
            go to 20
   10       continue
               b(1) = y(1)/x(1,1)
   20       continue
   30    continue
         if (cr) rsd(1) = 0.0d0
      go to 250
   40 continue
c
c        set up to compute qy or qty.
c
         if (cqy) call dcopy(n,y,1,qy,1)
         if (cqty) call dcopy(n,y,1,qty,1)
         if (.not.cqy) go to 70
c
c           compute qy.
c
            do 60 jj = 1, ju
               j = ju - jj + 1
               if (qraux(j) .eq. 0.0d0) go to 50
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  t = -ddot(n-j+1,x(j,j),1,qy(j),1)/x(j,j)
                  call daxpy(n-j+1,t,x(j,j),1,qy(j),1)
                  x(j,j) = temp
   50          continue
   60       continue
   70    continue
         if (.not.cqty) go to 100
c
c           compute trans(q)*y.
c
            do 90 j = 1, ju
               if (qraux(j) .eq. 0.0d0) go to 80
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  t = -ddot(n-j+1,x(j,j),1,qty(j),1)/x(j,j)
                  call daxpy(n-j+1,t,x(j,j),1,qty(j),1)
                  x(j,j) = temp
   80          continue
   90       continue
  100    continue
c
c        set up to compute b, rsd, or xb.
c
         if (cb) call dcopy(k,qty,1,b,1)
         kp1 = k + 1
         if (cxb) call dcopy(k,qty,1,xb,1)
         if (cr .and. k .lt. n) call dcopy(n-k,qty(kp1),1,rsd(kp1),1)
         if (.not.cxb .or. kp1 .gt. n) go to 120
            do 110 i = kp1, n
               xb(i) = 0.0d0
  110       continue
  120    continue
         if (.not.cr) go to 140
            do 130 i = 1, k
               rsd(i) = 0.0d0
  130       continue
  140    continue
         if (.not.cb) go to 190
c
c           compute b.
c
            do 170 jj = 1, k
               j = k - jj + 1
               if (x(j,j) .ne. 0.0d0) go to 150
                  info = j
c           ......exit
                  go to 180
  150          continue
               b(j) = b(j)/x(j,j)
               if (j .eq. 1) go to 160
                  t = -b(j)
                  call daxpy(j-1,t,x(1,j),1,b,1)
  160          continue
  170       continue
  180       continue
  190    continue
         if (.not.cr .and. .not.cxb) go to 240
c
c           compute rsd or xb as required.
c
            do 230 jj = 1, ju
               j = ju - jj + 1
               if (qraux(j) .eq. 0.0d0) go to 220
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  if (.not.cr) go to 200
                     t = -ddot(n-j+1,x(j,j),1,rsd(j),1)/x(j,j)
                     call daxpy(n-j+1,t,x(j,j),1,rsd(j),1)
  200             continue
                  if (.not.cxb) go to 210
                     t = -ddot(n-j+1,x(j,j),1,xb(j),1)/x(j,j)
                     call daxpy(n-j+1,t,x(j,j),1,xb(j),1)
  210             continue
                  x(j,j) = temp
  220          continue
  230       continue
  240    continue
  250 continue
      return
      end







