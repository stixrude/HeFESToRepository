        subroutine svdsub(m,n,s,msp,nsp,b,q1,q2,n1,nnulls)
C  Use LAPACK routines to perform the SVD decomposition of matrix s
        include 'P1'
        include 'const.inc'

	integer m,n,msp,nsp,nnulls,i,info,iprt,j,k,k1,k2,nc,nnull,nspec
        integer icol1(nsp),icol2(nsp)
        double precision s(msp,nsp)
        double precision b(msp)
        double precision q1(nsp,nsp),q2(nsp,nsp),n1(nsp)
        double precision u(msp,msp),w(nsp),v(nsp,nsp)
	double precision work(5*nsp),rw(nsp),y(nsp),swork(msp,nsp)
	double precision start,finish,gsvd
        double precision, parameter :: small=1.e-6
	data gsvd/0.0/
	call cpu_time(start)
        iprt = 0
C-->  Following assignments follow if s is the stoichiometric coefficient matrix
        nc = m
        nspec = n
C<---
        nnull = n - m
	do 121 i=1,msp
	 do 121 j=1,nsp
	  w(j) = 0.
	  swork(i,j) = s(i,j)
	  n1(j) = 0.
	  y(j) = 0.
121	continue

c       write(31,*) 's'
c       do 16 i=1,m
c        write(31,'(99f5.2)') (swork(i,j),j=1,n)
c16     continue

	call dgesvd('A','A',m,n,swork,msp,w,u,msp,v,nsp,work,5*nsp,info)
        k1 = 0
        k2 = 0
        do 13 k=1,n
c	 write(31,*) 'k,w,small',k,w(k),small
         if (w(k) .lt. small) then
          w(k) = 0.
	  rw(k) = 0.
          k2 = k2 + 1
          icol2(k2) = k
         else
	  rw(k) = 1./w(k)
          k1 = k1 + 1
          icol1(k1) = k
         end if
13      continue
clapack
	k = 0
	do 131 j=n,n-k2+1
	 k = k + 1
	 icol2(k) = j
131	continue
	do 141 j=1,n-k2
	 icol1(j) = j
141	continue
c solve linear problem: x = V*w*U^T*b
c y=U^T*b
	call dgemv('T',m,m,one,u,msp,b,ione,zero,y,ione)
c y_j*1/w_j with small values of 1/w_j having been set to zero
	do 151 i=1,n
	 y(i) = y(i)*rw(i)
151	continue
c  V*y note that dgesvd returned V^T
	call dgemv('T',n,n,one,v,nsp,y,ione,zero,n1,ione)
clapack
        if (k1 .ne. m) write(31,*) 'WARNING: matrix is ill-conditioned',k1,k2,m,n,nnull
        nnulls = k2

        do 14 j=1,k1
         do 14 i=1,n
clapack returns V^T
          q1(i,j) = v(icol1(j),i)
14      continue

        do 15 j=1,k2
         do 15 i=1,n
clapack returns V^T
          q2(i,j) = v(icol2(j),i)
15      continue

	call cpu_time(finish)
	gsvd = gsvd + finish - start
c	write(31,*) 'SVD time',gsvd

        return
        end
