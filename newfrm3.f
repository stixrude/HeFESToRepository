        subroutine newfrm3(s,dn,n,nspec,nc,ncs)

        include 'P1'
        include 'const.inc'

	integer nspec,nc,ncs,i,info,irow,ispec,j,nco,nnull,nnulls
	double precision s(ncompp,nspecp)
	double precision dn(nspecp),nnew3(nspecp),n(nspecp)
        double precision d(nspecp,nspecp),q11(nspecp,nspecp)
        double precision u(nspecp,nspecp),v(nspecp,nspecp),bb(nspecp)
	double precision btest(ncompp)
	double precision st(nspecp,ncompp),tau(nspecp),work(nspecp)
        double precision, parameter :: small = 1.e-5
	nco = nc - ncs
	nnull = nspec - nc

c        write(31,*) 'Enter newfrm3',nspec,nc,ncs,nco
c	write(31,'(a13,71f9.5)') 'n',(n(i),i=1,nspec)

        call transpose(nc,nspec,ncompp,nspecp,s,st)
        call DGEQRF( nspec, nspec, st, nspecp, TAU, WORK, nspecp, INFO )
        call DORGQR( nspec, nspec, nspec, st, nspecp, TAU, WORK, nspecp, INFO )

        do 2 i=1,nspecp
         nnew3(i) = 0.
	 bb(i) = 0.
         do 2 j=1,nspecp
          d(i,j) = 0.
2       continue
	irow = 0
        do 21 i=nco+1,nc
	 irow = irow + 1
	 do 23 ispec=1,nspec
	  d(irow,ispec) = s(i,ispec)
	  q11(ispec,irow) = st(ispec,i)
23	 continue
21      continue
c	do 24 i=1,ncs
c        write(31,'(a13,71f9.5)') 'dn',(dn(i),i=1,nspec)

!  Solve linear problem D*dn = D*q1*nnew3
!  where dn contains only the species to be added, and D contains non-zero elements corresponding all to absent species.
!  cf. Harvie et al. (1987) eq. 61.

        call dgemv('Normal d',ncs,nspec,one,d,nspecp,dn,ione,zero,bb,ione)
        call dgemm('Normal d','Normal q11',ncs,ncs,nspec,one,d,nspecp,q11,nspecp,zero,u,nspecp)
	call svdsub(ncs,ncs,u,nspecp,nspecp,bb,v,v,nnew3,nnulls)
	if (nnulls .gt. 0) write(31,*) 'WARNING: ill-conditioned matrix in newfrm3',nnulls

	do 22 ispec=1,nspec
	 do 22 j=1,ncs
	  n(ispec) = n(ispec) + q11(ispec,j)*nnew3(j)
22	 continue
	continue

	do 27 j=1,nco
	 btest(j) = 0.
	 do 27 ispec=1,nspec
	  btest(j) = btest(j) + s(j,ispec)*n(ispec)
27	continue
	  
	do 26 ispec=1,nspec
26	dn(ispec) = 0.

c	write(31,'(a13,71f9.5)') 'n',(n(i),i=1,nspec)

        return
        end
