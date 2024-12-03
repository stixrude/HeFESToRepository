        subroutine newfrm(q2,n,n1,nnew,nspec,nnull)

        include 'P1'
        include 'const.inc'

	integer nspec,nnull,i,info
	double precision vsum
	logical valid
	integer ipiv(nspecp)
        double precision n(nspecp),n1(nspecp)
        double precision q2(nspecp,nspecp)
        double precision nnew(nspecp),ndiff(nspecp)
        double precision u(nspecp,nspecp),bb(nspecp)
        double precision bbx(nspecp,nspecp)
	double precision, parameter :: fac=2., tol=1.e-5
	integer, parameter :: itermx=9

c	write(31,*) 'begin newfrm',nnull
c        write(31,'(71f9.5)') (n(i),i=1,nspec)
        do 1 i=1,nspec
         ndiff(i) = n(i) - n1(i)
1       continue
	do 11 i=1,nnull
	 nnew(i) = 0.
11	continue

!  Solve linear problem q2^T*q2*nnew = q2^T*(n-n1)
        call dgemv('Transpose q2',nspec,nnull,one,q2,nspecp,ndiff,ione,zero,bb,ione)
        call dgemm('Transpose q2','N',nnull,nnull,nspec,one,q2,nspecp,q2,nspecp,zero,u,nspecp)
        do 3 i=1,nnull
         bbx(i,1) = bb(i)
3       continue
        call dgesv(nnull,ione,u,nspecp,ipiv,bbx(1,1),nspecp,info)
        do 2 i=1,nnull
         nnew(i) = bbx(i,1)
2       continue

	if (.not. valid(vsum,nnew)) then
	 write(31,*) 'WARNING: newfrm ends with invalid n-vector',nnull,vsum
         write(31,'(71f9.5)') (n(i),i=1,nspec)
	end if

        call nform(nnew,n,n1,q2,nspec,nnull)

	return
	end
