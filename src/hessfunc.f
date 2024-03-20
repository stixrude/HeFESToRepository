	subroutine hessfunc(nnew,hespro)

        include 'P1'
        include 'chem.inc'
        include 'const.inc'

	integer ispec,i,j
        double precision nnew(nspecp)
	double precision hess(nspecp,nspecp),hespro(nspecp,nspecp),qh(nspecp,nspecp)

C  Compute projected Hessian

        do 92 ispec=1,nspec
         call hessian(ispec,n,hess)
92      continue
        call dgemm('Transpose Q2','Normal H',nnull,nspec,nspec,one,q2,nspecp,hess,nspecp,zero,qh,nspecp)
        call dgemm('Normal H','Normal Q2',nspec,nnull,nspec,one,qh,nspecp,q2,nspecp,zero,hespro,nspecp)
c	write(31,*) 'q2 in hessfunc'
c	write(31,*) ((q2(i,j),j=1,nnull),i=1,nspec)
c	write(31,*) 'projected Hessian'
c	do 1 i=1,nnull
c	 write(31,*) (hespro(i,j),j=1,nnull)
c1	continue

	return
	end

