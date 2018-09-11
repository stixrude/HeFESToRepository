	subroutine hessfunc(nnew,hespro)

        include 'P1'
        include 'chem.inc'
        include 'const.inc'

	integer ispec
        double precision nnew(nspecp)
	double precision hess(nspecp,nspecp),hespro(nspecp,nspecp),qh(nspecp,nspecp)

C  Compute projected Hessian

        do 92 ispec=1,nspec
         call hessian(ispec,n,hess)
92      continue
        call dgemm('Transpose Q2','Normal H',nnull,nspec,nspec,one,q2,nspecp,hess,nspecp,zero,qh,nspecp)
        call dgemm('Normal H','Normal Q2',nspec,nnull,nspec,one,qh,nspecp,q2,nspecp,zero,hespro,nspecp)

	return
	end

