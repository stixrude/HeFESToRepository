           subroutine myconstraint(val, ndim, x, grad, need_gradient, d)
	    include 'P1'
	    include 'chem.inc'
	    include 'const.inc'

	integer ndim,i,ispec,jconstr,nvet,nvep
	double precision dhdpmol,dhdtmol,dvdpmol,dsdtmol,ehugo,eps,phugo,vtarg,starg,tlast,vhugo,sfunc,vfunc,Tsmall,wmagg
            integer need_gradient
	    logical chcalc,adcalc,hucalc
            common /tfindc/ tlast,vtarg,starg,phugo,vhugo,ehugo,dvdpmol,dsdtmol,dhdpmol,dhdtmol,wmagg,chcalc,adcalc,hucalc,nvet,nvep
            double precision val, x(ndim), grad(ndim), d(nspecp)

           call nform(x,n,n1,q2,nspec,nnull)
	    eps = 1.e-8
	    Tsmall = 1.e-5
	    val = eps
	    do 1 i=1,nspec
1	    val = val + d(i)*n(i)
	    if (d(nspec+nvet) .ne. 0. .and. nvet .ne. 0) then
	     if (adcalc) val = -(exp(sfunc(x)/wmagg) - exp(starg))
c	     if (adcalc) val = -(sfunc(x)/wmagg - starg)
c	     if (adcalc) val = -(sfunc(x) - starg*wmagg)
	    end if
	    if (d(nspec+nvep) .ne. 0. .and. nvep .ne. 0) then
c	     if (chcalc) val = vfunc(x) - wmagg/vtarg
	     if (chcalc) val = vfunc(x)/wmagg - 1./vtarg
	    end if
            if (need_gradient.ne.0) then
             call dgemv('Transpose q2',nspec,nnull,one,q2,nspecp,d,ione,zero,grad,ione)
	     if (adcalc .and. nvet .ne. 0) grad(nnull+nvet) = 0.
	     if (chcalc .and. nvep .ne. 0) grad(nnull+nvep) = 0.
	     if (d(nspec+nvet) .ne. 0. .and. nvet .ne. 0) then
	      if (adcalc) call dsfunc(x,grad)
	     end if
	     if (d(nspec+nvep) .ne. 0. .and. nvep .ne. 0) then
	      if (chcalc) call dvfunc(x,grad)
	     end if
            endif
	    jconstr = 0
	    do 4 ispec=1,nspec+1
4	    if (d(ispec) .ne. 0.) jconstr = ispec
c	write(31,'(a8,i5,999e12.5)') 'myconstr',jconstr,x(ndim),val,grad,(n(i),i=1,nspec),starg,wmagg/vtarg
	    return
            end
