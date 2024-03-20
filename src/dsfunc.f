	subroutine dsfunc(x,sgrad)
	include 'P1'
        include 'chem.inc'
        include 'const.inc'
	logical chcalc,adcalc,hucalc
	integer i,nvet,nvep
	double precision x,dhdpmol,dhdtmol,dvdpmol,dsdtmol,ehugo,phugo,vtarg,starg,tlast,vhugo,wmagg,entagg,sfunc
	double precision sgrad(nspecp)
        common /tfindc/ tlast,vtarg,starg,phugo,vhugo,ehugo,dvdpmol,dsdtmol,dhdpmol,dhdtmol,wmagg,chcalc,adcalc,hucalc,nvet,nvep

	call dgemv('Transpose q2',nspec,nnull,one,q2,nspecp,sspeca,ione,zero,sgrad,ione)
	sgrad(nnull+nvet) = dsdtmol
	entagg = sfunc(x)

	do 1 i=1,nnull+nvet
	 sgrad(i) = -sgrad(i)/wmagg*exp(entagg/wmagg)
c	 sgrad(i) = -sgrad(i)/wmagg
c	 sgrad(i) = -sgrad(i)
1	continue

	return
	end
