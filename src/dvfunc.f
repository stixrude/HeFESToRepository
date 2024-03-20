	subroutine dvfunc(x,vgrad)
	include 'P1'
        include 'chem.inc'
        include 'const.inc'
	logical chcalc,adcalc,hucalc
	integer i,nvet,nvep
	double precision x,dhdpmol,dhdtmol,dvdpmol,dsdtmol,ehugo,phugo,vtarg,starg,tlast,vhugo,wmagg,volagg,vfunc
	double precision vgrad(nspecp)
        common /tfindc/ tlast,vtarg,starg,phugo,vhugo,ehugo,dvdpmol,dsdtmol,dhdpmol,dhdtmol,wmagg,chcalc,adcalc,hucalc,nvet,nvep

	call dgemv('Transpose q2',nspec,nnull,one,q2,nspecp,vspeca,ione,zero,vgrad,ione)
	vgrad(nnull+nvep) = dvdpmol
	volagg = vfunc(x)

	do 1 i=1,nnull+nvep
	 vgrad(i) = vgrad(i)/wmagg
1	continue

	return
	end
