        subroutine dfunc(nnew,xi)

        include 'P1'
        include 'chem.inc'
	include 'const.inc'
	integer nvet,nvep
	double precision dhdpmol,dhdtmol,dvdpmol,dsdtmol,ehugo,phugo,vtarg,starg,tlast,vhugo,wmagg
        double precision nnew(nspecp)
        double precision cpa(nspecp),xi(nspecp)
	logical chcalc,adcalc,hucalc
        common /chempot/ cpa
        common /tfindc/ tlast,vtarg,starg,phugo,vhugo,ehugo,dvdpmol,dsdtmol,dhdpmol,dhdtmol,wmagg,chcalc,adcalc,hucalc,nvet,nvep

C  The projected gradient
C
C  xi = q2^T*cp
C
	call dgemv('Transpose q2',nspec,nnull,one,q2,nspecp,cpa,ione,zero,xi,ione)

	if (adcalc) xi(nnull+nvet) = dhdtmol
	if (chcalc) xi(nnull+nvep) = dhdpmol

        return
        end
