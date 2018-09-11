        subroutine qcalc(nnew,qual)

        include 'P1'
        include 'chem.inc'

	integer i,nvet,nvep
	double precision qual,dhdpmol,dhdtmol,dvdpmol,dsdtmol,ehugo,fdumm,phugo,vtarg,starg,tlast,vhugo,func,wmagg
	logical chcalc,adcalc,hucalc,adcalcsave
        double precision nnew(nspecp),xi(nspecp)
        common /tfindc/ tlast,vtarg,starg,phugo,vhugo,ehugo,dvdpmol,dsdtmol,dhdpmol,dhdtmol,wmagg,chcalc,adcalc,hucalc,nvet,nvep
	adcalcsave = adcalc
c	adcalc = .false.
	
        fdumm = func(nnew)

        call dfunc(nnew,xi)
        qual = 0.0
        do 5 i=1,nnull
5       qual = qual + xi(i)**2
        if (nnull .gt. 0) qual = sqrt(qual/float(nnull))

	adcalc = adcalcsave

        return
        end
