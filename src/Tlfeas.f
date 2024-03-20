	subroutine Tlfeas(Tfeas,nnew)

	include 'P1'
	include 'chem.inc'
	logical chcalc,adcalc,hucalc
	integer iter,int,nvet,nvep
	double precision tlast,vtarg,starg,phugo,vhugo,ehugo,dvdpmol,dsdtmol,dhdpmol,dhdtmol,wmagg
	double precision Tfeas,val,nnew(nspecp),func,To,x(nspecp),sfunc
	double precision Tspin
	double precision, parameter :: Tsmall=1.e-5
	integer, parameter :: itermax=5
        common /tfindc/ tlast,vtarg,starg,phugo,vhugo,ehugo,dvdpmol,dsdtmol,dhdpmol,dhdtmol,wmagg,chcalc,adcalc,hucalc,nvet,nvep

	if (nnew(nnull+nvet) .le. Tsmall) then
	 write(31,*) 'Resetting To in Tlfeas',nnew(nnull+nvet),Tsmall
	 To = 1000.
	 nnew(nnull+nvet) = To
	end if
	iter = 0
10	continue
	To = nnew(nnull+nvet)
	val = func(nnew)

	Tfeas = log(To) + (starg*wmagg - sfunc(x))/(To*dsdtmol)
	Tfeas = min(exp(Tfeas),Tspin(int))
c	write(31,*) 'Tfeas = ',Tfeas,To,sfunc(x),starg,wmagg,starg*wmagg,dsdtmol
	iter = iter + 1
	if (iter .gt. itermax) then
	 write(31,*) 'WARNING: Tlfeas not converging',Tfeas,nnew(nnull+nvet)
	 return
	end if
	if (abs(Tfeas - To) .gt. Tsmall) then
	 nnew(nnull+nvet) = Tfeas
	 go to 10
	end if

	return
	end
