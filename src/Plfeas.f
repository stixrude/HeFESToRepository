	subroutine Plfeas(Pfeas,nnew)
C  Find a feasible pressure to initiate minimization
C  vtarg is the target density (g/cm^3)

	include 'P1'
	include 'chem.inc'
	logical chcalc,adcalc,hucalc
	integer iter,int,nvet,nvep,ispec
	double precision tlast,vtarg,starg,phugo,vhugo,ehugo,dvdpmol,dsdtmol,dhdpmol,dhdtmol,wmagg
	double precision Pfeas,val,nnew(nspecp),func,Po,x(nspecp),vfunc,pressure,pressurel,apar,Pi,Ti
	double precision, parameter :: Psmall=1.e-5
	integer, parameter :: itermax=5
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /tfindc/ tlast,vtarg,starg,phugo,vhugo,ehugo,dvdpmol,dsdtmol,dhdpmol,dhdtmol,wmagg,chcalc,adcalc,hucalc,nvet,nvep

	iter = 0
10	continue
C  If there is only one species, then find the feasible pressure by calculating it directly from pressure(l)
	Po = nnew(nnull+nvep)
	val = func(nnew)
	if (nspec .eq. 1) then
	 ispec = 1
	 Pi = 0.
c	 write(31,*) 'Calling pressure(l) from Plfeas',Pi,Ti,wmagg,vtarg,wmagg
	 if (apar(ispec,31) .eq. 0) then
	  Po = pressure(apar(ispec,3)/vtarg)
	 else 
	  Po = pressurel(apar(ispec,3)/vtarg)
	 end if
	 nnew(nnull+nvep) = Po
	 val = func(nnew)
	end if

C  If there are multiple species find feasible pressure iteratively using calcuated value of dvdpmol of the assemblage
	Pfeas = Po + (wmagg/vtarg - vfunc(x))/dvdpmol
c	write(31,*) 'Pfeas = ',Pfeas,Po,vfunc(x),vtarg,wmagg,wmagg/vtarg,dvdpmol
	iter = iter + 1
	if (iter .gt. itermax) then
	 write(31,*) 'WARNING: Plfeas not converging',Pfeas,nnew(nnull+nvep)
	 return
	end if
	if (abs(Pfeas - Po) .gt. Psmall) then
	 nnew(nnull+nvep) = Pfeas
	 go to 10
	end if

	return
	end
