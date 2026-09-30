	double precision function volumeh(ispec,x1)
	include 'P1'
	include 'hydrogen.inc'
	logical isochor
	integer, save :: ncall
	double precision logTi,x1,pressureh
	double precision p1,p2,plow,pupp
	double precision rholow,rhoupp,vlow,vupp,x2,zeroin
	integer ispec,jspec,ires,itemp,jpress,jpresslow,jpressupp
	double precision apar,Ti,Pi
	double precision, parameter :: logTimin = 2.d0, dlogTi = 0.05d0, ten=10.d0
	data ncall/0/
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /chor/ jspec,isochor
	external pressureh
	double precision, parameter :: tol=1.e-15
	ncall = ncall + 1

	logTi = log10(Ti)
C Set bounds on the volume
	rholow = -8.
	rhoupp = 2.*(logTi - 3.3)		!	chabrieretal_19 Eq. 3
	if (logTi .le. 3.15) then 
	 rhoupp = (logTi - 3.9166666)/2.55555	!	My fit to bound on the H2 melt curve (chabrieretal_19 Eq. 1) in rho-T space
	end if
	rhoupp = +6.				!       Full range of chabrieretal_19 table
	vupp = hmass/ten**rholow
	vlow = hmass/ten**rhoupp

	itemp = min(nt,max(1,nint((logTi - logTimin)/dlogTi) + 1))
	do 1 jpress=1,np
c	 if (pressh(itemp,jpress) .gt. log10(Pi)) go to 10
	 vlow = hmass/10**rhoh(jpress)
	 if (pressureh(vlow) .gt. 0.d0) go to 10
1	continue
10	continue
	jpresslow = max(1,jpress-1)
	jpressupp = min(jpress,np)
	rholow = rhoh(jpresslow)
	rhoupp = rhoh(jpressupp)
	vupp = hmass/10**rholow
	vlow = hmass/10**rhoupp
	if (vlow .gt. x1 .or. vupp .lt. x1) x1 = (vupp + vlow)/2.d0
c	print*, 'in volumeh', itemp,jpress,rhoh(jpresslow),rhoh(jpressupp),vlow,vupp

	x2 = x1*(1. + 1.e10*tol)
        p1 = pressureh(x1)
        p2 = pressureh(x2)
        plow = pressureh(vlow)
        pupp = pressureh(vupp)
c	print*, 'volumeh',x1,x2,vlow,vupp,p1,p2,plow,pupp
	call cage(pressureh,x1,x2,vlow,vupp,ires)
	if (ires .eq. 1) then
	 volumeh = zeroin(x1,x2,pressureh,tol)
	 return
	end if

	print*, 'volumeh failed to find V cage',ispec,Pi,Ti,x1,x2,rholow,rhoupp,vlow,vupp,plow,pupp,p1,p2

	return
	end
