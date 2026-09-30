	double precision function pressureh(Vi)
	include 'P1'
	include 'hydrogen.inc'
	logical isochor
	integer, save :: ncall
	double precision logTi,logRi,Ri,Vi
	double precision Fi,Pispline
	double precision y2at(np,nt),yt,dydxt,d2ydx2t
	integer jspec
	double precision apar,Ti,Pi
	double precision, parameter :: ten=10.d0
	data ncall/0/
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /chor/ jspec,isochor
	ncall = ncall + 1

	logTi = log10(Ti)
	Ri = hmass/Vi
	logRi = log10(Ri)

        if (ncall .eq. 1) call splie2(rhoh,temph,Fhtranspose,np,nt,y2at)
        call splin2(rhoh,temph,Fhtranspose,y2at,np,nt,logRi,logTi,yt,dydxt,d2ydx2t)
	Fi = ten**yt
        Pispline = Ri*Fi*dydxt

	pressureh = Pispline - Pi
c	print*, 'pressureh',Vi,Pi,Pispline

	return
	end
