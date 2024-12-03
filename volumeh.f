	double precision function volumeh(ispec,x1)
	include 'P1'
	include 'hydrogen.inc'
	logical isochor
	integer jlop,klop,jlot,klot,i,j,it,jp
	integer, save :: ncall
	integer, parameter :: mint = 2
	double precision logTi,logPi,y,dy,x1,rho,yspline,rhospline,yy
	double precision ya(mint,mint),y2a(nt,np)
	integer ispec,jspec
	real(8) apar,Ti,Pi
	data ncall/0/
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /chor/ jspec,isochor
	ncall = ncall + 1

	logTi = log10(Ti)
	logPi = log10(Pi)

	if (ncall .eq. 1) call splie2(temph,pressh,rhoh,nt,np,y2a)
	call splin2(temph,pressh,rhoh,y2a,nt,np,logTi,logPi,yspline)

	rhospline = 10**yspline
	volumeh = hmass/rhospline

	return
	end
