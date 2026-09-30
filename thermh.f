	subroutine thermh(ispec,Vi,volnl,Cp,Cv,gamma,K,Ks,alp,Ftot,ph,ent,deltas,tcal,
     &               zeta,Gsh,uth,uto,thet,q,etas,dGdT,pzp,Sel,Eel,Pel,Eig,Cvel,Pig,P,E)

	include 'P1'
	include 'hydrogen.inc'
        double precision Vi,volnl,Cp,Cv,gamma,K,Ks,alp,Ftot,ph,ent,deltas,
     &                   tcal,zeta,Gsh,uth,uto,thet,q,etas,dGdT,pzp,Sel,Eel,Pel,Eig,Cvel,Pig,P,E
        double precision apar,Ti,Pi,Fpv
	integer ispec
	integer, save :: ncall
        double precision logTi,logRi,y,x1,Ri,pressurecheck,yt,dydx,d2ydx2,dydxt,d2ydx2t
	double precision Fi,dlnsdlnpi
        double precision y2a(nt,np),y2at(np,nt)
	double precision yup,ydn,ytup,ytdn,dydxup,dydxdn,dydxtup,dydxtdn,dpdtt
	double precision, parameter :: deltaT = 1.d-3, deltar = 1.d-3, ten = 10.d0
	data ncall/0/
        common /state/ apar(nspecp,nparp),Ti,Pi
	ncall = ncall + 1
	ispec = 1
	x1 = 1.0

        logTi = log10(Ti)
        Ri = hmass/Vi
        logRi = log10(Ri)

        if (ncall .eq. 1) call splie2(temph,rhoh,Fh,nt,np,y2a)
        call splin2(temph,rhoh,Fh,y2a,nt,np,logTi,logRi+deltar,yup,dydxup,d2ydx2)
        call splin2(temph,rhoh,Fh,y2a,nt,np,logTi,logRi-deltar,ydn,dydxdn,d2ydx2)
        call splin2(temph,rhoh,Fh,y2a,nt,np,logTi,logRi,y,dydx,d2ydx2)

        if (ncall .eq. 1) call splie2(rhoh,temph,Fhtranspose,np,nt,y2at)
        call splin2(rhoh,temph,Fhtranspose,y2at,np,nt,logRi,logTi+deltaT,ytup,dydxtup,d2ydx2t)
        call splin2(rhoh,temph,Fhtranspose,y2at,np,nt,logRi,logTi-deltaT,ytdn,dydxtdn,d2ydx2t)
        call splin2(rhoh,temph,Fhtranspose,y2at,np,nt,logRi,logTi,yt,dydxt,d2ydx2t)

        dpdtt = (dydxtup - dydxtdn)/2.d0/deltaT/log(ten)
        dpdtt = dpdtt/dydxt + dydx

        Fi = ten**y
        pressurecheck = Ri*Fi*dydxt
        ent = -Fi/Ti*dydx
        E = (Fi + Ti*ent - fminceiling)*1000.d0*hmass
        ent = -Fi/Ti*dydx*1000.d0*hmass
        K = Pi*(d2ydx2t/dydxt/log(ten) + dydxt + 1.d0)
        alp = Pi/K*dpdtt/Ti
        dlnsdlnpi = dydxt/dydx*Pi/K*dpdtt
        Cp = ent*(d2ydx2/dydx/log(ten) + dydx - 1.d0 - dpdtt*dlnsdlnpi)

	gamma = alp*K/(0.001*Cp/Vi - alp**2*K*Ti)

	Cv = Cp/(1. + alp*gamma*Ti)
	KS = K*(1. + alp*gamma*Ti)

	P = Pi

	Ftot = E - ent*Ti 
	Fpv = 1000.*Pi*Vi
	Ftot = Ftot + Fpv + hG0

	return
	end
