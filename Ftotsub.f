        subroutine Ftotsub(ispec,Vi,Ftot)

        include 'P1'
        include 'const.inc'
        include 'theory.inc'

	integer ispec,ibv,ied,izp
	double precision vi,ftot,a3,a4,a5,alplan,anh,apar,be,beta,betlan,cp,cplan,cpve,cv,cvo,d,detasdv
	double precision entve,etas,f,fac,fbm,fel,fn,fo,fpv,fth,ftho,gam,gamma,gammo,ge,glan,go,gop,got,htl
	double precision ph,Pi,pzp,q,q2a2,qe1,qe2,qe3,qe4,qo,qorder,qp,slan,tc,theo,thet,Ti,to,uth,uto,vlan
	double precision vo,volve,vx,wd1,wd1o,wd2,wd2o,wd3,wd3o,we1,we1o,we2,we2o,we3,we3o,we4,we4o,wm,wol
	double precision wolo,wou,wouo,ws1,ws1o,ws2,ws2o,ws3,ws3o,xo,xt,zu
	double precision Etherm,Ctherm,Ftherm
        double precision Ko,Kop,Kopp,Kc,Kth,K,alp,bkve,volnl
	common /volent/ volve,entve,cpve,bkve
        common /state/ apar(nspecp,nparp),Ti,Pi
        d = 1.e-6

        call parset(ispec,apar,fn,zu,wm,To,Fo,Vo,Ko,Kop,Kopp,
     &                    wd1o,wd2o,wd3o,ws1o,ws2o,ws3o,
     &                    we1o,qe1,we2o,qe2,we3o,qe3,we4o,qe4,wouo,wolo,
     &                    gam,qo,be,ge,q2A2,
     &                    htl,ibv,ied,izp,
     &                    Go,Gop,Got)

        theo = wd1o
        gammo = gam
        call gamset(wd1o,wd2o,wd3o,ws1o,ws2o,ws3o,
     &                    we1o,we2o,we3o,we4o,wouo,wolo,
     &                    gammo,qo,Got,Vi,Vo,
     &                    wd1,wd2,wd3,ws1,ws2,ws3,
     &                    we1,we2,we3,we4,wou,wol,
     &                    gamma,q,etas,detasdv,qp)
        thet = wd1
        anh = q2A2
        xT = thet/Ti
        xo = thet/To
        qp = 0.0

C Birch-Murnaghan
        f = .5*((Vi/Vo)**(-2./3.) - 1.)
	a3 = 3.*(Kop - 4.)
        a4 = 9.*(Kopp + Kop*(Kop - 7.) + 143./9.)
	a5 = apar(ispec,43)
        if (Kopp .eq. 0.0) a4 = 0.0
        Fbm = 4500.*Ko*Vo*f*f*(1. + 1./3.*a3*f + 1./12.*a4*f*f + 1./60.*a5*f*f*f)
        Kc = Ko*(1. + 2.*f)**2.5*(1.
     &        + (7. + a3)*f
     &        + (9./2.*a3 + 1./2.*a4)*f*f
     &        + (11./6.*a4 + 1./6.*a5)*f*f*f
     &        + 13./24.*a5*f*f*f*f)

        uth = Etherm(Ti,fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,
     &                     we1,we2,we3,we4,qe1,qe2,qe3,qe4)
        uto = Etherm(To,fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,
     &                     we1,we2,we3,we4,qe1,qe2,qe3,qe4)
        Cv  = Ctherm(Ti,fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,
     &                     we1,we2,we3,we4,qe1,qe2,qe3,qe4)
        Cvo = Ctherm(To,fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,
     &                     we1,we2,we3,we4,qe1,qe2,qe3,qe4)
        Fth = Ftherm(Ti,fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,
     &                     we1,we2,we3,we4,qe1,qe2,qe3,qe4)
        Ftho = Ftherm(To,fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,
     &                     we1,we2,we3,we4,qe1,qe2,qe3,qe4)

        pzp = 0.0
        if (abs(izp) .eq. 1) pzp = izp*0.001*(9./8.)*fn*Rgas*thet*gamma/Vi
        ph = .001*(gamma/Vi)*(uth - uto)
        Kth = (gamma + 1. - q)*(ph + pzp) - .001*(gamma**2/Vi)*
     &        (Ti*Cv - To*Cvo)
        K = Kc + Kth
        fac = gamma*gamma*Cv*Ti/(1000.*Vi*K)
        Cp = Cv*(1. + fac)
        alp = .001*gamma*Cv/(Vi*K)

C  Fluid
        if (htl .ne. 0.) then
c         gamma = gammo + gammo*qo/Vo*(Vi-Vo)
c         q = Vi/gamma*gammo*qo/Vo
c        thet = theo*exp(-(gammo-gammo*qo)*log(Vi/Vo) - (gamma-gammo))
         Vx = 77.8
         Vx = Vo
         gamma = gammo + gammo*qo/Vx*(Vi-Vx)
         q = Vi/gamma*gammo*qo/Vx
         thet = theo*exp(-(gammo-gammo*qo)*log(Vi/Vx) - (gamma-gammo))
         Fth = fn*Rgas*Ti*3.*htl*(log(thet/Ti) - 1./3.)
         Ftho = fn*Rgas*To*3.*htl*(log(thet/To) - 1./3.)
        end if
        beta = be*(Vi/Vo)**(ge)
	entve = (uth - Fth)/Ti + beta*Ti

        Fpv = 1000.*Pi*Vi
        Fel = -(beta/2.)*(Ti*Ti - To*To)
        Ftot = 1000.*Fo + Fbm + Fth - Ftho + Fpv + Fel

c	print '(a8,99f16.5)', 'Ftotsub',Pi,Vi,Ti,thet,theo,gamma,gammo,q,Fth,Ftho,uth,uto,Fpv,Ftot

C Landau contributions
C Choose: landau for inv251010 and earlier, landauqr for later
	if (iltyp .eq. 1) then
	 call landauqr(ispec,Vi,qorder,Tc,glan,cplan,alplan,betlan,slan,vlan)
	else if (iltyp .eq. 2) then
	 call landau(ispec,Vi,qorder,Tc,glan,cplan,alplan,betlan,slan,vlan)
	else
	 print*, 'WARNING: Landau type not chosen.  Landau terms are not computed.'
	end if
	volnl = Vi
	Vi = Vi + vlan
        Ftot = Ftot + glan
	volve = Vi
	entve = entve + slan
	cpve = Cp + cplan
	bkve = 1./(1./K + betlan)
	bkve = Vi/volnl/(1./K + betlan)
c	print*, 'in Ftotsub',(uth-Fth)/Ti,beta*Ti,slan,entve
        
        return
        end
