        function pressure(Vi)

        include 'P1'
        include 'const.inc'

	integer ibv,ied,ispec,izp
	double precision Vi,a3,a4,a5,anh,apar,be,beta,d,detasdv,eta,etas,f,fn,fo,g,gam,gamma,gammo,ge,go,gop
	double precision got,htl,pa,pc,pel,ph,pi,pzp,q,q2a2,qe1,qe2,qe3,qe4,qo,qp,theo,thet,ti,to,uth,uto,vo
	double precision vx,wd1,wd1o,wd2,wd2o,wd3,wd3o,we1,we1o,we2,we2o,we3,we3o,we4,we4o,wm,wol,wolo,wou
	double precision wouo,ws1,ws1o,ws2,ws2o,ws3,ws3o,xv,zu,pressure,Etherm
        logical aniso
	logical isochor
        double precision Ko,Kop,Kopp
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /chor/ ispec,isochor
        d = 1.e-6
        aniso = .false.

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

        uth = Etherm(Ti,fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,
     &                     we1,we2,we3,we4,qe1,qe2,qe3,qe4)
        uto = Etherm(To,fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,
     &                     we1,we2,we3,we4,qe1,qe2,qe3,qe4)

        ph = .001*(gamma/Vi)*(uth - uto)
C  Fluid
        if (htl .ne. 0.) then
	 print*, 'Warning entering old fluid calculation in pressure.f'
	 write(31,*) 'Warning entering old fluid calculation in pressure.f'
c         gamma = gammo + gammo*qo/Vo*(Vi-Vo)
c         q = Vi/gamma*gammo*qo/Vo
         Vx = 77.8
         Vx = Vo
         gamma = gammo + gammo*qo/Vx*(Vi-Vx)
         q = Vi/gamma*gammo*qo/Vx
         ph = 0.001*(gamma/Vi)*3.*htl*fn*Rgas*(Ti-To)
        end if
        pa = .001*3.*fn*Rgas*anh*(Ti**2 - To**2)/Vi
        pzp = 0.0
        if (abs(izp) .eq. 1) pzp = izp*0.001*(9./8.)*fn*Rgas*thet*gamma/Vi
c       if (izp .eq. 1) pzp = 0.001*(3./2.)*
c     &    Wav(fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,
c     &              we1,we2,we3,we4,qe1,qe2,qe3,qe4)
        beta = be*(Vi/Vo)**(ge)
        pel = .001*.5*ge*beta*(Ti*Ti - To*To)/Vi

        if (ibv .eq. 2) then
C Lagrangian correct to 3rd order only
         g = -.5*((Vi/Vo)**(+2./3.) - 1.)
         pc = 3.*Ko*g*(1. - 2.*g)**(-1./2.)*(1. + 1.5*Kop*g)
        else if (ibv .eq. 1) then
C Vinet
         eta = 3.*(Kop - 1.)/2.
         xV = (Vi/Vo)**(1./3.)
         pc = 3.*Ko*(1. - xV)*exp(eta*(1. - xV))/xV**2
        else   
C Birch-Murnaghan
         f = .5*((Vi/Vo)**(-2./3.) - 1.)
         a3 = 3.*(Kop - 4.)
         a4 = 9.*(Kopp + Kop*(Kop - 7.) + 143./9.)
         a5 = apar(ispec,43)
         if (Kopp .eq. 0.0) a4 = 0.0
         pc = 3.*Ko*f*(1. + 2.*f)**2.5*(1. + 1./2.*a3*f + 1./6.*a4*f*f + 1./24.*a5*f*f*f)
        end if   

        pressure = pc + ph + pa + pel + pzp - Pi

c        write(31,*) 'pressure = ',ispec,Pi,Ti,Vi,pressure,pc,ph,thet,gamma,uth,uto
cc     &    ,3.*Ko*f*(1. + 2.*f)**2.5*(1. + 1./2.*a3*f),f

        return
        end
