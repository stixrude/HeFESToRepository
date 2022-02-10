        subroutine therm(ispec,Vi,volnl,Cp,Cv,gamma,K,Ks,alp,Ftot,ph,ent,deltas,
     &                   tcal,zeta,Gsh,uth,uto,thet,q,etas,dGdT,pzp,Vdeb,gamdeb)

        include 'P1'
        include 'const.inc'
        include 'theory.inc'

	character*80 phname(nphasep),sname(nspecp)
	integer ispec,ibv,ied,izp
	double precision Vi,volnl,Cp,Cv,gamma,alp,Ftot,ph,ent,deltas,tcal,zeta,Gsh,uth,uto,thet,q,etas,dGdT
        double precision pzp,a1,a2,a3,a4,a5,agt,alpk,alplan,anh,apar,b0,b1,b2,b3,be,beta,betlan,ca1,ca2,ca3
	double precision cplan,cvel,cvelo,cvn,cvo,cvon,cvp,d,del,deltat,detasdp,detasdv,dijkl,dkdt,dkdtp
        double precision dksdtp,duth,eta,ezp,f,fac,fbm,fel,flan,fn,fo,fpv,fth,ftho,g,gam,gammo,ge,gij,glan
	double precision go,gop,got,gp,gpp,gshp,h,helm,helmlan,htl,pa,pc,pi,q2a2,qe1,qe2,qe3,qe4,qo
	double precision qorder,qp,slan,tc,theo,ti,to,vlan,vo,vx,wd1,wd1o,wd2,wd2o,wd3,wd3o,we1,we1o,we2
        double precision we2o,we3,we3o,we4,we4o,wm,wol,wolo,wou,wouo,ws1,ws1o,ws2,ws2o,ws3,ws3o,xo,xt,xv
	double precision zeto,zu,pel,Kel,gel,qel
	double precision Etherm,Ctherm,Ztherm,Ftherm
        double precision Ko,Kop,Kopp,Kth,K,Kc,Ks
        double precision Kp,Ksp,Kpc,Kpth
	double precision alpKc,alpc,agTc,Ksc,DKDTc,deltaTc,Kspc,rho,Vpc,Vsc,Gshc,Gshpc,Vdeb,Vdeb3
	double precision dlnvpdlnv,dlnvsdlnv,dlnvdebdlnv,gamdeb

        common /state/ apar(nspecp,nparp),Ti,Pi
        common /names/ phname,sname
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
        if (thet .eq. 0.) thet = ws1
        if (thet .eq. 0.) thet = we1
        if (thet .eq. 0.) thet = (wou + wol)/2.
        anh = q2A2
        xT = thet/Ti
        xo = thet/To
c	print*, 'gamlat,qlat',gamma,q
        uth = Etherm(Ti,fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,
     &                     we1,we2,we3,we4,qe1,qe2,qe3,qe4)
        uto = Etherm(To,fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,
     &                     we1,we2,we3,we4,qe1,qe2,qe3,qe4)
        Cv  = Ctherm(Ti,fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,
     &                     we1,we2,we3,we4,qe1,qe2,qe3,qe4)
        Cvo = Ctherm(To,fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,
     &                     we1,we2,we3,we4,qe1,qe2,qe3,qe4)
C Landau contributions
C Choose: landau for inv251010 and earlier, landauqr for later
	if (iltyp .eq. 1) then
	 call landauqr(ispec,Vi,qorder,Tc,glan,cplan,alplan,betlan,slan,vlan)
	else if (iltyp .eq. 2) then
	 call landau(ispec,Vi,qorder,Tc,glan,cplan,alplan,betlan,slan,vlan)
        else
         print*, 'WARNING: Landau type not chosen.  Landau terms are not computed.'
	end if
C  Fluid
        if (htl .ne. 0.) then
	 print*, 'Warning, entering old fluid calculation in therm.f'
         Vx = 38.9
         Vx = Vo
         cvp = 0.0
         Cv = 3.*fn*Rgas*(htl + cvp*(Vi/Vx-1.))
         Cvo = 3.*fn*Rgas*(htl+ cvp*(Vi/Vx-1.))
         uth = Cv*Ti
         uto = Cvo*To
c        uth = 3.*fn*Rgas*Ti*htl
c        uto = 3.*fn*Rgas*To*htl
c        gamma = gammo + gammo*qo/Vo*(Vi-Vo)
c        q = Vi/gamma*gammo*qo/Vo
c         thet = theo*exp(-(gammo-gammo*qo)*log(Vi/Vo) - (gamma-gammo))
         gamma = gammo + gammo*qo/Vx*(Vi-Vx)
         q = Vi/gamma*gammo*qo/Vx
         thet = theo*exp(-(gammo-gammo*qo)*log(Vi/Vx) - (gamma-gammo))
        end if
        Cvn = Cv/(3.*fn*Rgas)
        Cvon = Cvo/(3.*fn*Rgas)
        call thetacal(Cvn,tcal)
        tcal = Ti*tcal
        zeta = Ztherm(Ti,fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,
     &                     we1,we2,we3,we4,qe1,qe2,qe3,qe4,gamma)/Cvn
        zeto = Ztherm(To,fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,
     &                     we1,we2,we3,we4,qe1,qe2,qe3,qe4,gamma)/Cvon
c       qp = 0.0
        a1 = 1.5*(Kop - 4.)
        a2 = 0.0

        if(ibv .eq. 1)then
C Vinet
           eta = 3.*(Kop - 1.)/2.
           xV = (Vi/Vo)**(1./3.)
           Kc = Ko*(2. + (eta - 1.)*xV -
     &          eta*xV*xV)*exp(eta*(1. - xV))/xV**2
        else if (ibv .eq. 2) then
C Lagrangian
         g = -.5*((Vi/Vo)**(+2./3.) - 1.)
         Kc = Ko*(1. - 2.*g)**(-1./2.)
     &      *(1. + g*(3.*Kop - 1.)
     &      - 4.5*Kop*g*g)
	else
C Birch-Murnaghan
           f = .5*((Vi/Vo)**(-2./3.) - 1.)
           a3 = 3.*(Kop - 4.)
           a4 = 9.*(Kopp + Kop*(Kop - 7.) + 143./9.)
           a5 = apar(ispec,43)
           if (Kopp .eq. 0.0) a4 = 0.0
           Kc = Ko*(1. + 2.*f)**2.5*(1. 
     &        + (7. + a3)*f
     &        + (9./2.*a3 + 1./2.*a4)*f*f
     &        + (11./6.*a4 + 1./6.*a5)*f*f*f
     &        + 13./24.*a5*f*f*f*f)
        endif

	ezp = 0.0
        pzp = 0.0
        pa = .001*3.*fn*Rgas*anh*(Ti**2 - To**2)/Vi
        if (abs(izp) .eq. 1) pzp = izp*0.001*(9./8.)*fn*Rgas*thet*gamma/Vi
        if (abs(izp) .eq. 1) ezp = izp*(9./8.)*fn*Rgas*thet
c        if (izp .eq. 1) pzp = 0.001*(3./2.)*
c     &    Wav(fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,
c     &              we1,we2,we3,we4,qe1,qe2,qe3,qe4)

        ph = .001*(gamma/Vi)*(uth - uto)
C  Fluid
        if (htl .ne. 0.) ph = 0.001*(gamma/Vi)*3.*htl*fn*Rgas*(Ti-To)
        a4 = 1.5*(Ko*Kopp + Kop*(Kop - 7.) + 143./9.)
        if (Kopp .eq. 0.0) a4 = 0.0
        pc = 3.*Ko*f*(1. + 2.*f)**2.5*
     &       (1. + (1.5*Kop - 6.)*f + a4*f*f)

        Fbm = 4500.*Ko*Vo*f*f*(1. + (Kop - 4.)*f + 1./2.*a4*f*f)
        Fth = Ftherm(Ti,fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,
     &                     we1,we2,we3,we4,qe1,qe2,qe3,qe4)
        Ftho = Ftherm(To,fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,
     &                     we1,we2,we3,we4,qe1,qe2,qe3,qe4)
C  Fluid
        if (htl .ne. 0.) then
c        Fth = 3.*fn*Rgas*Ti*htl*(log(thet/Ti) - 1./3.)
c        Ftho = 3.*fn*Rgas*To*htl*(log(thet/To) - 1./3.)
         Fth = Cv*Ti*(log(thet/Ti) - 1./3.)
         Ftho = Cvo*To*(log(thet/To) - 1./3.)
        end if
        Fpv = 1000.*Pi*Vi
        Kth = (gamma + 1. - q)*(ph + pzp) - .001*(gamma**2/Vi)*
     &        (Ti*Cv - To*Cvo)
        beta = be*(Vi/Vo)**(ge)
        Fel = -(beta/2.)*(Ti*Ti - To*To)
        Cvel = beta*Ti
        Cvelo = beta*To
        gamma = (gamma*Cv + ge*Cvel)/(Cv + Cvel)
        Cv = Cv + Cvel
        Cvo = Cvo + Cvelo
        Cvn = Cv/(3.*fn*Rgas)
        Cvon = Cvo/(3.*fn*Rgas)
        pel = 0.001*0.5*ge*beta*(Ti*Ti - To*To)/Vi
	qel = 0.0
	Kel = pel*(1. - ge - qel)
        Ftot = 1000.*Fo + Fbm + Fth - Ftho + Fpv + Fel
        ent = (uth - Fth)/Ti + beta*Ti

        K = Kc + Kth + Kel
        fac = gamma*gamma*Cv*Ti/(1000.*Vi*K)
        Cp = Cv*(1. + fac)
        alp = .001*gamma*Cv/(Vi*K)
        Del = .001*gamma**2*Cv*Ti/Vi
        Ks = K + Del
c       Ks = Kc + Kth + Del
        Ca1 = 7. + 2.*a1
        Ca2 = 9.*a1 + 3.*a2
        Ca3 = 11.*a2
        B1 = 7.*Ca1 + 2.*Ca2
        B2 = 9.*Ca2 + 3.*Ca3
        B3 = 11.*Ca3
        Kpc = (Ko/(3.*K))*(2.*f + 1)**2.5*(3.*Kop +
     &        B1*f + B2*f*f + B3*f*f*f)
c        zeta = 3.*gamma - (9.*fn*Rgas*gamma*xT*xT*exp(xT))/
c     &         (CV*(exp(xT) - 1.)**2)
c        zeto = 3.*gamma - (9.*fn*Rgas*gamma*xo*xo*exp(xo))/
c     &         (CVo*(exp(xo) - 1.)**2)
c        zeta = Ztherm(Ti,fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,
c     &                     we1,we2,we3,we4,qe1,qe2,qe3,qe4,gamma)/Cvn
c        zeto = Ztherm(To,fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,
c     &                     we1,we2,we3,we4,qe1,qe2,qe3,qe4,gamma)/Cvon
        Kpth = q*(qp - gamma)*ph + (gamma + 1. - q)*Kth +
     &         .001*gamma**2*(2.*q - 1.)/Vi*(Ti*CV - To*CVo) +
     &         .001*gamma**2/Vi*(Ti*CV*zeta - To*CVo*zeto)
        Kpth = Kpth/K
        Kp = Kpc + Kpth
        alpK = .001*gamma*CV/Vi
        DKDT = 1. - q - zeta
        deltaT = Kp - DKDT
        DKDTP = alpK*deltaT
        alp = alpK/K
        agT = alp*gamma*Ti
        Ks = K*(1. + agT)
        Ksp = (1. + agT)*Kp - agT*(deltaT + q)
        deltas = (Ksp + q - gamma - 1.)/(1. + agT)
        DKsDTP = deltas*alpK*(1. + agT)
c       print*, 'Ksp,Kp,Kpc,qp = ',Ksp,Kp,Kpc,qp

C Landau contributions
	Helm = Ftot - Fpv
        Volnl = Vi
        Cp = Cp + cplan
        ent = ent + slan
        Ftot = Ftot + glan 
        Vi = Vi + vlan
        Helmlan = Ftot - 1000.*Pi*Vi
	flan = Helmlan - Helm
        K = Vi/volnl/(1./K + betlan)
        alp = Volnl/Vi*(alp + alplan)
        Cv = Cp - 1000.*Ti*Vi*alp**2*K
        gamma = 1000.*Vi*alp*K/Cv
        Ks = K*(1. + alp*gamma*Ti)
c        write(71,'(i5,60f12.5)') ispec,Ti,glan/1000.,slan,vlan,qorder,Tc,Helm/1000.,Helmlan/1000.
c     &   ,flan/1000.,1000.*Pi*vlan,alplan,betlan,cplan,Ftot/1000.,gamma

        dijkl = -1
        gij = 0.
        Gp = Gop
c        Gpp = apar(ispec,38)
c       etas = Got
        if (Gop .eq. 0.) Gp = 7.*Go/(3.*Ko) - dijkl
        h = -2.*gammo*(Got - Gp) - 2.*gammo
c       b0 = Go + (2*h - gij*gij)/4.*0.001*(uth - uto)/Vo
        b0 = Go 
c       b1 = 3.*Ko*(Gp + dijkl) - 7.*Go
c       Gsh = (1. + 2.*f)**3.5*(b0 + b1*f) - pc*dijkl 
        b1 = 3.*Ko*Gp - 5.*Go
!       ivtyp = 1 : Volume dependence of G: Full expansion
!       ivtyp = 2 : Volume dependence of G: Truncated expansion
!       ivtyp = 3 : Cammarano et al. (2003)
!       ittyp = 1 : Temperature dependence of G: Full theory
!       ittyp = 2 : Temperature dependence of G: eta prop. V
!       ittyp = 3 : Temperature dependence of G: eta constant
!       ittyp = 4 : Temperature dependence of G: eta prop. gq
!       htl .ne. 0 : Fluid thermal eos
        if (ivtyp .eq. 1) then
         b2 = 6.*Ko*Gp - 24.*Ko - 14.*Go + 4.5*Ko*Kop
         Gsh = (1. + 2.*f)**2.5*(b0 + b1*f + b2*f*f)                    ! Full finite strain theory
         Gshp = 1./(3.*K)*(5.*Gsh + (1. + 2.*f)**(7./2.)*(b1 + 2.*b2*f))
	 Gshc = Gsh
	 Gshpc = Gshp
c        print*, 'Gshp = ',Gshp
        end if
        if (ivtyp .eq. 2) then
         Gsh = (1. + 2.*f)**2.5*(b0 + b1*f)                             ! As submitted to JGR (LVZ)
        end if
        if (ivtyp .eq. 3) then
         b2 = 9.*Ko**2*Gpp - 3.*(5.*Go - 3.*Ko*Gop)*(Kop - 4.) + 5.*Go*(3.*Kop - 5.)
c        print*, Gpp,b0,b1,b2
         Gsh = (1. + 2.*f)**2.5*(b0 + b1*f + b2*f*f)                    ! Cammarano et al. (2003)
        end if
        if (ibv .eq. 2) then
         b0 = Go
         b1 = 3.*Ko*Gop - Go
         b2 = -(2.*Go + 6.*Ko*Gop + 1.5*Ko*(3.*Kop + 4.))
         Gsh = (1. - 2.*g)**(-1./2.)*(b0 + b1*g + b2*g*g)               ! Lagrangian
        end if
C  Softening in stishovite
        if (sname(ispec)(1:2) .eq. "st") call stishtran(Pi,Ti,Gsh)
        if (ittyp .eq. 1) then
         Gsh = Gsh - 0.001*etas/Vi*(uth - uto + ezp)                       ! Full finite strain theory
c        Gshp = Gshp + 0.001/K*(detasdv*(uth-uto) - etas/Vi*(uth-uto) + etas*gamma/Vi*(CV*Ti - CVo*To + uth - uto))
         duth = uth - uto
         detasdp = -Vi/K*detasdv
         Gshp = Gshp - 0.001/Vi*(duth*detasdp + etas/K*duth + etas*gamma/K*(duth - CV*Ti + CVo*To))
         dGdT = -0.001*etas*CV/Vi - alp*K*Gshp
c        print*, 'Gshp,dGdT = ',Gshp,1000.*dGdT,Ksp,Gsh,Ks
        end if
        if (ittyp .eq. 2) then
         Gsh = Gsh - 0.001*Got/Vo*(uth-uto)                             ! As submitted to JGR (LVZ): prop. to V
        end if
        if (ittyp .eq. 3) then
         Gsh = Gsh - 0.001*Got/Vi*(uth-uto)                             ! constant
        end if
        if (ittyp .eq. 4) then
         Gsh = Gsh - 0.001*etas/Vi*(gamma/gammo)*(q/qo)*(uth-uto)       ! proportional to gq
        end if
c       print*, Pi,Ti,b0,b1,Gsh

	alpKc = .001*gamma*Cvo/Vi
	alpc = alpKc/Kc
        agTc = alpc*gamma*To
	Ksc = Kc*(1. + agTc)
        DKDTc = 1. - q - zeto
        deltaTc = Kpc - DKDTc
        Kspc = (1. + agTc)*Kpc - agTc*(deltaTc + q)
	rho = wm/Vi
	Vpc = sqrt((Ksc + 4./3.*Gshc)/rho)
	Vsc = sqrt(Gshc/rho)

        Vdeb3 = 2./3./Vsc**3 + 1./3./Vpc**3
        Vdeb = 1./Vdeb3**(1./3.)
        if (Vdeb3 .lt. 0.) Vdeb = 1./abs(Vdeb3)**(1./3.)
        dlnvsdlnv = 0.5*(1. - Ksc/Gshc*Gshpc)
        dlnvpdlnv = 0.5*(1. - Ksc/(Ksc + 4./3.*Gshc)*(Kspc + 4./3.*Gshpc))
        dlnvdebdlnv = 2./3.*(Vdeb/Vsc)**3*dlnvsdlnv + 1./3.*(Vdeb/Vpc)**3*dlnvpdlnv
        gamdeb = -(dlnvdebdlnv - 1./3.)

c	write(31,*) 'in therm',ispec,Pi,Ti,Vi,Ks

        return
        end
