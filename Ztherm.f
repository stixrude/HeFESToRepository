        double precision function Ztherm(Ti,fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,
     &                     we1,we2,we3,we4,qe1,qe2,qe3,qe4,gamma)

	double precision Ti,fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,we1,we2,we3,we4,qe1,qe2,qe3,qe4,gamma
	double precision ce1,ce2,ce3,ce4,cel,ceu,co,cs1,cs2,cvd1,cvd2,cvd3,cve1,cve2,cve3,do,dx,qe,qo,su
	double precision wdav,wo,wsav,x1,x2,xe1,xe2,xe3,xe4,xl,xs,xu,zeta1,zeta2,zeta3,zetad
	double precision zetad1,zetad2,zetad3,zetae,zetae1,zetae2,zetae3,zetae4,zetao,zetas,heat
        double precision, parameter :: eps = 1.e-3
        logical aniso
        aniso = .false.

        su = fn*zu
        wo = (wou + wol)/2.
        do = 0.
        if (wo .ne. 0.) do = (wou - wol)/wo
        wdav = (wd1 + wd2 + wd3)/3.
        wsav = (ws1 + ws2 + ws3)/3.
        if (wdav .eq. 0. .and. wsav .eq. 0.) su = 1.e15
        qe = qe1 + qe2 + qe3 + qe4
        qo = 1. - 1./su - qe
        if (wo .eq. 0.) qo = 0.
        if (wdav .ne. wd1/3.) aniso = .true.
        if (wsav .ne. ws1/3.) aniso = .true.
        if (qe .eq. 0. .and. qo .eq. 0.) su = 1.
C  If qo = 0, then assign all non-Einstein modes to the acoustic band.  This may not work for Ztherm.
        if (qo .eq. 0.) then
         su = 1./(1. - qe)
        end if

        zetad = 0.
        zetad1 = 0.
        zetad2 = 0.
        zetad3 = 0.
        zetas = 0.
        zetae = 0.
        zetae1 = 0.
        zetae2 = 0.
        zetae3 = 0.
        zetae4 = 0.
        zetao = 0.

C  Debye
        
        Cvd1 = Heat(wd1/Ti,do,1)
        Cve1 = Heat(wd1/Ti,do,2)
        zetad = 1./su*3.*gamma*(Cvd1 - Cve1)
        if (aniso) then
         Cvd2 = Heat(wd2/Ti,do,1)
         Cve2 = Heat(wd2/Ti,do,2)
         Cvd3 = Heat(wd3/Ti,do,1)
         Cve3 = Heat(wd3/Ti,do,2)
         zeta1 = zetad/3.
         zeta2 = 1./(3.*su)*3.*gamma*(Cvd2 - Cve2)
         zeta3 = 1./(3.*su)*3.*gamma*(Cvd3 - Cve3)
         zetad = zeta1 + zeta2 + zeta3
        end if

C  Sin
        xs = ws1/Ti
        x1 = xs - eps
        x2 = xs + eps
        Cs1 = Heat(x1,do,3)
        Cs2 = Heat(x2,do,3) 
        zetas = -xs*gamma*(Cs2 - Cs1)/(2.*eps)/su
        if (aniso) then
         zeta1 = zetas/3.
         xs = ws2/Ti
         x1 = xs - eps
         x2 = xs + eps
         Cs1 = Heat(x1,do,3)
         Cs2 = Heat(x2,do,3)
         zeta2 = -xs*gamma*(Cs2 - Cs1)/(2.*eps)/(3.*su)
         xs = ws3/Ti
         x1 = xs - eps
         x2 = xs + eps
         Cs1 = Heat(x1,do,3)
         Cs2 = Heat(x2,do,3)
         zeta3 = -xs*gamma*(Cs2 - Cs1)/(2.*eps)/(3.*su)
         zetas = zeta1 + zeta2 + zeta3
        end if

C  Einstein
        xe1 = we1/Ti
        xe2 = we2/Ti
        xe3 = we3/Ti
        xe4 = we4/Ti
        Ce1 = Heat(xe1,do,2)
        Ce2 = Heat(xe2,do,2)
        Ce3 = Heat(xe3,do,2)
        Ce4 = Heat(xe4,do,2)
        if (xe1 .ne. 0.) 
     &  zetae1 = qe1*Ce1*(xe1*exp(xe1)/(exp(xe1) - 1.) - 1. - xe1/2.)
        if (xe2 .ne. 0.) 
     &  zetae2 = qe2*Ce2*(xe2*exp(xe2)/(exp(xe2) - 1.) - 1. - xe2/2.)
        if (xe3 .ne. 0.) 
     &  zetae3 = qe3*Ce3*(xe3*exp(xe3)/(exp(xe3) - 1.) - 1. - xe3/2.)
        if (xe4 .ne. 0.) 
     &  zetae4 = qe4*Ce4*(xe4*exp(xe4)/(exp(xe4) - 1.) - 1. - xe4/2.)
        zetae = 2.*gamma*(zetae1 + zetae2 + zetae3 + zetae4)

C  Optic Continuum
        xu = wou/Ti
        xl = wol/Ti
        dx = xu - xl
        Ceu = Heat(xu,do,2)
        Cel = Heat(xl,do,2)
        Co  = Heat(wo/Ti,do,4)
        if (dx .ne. 0.) zetao = qo*gamma*(Co - xu*Ceu/dx + xl*Cel/dx)

        Ztherm = zetad + zetas + zetae + zetao

        return
        end
