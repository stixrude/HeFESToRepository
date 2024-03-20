        subroutine gamset(wd1o,wd2o,wd3o,ws1o,ws2o,ws3o,
     &                    we1o,we2o,we3o,we4o,wouo,wolo,
     &                    gammo,qo,Got,Vi,Vo,
     &                    wd1,wd2,wd3,ws1,ws2,ws3,
     &                    we1,we2,we3,we4,wou,wol,
     &                    gamma,q,etas,detasdv,qp)
        include 'theory.inc'
	double precision wd1o,wd2o,wd3o,ws1o,ws2o,ws3o,we1o,we2o,we3o,we4o,wouo,wolo,gammo,qo,got,vi,vo
	double precision wd1,wd2,wd3,ws1,ws2,ws3,we1,we2,we3,we4,wou,wol,gamma,q,etas,detasdv,qp,a,b,bs
	double precision e,f,g,h,ratio

!  ityp = 1:
!  Full Finite-Strain Expansion in theta^2
!  Davies (1974) J. Phys. Chem. Solids, 35, 1513.  Stixrude and Lithgow-Bertelloni, GJI, 2005
!  ityp = 2:
!  Full Finite-Strain Expansion in theta
!  Davies (1974) J. Phys. Chem. Solids, 35, 1513.
!  ityp = 3:
!  qo = constant
!  as in LVZ paper

        E = -0.5*((Vi/Vo)**(-2./3.) - 1.)
        g = -6.*gammo
        h = g*(2. + 3.*qo + g)

        wd1 = wd1o*sqrt(1. + g*E + 0.5*h*E**2)
        wd2 = wd2o*sqrt(1. + g*E + 0.5*h*E**2)
        wd3 = wd3o*sqrt(1. + g*E + 0.5*h*E**2)
        ws1 = ws1o*sqrt(1. + g*E + 0.5*h*E**2)
        ws2 = ws2o*sqrt(1. + g*E + 0.5*h*E**2)
        ws3 = ws3o*sqrt(1. + g*E + 0.5*h*E**2)
        we1 = we1o*sqrt(1. + g*E + 0.5*h*E**2)
        we2 = we2o*sqrt(1. + g*E + 0.5*h*E**2)
        we3 = we3o*sqrt(1. + g*E + 0.5*h*E**2)
        we4 = we4o*sqrt(1. + g*E + 0.5*h*E**2)
        wou = wouo*sqrt(1. + g*E + 0.5*h*E**2)
        wol = wolo*sqrt(1. + g*E + 0.5*h*E**2)
        gamma = -(1. - 2.*E)*(g + h*E)/(6.*(1. + g*E + 0.5*h*E**2))
        q = 2.*gamma - 2./3. + h/3.*(1. - 2.*E)/(g + h*E)

        f = 0.5*((Vi/Vo)**(-2./3.) - 1.)
        a = 6.*gammo
        b = gammo*(36.*gammo - 18.*qo - 12.)
        bs = -2.*gammo - 2.*Got

        wd1 = wd1o*sqrt(1. + a*f + 0.5*b*f**2)
        wd2 = wd2o*sqrt(1. + a*f + 0.5*b*f**2)
        wd3 = wd3o*sqrt(1. + a*f + 0.5*b*f**2)
        ws1 = ws1o*sqrt(1. + a*f + 0.5*b*f**2)
        ws2 = ws2o*sqrt(1. + a*f + 0.5*b*f**2)
        ws3 = ws3o*sqrt(1. + a*f + 0.5*b*f**2)
        we1 = we1o*sqrt(1. + a*f + 0.5*b*f**2)
        we2 = we2o*sqrt(1. + a*f + 0.5*b*f**2)
        we3 = we3o*sqrt(1. + a*f + 0.5*b*f**2)
        we4 = we4o*sqrt(1. + a*f + 0.5*b*f**2)
        wou = wouo*sqrt(1. + a*f + 0.5*b*f**2)
        wol = wolo*sqrt(1. + a*f + 0.5*b*f**2)

        ratio = 1./(1. + a*f + 0.5*b*f**2)
        gamma = (1./3.)*0.5*ratio*(2.*f + 1.)*(a + b*f)
        q = (1./9.)*(18.*gamma*gamma - 6.*gamma 
     &    - 0.5*ratio*(2.*f + 1.)**2*b)/gamma
        qp = (1./9.)*(36.*gamma*gamma -6.*gamma -b*ratio**2/(6.*q)*(a+b*f)*(2*f+1.)**3 + 2.*b*ratio/(3.*q)*(2.*f+1.)**2)/gamma - q
        etas = Got*gamma/gammo
        etas = -gamma - 0.5*ratio*(2.*f + 1.)**2*bs
        detasdv = -gamma*q/Vi + 2.*bs*(1. + 2.*f)**(7./2.)*ratio/(3.*Vo) - bs*(1. + 2.*f)**(9./2.)*ratio**2/(6.*Vo)*(a + b*f)
c        write(31,*) 'Vi,gamma,q,qp =',Vi,gamma,q,qp,1./ratio,f,a,0.5*b
        if (g .eq. 0. .and. qo .eq. 0.) q = qo

        if (ityp .eq. 1) then
         return
        end if

        f = 0.5*((Vi/Vo)**(-2./3.) - 1.)
        a = 3.*gammo
        b = gammo*(9.*gammo - 9.*qo - 6.)
        bs = -gammo - Got

        wd1 = wd1o*(1. + a*f + 0.5*b*f**2)
        wd2 = wd2o*(1. + a*f + 0.5*b*f**2)
        wd3 = wd3o*(1. + a*f + 0.5*b*f**2)
        ws1 = ws1o*(1. + a*f + 0.5*b*f**2)
        ws2 = ws2o*(1. + a*f + 0.5*b*f**2)
        ws3 = ws3o*(1. + a*f + 0.5*b*f**2)
        we1 = we1o*(1. + a*f + 0.5*b*f**2)
        we2 = we2o*(1. + a*f + 0.5*b*f**2)
        we3 = we3o*(1. + a*f + 0.5*b*f**2)
        we4 = we4o*(1. + a*f + 0.5*b*f**2)
        wou = wouo*(1. + a*f + 0.5*b*f**2)
        wol = wolo*(1. + a*f + 0.5*b*f**2)

        gamma = (1./3.)*(wd1o/wd1)*(2.*f + 1.)*(a + b*f)
        q = (1./9.)*(9.*gamma*gamma - 6.*gamma 
     &    - (wd1o/wd1)*(2.*f + 1.)**2*b)/gamma
        etas = -gamma - (wd1o/wd1)*(2.*f + 1.)**2*bs
        qp = 0.

        if (ityp .eq. 2) then
         return
        end if

C  q=constant                                   ! As submitted to JGR (LVZ)
        q = qo
        gamma = gammo*(Vi/Vo)**q
        wd1 = wd1o*exp((gammo - gamma)/q)
        wd2 = wd2o*exp((gammo - gamma)/q)
        wd3 = wd3o*exp((gammo - gamma)/q)
        ws1 = ws1o*exp((gammo - gamma)/q)
        ws2 = ws2o*exp((gammo - gamma)/q)
        ws3 = ws3o*exp((gammo - gamma)/q)
        we1 = we1o*exp((gammo - gamma)/q)
        we2 = we2o*exp((gammo - gamma)/q)
        we3 = we3o*exp((gammo - gamma)/q)
        we4 = we4o*exp((gammo - gamma)/q)
        wou = wouo*exp((gammo - gamma)/q)
        wol = wolo*exp((gammo - gamma)/q)

        if (ityp .eq. 3) then
         return
        end if
        
        return
        end
