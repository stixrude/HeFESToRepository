        double precision function Etherm(Ti,fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,
     &                     we1,we2,we3,we4,qe1,qe2,qe3,qe4)

        include 'P1'
        include 'const.inc'

	integer ifour
	double precision Ti,fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,we1,we2,we3,we4,qe1,qe2,qe3,qe4,do,qe,qo
	double precision su,ud,ue,uo,us,wdav,wo,wsav,Ener
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
C  If qo = 0, then assign all non-Einstein modes to the acoustic band
        if (qo .eq. 0.) then
         su = 1./(1. - qe)
        end if

C  Debye
        
        ud = 1./su*Ener(wd1/Ti,do,ione)
        if (aniso) ud = ud/3. + 
     &             (Ener(wd2/Ti,do,ione) + Ener(wd3/Ti,do,ione))/(3.*su)

C  Sin
        us = 1./su*Ener(ws1/Ti,do,ithree)
        if (aniso) us = us/3. + 
     &             (Ener(ws2/Ti,do,ithree) + Ener(ws3/Ti,do,ithree))/(3.*su)

C  Einstein
        ue = qe1*Ener(we1/Ti,do,itwo) + qe2*Ener(we2/Ti,do,itwo) +
     &       qe3*Ener(we3/Ti,do,itwo) + qe4*Ener(we4/Ti,do,itwo)

C  Optic Continuum
	ifour = 4
        uo = qo*Ener(wo/Ti,do,ifour)

        Etherm = 3.*fn*Rgas*Ti*(ud + us + ue + uo)

        return
        end
