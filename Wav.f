        double precision function Wav(fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,
     &                     we1,we2,we3,we4,qe1,qe2,qe3,qe4)

	double precision fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,we1,we2,we3,we4,qe1,qe2,qe3,qe4,fac,qe,qo
	double precision su,ud,ue,uo,us,wd,wdav,wo,wsav
        logical aniso
        aniso = .false.

        su = fn*zu
        wo = (wou + wol)/2.
        wdav = (wd1 + wd2 + wd3)/3.
        wsav = (ws1 + ws2 + ws3)/3.
        if (wdav .eq. 0. .and. wsav .eq. 0.) su = 1.e15
        qe = qe1 + qe2 + qe3 + qe4
        qo = 1. - 1./su - qe
        if (wo .eq. 0.) qo = 0.
        if (wdav .ne. wd1/3.) aniso = .true.
        if (wsav .ne. ws1/3.) aniso = .true.
        if (qe .eq. 0. .and. qo .eq. 0.) su = 1.

C  Debye
        
        wd = 1./su*0.75*wd1
        if (aniso) ud = ud/3. + (0.75*wd2 + 0.75*wd3)/(3.*su)

C  Sin
        fac = 0.886
        us = 1./su*fac*ws1
        if (aniso) ud = us/3. + (fac*ws2 + fac*ws3)/(3.*su)

C  Einstein
        ue = qe1*we1 + qe2*we2 + qe3*we3 + qe4*we4

C  Set contribution from OH stretch to zero
        ue = 0.

C  Optic Continuum
        uo = qo*wo

        Wav = ud + us + ue + uo

        return
        end
