        double precision function Ftherm(Ti,fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,
     &                     we1,we2,we3,we4,qe1,qe2,qe3,qe4)

        include 'P1'
        include 'const.inc'

        double precision Ti,fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,we1,we2,we3,we4,qe1,qe2,qe3,qe4,do,qe,qo
        double precision su,wdav,wo,wsav,Helm,Fd,Fe,Fo,Fs
        logical aniso
        aniso = .false.

	Ftherm = 0.
	if (Ti .le. 0.) return

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
	Fd = 0.
	Fs = 0.
	Fe = 0.
	Fo = 0.

C  Debye
        Fd = 1./su*Helm(wd1/Ti,do,1)
        if (aniso) Fd = Fd/3. + 
     &             (Helm(wd2/Ti,do,1) + Helm(wd3/Ti,do,1))/(3.*su)
        Fd = 3.*fn*Rgas*Ti*Fd + 9./8.*fn*Rgas*wd1

C  Sin
        Fs = 1./su*Helm(ws1/Ti,do,3)
        if (aniso) Fs = Fs/3. + 
     &             (Helm(ws2/Ti,do,3) + Helm(ws3/Ti,do,3))/(3.*su)
	Fs = 3.*fn*Rgas*Ti*Fs

C  Einstein
        Fe = qe1*Helm(we1/Ti,do,2) + qe2*Helm(we2/Ti,do,2) +
     &        qe3*Helm(we3/Ti,do,2) + qe4*Helm(we4/Ti,do,2)
	Fe = 3.*fn*Rgas*Ti*Fe

C  Optic Continuum
        Fo = qo*Helm(wo/Ti,do,4)
	Fo = 3.*fn*Rgas*Ti*Fo

        Ftherm = Fd + Fs + Fe + Fo
c        Ftherm = 3.*fn*Rgas*Ti*(Fd + Fs + Fe + Fo)
c        Ftherm = 3.*fn*Rgas*Ti*(Fd + Fs + Fe + Fo) + 9./8.*fn*Rgas*wd1

        return
        end
