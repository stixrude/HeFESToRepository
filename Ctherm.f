        double precision function Ctherm(Ti,fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,
     &                     we1,we2,we3,we4,qe1,qe2,qe3,qe4)

        include 'P1'
        include 'const.inc'

        double precision Ti,fn,zu,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol,we1,we2,we3,we4,qe1,qe2,qe3,qe4,do,qe,qo
        double precision su,wdav,wo,wsav,Heat,Cd,Ce,Co,Cs,hs1,hs2,hs3,he1,he2,ho
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
        Cd  = 1./su*Heat(wd1/Ti,do,1)
        if (aniso) Cd = Cd/3. + 
     &             (Heat(wd2/Ti,do,1) + Heat(wd3/Ti,do,1))/(3.*su)

C  Sin
        Cs  = 1./su*Heat(ws1/Ti,do,3)
        if (aniso) Cs = Cs/3. + 
     &             (Heat(ws2/Ti,do,3) + Heat(ws3/Ti,do,3))/(3.*su)

C  Einstein
        Ce  = qe1*Heat(we1/Ti,do,2) + qe2*Heat(we2/Ti,do,2) +
     &        qe3*Heat(we3/Ti,do,2) + qe4*Heat(we4/Ti,do,2)

C  Optic Continuum
        Co  = qo*Heat(wo/Ti,do,4)

        Ctherm = 3.*fn*Rgas*(Cd + Cs + Ce + Co)

        return
        end
