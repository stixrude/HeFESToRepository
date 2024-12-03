        double precision function vdos(ispec,apar,w)

        include 'P1'
        include 'const.inc'
        logical aniso

	integer ispec,ibv,ied,izp
	double precision w,be,detasdv,do,etas,fn,fo,gam,gamma,gammo,ge,go,gop,got,gplt,htl,qe1,qe2,qe3,qe4
	double precision qo,qp,su,to,ud,ue,ueplt,uo,us,vi,vo,wd1,wd1o,wd2,wd2o,wd3,wd3o,wdav,we1,we1o,we2,we2o
	double precision we3,we3o,we4,we4o,wm,wo,wol,wolo,wou,wouo,wplt,ws1,ws1o,ws2,ws2o,ws3,ws3o,wsav,zu
	double precision q,q2a2,qe
	double precision vcon
	double precision Ko,Kop,Kopp
        double precision apar(nspecp,nparp)
        double precision, parameter :: feps=0.1, dw = 10.0
        aniso = .false.

        call parset(ispec,apar,fn,zu,wm,To,Fo,Vo,Ko,Kop,Kopp,
     &                    wd1o,wd2o,wd3o,ws1o,ws2o,ws3o,
     &                    we1o,qe1,we2o,qe2,we3o,qe3,we4o,qe4,wouo,wolo,
     &                    gam,qo,be,ge,q2A2,
     &                    htl,ibv,ied,izp,
     &                    Go,Gop,Got)
        Vi = Vo
        gammo = gam
        call gamset(wd1o,wd2o,wd3o,ws1o,ws2o,ws3o,
     &                    we1o,we2o,we3o,we4o,wouo,wolo,
     &                    gammo,qo,Got,Vi,Vo,
     &                    wd1,wd2,wd3,ws1,ws2,ws3,
     &                    we1,we2,we3,we4,wou,wol,
     &                    gamma,q,etas,detasdv,qp)

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
        gplt = do*wo/2.
        wplt = 1./hcok

        ud = 0.
        us = 0.
        ue = 0.
        ueplt = 0.
        uo = 0.

C  Debye
        
        if (wdav .ne. 0.) then
         if (w .le. wd1) ud = 1./su*vcon(w/wd1,do,1)/wd1
         if (aniso) then
          ud = ud/3. 
          if (w .le. wd2) ud = ud + vcon(w/wd2,do,1)/wd2/(3.*su) 
          if (w .le. wd3) ud = ud + vcon(w/wd3,do,1)/wd3/(3.*su)
         end if
        end if

C  Sin
        if (wsav .ne. 0.) then
         if (w .lt. ws1) us = 1./su*vcon(w/ws1,do,3)/ws1
          if (aniso) then
          us = us/3.
          if (w .lt. ws2) us = us + vcon(w/ws2,do,3)/ws2/(3.*su) 
          if (w .lt. ws3) us = us + vcon(w/ws3,do,3)/ws3/(3.*su)
         end if
        end if

C  Einstein
        if ((w-we1) .ge. -dw/2. .and. (w-we1) .lt. dw/2.) then
         ue = ue + qe1*vcon(w,dw,2)
         ueplt = ueplt + qe1/gplt
        end if
        if ((w-we2) .ge. -dw/2. .and. (w-we2) .lt. dw/2.) then
         ue = ue + qe2*vcon(w,dw,2)
         ueplt = ueplt + qe2/gplt
        end if
        if ((w-we3) .ge. -dw/2. .and. (w-we3) .lt. dw/2.) then
         ue = ue + qe3*vcon(w,dw,2)
         ueplt = ueplt + qe3/gplt
        end if
        if ((w-we4) .ge. -dw/2. .and. (w-we4) .lt. dw/2.) then
         ue = ue + qe4*vcon(w,dw,2)
         ueplt = ueplt + qe4/gplt
        end if

C  Optic Continuum
        if (w .ge. wol .and. w .lt. wou) uo = qo*vcon(wo,wo*do,4)

        vdos = ud + us + ue + uo
        return

        return
        end
