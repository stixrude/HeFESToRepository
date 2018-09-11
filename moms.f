        subroutine moms(ispec,apar,fmom,thetaV)

        include 'P1'
        include 'const.inc'
        logical aniso

	integer ispec,ib,ied,ier,im,izp,jm,jspec,last,neval,ibv
	double precision thetav,a,abserr,b,be,detasdv,do,etas,faci2,fm,fn,fo,gam,gamma,gammo,ge,go,gop,got
	double precision gplt,htl,q,q2a2,qe,qe1,qe2,qe3,qe4,qo,qp,rho,rkmax,rwav,ss,ssi,su,to,va,vavm3,vd
	double precision vi,vo,vp,vs,wd1,wd1o,wd2,wd2o,wd3,wd3o,wdav,we1,we1o,we2,we2o,we3,we3o,we4,we4o,wm
	double precision wo,wol,wolo,wou,wouo,wplt,ws1,ws1o,ws2,ws2o,ws3,ws3o,wsav,wsavm3,zu
        double precision wb(14)
	double precision fmom(24)
        double precision apar(nspecp,nparp)
	double precision Ko,Kop,Kopp
	integer iord(10)
	double precision alist(10),blist(10),elist(10),rlist(10)
        double precision, parameter :: feps=0.1, dw = 10.0, epsabs=1.e-12, epsrel=1.e-12
	double precision, parameter :: limit=10
        external midpnt,vdosm,midsqu
        common /vdosc/ jspec,im
	jspec = ispec
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

	rho = wm/Vo
	VS = sqrt(Go/rho)
	VP = sqrt((Ko + 4./3.*Go)/rho)
	VD = (1./3./VP**3 + 2./3./VS**3)**(-1./3.)
	Va = Vo/fn
	thetaV = 251.2/Va**(1./3.)*VD
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

        Va = Vi/fn
c       fac = 251.2/Va**(1./3.)
c       faci = 4.*fac/su**(1./3.)/(2.*pirad)
        rwav = (4.*pirad*Va*su/(3.*avn))**(1./3.)
        rkmax = 2.*pirad/rwav
        faci2 = rkmax*2./pirad*hbok*cmkm
        if (aniso) then
         wsavm3 = ((1./ws1**3 + 1./ws2**3 + 1./ws3**3)/3.)**(-1./3.)
        else
         wsavm3 = ws1
        end if
        vavm3 = wsavm3/faci2
        ss = hbok*(Va/(avn*6.*pirad**2))**(-1./3.)*vavm3*cmkm
        ss = wsavm3*pirad/(2.*su**(-1./3.))
        do 4 jm=1,24
4       fmom(jm) = 0.
        if (ws1 .eq. 0.) then
	 fmom(4) = wd1
	 return
	end if
        im = -3
        print*, im,ss
        fmom(1) = ss

        do 3 ib=1,14
3       wb(ib) = 0.
        wb(1) = 0.
        wb(2) = ws1
        wb(3) = ws2
        wb(4) = ws3
        wb(5) = wol
        wb(6) = wou
        if (we1 .ne. 0.) then
         wb(7) = we1 - dw/2.
         wb(8) = we1 + dw/2.
        end if
        if (we2 .ne. 0.) then
         wb(9) = we2 - dw/2.
         wb(10) = we2 + dw/2.
        end if
        if (we3 .ne. 0.) then
         wb(11) = we3 - dw/2.
         wb(12) = we3 + dw/2.
        end if
        if (we4 .ne. 0.) then
         wb(13) = we4 - dw/2.
         wb(14) = we4 + dw/2.
        end if
	call dsort(wb,wb,14,1)
        do 1 jm=2,10
         im = jm - 4
         fm = im
         ss = 0.
         do 2 ib=1,13
          a = wb(ib)
          b = wb(ib+1)
          if (b .eq. 0.) go to 2
          if (b .eq. ws1 .or.
     &        b .eq. ws2 .or.
     &        b .eq. ws3) then
	   call dqagse(vdosm,a,b,epsabs,epsrel,limit,ssi
     &      ,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
          else
	   call dqagse(vdosm,a,b,epsabs,epsrel,limit,ssi
     &      ,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
          end if
          ss = ss + ssi
c         print*, ib,a,b,ssi,ss
2        continue
         if (im .ne. 0) then
          ss = ((fm + 3.)/3.*ss)**(1./fm)
         else
          ss = exp(1./3.)*exp(ss)
         end if
         print*, im,ss,last,neval,abserr
         fmom(jm) = ss
1       continue
	print*, thetaV,VD,VP,VS,Va

        return
        end
