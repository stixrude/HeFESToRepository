        function gspec(ispec)

        include 'P1'

	double precision alp,apar,be,Cp,cpve,Cv,deltas,dgdt,ent,entve,etas,fn,Fo,Ftot
	double precision gam,gamma,ge,Go,Got,Gop,Gsh,htl,ph,Pi,pzp,q2a2,qo,qq,tcal,thet
	double precision Ti,To,uth,uto,Vo,vol,volnl,volve
	double precision qe1,qe2,qe3,qe4,we1,we2,we3,we4,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol
	double precision wm,x1,zeta,zu,gspec,volume,volumel,bkve
	integer ispec,i,ibv,ied,iphyflag,izp,ncall
        double precision Ko,Kop,Kopp,K,Ks,vdeb,gamdeb
	double precision E,Eel,Eig,P,Pel,Pig,Sel,Cvel
	double precision x1a(nspecp)
	logical spinod(nspecp),spinph(nphasep)
	common /volent/ volve,entve,cpve,bkve
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /prop/ vol,Cp,Cv,gamma,K,Ks,alp,Ftot,ph,ent,deltas,
     &                tcal,zeta,Gsh,uth,uto,thet,qq,etas,dGdT,pzp,Vdeb,gamdeb
        common /phycom/ iphyflag
	common /spinc/ spinod,spinph
	save x1a
        double precision, parameter :: Tsmall = 1.e-5,xsmall=0.1,tol=1.e-5,vsmall=0.5
        external pressure
        data ncall/0/

        ncall = ncall + 1
	spinod(ispec) = .false.

        if (ncall .eq. 1) then
         do 1 i=1,nspecp
          x1a(i) = 0.0
1        continue
        end if

        call parset(ispec,apar,fn,zu,wm,To,Fo,Vo,Ko,Kop,Kopp,
     &                    wd1,wd2,wd3,ws1,ws2,ws3,
     &                    we1,qe1,we2,qe2,we3,qe3,we4,qe4,wou,wol,
     &                    gam,qo,be,ge,q2A2,
     &                    htl,ibv,ied,izp,
     &                    Go,Gop,Got)

c	print '(a11,i5,99f12.5)', 'x1a volumes',ispec,(x1a(jj),jj=1,47)
        if (Ti .lt. Tsmall) then
	 print*, 'WARNING: Temperature is less than Tsmall',Ti
	 write(31,*) 'WARNING: Temperature is less than Tsmall',Ti
	 Ti = Tsmall
	end if
        x1 = Vo
        if (x1a(ispec) .gt. xsmall) x1 = x1a(ispec)
c	print*, 'gspec',htl
	if (htl .eq. 0.) then
         vol = volume(ispec,x1)
	else
	 vol = volumel(ispec,x1)
	end if
        if (vol .gt. xsmall .and. vol .lt. 10.*Vo) x1a(ispec) = vol
c	print '(a11,i5,99f12.5)', 'x1a vol    ',ispec,vol
c	print '(a11,i5,99f12.5)', 'x1a after v',ispec,(x1a(jj),jj=1,47)
C  Volume failed for physical reasons: species is spinodally unstable.
C  Reset the volume to a very small value.  This will cause
C  the Gibbs free energy to be very large and the species to be unstable.
        if (vol .eq. -1.) then
         print '(a33,i5,12f12.5)', 'WARNING: Spinodal instability in species',ispec,Pi,Ti,vol,Vo,x1
         write(31,'(a33,i5,12f12.5)') 'WARNING: Spinodal instability in species',ispec,Pi,Ti,vol,Vo,x1
         vol = Vo
	 x1a(ispec) = vol
	 spinod(ispec) = .true.
        end if
C  Volume failed for numerical or logical reasons
        if (vol .eq. -2.) then
         write(31,*) 'WARNING: Volume failed for logical reasons',Pi,Ti,ispec,vol,x1a(ispec)
         vol = Vo
	 x1a(ispec) = vol
        end if
c       print*, 'gspec',Pi,Ti,vol
        volnl = vol
        if (iphyflag .eq. 1) then
	 if (htl .eq. 0.) then
          call therm(ispec,vol,volnl,Cp,Cv,gamma,K,Ks,alp,Ftot,ph,ent,deltas,tcal,
     &               zeta,Gsh,uth,uto,thet,qq,etas,dGdT,pzp,Vdeb,gamdeb)
	 else
          call therml(ispec,vol,volnl,Cp,Cv,gamma,K,Ks,alp,Ftot,ph,ent,deltas,tcal,
     &               zeta,Gsh,uth,uto,thet,qq,etas,dGdT,pzp,Sel,Eel,Pel,Cvel,Eig,Pig,P,E)
	 end if
        end if
	if (htl .eq. 0.) then
         call Ftotsub(ispec,volnl,Ftot)
	else
	 call Ftotsubl(ispec,volnl,Ftot)
	end if

c	volve = vol
        gspec = Ftot

        return
        end
