        function gspec(ispec)

        include 'P1'

	double precision alp,apar,be,Cp,cpve,Cv,deltas,dgdt,ent,entve,etas,fn,Fo,Ftot
	double precision gam,gamma,ge,Go,Got,Gop,Gsh,htl,ph,Pi,pzp,q2a2,qo,qq,tcal,thet
	double precision Ti,To,uth,uto,Vo,vol,volnl,volve
	double precision qe1,qe2,qe3,qe4,we1,we2,we3,we4,wd1,wd2,wd3,ws1,ws2,ws3,wou,wol
	double precision wm,x1,zeta,zu,gspec,volume,volumel,volumew,bkve,volumeh
	integer ispec,i,ibv,ied,iphyflag,izp,ncall
        double precision Ko,Kop,Kopp,K,Ks,vdeb,gamdeb
	double precision E,Eel,Eig,P,Pel,Pig,Sel,Cvel,videal
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
c	 print*, 'WARNING: Temperature is less than Tsmall',Ti
c	 write(31,*) 'WARNING: Temperature is less than Tsmall',Ti
c	 Ti = Tsmall
	end if
        x1 = Vo
        if (htl .eq. 4) x1 = videal(fn,Pi,Ti)
        if (x1a(ispec) .gt. xsmall) x1 = x1a(ispec)
	if (htl .eq. 0.) vol = volume(ispec,x1)
	if (htl .eq. 1. .or. htl .eq. 4) vol = volumel(ispec,x1)
	if (htl .eq. 3.) vol = volumew(ispec,x1)
	if (htl .eq. 5.) vol = volumeh(ispec,x1)
c	print '(a11,i5,99f12.5)', 'x1a vol    ',ispec,vol
c	print '(a11,i5,99f12.5)', 'x1a after v',ispec,(x1a(jj),jj=1,47)
C  Volume failed for physical reasons: species is spinodally unstable.
C  In this case, we can try resetting the volume to a very small value, which causes
C  the Gibbs free energy to be very large and the species to be unstable, or
C  we could remove the species by setting spinod(ispec) = T and reset vol = Vo, or
C  remove the species AND set vol=vlow.
c	print*, 'In gspec volume =',vol
        if (vol .eq. -1.) then
c         print '(a33,i5,12f12.5)', 'WARNING: Spinodal instability in species',ispec,Pi,Ti,vol,Vo,x1
c         write(31,'(a33,i5,12f12.5)') 'WARNING: Spinodal instability in species',ispec,Pi,Ti,vol,Vo,x1
	 vol = x1
	 spinod(ispec) = .true.
        end if
C  Volume failed for numerical or logical reasons.  This can occur e.g. if the volume exceeds the bounds set by vibrational limits.
C  Again, we must either cause the Gibbs free energy of the species to be arbitrarily large, or remove it, or both.
        if (vol .eq. -2.) then
c         print '(a39,i5,12f12.5)', 'WARNING: Volume encountered instability',ispec,Pi,Ti,vol,Vo,x1
c         write(31,'(a39,i5,12f12.5)') 'WARNING: Volume encountered instability',ispec,Pi,Ti,vol,Vo,x1
	 vol = x1
	 spinod(ispec) = .true.
        end if
        volnl = vol
	x1a(ispec) = vol
C  Now that we have found the volume, we can compute the physical properties and free energy of the species.
        if (iphyflag .eq. 1) then
	 if (htl .eq. 0.) then
          call therm(ispec,vol,volnl,Cp,Cv,gamma,K,Ks,alp,Ftot,ph,ent,deltas,tcal,
     &               zeta,Gsh,uth,uto,thet,qq,etas,dGdT,pzp,Vdeb,gamdeb)
	 end if
	 if (htl .eq. 1. .or. htl .eq. 4) then
          call therml(ispec,vol,volnl,Cp,Cv,gamma,K,Ks,alp,Ftot,ph,ent,deltas,tcal,
     &               zeta,Gsh,uth,uto,thet,qq,etas,dGdT,pzp,Sel,Eel,Pel,Cvel,Eig,Pig,P,E)
	 end if
	 if (htl .eq. 2.) then
          call thermg(ispec,vol,volnl,Cp,Cv,gamma,K,Ks,alp,Ftot,ph,ent,deltas,tcal,
     &               zeta,Gsh,uth,uto,thet,qq,etas,dGdT,pzp,Sel,Eel,Pel,Cvel,Eig,Pig,P,E)
	 end if
	 if (htl .eq. 3.) then
          call thermw(ispec,vol,volnl,Cp,Cv,gamma,K,Ks,alp,Ftot,ph,ent,deltas,
     &                   tcal,zeta,Gsh,uth,uto,thet,qq,etas,dGdT,pzp,Sel,Eel,Pel,Cvel,Eig,Pig,P,E)
	 end if
	 if (htl .eq. 5.) then
          call thermh(ispec,vol,volnl,Cp,Cv,gamma,K,Ks,alp,Ftot,ph,ent,deltas,
     &                   tcal,zeta,Gsh,uth,uto,thet,qq,etas,dGdT,pzp,Sel,Eel,Pel,Cvel,Eig,Pig,P,E)
	 end if
        end if
	if (htl .eq. 0) call Ftotsub(ispec,volnl,Ftot)
	if (htl .eq. 1 .or. htl .eq. 4) call Ftotsubl(ispec,volnl,Ftot)
	if (htl .eq. 2)  call Ftotsubg(ispec,volnl,Ftot)
	if (htl .eq. 3)  call Ftotsubw(ispec,volnl,Ftot)
	if (htl .eq. 5)  call Ftotsubh(ispec,volnl,Ftot)

c	volve = vol
        gspec = Ftot
c	print*, 'in gspec',Ftot,Cp,Cv,vol,volnl,volve,ispec,htl
c 	write(31,*) 'in gspec',Ftot,Cp,Cv,vol,volnl,volve,ispec,htl

        return
        end
