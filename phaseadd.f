        subroutine phaseadd(add,allow)

        include 'P1'
        include 'chem.inc'
        include 'lag.inc'
        include 'absent.inc'
	include 'numpar.inc'
	include 'const.inc'

	integer i,ibv,ic,ied,iguess,iph,ispec,iter,itermin,izp,jspec,m,ncoph,ncph,ncsph,ndim,nnullsv,nvet,nvep,k,ispec1,ispec2
	integer itersum,inull,igmin,imin,ires,iguessmax
	double precision affmin,apar,be,dhdpmol,dhdtmol,dvdpmol,dsdtmol,ehugo,fn,fo,fret,fretmin,fretmax,gam,ge,go,gop,got,htl,phugo
	double precision Pi,q2a2,qe1,qe2,qe3,qe4,qo,qual,qualmin,vtarg,starg,Ti,tlast,to,vhugo,vo,vsum,wd1,wd2,wd3
	double precision we1,we2,we3,we4,wm,wol,wou,ws1,ws2,ws3,zu,wmagg,vec1(nspecp),vec2(nspecp),test,dnrm2,ddot
	double precision lb(nspecp),ub(nspecp),fretguess,fretgmin,func
        logical add,allow(nphasep)
        logical spinod(nspecp),spinph(nphasep)
	logical valid,exsolve
	logical absentssv(nspecp),pabounds
	logical chcalc,adcalc,hucalc,adcalcsave,chcalcsave,tfix,pfix
	double precision nnsv(nspecp),n1sv(nspecp),q2sv(nspecp,nspecp)
	double precision sph(ncompp,nspecp),bph(ncompp),q1ph(nspecp,nspecp)
	integer lafmin
        double precision Ko,Kop,Kopp
        double precision nnew(nspecp)
	double precision affin(nphasep)
        double precision nnewmin(nspecp),nafmin(nspecp)
        character*80 phname(nphasep),sname(nspecp)
        double precision, parameter :: ftol=1.e-12, ntol=1.e-2, alarge=10.
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /names/ phname,sname
        common /spinc/ spinod,spinph
        common /tfindc/ tlast,vtarg,starg,phugo,vhugo,ehugo,dvdpmol,dsdtmol,dhdpmol,dhdtmol,wmagg,chcalc,adcalc,hucalc,nvet,nvep
        common /tfixc/ tfix,pfix
	common /nlbounds/ lb,ub,pabounds
	pabounds = .true.
	add = .false.
	adcalcsave = adcalc
	chcalcsave = chcalc

C  nphase(lph) = number of moles of phase lph
C  mphase(lph) = number of species in phase lph
C  iphase(lph) = species number of the initial species in phase lph
C  absent(lph) = T if phase lph is constrained to be absent

        do 7 iph=1,nph
         if (absent(iph)) go to 10
7       continue
	pabounds = .false.
        return
10      continue

        nphpres = 0
        do 91 iph=1,nph
	 affin(iph) = +1.e15
         if (absent(iph)) go to 91
         nphpres = nphpres + 1
91      continue

C--->  Save chem parameters and matrices needed globally (e.g. by func)
	nnullsv = nnull
	do 41 ispec=1,nspec
	 absentssv(ispec) = absents(ispec)
	 n1sv(ispec) = n1(ispec)
	 nnsv(ispec) = n(ispec)
	 do 41 jspec=1,nspec
	  q2sv(ispec,jspec) = q2(ispec,jspec)
41	continue
C<---

c	if (adcalc) then
c	 adcalc = .false.
c	 val = func(nnew)
c	end if
        call lagcomp(s,cpcomp,lagc,nspec,nco,absents)

        write(31,'(a)') 'Check for phase addition or subtraction'
        write(31,'(a7,2x,9a12)') 'phase','affinity','best_guess','total_guess','iterations','total_iter','quality','composition'
        affmin = +1.e15

	lafmin = nphasep
C  Loop over phases, excluding those that are 1) already present 2) not allowed 3) spinodally unstable 4) redundant
        do 3 lph=1,nph
	 exsolve = .false.
         if (.not. absent(lph)) go to 3
	 if (.not. allow(lph)) go to 3
	 if (spinph(lph)) go to 3
C  Avoid redundant affinity calculation
	 if (lph .gt. 1) then
	  if (phname(lph) .eq. phname(lph-1)) then
	   if (absent(lph-1) .and. allow(lph-1)) go to 3
	   if (.not. absent(lph-1)) exsolve = .true.
	  end if
	 end if
         m = mphase(lph)
         fretmin = +1.e15
         fretmax = -1.e15
         ispec = iphase(lph)
         call parset(ispec,apar,fn,zu,wm,To,Fo,Vo,Ko,Kop,Kopp,
     &                   wd1,wd2,wd3,ws1,ws2,ws3,
     &                   we1,qe1,we2,qe2,we3,qe3,we4,qe4,wou,wol,
     &                   gam,qo,be,ge,q2A2,
     &                   htl,ibv,ied,izp,
     &                   Go,Gop,Got)

C  Set up the linear problem sum_i n_i = 1 augmented by n_i = 0 for spinodally unstable species
	 ncoph = 1
	 do 33 ic=1,ncompp
	  bph(ic) = 0.
	  if (ic .eq. 1) bph(ic) = 1.
	  do 33 ispec=1,nspec
	   sph(ic,ispec) = 0.
	   if (ic .eq. 1 .and. f(lph,ispec) .ne. 0.) then
c	    if (.not. spinod(ispec)) sph(ic,ispec) = 1.
	    sph(ic,ispec) = 1.
	   end if
33	 continue
	 do 31 jspec=1,nspec
	  absents(jspec) = .true.
	  nnew(jspec) = 0.
31	 continue
	 do 32 jspec=iphase(lph),iphase(lph)+mphase(lph)-1
	  if (.not. spinod(jspec)) absents(jspec) = .false.
32	 continue
	 call sform(sph,bph,n1,q1ph,q2,nspec,ncoph,ncph,ncsph,nnull,nnulls,absents)

C  Find the phase composition that minimizes SLB11 Eq. B10.  Consider the series of initial guesses:
C  for m>1 guesses i=1,m where m is the number of species in the phase
C  guess i has all n_j=0 except n_i
C  guess m+1 has n_j=1/m all j
C  Evaluate the affinity for each of these initial guesses.
C  Use these initial guesses to find bounds for the minimization problem.
	 itersum = 0
	 igmin = 0
	 do 17 inull=1,nspecp
	  lb(inull) = +1.e15
	  ub(inull) = -1.e15
17	 continue
	 fretgmin = +1.e15
	 iguessmax = m+1
	 if (m .eq. 1) iguessmax = 1
         do 19 iguess=1,iguessmax
          do 18 i=1,m
           n(iphase(lph)+i-1) = ssmall
           if (i .eq. iguess) n(iphase(lph)+i-1) = 1.0 - float(m-1)*ssmall
18         continue
          if (iguess .eq. m+1) then
           do 161 i=1,m
161          n(iphase(lph)+i-1) = 1./float(m)
          end if
c	  write(31,*) iguess,(n(ispec),ispec=iphase(lph),iphase(lph)+mphase(lph)-1)
          if (m .eq. 1) n(ispec) = 1.0
          iter = 0
	  call newfrm(q2,n,n1,nnew,nspec,nnull,absents)
	  fretguess = func(nnew)
c	  write(31,*) iguess,(nnew(inull),inull=1,nnull),(n(ispec),ispec=iphase(lph),iphase(lph)+mphase(lph)-1)
c     &     ,(n1(ispec),ispec=iphase(lph),iphase(lph)+mphase(lph)-1),fretguess/fn
	  fretgmin = min(fretgmin,fretguess)
	  if (fretgmin .eq. fretguess) igmin = iguess
	  do 62 inull=1,nnull
	   lb(inull) = min(lb(inull),nnew(inull))
	   ub(inull) = max(ub(inull),nnew(inull))
62	  continue
19	continue
c	write(31,*) 'Bounds',(lb(inull),inull=1,nnull),(ub(inull),inull=1,nnull),fretgmin/fn
C  Skip phase because the best guess at its affinity it very unfavorable
	if (fretgmin/fn .gt. alarge) then
	 itermin = 0
	 itersum = 0
	 qualmin = 0
         write(31,'(i3,1x,a5,f12.4,4i12,99f12.4)') lph,phname(lph),fretgmin/fn,igmin,iguessmax,itermin,itersum,qualmin
c     &    ,(n(i),i=iphase(lph),iphase(lph)+mphase(lph)-1)
	 go to 3
	end if

C  Now do the minimization starting with each of the above initial guesses.  Save the one that produces the minimum affinity.
	 imin = 0
         do 9 iguess=1,iguessmax
          do 8 i=1,m
           n(iphase(lph)+i-1) = ssmall
           if (i .eq. iguess) n(iphase(lph)+i-1) = 1.0 - float(m-1)*ssmall
8         continue
          if (iguess .eq. m+1) then
           do 61 i=1,m-1
61          n(iphase(lph)+i-1) = 1./float(m)
          end if
          if (m .eq. 1) n(ispec) = 1.0
          iter = 0
	  call newfrm(q2,n,n1,nnew,nspec,nnull,absents)
	  ndim = nnull
c	  if (adcalc) then
c	   ndim = nnull + 1
c	   nnew(ndim) = Ti
c	  end if
	  if (adcalc) tfix = .true.
	  if (chcalc) pfix = .true.
          call nlmin_L(nnew,ndim,fret,iter,ires)
C  do not consider further failed minimizations.  Allow ires=-4 (Halted because roundoff errors limited progress).
	  if (ires .lt. 0 .and. ires .ne. -4) then
	   write(31,*) 'Minimization starting from guess',iguess,' out of',iguessmax,' failed',fret,ires
	   go to 9
	  end if
c	  if (adcalc) Ti = nnew(ndim)
          call nform(nnew,n,n1,q2,nspec,nnull)
	  call qcalc(nnew,qual)
	  if (.not. valid(vsum,nnew)) write(31,*) 'WARNING: invalid solution in phaseadd',vsum,lph,iguess

	  itersum = itersum + iter
          fretmin = min(fretmin,fret)
          fretmax = max(fretmax,fret)
c	  if (iguess .eq. igmin) fretgmin = fret
          if (fret .eq. fretmin) then
           do 11 i=1,m-1
11         nnewmin(i) = nnew(i)
	   itermin = iter
	   qualmin = qual
	   imin = iguess
          end if
c	  write(31,*) lph,iguess,iter,fret,fretmin,qual,qualmin,ndim,(n(i),i=iphase(lph),iphase(lph)+mphase(lph)-2)
c	  if (abs(fret-fretmin) .gt. 1.e-6) 
c     &     write(31,*) 'WARNING: phaseadd solution depends on initial conditions',fret/fn,fretmin/fn,(fret-fretmin)/fn,lph,iguess
9	 continue
15       fret = fretmin
         do 12 i=1,m-1
12       nnew(i) = nnewmin(i)
         call nform(nnew,n,n1,q2,nspec,nnull)
	 qual = qualmin
         affin(lph) = fret/fn
	 if (fretmax-fretmin .gt. 1.e-6) write(31,*) 'WARNING: phaseadd solution depends on initial conditions '
     &     ,phname(lph)(1:6),imin,fretmax/fn,fretmin/fn,(fretmax-fretmin)/fn
	 if (fretmin - fretgmin .gt. 0.) write(31,*) 'Minimizations did not improve on initial guesses',lph,fretmin/fn,fretgmin/fn
     &    ,(fretgmin - fretmin)/fn,igmin,imin
	 if (fretgmin*fretmin .lt. 0. .and. fretgmin/fn .gt. alarge) write(31,*) 
     &    'Best initial guess much worse than best minimization',lph,fretmin/fn,fretgmin/fn,(fretgmin - fretmin)/fn,igmin,imin
         write(31,'(i3,1x,a5,f12.4,4i12,99f12.4)') lph,phname(lph),affin(lph),imin,iguessmax,itermin,itersum,qualmin,
     &    (n(i),i=iphase(lph),iphase(lph)+mphase(lph)-1)
C  Avoid adding phases of redundant composition
	 if (exsolve) then
	  k = 0
          do 24 ispec1=iphase(lph-1),iphase(lph-1)+mphase(lph-1)-1
           ispec2 = iphase(lph) + k
           vec1(k+1) = nnsv(ispec1)
           vec2(k+1) = n(ispec2)
           k = k + 1
24        continue
          test = ddot(k,vec1,ione,vec2,ione)/(dnrm2(k,vec1,ione)*dnrm2(k,vec2,ione)) - 1.
          write(31,'(a29,1x,a4,99e16.5)') 'Attempting to exsolve a phase',phname(lph-1),(vec1(k),k=1,mphase(lph-1))
     &     ,(vec2(k),k=1,mphase(lph)),test
          if (abs(test) .lt. ntol) then
	   go to 3
	  end if
	 end if
         affmin = min(affmin,affin(lph))
         if (affin(lph) .eq. affmin) then
          lafmin = lph
          do 21 i=1,nspec
21        nafmin(i) = n(i)
         end if
3       continue

	if (lafmin .ne. nphasep) then
         write(31,'(i3,1x,a5,90f12.4)') lafmin,phname(lafmin),affmin,
     &    (nafmin(i),i=iphase(lafmin),iphase(lafmin)+mphase(lafmin)-1)
 	 write(31,*) 'phaseadd Ti',Ti
	end if

C  If a phase with negative affinity is not found, set lafmin to an arbitrary value
	if (lafmin .eq. nphasep) lafmin = 1

        do 51 i=1,nspec
51      dn(i) = 0.
	lph = lafmin
        if (affin(lph) .lt. 0.0 .and. allow(lph)) then
C  -> Test only!  
c        if (affin(lph)-100. .lt. 0.0 .and. allow(lph)) then
C  <-
         if (absent(lph)) then
          add = .true.
          absent(lph) = .false.
          allow(lph) = .false.
          print*, 'Adding Phase',lph,phname(lph)
          write(31,*) 'Adding Phase',lph,phname(lph)
          nphpres = nphpres + 1
          do 23 i=1,nspec
           if (f(lph,i) .ne. 0) then
            if (.not. spinod(i)) then
	     absentssv(i) = .false.
	     dn(i) = mphase(lph)*tradd*nafmin(i)
c	     dn(i) = tradd*nafmin(i)
	    end if
           end if
23        continue
         end if
        end if

	do 53 ispec=1,nspec
53	cpcomp(ispec) = 0.

	adcalc = adcalcsave
	chcalc = chcalcsave
	  if (adcalc) tfix = .false.
	  if (chcalc) pfix = .false.

C--->  Restore chem parameters and matrices
        nnull = nnullsv
        do 47 ispec=1,nspec
         absents(ispec) = absentssv(ispec)
         n1(ispec) = n1sv(ispec)
         n(ispec) = nnsv(ispec)
         do 47 jspec=1,nspec
          q2(ispec,jspec) = q2sv(ispec,jspec)
47      continue
C<---
	write(31,*) 'Ti at end of phaseadd',Ti
	write(31,*) 'absents at end of phaseadd',(absents(ispec),ispec=1,nspec)
	pabounds = .false.
        return
        end
