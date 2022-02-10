        double precision function func(nnew)

        include 'P1'
        include 'chem.inc'
        include 'absent.inc'
        include 'lag.inc'

	double precision apar,chempot,volve,entve,cpve,dhdpmol,dhdtmol,dvdpmol,dsdtmol,ehugo,vhugo,phugo
	double precision ersum,fnagg,hcalc,Pi,Ti,Pold,Told,rsum,smag,smixi,vtarg,starg,Tlast
	double precision vsum,gspec,depth,dgdtmol,wm,wmagg,fcalc,bkve,volagg,volsum
	double precision start,finish,gloop
	integer iph,ispec,iter,npres,i,nvet,nvep
        double precision nnew(nspecp)
        double precision cpa(nspecp)
	double precision fnpha(nspecp)
	double precision sspeco(nspecp),vspeco(nspecp)
	logical valid
	logical chcalc,adcalc,hucalc,tfix,pfix
        common /chempot/ cpa
        common /state/ apar(nspecp,nparp),Ti,Pi
	common /volent/ volve,entve,cpve,bkve
        common /tfindc/ tlast,vtarg,starg,phugo,vhugo,ehugo,dvdpmol,dsdtmol,dhdpmol,dhdtmol,wmagg,chcalc,adcalc,hucalc,nvet,nvep
	common /tfixc/ tfix,pfix
	save Pold,Told
	data iter/0/
	data gloop/0.0/
	data Pold/-1.e15/,Told/-1.e15/

	iter = iter + 1
	hcalc = 0.
	fcalc = 0.
c	if (tfix) write(31,*) 'tfix is true in func'
c	if (pfix) write(31,*) 'pfix is true in func'
	if (adcalc) then
	 hcalc = 1.
	 Ti = nnew(nnull+nvet)
c	 if (.not. tfix) Ti = nnew(nnull+nvet)
c	 if (tfix) hcalc = 0.
	end if
	if (chcalc) then
	 fcalc = 1.
	 Pi = nnew(nnull+nvep)
c	 if (.not. pfix) Pi = nnew(nnull+nvep)
c	 if (pfix) fcalc = 0.
	end if
	call cpu_time(start)
	if (Pi .ne. Pold .or. Ti .ne. Told) then
	 do 5 ispec=1,nspec
c	  if (absents(ispec)) go to 5
	  gspeca(ispec) = gspec(ispec)
	  vspeco(ispec) = volve
	  sspeco(ispec) = entve
	  cspeca(ispec) = cpve
	  bspeca(ispec) = bkve
c	  write(31,*) 'func 0 props',ispec,Pi,Ti,gspeca(ispec),vspeco(ispec),sspeco(ispec),cspeca(ispec),Pold,Told
5	 continue
	 Pold = Pi
	 Told = Ti
	end if
	call cpu_time(finish)
	gloop = gloop + finish - start

        call nform(nnew,n,n1,q2,nspec,nnull)
c	write(31,*) 'func',(n(i),i=1,nspec),(gspeca(i),i=1,nspec),Pi,Pold,Ti,Told
	
	if (.not. valid(vsum,nnew)) then
c	 write(31,*) 'WARNING: invalid solution in func',vsum
c	 write(31,*) (n(j),j=1,nspec)
c	 write(31,*) (nnew(j),j=1,nnull)
	end if
        ersum = 0.0
	dgdtmol = 0.0
	dhdtmol = 0.0
	dsdtmol = 0.0
	dvdpmol = 0.0
	wmagg = 0.
	volagg = 0.
        do 1 ispec=1,nspec
         cpa(ispec) = 0.0
         if (absents(ispec)) go to 1
	 wm = apar(ispec,3)
	 wmagg = wmagg + n(ispec)*wm
         call cp(ispec,n,chempot,rsum,volsum,smixi,smag)
         cpa(ispec) = gspeca(ispec)/1000. + chempot/1000.
	 sspeca(ispec) = sspeco(ispec) + smixi + smag
c	print*, 'in func',ispec,sspeco(ispec),smixi,smag,sspeca(ispec)
	 vspeca(ispec) = vspeco(ispec) + volsum
	 volagg = volagg + n(ispec)*vspeca(ispec)
	 dhdtmol = dhdtmol + n(ispec)*cspeca(ispec)/1000.
	 dgdtmol = dgdtmol - n(ispec)*sspeca(ispec)
	 dvdpmol = dvdpmol - n(ispec)*vspeca(ispec)/bspeca(ispec)
	 cpa(ispec) = cpa(ispec) - cpcomp(ispec)
c	 cpa(ispec) = cpa(ispec) - cpcomp(ispec) + hcalc*Ti*sspeca(ispec)/1000. - fcalc*Pi*vspeca(ispec)
         ersum = ersum + n(ispec)*rsum
c	 write(31,*) 'func final props',Pi,Ti,ispec,n(ispec),gspeca(ispec),vspeca(ispec)
c     &    ,sspeca(ispec),sspeco(ispec),smixi,smag,vspeco(ispec),volsum
c	 write(31,'(a13,i5,71f29.15)') 'in func',ispec,n(ispec),chempot/1000.,rsum,gspeca(ispec)/1000.,cpcomp(ispec),cpa(ispec)
c     &    ,vspeca(ispec),bspeca(ispec),dvdpmol
1       continue
	dsdtmol = 1000.*dhdtmol/Ti
	dhdtmol = 0.
	dhdpmol = 0.
	if (adcalc) dhdtmol = (dgdtmol + hcalc*starg*wmagg)/1000.
	if (chcalc) dhdpmol = volagg - fcalc*wmagg/vtarg
c	if (chcalc) dhdpmol = -fcalc*Pi*dvdpmol

        func = 0.0
        do 4 ispec=1,nspec
         if (absents(ispec)) go to 4
c         func = func + n(ispec)*(cpa(ispec) + hcalc*Ti*sspeca(ispec)/1000.)
         func = func + n(ispec)*cpa(ispec)
c	 write(31,*) 'Chemical potentials',ispec,func,n(ispec),cpa(ispec)
4       continue
	func = func + hcalc*Ti*starg/1000.*wmagg - fcalc*Pi*wmagg/vtarg
c	write(31,*) 'func = ',Pi,Ti,func,starg,vtarg,wmagg,hcalc,fcalc,dgdtmol,dhdtmol,hcalc*starg*wmagg,-fcalc*Pi*wmagg/vtarg
c     &   ,nnew(nnull+nvet),nnew(nnull+nvep)
c	write(31,*) 'func = ',Pi,Ti,func
c  ,vtarg,wmagg,wmagg/vtarg,volagg,-fcalc*Pi*wmagg/vtarg,nnew(nnull+nvep)

	fnagg = 0.
        do 11 iph=1,nph
	 fnpha(iph) = 0.
         do 12 ispec=1,nspec
          if (f(iph,ispec) .eq. 0) go to 12
	  fnpha(iph) = fnpha(iph) + apar(ispec,1)*n(ispec)
12	 continue
	 fnagg = fnagg + fnpha(iph)
11	continue
	npres = 0
	do 13 iph=1,nph
	 fnpha(iph) = fnpha(iph)/fnagg
	 if (fnpha(iph) .gt. 1.e-6) npres = npres + 1
13	continue

c	write(31,*) 'calls to func',iter,gloop

c	if (npres .gt. 1) write(95,700) Pi,depth(Pi),Ti,(fnpha(iph),iph=1,nph),func

        return
700     format(f6.2,f8.2,f8.2,105f22.12)
        end
