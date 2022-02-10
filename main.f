        include 'P1'
        include 'dos.inc'
        include 'const.inc'
        include 'chem.inc'
        include 'lag.inc'
        include 'absent.inc'
        include 'elem.inc'

	double precision alp,apar,Cp,Cv,cpve,deltas,depthp,dgdt,dp,dt,dhdpmol,dhdtmol,dvdpmol,dsdtmol
	double precision ehugo,ent,entve,etas,freeagg,ftot,gamma,gsh,P1,ph,phugo,pi,pzp
        double precision qq,rho,vtarg,starg,superad,T1,tcal,tfreeze,thet,Ti,Tiphy,tlast
	double precision uth,uto,vhugo,vol,volve,wmagg,zeta,depth,vdeb,gamdeb,bkve
	double precision start,finish,gmainloop,gibbs,func,Pisave,Tisave,Pfrozen,Tfrozen
	integer i,ibulk,ic,icfe,ip,iphyflag,ispec,itersum,itertot,jt,icalc,noln
	integer nbulk,nt,np,nvet,nvep,nobm,noth,maxij,lineart,nbm,mfit
        character*2 atom(natomp),comp(natomp)
        character*80 phname(nphasep),sname(nspecp)
        double precision K,Ks
        double precision nnew(nspecp),nnewold(nspecp)
        double precision wreg(nphasep,nsitep,nspecp,nspecp),vreg(nphasep,nsitep,nspecp,nspecp)
        double precision binit(ncompp),dbulk(ncompp)
        double precision cpa(nspecp)
        double precision wox(natomp),stox(natomp),wcomp(natomp),stcomp(natomp)
	double precision aliqc
        logical adiabat,chcalc,adcalc,hucalc,tfix,pfix,frozen
        logical isochor,sfix,vofix,cvfix
        logical spinod(nspecp),spinph(nphasep)
        common /names/ phname,sname
        common /chor/ ispec,isochor
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /volent/ volve,entve,cpve,bkve
        common /prop/ vol,Cp,Cv,gamma,K,Ks,alp,Ftot,ph,ent,deltas,
     &                tcal,zeta,Gsh,uth,uto,thet,qq,etas,dGdT,pzp,Vdeb,gamdeb
        common /phycom/ iphyflag
        common /regcom/ wreg,vreg
        common /atomc/ stox,wox,wcomp,stcomp,atom,comp
c        common /earth/ Ps(2000),Ts(2000),dPREM(2000),PPREM(2000),
c     &                 superad,nPREM,ns,adiabat,chcalc,adcalc,hucalc
        common /mag/ icfe
        common /tfindc/ tlast,vtarg,starg,phugo,vhugo,ehugo,dvdpmol,dsdtmol,dhdpmol,dhdtmol,wmagg,chcalc,adcalc,hucalc,nvet,nvep
        common /chempot/ cpa
        common /spinc/ spinod,spinph
        common /tfixc/ tfix,pfix
        common /liqc/ aliqc(nparp,nparp),mfit,nobm,noth,sfix,vofix,cvfix,maxij,lineart,nbm,noln
        tfix = .false.
	pfix = .false.
        starg = 0.
        tlast = 0.
        phugo = 0.
        vhugo = 0.
        ehugo = 0.
        nbm = 2
        nobm = 3
        noth = 1
	noln = 1
        maxij = nobm + noth
        lineart = 9

        na = 0
        call setup(binit,dbulk,P1,dP,nP,T1,dT,nT,superad,tfreeze,nbulk,adiabat,chcalc,adcalc,hucalc,frozen,nvet,nvep)

	if (frozen) Pfrozen = superad
	if (frozen) Tfrozen = tfreeze
        iphyflag = 0
        itertot = 0
	Pi = 0.
	if (adcalc) Ti = superad
	if (chcalc) Pi = superad
	starg = 0.0
	vtarg = 1.0
	icalc = 0
	call cpu_time(start)
        do 3 ibulk = 1,nbulk+1
         do 31 ic = 1,nco
31       b(ic) = binit(ic) + dbulk(ic)*float(ibulk - 1)
         if (ibulk .gt. 1) then
          do 15 i=1,nspec
           absents(i) = absentso(i)
15        continue
          do 112 i=1,nph
           absent(i) = absento(i)
112       continue
          call sform(s,b,n1,q1,q2,nspec,nco,nc,ncs,nnull,nnulls,absents)
         end if
         print*, 'Bulk Composition',nnull,nnulls
         print*, (b(ic),ic=1,nco)
         write(31,*) 'Bulk Composition',nnull,nnulls
         write(31,*) (b(ic),ic=1,nco)
        do 1 jt=1,nT+1
         if (jt .gt. 1) then
          do 11 i=1,nspec
           absents(i) = absentso(i)
           nnew(i) = nnewold(i)
11        continue
          do 111 i=1,nph
            absent(i) = absento(i)
111       continue
          call sform(s,b,n1,q1,q2,nspec,nco,nc,ncs,nnull,nnulls,absents)
	  if (adcalc) then
	   Ti = nnew(nnull+nvet)
	   print*, 'assign Ti',jT,(nnew(i),i=1,nnull+nvet),Ti
	  end if
         end if
         do 2 ip=1,nP+1
          if (nbulk .eq. 0) print*, jt,ip
          if (.not. chcalc) Pi = P1 + dP*float(ip-1)
	  if (chcalc) vtarg = P1 + dP*float(ip-1)
          if (.not. adcalc) Ti = T1 + dT*float(jt-1)
	  if (adcalc) starg = T1 + dT*float(jt-1)
          iphyflag = 1
          call PTfind(Pi,Ti,adiabat,adcalc,hucalc,superad)
          iphyflag = 0
          Tiphy = Ti
c          Ti = max(Tiphy,tfreeze)
c          if (Ti .ne. Tiphy) print*, 'Frozen calculation: T_pet=',Ti,' T_phy=',Tiphy
          write(31,'(//,a)') '------------------------- Pressure (GPa), Depth (km), Temperature (K) -------------------------'
          depthp = depth(Pi)
          write(31,700)  Pi,depthp,Ti
c          if (adcalc .and. tlast .ne. 0.) Ti = tlast
	  icalc = icalc + 1
	  if (frozen) then
	   if (icalc .eq. 1) then
	    Pisave = Pi
	    Tisave = Ti
	    Pi = Pfrozen
	    Ti = Tfrozen
	    print*, 'frozen',Pi,Ti
            call petsub(nnew,itersum,ione)
	    Pi = Pisave
	    Ti = Tisave
	   end if
	   gibbs = func(nnew)
	   go to 10
	  end if
          call petsub(nnew,itersum,ione)
10        tlast = Ti
          itertot = itertot + itersum
          iphyflag = 1
          if (.not. adcalc) Ti = Tiphy
          call physub(nnew,rho,wmagg,freeagg,1)
          if (adcalc) write(75,'(4f16.5)') depth(Pi),Pi,Ti,ent
          iphyflag = 0
c         if (jt .eq. 1) then
          if (ip .eq. 1) then
           do 12 i=1,nspec
            absentso(i) = absents(i)
            nnewold(i) = nnew(i)
12         continue
           do 121 i=1,nph
            absento(i) = absent(i)
121        continue
          end if
2        continue
1       continue
3       continue
	call cpu_time(finish)
	gmainloop = finish - start
	print*, 'time for main loop = ',gmainloop
        print '(a34,5i12)', 'Number of iterations, per point = '
     &    ,itertot,nint(float(itertot)/float((nbulk+1)*(nP+1)*(nT+1))),nbulk+1,nP+1,nT+1

        stop 
700     format(f6.2,f8.2,f8.2,105f9.5)
        end
