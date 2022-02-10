        subroutine physub(nnew,rho,wmagg,freeagg,iprint)

        include 'P1'
        include 'chem.inc'
        include 'const.inc'
        include 'lag.inc'
        include 'absent.inc'
	
	integer iprint,i,ibv,ic,icfe,ied,ii,iimax,iph,iphyflag,ispec,izp,jspec,ncall,nvet,nvep
	integer ipiv(nspecp,nspecp),info,iter,isign,nfast,nnullfast,j,iiron,img,ife
	double precision swork(nspecp,nspecp)
	double precision alpmet,alptot,cpmet,cptot,cvagg,dhdpmol,dhdtmol,dvdpmol,dsdtmol,phugo,xfe
	double precision ehugo,vtarg,starg,tlast,tmelt,vhugo,asqrt,tlindeman
        double precision rho,wmagg,freeagg,alp,alpagg,alpph,apar,cpph,cpphtot,cvphtot,gamphtot,alpphtot
	double precision buktphtot,bukphreuss
        double precision baggv,baggr,baggh,be,bhashm,bhashp,bmax,bmin
        double precision btaggh,btaggr,btaggv,bukph,buktph,cap,chempot,cpagg,cv
        double precision delagg,delph,deltas,dgdt,dgdtagg,dgdtaggr,dgdtaggv
        double precision dgdtph,dlnvbdt,dlnvbdtph,dlnvsdt,ent,entagg
        double precision enth,enthagg,entoph,entph,etas,etot,fdumm,fgammam,eintagg
        double precision dlnvpdt,fgammap,flambdam,flambdap,fn,fnagg,fnph
        double precision fnphamax,fo,ftot,gaggh,gaggv,gam,gamma,ge,ghashm,ghashp
        double precision giboph,gibph,gmax,gmin,go,gop,got,gaggr,gruagg,gsh
        double precision gshph,gsolph,gtot,henoph,henph,hsolph,htl,htot
        double precision ph,q2a2,qe1,qe2,qe3,qe4,qo,qp,qq,qs,rhoo,rhoph
        double precision pi,pzp,rsum,sconf,smag,smix,smixi,ssum,stsum,tcal,thet
        double precision ti,to,uth,uto,vabsagg,vabsoagg,vbh,vbr,vbv,vdeb,vdeb3
        double precision vol,volagg,vo,volph,volpho,vph,vphasha,vphashm,vphashp
        double precision vpmax,vpmin,vpr,vpred,vpv,vshasha,vshashm,vshashp
        double precision vsv,vsmax,wmaggc,volsum,volxs
        double precision vsmin,vsr,vsh,vsred,wd1,wd2,wd3
        double precision we1,we2,we3,we4,wm,wmph,wol,wou,ws1,ws2,ws3,wsum
        double precision ximax,ximin,zeta,zu,gspec,depth,epsilon1,epsilon2,epsilon,wmav
	double precision ksp,gshp,dlnvsdlnv,dlnvpdlnv,dlnvdebdlnv,gamdeb,gamp,tcalagg,vdebagg,vdebold
	double precision start,finish,gphysub,Gqcalc
	double precision Gwu,Giron(nspecp),Gironmin,fug,fugiw,fnaggtransform
        character*2 atom(natomp),comp(natomp)
        character*80 phname(nphasep),sname(nspecp)
	double precision nnew(nspecp)
	double precision st(ncompp),wco(nphasep,ncompp)
        double precision cpa(nspecp)
        double precision K,Ks,Ko,Kop,Kopp
        double precision vbpha(nphasep),vspha(nphasep),vppha(nphasep),rhopha(nphasep)
        double precision volpha(nphasep),fnpha(nphasep),wmpha(nphasep)
        double precision entpha(nphasep),henpha(nphasep),gibpha(nphasep)
	double precision entopha(nphasep),henopha(nphasep),gibopha(nphasep)
        double precision vabs(nphasep),vabso(nphasep)
        double precision gpha(nphasep),bpha(nphasep),btpha(nphasep)
        double precision ntemp(nspecp)
        double precision wox(natomp),stox(natomp),wcomp(natomp),stcomp(natomp)
	double precision hespro(nspecp,nspecp),dmdt(nspecp),dndt(nspecp),dmdtp(nspecp),dndtp(nspecp)
	double precision dfdp(nphasep),dnatomdp(nspecp),dfdptotal,dfdpabstotal
	double precision dmdp(nspecp),dmdpp(nspecp),dndpp(nspecp),dndp(nspecp),bmet,btot,bstot,gamtot,cvtot
	double precision dndpfast(nspecp),dndtfast(nspecp),hess(nspecp,nspecp),dndptemp(nspecp,ione),q2save(nspecp,nspecp)
	double precision vwork(nspecp,nspecp),checkt,checkp,ddot,wwork(nspecp),dndttest(nspecp)
	double precision btiso,bsiso,cpiso,cviso,alpiso,gamiso,phasebuoyancyparameter,deltaent,deltavol,ClapeyronSlope
	double precision Tro2,ao2,bo2,co2,do2,ho2gas,so2gas,go2gas,STro2
	double precision, parameter :: rhomantle = 4422.76, gmantle = 10.0, hmantle = 2891000.
	double precision, parameter :: pmantle = rhomantle*gmantle*hmantle/1.e9
	logical chcalc,adcalc,hucalc
	integer, parameter :: iaddmax=5
        double precision, parameter :: Tsmall = 1.e-5, asmall=1.e-15, nsmall=1.e-15, dfdpsmall=1.e-13
        common /names/ phname,sname
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /prop/ vol,Cap,Cv,gamma,K,Ks,alp,Ftot,ph,ent,deltas,
     &                tcal,zeta,Gsh,uth,uto,thet,qq,etas,dGdT,pzp,Vdeb,gamdeb
        common /chempot/ cpa
        common /mag/ icfe
        common /atomc/ stox,wox,wcomp,stcomp,atom,comp
        common /tfindc/ tlast,vtarg,starg,phugo,vhugo,ehugo,dvdpmol,dsdtmol,dhdpmol,dhdtmol,wmaggc,chcalc,adcalc,hucalc,nvet,nvep
        common /phycom/ iphyflag
        data ncall/0/
	data gphysub/0.0/
	call cpu_time(start)
        ncall = ncall + 1
        ent = 0.0
	iphyflag = 1

C Headers
        if (ncall .eq. 1) then
        phname(nphasep) = '                    '
        write(89,'(a4,a7,a8,12a11)') 'spec','Pi','Ti','Enthalpy','Entropy','Cp','Gibbs','Volume','Smix','Sconf','CV'
        write(66,'(99a12)') 'Pi','depth','Ti',(phname(iph)(1:8),iph=1,nph)
        write(661,'(99a12)') 'Pi','depth','Ti',(phname(iph)(1:8),iph=1,nph)
        write(67,*) '  Pi    depth   Ti',
     &   phname(nphasep)(1:7),(phname(iph)(1:8),iph=1,nph)
        write(58,'(99a12)') 'Pi','depth','Ti','rho','KS','G','VBh','VSh','VPh','KSr','Gr','VBr','VSr','VPr'
     &                                                                        ,'KSv','Gv','VBv','VSv','VPv'
        write(59,'(99a14)') 'Pi','depth','Ti','Vol','KS','KT','alpha','Heat C','thet','g','q','Velocity','Ppart1','Ppart2'
        write(61,'(a6,a9,a9,5x,99a12)') 'Pi','depth','Ti',
     &   ('rh'//phname(iph)(1:11),iph=1,nph)
        write(62,*) '  Pi    depth   Ti',
     &   phname(nphasep)(1:7),('vb'//phname(iph)(1:11),iph=1,nph)
        write(63,*) '  Pi    depth   Ti',
     &   phname(nphasep)(1:7),('vs'//phname(iph)(1:11),iph=1,nph)
        write(64,*) '  Pi    depth   Ti',
     &   phname(nphasep)(1:7),('vp'//phname(iph)(1:11),iph=1,nph)
        write(65,*) '  Pi    depth   Ti',
     &   phname(nphasep)(1:7),('wm'//phname(iph)(1:11),iph=1,nph)
        write(68,*) '  Pi    depth   Ti',
     &   phname(nphasep)(1:7),('vol'//phname(iph)(1:11),iph=1,nph)
        end if

C  Aggregate Properties
        entagg = 0.0
	freeagg = 0.0
        volagg = 0.0
        wmagg = 0.0
        cpagg = 0.0
	cvagg = 0.0
        baggr = 0.0
        btaggr = 0.0
        gaggr = 0.0
        baggv = 0.0
        btaggv = 0.0
        gaggv = 0.0
        alpagg = 0.0
        delagg = 0.0
        fnagg = 0.0
        vabsagg = 0.0
        vabsoagg = 0.
        dgdtaggv = 0.
        dgdtaggr = 0.
        smix = 0.0
	epsilon1 = 0.0
	epsilon2 = 0.0
	vdebagg = 0.0
	Gwu = 0.
	Gironmin = +1.e15
	do 222 iiron=1,3
222	Giron(iiron) = 0.

	call nform(nnew,n,n1,q2,nspec,nnull)

	do 152 ispec=1,nspec
C  Following line is added to aid computation of moduli which diverge if n of any constituent species is exactly zero.
	 if (absents(ispec)) n(ispec) = nsmall/10.
C  reverse the signs of dmdt and dmdp here since e.g. H (dn/dP) = -dmu/dP by the triple product rule
	 dmdt(ispec) = sspeca(ispec)/1000.
	 dmdp(ispec) = -vspeca(ispec)
152	continue
C  Project temperature derivative of the chemical potential dm/dt
	call dgemv('Transpose q2',nspec,nnull,one,q2,nspecp,dmdt,ione,zero,dmdtp,ione)
C  Project pressure derivative of the chemical potential dm/dp
	call dgemv('Transpose q2',nspec,nnull,one,q2,nspecp,dmdp,ione,zero,dmdpp,ione)
C  Compute projected hessian
	call hessfunc(nnew,hespro)
C  Solve linear problem: H^P dn^P/dT = dm^P/dT
	call svdsub(nnull,nnull,hespro,nspecp,nspecp,dmdtp,vwork,vwork,dndtp,nnulls)
C  Solve linear problem: H^P dn^P/dP = dm^P/dP
	call svdsub(nnull,nnull,hespro,nspecp,nspecp,dmdpp,vwork,vwork,dndpp,nnulls)
C  Form dn/dT from dn^P/dT
	call nform(dndtp,dndt,zervec,q2,nspec,nnull)
C  Form dn/dP from dn^P/dP
	call nform(dndpp,dndp,zervec,q2,nspec,nnull)
C  Check results by computing mu_i dn_i/dT and mu_i dn_i/dP, which should be zero.
	checkt = ddot(nspec,cpa,ione,dndt,ione)
	checkp = ddot(nspec,cpa,ione,dndp,ione)
C  Compute dndt using alternative procedure as detailed in thermal expansivity paper
C  This may not be the most efficient way to compute dndt, because it involves an nspecxnspec matrix, rather than an nnullxnnull matrix (the projected Hessian).
C  However, the following method looks better on paper, and permits a more direct connection to V_i and S_i
C  Inverse of projected Hessian
	call dlacpy('All',nnull,nnull,hespro,nspecp,swork,nspecp)
	call dgetrf(nnull,nnull,swork,nspecp,ipiv,info)
c	write(31,*) 'dgetrf info',info
	call dgetri(nnull,swork,nspecp,ipiv,wwork,nspecp,info)
c	write(31,*) 'dgetri info',info
c	write(31,*) 'q2'
c	do 1542 i=1,nspec
c	 write(31,*) (q2(i,j),j=1,nnull)
c1542	continue
c	write(31,*) 'Inverse of projected Hessian'
c	do 154 i=1,nnull
c	 write(31,*) (swork(i,j),j=1,nnull)
c154	continue
	call dgemm('Normal','Transpose q2',nnull,nspec,nnull,one,swork,nspecp,q2,nspecp,zero,vwork,nspecp)
c	do 1541 i=1,nspec
c	 write(31,*) (vwork(i,j),j=1,nspec)
c1541	continue
	call dgemm('Normal q2','Normal',nspec,nspec,nnull,one,q2,nspecp,vwork,nspecp,zero,swork,nspecp)
c	write(31,*) 'New projected Hessian'
c	do 155 i=1,nspec
c	 write(31,*) (swork(i,j),j=1,nspec)
c155	continue
	call dgemv('Normal',nspec,nspec,one,swork,nspecp,dmdt,ione,zero,dndttest,ione)
	do 1561 ispec=1,nspec
	 dnatomdp(ispec) = apar(ispec,1)*dndp(ispec)
1561	continue
	call dgemv('Normal',nph,nspec,one,f,nphasep,dnatomdp,ione,zero,dfdp,ione)
c	write(31,*) 'dndt',(dndt(ispec),ispec=1,nspec)
c	write(31,*) 'dndttest',(dndttest(ispec),ispec=1,nspec)
c	write(31,*) 'dndp',(dndp(ispec),ispec=1,nspec)
c	write(31,*) 'dmdt',(dmdt(ispec),ispec=1,nspec)
c	write(31,*) 'dmdp',(dmdp(ispec),ispec=1,nspec)
c	write(31,*) 'dndt projected',(dndtp(i),i=1,nnull)
c	write(31,*) 'dmdt projected',(dmdtp(i),i=1,nnull)
c	write(31,*) 'dfdp',(dfdp(i),i=1,nph)
c	write(31,*) 'm',(cpa(ispec),ispec=1,nspec)
c	write(31,*) 'checkt mu_i dn_i/dT',checkt
c	write(31,*) 'checkp mu_i dn_i/dP',checkp
c	write(31,*) 'vspeca',(vspeca(ispec),ispec=1,nspec)
c	write(31,*) 'sspeca',(sspeca(ispec),ispec=1,nspec)
C  Compute dndpfast and dndtfast which include only rapidly transforming phases
	isign = -1.0
	nfast = 0
	do 158 i=1,nspecp
	 do 158 j=1,nspecp
	  q2save(i,j) = q2(i,j)
158	continue
	do 156 iph=1,nph
	 do 157 ispec=1,nspec
	  if (f(iph,ispec) .eq. 0.) go to 157
	  q2(ispec,ione) = 0.
          if (ispec .ge. iophase(iph) .and. ispec .le. iophase(iph)+mophase(iph) - 1) then
	   nfast = nfast + 1
	   isign = -1.0*isign
	   q2(ispec,ione) = float(isign)/sqrt(2.)
	  end if
157	 continue
156	continue
	nnullfast = 1
	if (nfast .gt. 2) write(31,*) 'WARNING: nfast greater than two',nfast
C  Project temperature derivative of the chemical potential dm/dt
	call dgemv('Transpose q2fast',nspec,nnullfast,one,q2,nspecp,dmdt,ione,zero,dmdtp,ione)
C  Project pressure derivative of the chemical potential dm/dp
	call dgemv('Transpose q2fast',nspec,nnullfast,one,q2,nspecp,dmdp,ione,zero,dmdpp,ione)
C  Compute projected hessian
	call hessfunc(nnew,hespro)
C  Solve linear problem: H^P dn^P/dT = dm^P/dT
	call svdsub(nnullfast,nnullfast,hespro,nspecp,nspecp,dmdtp,vwork,vwork,dndtp,nnulls)
C  Solve linear problem: H^P dn^P/dP = dm^P/dP
	call svdsub(nnullfast,nnullfast,hespro,nspecp,nspecp,dmdpp,vwork,vwork,dndpp,nnulls)
C  Form dn/dT from dn^P/dT
	call nform(dndtp,dndtfast,zervec,q2,nspec,nnull)
C  Form dn/dP from dn^P/dP
	call nform(dndpp,dndpfast,zervec,q2,nspec,nnull)
c	write(31,*) 'dndtfast',(dndtfast(ispec),ispec=1,nspec)
c	write(31,*) 'dndpfast',(dndpfast(ispec),ispec=1,nspec)

	do 159 i=1,nspecp
	 do 159 j=1,nspecp
	  q2(i,j) = q2save(i,j)
159	continue

        do 1 iph=1,nph

C  Phase Properties
         volph = 0.0
         volpho = 0.0
         rhoph = 0.0
         bukph = 0.0
         buktph = 0.0
         gshph = 0.0
         alpph = 0.0
	 cpph = 0.0
         delph = 0.0
         wmph = 0.0
         fnph = 0.0
         entph = 0.0
         henph = 0.0
         gibph = 0.0
	 entoph = 0.0
	 henoph = 0.0
	 giboph = 0.0
         gsolph = 0.
         hsolph = 0.
         dgdtph = 0.
         rhopha(iph) = 0.
         vbpha(iph) = 0.
         vspha(iph) = 0.
         vppha(iph) = 0.
         volpha(iph) = 0.
         fnpha(iph) = 0.
         entpha(iph) = 0.
         henpha(iph) = 0.
         gibpha(iph) = 0.
         wmpha(iph) = 0.
         vabs(iph) = 0.
         vabso(iph) = 0.
         gpha(iph) = 0.
         bpha(iph) = 0.
         do 2 ispec=1,nspec
          if (f(iph,ispec) .eq. 0) go to 2
          call parset(ispec,apar,fn,zu,wm,To,Fo,Vo,Ko,Kop,Kopp,
     &                    wd1,wd2,wd3,ws1,ws2,ws3,
     &                    we1,qe1,we2,qe2,we3,qe3,we4,qe4,wou,wol,
     &                    gam,qo,be,ge,q2A2,
     &                    htl,ibv,ied,izp,
     &                    Go,Gop,Got)
          do 13 jspec=1,nspec
13        ntemp(jspec) = 0.
          ntemp(ispec) = 1.
          call cp(ispec,ntemp,chempot,rsum,volsum,sconf,smag)
          sconf = -chempot/Ti
          fdumm = gspec(ispec)
          call cp(ispec,n,chempot,rsum,volsum,smixi,smag)
c         print*, 'Regular solution term',phname(iph)(1:5),sname(ispec)(1:4),rsum,ent
          smix = smix + n(ispec)*smixi
c  with electronic contribution
c         if (icfe .ne. 0) smag = s(icfe,ispec)*Rgas*(log(2.*2. + 1.) + 3.*log(3.))
          Ftot = Ftot - Ti*smag
          ent = ent + smag
c          entagg = entagg + n(ispec)*ent
	  entagg = entagg + n(ispec)*sspeca(ispec)
c	print*, 'in physub',ispec,n(ispec),sspeca(ispec),entagg
	  freeagg = freeagg + n(ispec)*(fdumm + chempot)
	  volxs = vol + volsum
          volagg = volagg + n(ispec)*vspeca(ispec)
c	write(31,*) 'volumes in physub',ispec,n(ispec),vspeca(ispec),vol,volsum,volxs
          wmagg = wmagg + n(ispec)*wm
          wmph = wmph + n(ispec)*wm
          if (Ks .gt. 0.)  bukph = bukph + n(ispec)*vol/Ks
          if (K .gt. 0.)   buktph = buktph + n(ispec)*vol/K
          gshph = gshph + n(ispec)*vol/Gsh
c	write(31,*) 'gshph calc',iph,ispec,n(ispec),vol,Gsh,gshph
          volph = volph + n(ispec)*volxs
          volpho = volpho + n(ispec)*Vo
          fnph = fnph + n(ispec)
          if (.not. absents(ispec)) cpagg = cpagg + n(ispec)*Cap
          cvagg = cvagg + n(ispec)*Cv
          if (.not. absents(ispec)) alpagg = alpagg + n(ispec)*vol*alp
          delagg = delagg + n(ispec)*vol*alp*(1. + deltas)/Ks
	  vdebagg = vdebagg + n(ispec)*vol*vdeb
          if (.not. absents(ispec)) alpph = alpph + n(ispec)*vol*alp
	  if (.not. absents(ispec)) cpph = cpph + n(ispec)*Cap
          delph = delph + n(ispec)*vol*alp*(1. + deltas)/Ks
          Gtot = Ftot
          Htot = Gtot + Ti*ent
          Etot = Htot - 1000.*Pi*vol
          entph = entph + n(ispec)*(ent + smixi)
          henph = henph + n(ispec)*cpa(ispec) + Ti*n(ispec)*(ent + smixi)/1000.
          gibph = gibph + n(ispec)*cpa(ispec)
	  entoph = entoph + n(ispec)*ent
          henoph = henoph + n(ispec)*Htot/1000.
	  giboph = giboph + n(ispec)*Gtot/1000.
          gsolph = gsolph + n(ispec)*(Gtot + chempot)
          dgdtph = dgdtph + n(ispec)*vol/Gsh**2*dGdT
	  epsilon1 = epsilon1 + n(ispec)*wm
	  epsilon2 = epsilon2 + n(ispec)*wm**2
	  if (sname(ispec) .eq. 'wu  ') Gwu = Gtot/2.
	  if (sname(ispec) .eq. 'fea ') Giron(1) = Gtot
	  if (sname(ispec) .eq. 'feg ') Giron(2) = Gtot
	  if (sname(ispec) .eq. 'fee ') Giron(3) = Gtot
	  do 223 iiron=1,3
c	   write(31,*) 'Find minimum iron Gibbs free energy'
c	   write(31,*) iiron,Giron(iiron),Gironmin,Gwu
	   if (Giron(iiron) .eq. 0.) go to 223
	   Gironmin = min(Gironmin,Giron(iiron))
223	  continue
          if (iprint .eq. 1) then
                write(89,'(i3,a4,f7.2,f8.2,12f11.4)') ispec,sname(ispec),Pi,Ti,Etot/1000.,ent,Cap
     ;     ,Gtot/1000.,vol,smixi,sconf,Cv,thet
c                write(89,'(i3,a4,f7.2,f8.2,12f11.4)') ispec,sname(ispec),Pi,Ti,Gtot/1000.-Pi*vol,ent,Cap
c     ;     ,Gtot/1000.,vol,smixi,sconf,Cv,thet
c     ;    write(89,'(i4,f7.2,f8.2,22f11.4)') ispec,Pi,Ti,Etot/1000.,ent,Cap
c     ;     ,Gtot/1000.,vol,smixi,sconf,Cv,K
	  end if
2        continue

C Include only rapidly transforming phases in metamorphic terms
C Phase-wise metamorphic component is used only for the calculation of the adiabatic bulk modulus
	 alpmet = 0.
	 cpmet = 0.
	 bmet = 0.
         do 151 ispec=1,nspec
          if (ispec .ge. iophase(iph) .and. ispec .le. iophase(iph)+mophase(iph) - 1) then
           alpmet = alpmet + dndtfast(ispec)*vspeca(ispec)
           cpmet = cpmet + Ti*dndtfast(ispec)*sspeca(ispec)
           bmet = bmet + dndpfast(ispec)*dmdp(ispec)
          end if
151      continue

C  Previous calculation uses V averaging for K which is inconsistent with
C  calculation of delta as coded.  Changed to R averaging 8/31/01.
         bukphreuss = volph/bukph
         buktphtot = volph/(buktph + bmet)
         buktph = volph/buktph
         gshph = volph/gshph
         hsolph = gsolph + Ti*entph
         alpphtot = alpph/volph + alpmet/volph
         alpph = alpph/volph
	 cpphtot = cpph/wmph + cpmet/wmph
	 cvphtot = cpphtot - Ti*volph*alpphtot**2*buktphtot/wmagg*1000
	 if (cvphtot .gt. 0.) gamphtot = volph*alpphtot*buktphtot/(cvphtot*wmagg)*1000.
	 bukph = buktphtot*(1. + alpphtot*gamphtot*Ti)
c	 write(31,*) 'Phase properties',iph,cpphtot,cvphtot,gamphtot,buktph,volph,alpph,bukph,bukphreuss,wmph
	 if (alpmet .ne. 0.) then
          write(31,*) 'Fast alpmet,alpagg,alptot',alpmet/volagg,alpph,alpphtot
          write(31,*) 'Fast Pi,X_Fe,bmet,btaggh,btot,bstot',Pi,b(2),bmet,buktph,buktphtot,bukph
          write(31,*) 'Fast cpmet,cpagg,cptot,cvtot,gamtot',cpmet/wmph,cpph/wmph,cpphtot,cvphtot,gamphtot
	 end if
         delph = delph*bukph/(volph*alpph) - 1.
         rhoph = wmph/volph
         dlnvbdtph = alpph/2.*(delph - 1.)
         dgdtph = gshph**2*dgdtph/volph
         rhopha(iph) = rhoph
         vbpha(iph) = asqrt(bukph/rhoph)
         vspha(iph) = asqrt(gshph/rhoph)
         vppha(iph) = asqrt((bukph + 4./3.*gshph)/rhoph)
         volpha(iph) = volph/fnph
         if (fnph .gt. Tsmall) then
          entpha(iph) = entph/fnph
          henpha(iph) = henph/fnph/1000.
          gibpha(iph) = gibph/fnph/1000.
	  entopha(iph) = entoph/fnph
	  henopha(iph) = henoph/fnph/1000.
	  gibopha(iph) = giboph/fnph/1000.
          fnpha(iph) = fn*fnph
         end if
         wmpha(iph) = wmph
         vabs(iph) = volph
         vabso(iph) = volpho
         gpha(iph) = gshph
         bpha(iph) = bukph
         btpha(iph) = buktph
         dgdtph = dgdtph*gshph**2
C  Accumulate aggregate properties
         fnagg = fnagg + fnpha(iph)
         vabsagg = vabsagg + vabs(iph)
         vabsoagg = vabsoagg + vabso(iph)
         baggv = baggv + volph*bukph
         baggr = baggr + volph/bukph
         btaggv = btaggv + volph*buktph
         btaggr = btaggr + volph/buktph
         gaggv = gaggv + volph*gshph
         gaggr = gaggr + volph/gshph
c	 write(31,*) 'calculate gaggr from gshph',iph,volph,gshph,gaggr
         dgdtaggv = dgdtaggv + volph*dgdtph
         dgdtaggr = dgdtaggr + volph/gshph**2*dgdtph
1       continue
        gmin = +1.e15
        gmax = -1.e15
        bmin = +1.e15
        bmax = -1.e15
        vsmin = +1.e15
        vsmax = -1.e15
        vpmin = +1.e15
        vpmax = -1.e15
        do 11 iph=1,nph
         if (vabs(iph) .eq. 0.) go to 11
         gmin = min(gmin,gpha(iph))
         gmax = max(gmax,gpha(iph))
         bmin = min(bmin,bpha(iph))
         bmax = max(bmax,bpha(iph))
         vsmin = min(vsmin,vspha(iph))
         vsmax = max(vsmax,vspha(iph))
         vpmin = min(vpmin,vppha(iph))
         vpmax = max(vpmax,vppha(iph))
11      continue
        flambdap = 0.
        flambdam = 0.
        fgammap = 0.
        fgammam = 0.
        ximin = gmin/6.*((9.*bmin + 8.*gmin)/(bmin + 2.*gmin))
        ximax = gmax/6.*((9.*bmax + 8.*gmax)/(bmax + 2.*gmax))
        do 12 iph=1,nph
         if (vabs(iph) .eq. 0.) go to 12
         flambdap = flambdap + vabs(iph)/volagg/(bpha(iph) + 4./3.*gmax)
         flambdam = flambdam + vabs(iph)/volagg/(bpha(iph) + 4./3.*gmin)
         fgammap = fgammap + vabs(iph)/volagg/(ximax + gpha(iph))
         fgammam = fgammam + vabs(iph)/volagg/(ximin + gpha(iph))
12      continue
        bhashp = 1./flambdap - 4./3.*gmax
        bhashm = 1./flambdam - 4./3.*gmin
        ghashp = 1./fgammap - ximax
        ghashm = 1./fgammam - ximin
	wmaggc = wmagg
        rho = wmagg/volagg
        baggv = baggv/volagg
        baggr = volagg/baggr
        btaggv = btaggv/volagg
        btaggr = volagg/btaggr
        gaggv = gaggv/volagg
        gaggr = volagg/gaggr
        baggh = (baggv + baggr)/2.
c	write(31,*) 'baggh calc from baggv,baggr volagg',baggv,baggr,volagg
        btaggh = (btaggv + btaggr)/2.
        gaggh = (gaggv + gaggr)/2.
        Vbr = asqrt(baggr/rho)
        Vsr = asqrt(gaggr/rho)
        Vpr = asqrt((baggr + 4./3.*gaggr)/rho)
        Vbv = asqrt(baggv/rho)
        Vsv = asqrt(gaggv/rho)
        Vpv = asqrt((baggv + 4./3.*gaggv)/rho)
        Vbh = (Vbv + Vbr)/2.
        Vsh = (Vsv + Vsr)/2.
        Vph = (Vpv + Vpr)/2.
c	Vdeb = (2./3./Vsh**3 + 1./3./Vph**3)**(-1./3.)
	Vdeb3 = 2./3./Vsh**3 + 1./3./Vph**3
	Vdebold = 1./Vdeb3**(1./3.)
	if (Vdeb3 .lt. 0.) Vdebold = 1./abs(Vdeb3)**(1./3.)
	dlnvsdlnv = 0.5*(1. - btaggh/gaggh*Gshp) 
	dlnvpdlnv = 0.5*(1. - btaggh/(baggh + 4./3.*gaggh)*(Ksp + 4./3.*Gshp))
	dlnvdebdlnv = 2./3.*(Vdeb/Vsh)**3*dlnvsdlnv + 1./3.*(Vdeb/Vph)**3*dlnvpdlnv
	gamdeb = -(dlnvdebdlnv - 1./3.)
	gamp= -(dlnvpdlnv - 1./3.)
c	print*, Ksp,Gshp,dlnvsdlnv,dlnvpdlnv,dlnvdebdlnv,gamdeb
        Vshashp = asqrt(ghashp/rho)
        Vshashm = asqrt(ghashm/rho)
        Vshasha = asqrt(0.5*(ghashp + ghashm)/rho) 
        Vphashp = asqrt((bhashp + 4./3.*ghashp)/rho)
        Vphashm = asqrt((bhashm + 4./3.*ghashm)/rho)
        Vphasha = asqrt(0.5*(bhashp + 4./3.*ghashp + bhashm + 4./3.*ghashm)/rho) 
        alpagg = alpagg/volagg
        delagg = delagg*baggr/(volagg*alpagg) - 1.
        dlnvbdt=alpagg/2.*(delagg - 1.)
        cpagg = cpagg/wmagg
        cvagg = cvagg/wmagg
        gruagg = 1000.*baggh*alpagg/(rho*cpagg)
        dgdtaggv = dgdtaggv/volagg
        dgdtaggr = dgdtaggr/volagg*gaggr**2
        dgdtagg = 0.5*(dgdtaggv + dgdtaggr)
        dlnvsdt = 0.5*(dgdtagg/gaggh + alpagg)
        dlnvpdt = 0.5*((-delagg*alpagg*baggh + 4./3.*dgdtagg)/(baggh+4./3.*gaggh) + alpagg)
	rhoo = wmagg/vabsoagg
	vdebagg = vdebagg/volagg

        if (Ti .eq. Tsmall) delagg = 0.
        if (Ti .eq. Tsmall) dlnvbdt = 0.
        ent = entagg 
c        entagg = entagg + smix
	enthagg = freeagg + Ti*entagg
        eintagg = freeagg + Ti*entagg - 1000.*Pi*volagg
c	print*, rho,ent,entagg,smix,wmagg,entagg/wmagg,enthagg/wmagg/1000.,freeagg/wmagg/1000.,freeagg/1000.
c-> Include configurational entropy in computation of self-consistent isentrope
	ent = entagg
	enth = enthagg
c<-

	alpmet = 0.
	cpmet = 0.
	bmet = 0.
c	write(31,*) 'volagg',volagg
c	write(31,*) 'wmagg',wmagg
	do 15 ispec=1,nspec
	 alpmet = alpmet + dndt(ispec)*vspeca(ispec)/volagg
	 cpmet = cpmet + Ti*dndt(ispec)*sspeca(ispec)
	 bmet = bmet + dndp(ispec)*dmdp(ispec)
15	continue
C  Compute effective derivative properties including metamorphic contribution
	btot = volagg/(volagg/btaggr + bmet)
	alptot = alpagg + alpmet
	cptot = cpagg + cpmet/wmagg
	cvtot = cptot - Ti*volagg*alptot**2*btot/wmagg*1000.
	gamtot = volagg*alptot*btot/(cvtot*wmagg)*1000.
	bstot = btot*(1. + alptot*gamtot*Ti)
C  Compute isomorphic contributions to derivative properties
	btiso = btaggr
	cpiso = cpagg
	alpiso = alpagg
	cviso = cpiso - Ti*volagg*alpiso**2*btiso/wmagg*1000.
	gamiso = volagg*alpiso*btiso/(cviso*wmagg)*1000.
	bsiso = btiso*(1. + alpiso*gamiso*Ti)
	deltavol = ddot(nspec,dndp,ione,vspeca,ione)
	deltaent = ddot(nspec,dndp,ione,sspeca,ione)
	dfdpabstotal = 0.
	dfdptotal = 0.
	fnaggtransform = 0.
	do 1562 iph=1,nph
	 if (dfdp(iph) .gt. dfdpsmall) dfdptotal = dfdptotal + dfdp(iph)
	 dfdpabstotal = dfdpabstotal + abs(dfdp(iph))
	 fnaggtransform = fnaggtransform + 2.*fnpha(iph)*abs(dfdp(iph))
c	 write(31,*) 'do 1562',iph,fnpha(iph),dfdp(iph),dfdptotal,dfdpabstotal,fnaggtransform
1562	continue
	fnaggtransform = fnaggtransform/dfdpabstotal
	if (fnaggtransform .ne. 0.) dfdptotal = dfdptotal/fnaggtransform
c	write(31,*) 'After 1562',dfdptotal,dfdpabstotal,fnagg,fnaggtransform
	phasebuoyancyparameter = 0.
	if (dfdptotal .gt. 0.) phasebuoyancyparameter = alpmet/alpiso/dfdptotal/pmantle
	ClapeyronSlope = 0.
	if (abs(deltavol) .gt. asmall .and. phasebuoyancyparameter .ne. 0.) ClapeyronSlope = deltaent/deltavol
	write(69,'(99f25.16)') Pi,depth(Pi),Ti,phasebuoyancyparameter,ClapeyronSlope
     &                         ,dfdptotal,deltaent,deltavol,alpagg,alpmet,alpmet/alpagg,alptot,alptot/alpagg,cvtot,bstot
	write(31,*) 'alpmet,alpagg,alptot',alpmet,alpagg,alptot,phasebuoyancyparameter,dfdptotal,pmantle
	write(31,*) 'Delta S, Delta V, Gamma = ',deltaent,deltavol,ClapeyronSlope
	write(31,*) 'Pi,X_Fe,bmet,btaggh,btot,bstot',Pi,b(2),bmet,btaggh,btot,bstot
	write(31,*) 'cpmet,cpagg,cptot,cvtot,gamtot',cpmet/wmagg,cpagg,cptot,cvtot,gamtot
	dhdtmol = cptot*wmagg/1000.
	dsdtmol = cptot*wmagg/Ti
	dhdtmol = cpagg*wmagg/1000.
	dsdtmol = cpagg*wmagg/Ti

c	print*, 'in physub',ent,smix
        call qr19(depth(Pi),Ti,qs,qp)
        call vred(qs,qp,vsred,vpred)
        if (iprint .eq. 1) then
	fnphamax = -1.e15
	iimax = nphasep
	do 77 ii=1,nph
	 fnphamax=max(fnphamax,fnpha(ii))
	 if (fnphamax .eq. fnpha(ii)) iimax = ii
77	continue
	tmelt = tlindeman(vol,wmagg,fnagg,thet)
	wmav = wmagg/fnph
	epsilon = (epsilon2 - 2.*wmav*epsilon1 + fnph*wmav**2)/(fnph*wmav**2)
        call thetacal(cvagg*wmagg/(3.*fnagg*Rgas),tcalagg)
	tcalagg = tcalagg*Ti
C  Units in fort.56:
C  P(GPa) depth(km) T(K) rho(g/cm^3) VB(km/s) VS(km/s) VP(km/s) VSQ(km/s) VPQ(km/s) H(kJ/g) S(J/g/K) alpha(1e5 K^-1) cp(J/g/K) KT(GPa) Qs(-) Qp(-) rho_0(g/cm^3) dominant_phase
c        write(56,'(17f25.16,2x,a5)') Pi,depth(Pi),Ti,rho,Vbh,Vsh,Vph,Vsh*vsred,
c     &   Vph*vpred,enthagg/wmagg/1000.,entagg/wmagg,1.e5*alptot,cptot,btot,qs,qp,rhoo,phname(iimax)
        write(56,'(17f25.16,2x,a5)') Pi,depth(Pi),Ti,rho,Vbh,Vsh,Vph,Vsh*vsred,
     &   Vph*vpred,enthagg/wmagg/1000.,entagg/wmagg,1.e5*alptot,cptot,btot,qs,qp,rhoo,phname(iimax)
c       write(59,500) Pi,depth(Pi),Ti,vol,baggh,btaggh,1.e5*alpagg,cpagg,
c     &   gruagg,n(1),hsolph/1000.
c        write(59,500) Pi,depth(Pi),Ti,vol/81.8,baggh,btaggr,1.e5*alpagg,Cv/(fn*Rgas),cpagg*wmagg,
c     &   gruagg,qq,ph
c       write(599,500) Pi,depth(Pi),Ti,vol,baggh,btaggh,1.e5*alpagg,cpagg*wmagg,
c     &   -1000.*alpagg*baggh*delagg,1000.*dgdtagg,-dlnvsdt*1.e5,-dlnvpdt*1.e5,deltas
c       write(59,'(99f14.5)') Pi,depth(Pi),Ti,volagg,baggh,btaggh,1.e5*alpagg,cvagg*wmagg,
c     &   thet,gruagg,qq,Vdeb,ph,pzp,tmelt
       write(59,'(99f14.5)') Pi,depth(Pi),Ti,volagg,bsiso,btiso,1.e5*alpiso,cpiso,
     &   thet,gamiso,qq,Vdeb,ph,pzp,tmelt
        write(57,500) Pi,depth(Pi),Ti,volagg,fn,zu,1000.*Vdebold,gruagg,qq,wmav,epsilon,cvagg*wmagg/(3.*fnagg*Rgas)
     &    ,tcalagg,gamdeb,fn,fnph,fnagg,Cv/(3*fn*Rgas),thet
        write(58,'(99f12.5)') Pi,depth(Pi),Ti,rho,baggh,gaggh,Vbh,Vsh,Vph,baggr,gaggr,Vbr,Vsr,Vpr,baggv,gaggv,Vbv,Vsv,Vpv
c     &   ,Vshasha,Vphasha,Vshashp,Vshashm,Vphashp,Vphashm,Vdeb
c     &   ,vsmin,vsmax,vpmin,vpmax
        write(61,500) Pi,depth(Pi),Ti,(exp(log(rhopha(iph)-asmall)),iph=1,nph)
        write(62,500) Pi,depth(Pi),Ti,(exp(log(vbpha(iph)-asmall)),iph=1,nph)
        write(63,500) Pi,depth(Pi),Ti,(exp(log(vspha(iph)-asmall)),iph=1,nph)
        write(64,500) Pi,depth(Pi),Ti,(exp(log(vppha(iph)-asmall)),iph=1,nph)
        write(65,500) Pi,depth(Pi),Ti,(exp(log((wmpha(iph)-asmall)/fnpha(iph))),iph=1,nph)
        write(68,500) Pi,depth(Pi),Ti,(exp(log(volpha(iph)-asmall)),iph=1,nph)
        write(66,700) Pi,depth(Pi),Ti,(fnpha(iph)/fnagg,iph=1,nph)
        write(661,700) Pi,depth(Pi),Ti,(dfdp(iph)/fnagg,iph=1,nph)
        write(67,700) Pi,depth(Pi),Ti,(vabs(iph)/vabsagg,iph=1,nph)
        write(31,'(/,a)') 'Phase proportions: atomic, volume, mass'
        write(31,'(99a12)') 'Pi','depth','Ti',(phname(iph)(1:8),iph=1,nph)
        write(31,700) Pi,depth(Pi),Ti,(fnpha(iph)/fnagg,iph=1,nph)
        write(31,700) Pi,depth(Pi),Ti,(vabs(iph)/vabsagg,iph=1,nph)
        write(31,700) Pi,depth(Pi),Ti,(wmpha(iph)/wmagg,iph=1,nph)
        write(31,'(/,a)') 'Moduli and Velocities'
        write(31,'(a6,a8,a8,20a12)') 'Pi','depth','Ti','rho','KS','G ','VBh','VSh','VPh','VBr','VSr','VPr','VBv','VSv','VPv'
        write(31,500) Pi,depth(Pi),Ti,rho,baggh,gaggh,Vbh,Vsh,Vph
     &                                     ,Vbr,Vsr,Vpr,Vbv,Vsv,Vpv
        write(31,'(/,a)') 'Physical Properties'
        write(31,'(a6,a8,a8,20a12)') 'Pi','depth','Ti','rho','KS','KT','alpha','cP','g','S(J/g)','H(kJ/g)'
c        write(31,500) Pi,depth(Pi),Ti,rho,baggh,btaggh,1.e5*alpagg,cpagg,
c     &   gruagg,entagg/wmagg,enthagg/wmagg/1000.,wmagg
        write(31,500) Pi,depth(Pi),Ti,rho,bstot,btot,1.e5*alptot,cptot,
     &   gamtot,entagg/wmagg,enthagg/wmagg/1000.,wmagg
        write(31,'(/,a)') 'Physical properties of individual phases'
        write(31,'(a5,24a13)') 'phase','rho','VB','VS','VP','G','vol','ent','H ','G','S0','H0','G0','volo','KT','KS'
        do 25 iph=1,nph
	 if (fnpha(iph) .lt. 1.e-12) go to 25
         write(31,800) phname(iph),rhopha(iph),vbpha(iph),vspha(iph),
     &    vppha(iph),gpha(iph),volpha(iph),entpha(iph),henpha(iph),gibpha(iph),
     &    entopha(iph),henopha(iph),gibopha(iph),vabso(iph)/fnpha(iph)*5,btpha(iph),bpha(iph)
25      continue
	do 251 ic=1,nc
	 if (comp(ic) .eq. 'O ') then
	  write(31,*) 'oxygen chemical potential = ',lagc(ic)
	  ao2 = 47.255
	  bo2 = -4.550e-4
	  co2 = 4.402e+5
	  do2 = -393.5
	  Tro2 = 298.15
	  STro2 = 205.147
	  ho2gas = 0. + ao2*(Ti - Tro2) + 0.5*bo2*(Ti**2 - Tro2**2) - co2*(1./Ti - 1./Tro2) + 2.*do2*(sqrt(Ti) - sqrt(Tro2))
	  so2gas = 205.15 + ao2*log(Ti/Tro2) + bo2*(Ti - Tro2) - 0.5*co2*(1./Ti**2 - 1./Tro2**2) 
     &            - 2.*do2*(1./sqrt(Ti) - 1./sqrt(Tro2))
	  go2gas = (ho2gas - Ti*so2gas + STro2*Tro2)/1000.
	  write(31,*) 'oxygen gas enthalpy, entropy, gibbs (kJ/mol) = ',ho2gas,so2gas,go2gas
	  write(31,*) '2.*lagc(ic) - go2gas = ',2.*lagc(ic) - go2gas
	  fug = log10(exp(1000.*(2.*lagc(ic) - go2gas)/(Rgas*Ti)))
	  fug = (1000.*(2.*lagc(ic) - go2gas)/(Rgas*Ti))/log(10.)
	  write(31,*) 'log_10 oxygen fugacity = ',fug
	 end if
251	continue
        write(31,'(/,a)') 'Phase Compositions - Cations'
        write(31,'(12a12)') (comp(i),i=1,nco),'XFe'
	do 26 iph=1,nph
	 do 29 ic=1,nco
	  st(ic) = 0.
	  wco(iph,ic) = 0.
	  if (comp(ic) .eq. 'Mg') img = ic
	  if (comp(ic) .eq. 'Fe') ife = ic
29	 continue
	 if (fnpha(iph) .lt. 1.e-12) go to 26
	 ssum = 0.
	 do 27 ispec=1,nspec
          if (f(iph,ispec) .eq. 0) go to 27
	  wsum = 0.
	  do 28 ic=1,nco
	   st(ic) = st(ic) + s(ic,ispec)*n(ispec)
	   if (ispec .eq. iphase(iph)) ssum = ssum + s(ic,ispec)
	   if (comp(ic) .ne. 'O ') wco(iph,ic) = wco(iph,ic) + s(ic,ispec)*n(ispec)*wcomp(ic)/stcomp(ic)
c	   if (iph .eq. 1) print '(3i5,22f12.5)', iph,ispec,ic,wco(iph,ic),s(ic,ispec),n(ispec),wcomp(ic),stcomp(ic)
c	   if (iph .eq. 1) print '(3i5,22f12.5)', iph,ispec,ic,s(ic,ispec),n(ispec),st(ic),ssum
28	  continue
27	 continue
	 stsum = 0.
	 wsum = 0.
	 do 281 ic=1,nco
	  stsum = stsum + st(ic)
	  if (comp(ic) .ne. 'O ') wsum = wsum + wco(iph,ic)
281	 continue
	 if (stsum .ne. 0.) then
	  do 282 ic=1,nco
	   st(ic) = st(ic)/stsum*ssum
	   if (comp(ic) .ne. 'O ') wco(iph,ic) = wco(iph,ic)/wsum
282	  continue
	 end if
	 xfe = 0.0
	 if (ife .gt. 0 .and. img .gt. 0) then
	  if ((st(ife)+st(img)) .gt. 0.) xfe = st(ife)/(st(ife)+st(img))
	 end if
	 write(31,'(a5,12f12.5)') phname(iph),(st(ic),ic=1,nco),xfe
26	continue

        write(31,'(/,a)') 'Phase Compositions - Mass Standard Oxides (%)'
        write(31,'(12a12)') (comp(i),i=1,nco)
	do 283 iph=1,nph
	 if (fnpha(iph) .lt. 1.e-12) go to 283
	 write(31,'(a5,22f12.5)') phname(iph),(100.*wco(iph,ic),ic=1,nco)
283	continue
c	print*, (wcomp(ic),ic=1,nco)
c	print*, (stcomp(ic),ic=1,nco)

        end if

	iphyflag = 0

	call cpu_time(finish)
	gphysub = gphysub + finish - start
	write(31,*) 'time in physub',gphysub,ncall

        return
500     format(f6.2,f9.2,f9.2,35f12.5)
c700     format(f7.3,f9.3,f9.2,105f8.4)
700     format(99f12.5)
800     format(a5,24f13.7)
        end
