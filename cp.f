        subroutine cp(ispec,ncp,chempot,rsum,volsum,smixi,smag)

C  Calculate Contribution to Chemical Potential of Species ispec

        include 'P1'
        include 'chem.inc'
        include 'const.inc'
        include 'lag.inc'
        include 'absent.inc'
        include 'theory.inc'

	integer ispec,ia,ib,ic,icfe,iph,jsp,kst,nsitecp
	double precision chempot,rsum,smixi,smag,apar,Pi,qa,qb,siza,sizb,sizi,sum1,sum2,Ti,wregsz,oregsz,wregsave
	double precision oregsave
	double precision ourlog,dkron,hev,vregsz,volsum,ntot,cplandau
	double precision Vi,qorder,Tc,glan,cplan,alplan,betlan,slan,vlan,x,bx
        character*2 atom(natomp),comp(natomp)
        character*80 phname(nphasep),sname(nspecp)
        double precision nkp,nkpr,nikp(ncompp),ncp(nspecp)
        double precision wreg(nphasep,nsitep,nspecp,nspecp),vreg(nphasep,nsitep,nspecp,nspecp)
	double precision wox(natomp),stox(natomp),wcomp(natomp),stcomp(natomp)
	double precision nsum,msum,ratio,cideal,ssum
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /regcom/ wreg,vreg
        common /mag/ icfe
        common /names/ phname,sname
        common /atomc/ stox,wox,wcomp,stcomp,atom,comp

        iph = lphase(ispec)
        chempot = 0.0
        smixi = 0.0
        rsum = 0.0
	volsum = 0.0
        nsitecp = nsite(iph)
        do 1 kst=1,nsitecp
         nkp = 0.0
	 nkpr = 0.0
         sum1 = 0.0
         sum2 = 0.0
         do 3 ic=1,nco
          nikp(ic) = 0.0
          do 4 jsp=1,nspec
c          if (absents(jsp)) go to 4
           nkp = nkp + f(iph,jsp)*r(ic,jsp,kst)*ncp(jsp)
C  --> Comment out the following two lines to eliminate the entropy of mixing
c	   if (jsp   .ge. iophase(iph) .and. jsp   .le. iophase(iph)+mophase(iph)-1 .and.
c     &         ispec .ge. iophase(iph) .and. ispec .le. iophase(iph)+mophase(iph)-1 .and. jsp .ne. ispec) then
c	    write(31,*) 'in cp Found order-disorder',ispec,jsp,iph,kst,iophase(iph),mophase(iph),iastate(ispec,jsp,kst)
c	    go to 4
c	   end if
C <--
	   if (iastate(ispec,jsp,kst)) then	
c	    write(31,*) 'in cp Found iastate',ispec,jsp,kst,iph,iophase(iph),mophase(iph),iastate(ispec,jsp,kst)
            go to 4
	   end if
           nikp(ic) = nikp(ic) + f(iph,jsp)*r(ic,jsp,kst)*ncp(jsp)
4         continue
          nkpr = nkpr + nikp(ic)
          sum1 = sum1 + r(ic,ispec,kst)
          sum2 = sum2 + r(ic,ispec,kst)*ourlog(nikp(ic))
3        continue
         chempot = chempot - sum2
         chempot = chempot + sum1*ourlog(nkp)

	 if (kst .eq. 1) then
	  msum = 0.
	  nsum = 0.
	  ssum = 0.
	  do 51 ia=1,nspec
	   msum = msum + f(iph,ia)*apar(ia,41)*ncp(ia)
	   nsum = nsum + f(iph,ia)*abs(apar(ia,41))*ncp(ia)
	   ssum = ssum + f(iph,ia)*ncp(ia)
51        continue
	  ratio = msum/nsum
	  if (nsum .le. 0.) go to 52
	  do 5 ia=1,nspec-1
	   do 5 ib=ia+1,nspec
	    siza = apar(ia,41)
	    sizb = apar(ib,41)
	    sizi = apar(ispec,41)
	    qa = abs(siza)*ncp(ia)/nsum
	    qb = abs(sizb)*ncp(ib)/nsum
	    wregsave = wreg(iph,kst,ia,ib)
	    oregsave = wreg(iph,kst,ib,ia)
	    wreg(iph,kst,ia,ib) = wregsave + Pi*vreg(iph,kst,ia,ib)
	    wreg(iph,kst,ib,ia) = oregsave + Pi*vreg(iph,kst,ib,ia)
	    wregsz = (wreg(iph,kst,ia,ib) + wreg(iph,kst,ib,ia)*(qb - qa))*2.*abs(sizi)/(abs(siza) + abs(sizb))
	    oregsz = wreg(iph,kst,ib,ia)*2.*abs(sizi)/(abs(siza) + abs(sizb))
	    vregsz = (vreg(iph,kst,ia,ib) + vreg(iph,kst,ib,ia)*(qb - qa))*2.*abs(sizi)/(abs(siza) + abs(sizb))
c	    vregsz = vreg(iph,kst,ia,ib)*2.*abs(sizi)/(abs(siza) + abs(sizb))
	    wreg(iph,kst,ia,ib) = wregsave
	    wreg(iph,kst,ib,ia) = oregsave
	    rsum = rsum - f(iph,ia)*f(iph,ib)*(dkron(ispec,ia) - qa)*(dkron(ispec,ib) - qb)*wregsz*ratio
	    rsum = rsum + f(iph,ia)*f(iph,ib)*qa*qb*(dkron(ispec,ib) - qb - dkron(ispec,ia) + qa)*oregsz
	    volsum = volsum - f(iph,ia)*f(iph,ib)*(dkron(ispec,ia) - qa)*(dkron(ispec,ib) - qb)*vregsz*ratio
c	    rsum = rsum + f(iph,ia)*f(iph,ib)*qa*qb*(hev(sizi) - ratio)*wregsz
c	if (oregsz .ne. 0.) write(31,*) 'Asymmetric',ispec,ia,ib,sizi,siza,sizb,qa,qb,wregsz,msum,nsum,ratio,rsum
5	  continue
52	  continue
	 end if
1       continue

        smag = 0.0
        smixi = Rgas*chempot
	cideal = -Ti*Rgas*chempot
        chempot = -Ti*Rgas*chempot + 1000.*rsum - Ti*smag

C  Composition dependent Landau contribution

c	if (absent(iph)) return
	if (mphase(iph) .le. 1) return
	ntot = 0.
	do 31 jsp=1,nspec
	 ntot = ntot + f(iph,jsp)*n(jsp)
31	continue 

	cplandau = 0.
	Vi = 7.0
	do 32 jsp=1,nspec
	 if (f(iph,jsp) .eq. 0.) go to 32
	 x = n(jsp)/ntot
	 bx = -apar(jsp,38)
C Landau contributions
C Choose: landau for inv251010 and earlier, landauqr for later
	 if (iltyp .eq. 1) then
	  call landauqr(jsp,Vi,qorder,Tc,glan,cplan,alplan,betlan,slan,vlan)
	 else if (iltyp .eq. 2) then
	  call landau(jsp,Vi,qorder,Tc,glan,cplan,alplan,betlan,slan,vlan)
         else
          print*, 'WARNING: Landau type not chosen.  Landau terms are not computed.'
	 end if
	 cplandau = cplandau - x*apar(jsp,39)*bx*(dkron(ispec,jsp) - x)*(1. - qorder**2)
c	 write(31,*) 'in cplandau',ispec,jsp,cplandau,ntot,x,qorder
32	continue
c	write(31,*) 'in cplandau',ispec,cplandau,ntot,nphase(iph)

c	chempot = chempot + cplandau

        return
        end
