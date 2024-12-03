        subroutine hessian(ispec,ncp,hess)

C  Calculate ispec portion of hessian

        include 'P1'
        include 'chem.inc'
        include 'const.inc'
        include 'lag.inc'
        include 'absent.inc'

	integer ispec,ia,ib,ic,icfe,iph,j,jsp,jspec,kst,nsitecp
	double precision apar,chempot,dqadnj,dqbdnj,Pi,qa,qb,rsum,siza,sizb,sizi,smixi,sum1,sum2,Ti,wregsz,oregsave,wregsave
	double precision dkron,oregsz
        character*2 atom(natomp),comp(natomp)
        character*80 phname(nphasep),sname(nspecp)
        double precision nkp,nkpr,nikp(ncompp),ncp(nspecp)
        double precision wreg(nphasep,nsitep,nspecp,nspecp),vreg(nphasep,nsitep,nspecp,nspecp)
	double precision wox(natomp),stox(natomp),wcomp(natomp),stcomp(natomp)
	double precision nsum
	double precision sum1a(nspecp),sum2a(nspecp),rsuma(nspecp),hess(nspecp,nspecp)
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /regcom/ wreg,vreg
        common /mag/ icfe
        common /names/ phname,sname
        common /atomc/ stox,wox,wcomp,stcomp,atom,comp

        iph = lphase(ispec)
        chempot = 0.0
	do 9 jspec=1,nspec
	 hess(ispec,jspec) = 0.
	 rsuma(jspec) = 0.
9	continue
        smixi = 0.0
        rsum = 0.0
        nsitecp = nsite(iph)
        do 1 kst=1,nsitecp
         nkp = 0.0
	 nkpr = 0.0
         sum1 = 0.0
         sum2 = 0.0
	 do 11 j=1,nspec
	  sum1a(j) = 0.
	  sum2a(j) = 0.
11	 continue
         do 3 ic=1,nco
          nikp(ic) = 0.0
          do 4 jsp=1,nspec
           nkp = nkp + f(iph,jsp)*r(ic,jsp,kst)*ncp(jsp)
c	   if (jsp   .ge. iophase(iph) .and. jsp   .le. iophase(iph)+mophase(iph)-1 .and.
c     &         ispec .ge. iophase(iph) .and. ispec .le. iophase(iph)+mophase(iph)-1 .and. jsp .ne. ispec) go to 4
	   if (iastate(ispec,jsp,kst)) go to 4
           nikp(ic) = nikp(ic) + f(iph,jsp)*r(ic,jsp,kst)*ncp(jsp)
4         continue
          nkpr = nkpr + nikp(ic)
	  sum1 = sum1 + r(ic,ispec,kst)
	  do 6 jspec=1,nspec
	   sum1a(jspec) = sum1a(jspec) + r(ic,jspec,kst)*f(iph,jspec)
c	   if (jspec .ge. iophase(iph) .and. jspec .le. iophase(iph)+mophase(iph)-1 .and.
c     &         ispec .ge. iophase(iph) .and. ispec .le. iophase(iph)+mophase(iph)-1 .and. jspec .ne. ispec) go to 6
           if (iastate(ispec,jspec,kst)) go to 6
	   if (nikp(ic) .gt. 0.) sum2a(jspec) = sum2a(jspec) + r(ic,ispec,kst)*r(ic,jspec,kst)*f(iph,jspec)/nikp(ic)
6	  continue
3        continue
	 do 2 jspec=1,nspec
	  hess(ispec,jspec) = hess(ispec,jspec) - sum2a(jspec)
	  if (nkp .gt. 0.) hess(ispec,jspec) = hess(ispec,jspec) + sum1*sum1a(jspec)/nkp
2	 continue

	 if (kst .eq. 1) then
	  nsum = 0.
	  do 51 ia=1,nspec
	   nsum = nsum + f(iph,ia)*apar(ia,41)*ncp(ia)
51        continue
	  if (nsum .le. 0.) go to 52
	  do 5 ia=1,nspec-1
	   do 5 ib=ia+1,nspec
	    siza = apar(ia,41)
	    sizb = apar(ib,41)
	    sizi = apar(ispec,41)
	    qa = siza*ncp(ia)/nsum
	    qb = sizb*ncp(ib)/nsum
	    wregsave = wreg(iph,kst,ia,ib)
	    oregsave = wreg(iph,kst,ib,ia)
            wreg(iph,kst,ia,ib) = wregsave + Pi*vreg(iph,kst,ia,ib)
            wreg(iph,kst,ib,ia) = oregsave + Pi*vreg(iph,kst,ib,ia)
            wregsz = (wreg(iph,kst,ia,ib) + wreg(iph,kst,ib,ia)*(qb - qa))*2.*abs(sizi)/(abs(siza) + abs(sizb))
            oregsz = wreg(iph,kst,ib,ia)*2.*abs(sizi)/(abs(siza) + abs(sizb))
            wreg(iph,kst,ia,ib) = wregsave
            wreg(iph,kst,ib,ia) = oregsave
	    do 7 jspec=1,nspec
	     dqadnj = f(iph,jspec)*apar(ia,41)*(dkron(jspec,ia)*nsum - ncp(ia)*apar(jspec,41))/(nsum*nsum)
	     dqbdnj = f(iph,jspec)*apar(ib,41)*(dkron(jspec,ib)*nsum - ncp(ib)*apar(jspec,41))/(nsum*nsum)
	     rsuma(jspec) = rsuma(jspec) - f(iph,ia)*f(iph,ib)*((-dqbdnj*(dkron(ispec,ia) - qa) - dqadnj*(dkron(ispec,ib) - qb))
     &        *wregsz
     &        + (dkron(ispec,ia) - qa)*(dkron(ispec,ib) - qb)*wreg(iph,kst,ib,ia)
     &        *(dqbdnj - dqadnj)*2.*abs(sizi)/(abs(siza) + abs(sizb)))
	     rsuma(jspec) = rsuma(jspec) + f(iph,ia)*f(iph,ib)*oregsz*((dqadnj*qb + dqbdnj*qa)
     &        *(dkron(ispec,ib) - qb - dkron(ispec,ia) + qa) + qa*qb*(dqadnj - dqbdnj))
c	      write(31,*) ia,ib,ispec,jspec,rsuma(jspec),dqbdnj,dqadnj,dkron(ispec,ia),dkron(ispec,ib),qa,qb
7	    continue
5	  continue
52	  continue
	 end if
1       continue

	do 8 jspec=1,nspec
	 hess(ispec,jspec) = -Ti*Rgas*hess(ispec,jspec)/1000. + rsuma(jspec) 
c	 write(31,*) 'Hessian',ispec,jspec,hess(ispec,jspec),wregsz,qa,qb,f(iph,ispec),f(iph,jspec),dqadnj,dqbdnj
c     &     ,apar(ispec,41),apar(jspec,41),ncp(1),ncp(2),nsum
8	continue
c	write(31,*) 'Hessian',ispec,(-1000./(Rgas*Ti)*hess(ispec,jspec),jspec=1,nspec)
c	write(31,*) 'Hessian',ispec,(hess(ispec,jspec),jspec=1,nspec)

        return
        end
