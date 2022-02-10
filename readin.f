        subroutine readin(binit,bfinal,nbulk)

        include 'P1'
        include 'chem.inc'
        include 'const.inc'
        include 'lag.inc'
        include 'absent.inc'
        include 'theory.inc'
	include 'numpar.inc'

        integer nbulk,i,ia,ib,icfe,im,jm,iph,iphs,ispec
        integer iw,j,jco,jspec,k,ipar,natom,nbulk1,ncall,ndname,nmin,nw
	double precision apar,dw,Pi,Pinitial,Telastic,Ti,vdoswr,w,w1,w2,zeror,vdos
        logical lox
        character*80 line
        character*1 blank
        character*2 atom(natomp),comp(natomp),xatom
        character*80 dirname,fname,newfname,coxide
        character*80 phname(nphasep),sname(nspecp)
        double precision binit(ncompp),bfinal(ncompp)
        double precision wreg(nphasep,nsitep,nspecp,nspecp),vreg(nphasep,nsitep,nspecp,nspecp)
	double precision fmom(24)
	double precision wox(natomp),stox(natomp),wcomp(natomp),stcomp(natomp)
        common /names/ phname,sname
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /regcom/ wreg,vreg
        common /atomc/ stox,wox,wcomp,stcomp,atom,comp
        common /mag/ icfe

        integer, parameter :: npar=43
        integer, parameter :: lmax = 1000
        data ncall/0/
        blank = ' '
        lox = .false.
        coxide = 'oxides'
        write(28,'(a5,12a12)') 'spec','Telastic','thetam3','thetam2','thetam1','thetap0',
     &  'thetap1','thetap2','thetap3','thetap4','thetap5','thetap6'

        open(51,file='control',status='old')
        open(7,file='log',status='unknown')
        open(31,file='qout',status='unknown')
        open(71,file='landau',status='unknown')
c        write(71,'(a5,7a12)') 'Species','Temperature','Gconf','Sconf','Vconf','Qorder','Tc'
        write(31,*) 'Echo input'
        write(31,*) '-----------------Begin control file-----------------'
        do 99 i=1,lmax
         read(51,'(a)',end=98) line
         write(31,'(a80)') line
99      continue
98      continue
        write(31,*) '-----------------End control file-------------------'
        write(31,*) ' '
        rewind (51)

        read(51,*) Pinitial
        ncall = ncall + 1
        if (ncall .eq. 1) call atomset(natom,atom,stox,wox)
        read(51,*) nc
        read(51,'(a)') fname
        if (fname(1:6) .eq. coxide) lox = .true.
        icfe = 0
        do 1 i=1,nc
c        read(51,'(a2,6x,f12.5)') xatom,b(i)
         read(51,*) xatom,b(i)
         do 2 j=1,natom
          if (xatom .eq. atom(j)) then
           comp(i) = atom(j)
	   wcomp(i) = wox(j)
  	   stcomp(i) = stox(j)
           if (j .eq. 26) icfe = i
          end if
2        continue
1       continue
c       print*, (comp(i),i=1,nc)
c       print*, (b(i),i=1,nc)
        call back(51,nc)
        do 12 i=1,nc
c        read(51,'(a2,6x,2f12.5,i5)') xatom,binit(i),bfinal(i),nbulk
         read(51,*) xatom,binit(i),bfinal(i),nbulk
         if (i .eq. 1) nbulk1 = nbulk
         nbulk = nbulk1
         if (nbulk .eq. 0) bfinal(i) = binit(i)
12      continue

        ispec = 0
        iph = 0
        ncs = 0
        do 42 i=1,nphasep
         do 42 j=1,nspecp
          f(i,j) = 0.0
42      continue
        do 41 i=1,nphasep
41      absent(i) = .false.
        read(51,*) ityp,ivtyp,ittyp,iltyp
        print*, ityp,ivtyp,ittyp,iltyp
        write(31,'(/,a32)') ' ************ THEORY ************'
        if (ityp .eq. 1) write(31,*) 'Volume dependence of VDOS:  full theory quadratic in theta, ityp=',ityp
        if (ityp .eq. 2) write(31,*) 'Volume dependence of VDOS:  Full theory linear in theta, ityp=',ityp
        if (ityp .eq. 3) write(31,*) 'Volume dependence of VDOS:  q=constant, ityp=',ityp
        if (ivtyp .eq. 1) write(31,*) 'Volume dependence of G: Full expansion, ivtyp = ',ivtyp
        if (ivtyp .eq. 2) write(31,*) 'Volume dependence of G: Truncated expansion, ivtyp = ',ivtyp
        if (ittyp .eq. 1) write(31,*) 'Temperature dependence of G: Full theory, ittyp = ',ittyp
        if (ittyp .eq. 2) write(31,*) 'Temperature dependence of G: eta prop. V, ittyp = ',ittyp
        if (ittyp .eq. 3) write(31,*) 'Temperature dependence of G: eta constant, ittyp = ',ittyp
        if (ittyp .eq. 4) write(31,*) 'Temperature dependence of G: eta prop. gq, ittyp = ',ittyp
        if (iltyp .eq. 1) write(31,*) 'New Landau theory, post SLB11, iltyp = ',iltyp
        if (iltyp .eq. 2) write(31,*) 'Landau theory as in SLB11, iltyp = ',iltyp
        write(31,'(a32,/)') ' ********************************'
        write(31,*) '-----------------Begin parameter files-----------------'
        read(51,'(a80)') dirname
        ndname = len_trim(dirname)
        print*, dirname,ndname
        read(51,*) nspec
        do 43 iphs=1,nphasep
43      phname(iphs) = blank
        do 4 k=1,nspecp
         if (ispec .eq. nspec) go to 30
         read(51,'(a)') fname
         if (fname(1:5) .eq. 'phase') then
          iph = iph + 1
          call back(51,1)
          read(51,*) fname,phname(iph)
          read(51,*) zeror
          if (zeror .eq. 0.0) absent(iph) = .true.
          go to 4
         end if
         ispec = ispec + 1
         sname(ispec) = fname
         newfname = dirname(1:ndname)//'/'//fname          ! TEST
         open(1,file=newfname,status='old',err=20)
         go to 10
20       print*, 'Filename',fname,newfname,'Not found'
10       continue
         write(31,*) newfname
         do 97 i=1,lmax
          read(1,'(a80)',end=96) line
          write(31,'(a80)') line
97       continue
96       rewind 1
c         call formula(ispec,nc,nsitsp(ispec),lox,comp,s,r)
c	 rewind 1
         call cformula(ispec,nc,nsitsp(ispec),lox,comp,s,r,rep)
         do 31 i=1,nparp
31       apar(ispec,i) = 0.
	 apar(ispec,42) = 2.
         do 3 i=1,npar
          read(1,*,end=3,err=3) apar(ispec,i)
c         read(1,*,end=40,err=40) apar(ispec,i)
c40       write(2,*) i,apar(ispec,i)
3        continue
         f(iph,ispec) = 1.0
         close (1)
         w1 = 0.
         w2 = 1500.
         nw = 1500
         dw = (w2 - w1)/float(nw)
         do 32 iw=1,nw+1
          w = w1 + dw*float(iw-1)
          vdoswr = vdos(ispec,apar,w)
          write(29,*) w/hcok,vdoswr
32       continue
         call moms(ispec,apar,fmom,Telastic)
         if (fmom(1) .ne. 0.) then
          write(31,*) 'VDOS Moments'
          do 33 jm=1,10
           im = jm - 4
           write(31,*) im,fmom(jm)
33        continue
         end if
         write(28,'(a5,12f12.5)') sname(ispec),Telastic,(fmom(jm),jm=1,10)
4       continue
30      nph = iph
c	 write(31,*) 'After do 4',nspec,nc,nsite(lphase(1))
c         do 266 ispec=1,nspec
c          write(31,'(a5,240i3)') sname(ispec),(int(s(j,ispec)),j=1,nc),((int(r(j,ispec,k)),j=1,nc),k=1,nsite(lphase(ispec)))
c          write(*,'(a5,240i3)') sname(ispec),(int(s(j,ispec)),j=1,nc),((int(r(j,ispec,k)),j=1,nc),k=1,1),lphase(ispec)
c266      continue
c	print*, 'after 266'
        close (28)
        close (29)
        write(31,*) '-----------------End parameter files-----------------'
        write(31,*) 'End input echo'
        write(31,*) ' '
        if (lox) write(31,*) 'Oxides'
        write(31,*) 'Bulk Composition'
        do 95 i=1,nc
         write(31,*) comp(i),binit(i),bfinal(i),nbulk
95      continue
        if (icfe .ne. 0) 
     &   write(31,*) 'Magnetic atom present in component: ',comp(icfe)

        do 21 iph=1,nph
         if (.not. absent(iph)) go to 21
         do 22 i=1,nspec
          if (f(iph,i) .eq. 0) go to 22
          absents(i) = .true.
22       continue
21      continue
        do 23 lph=1,nph
         do 23 ispec=1,nspec
          if (f(lph,ispec) .eq. 0.0) go to 23
          lphase(ispec) = lph
          mphase(lph) = mphase(lph) + 1
          if (mphase(lph) .eq. 1 .and. iphase(lph) .eq. 0)
     &        iphase(lph) = ispec
23      continue

	do 38 lph=1,nph
	 iophase(lph) = 0
	 mophase(lph) = 0
38	continue
	do 34 lph=1,nph
	 do 35 ispec=iphase(lph),iphase(lph)+mphase(lph)-2
	  do 36 jspec=ispec+1,iphase(lph)+mphase(lph)-1
	   do 37 jco=1,nc
	    if (s(jco,ispec) .ne. s(jco,jspec)) go to 36
37	   continue
	   write(31,*) 'Found two species of the same phase with identical stoichiometries',ispec,jspec
	   write(31,*) 'Assume species of identical stoichiometry appear in sequence in control file'
	   go to 39
	   write(31,*) 'Assume species of identical stoichiometry do NOT model order-disorder'
	   go to 36
39	   write(31,*) 'Assume species of identical stoichiometry model order-disorder'
	   if (mophase(lph) .eq. 0) then
	    mophase(lph) = 1
	    iophase(lph) = ispec
	   end if
	   mophase(lph) = mophase(lph) + 1
36	  continue
35	 continue
34	continue

	do 44 lph=1,nph
	 if (mophase(lph) .gt. 0) then
	  write(31,'(a23,2x,71a6)') 'Order-disorder species:',(sname(ispec),ispec=iophase(lph),iophase(lph)+mophase(lph)-1)
	 end if
44	continue

	do 59 ispec=1,nspecp
	 do 59 jspec=1,nspecp
	  do 59 k=1,nsitep
	   iastate(ispec,jspec,k) = .false.
59	continue
	if (icfe .eq. 0) go to 541
	do 54 lph=1,nph
	 do 55 ispec=iphase(lph),iphase(lph)+mphase(lph)-2
	  do 56 jspec=ispec+1,iphase(lph)+mphase(lph)-1
	   do 57 k=1,nsitsp(lph)
	    if (r(icfe,ispec,k) .eq. 0.) go to 57
C  Ferric iron on the first site of hepv and hlpv are both high spin (the large, dodecahedral, 'Mg' site on which Fe is high spin throughout the mantle in FeSiO3 and FeFeO3).  Following assumes that hlpv follows hepv in the control file
C  Ferric iron on the first site of mag and on the third site of mag.  Assume no disorder on the first (tetrahedral) ferric site.  Following assumes that mag follows wu in the control file
	    if (sname(ispec) .eq. 'hepv' .and. sname(jspec) .eq. 'hlpv' .and. k .eq. 1) go to 57
	    if (sname(ispec) .eq. 'hppv' .and. sname(jspec) .eq. 'lppv' .and. k .eq. 1) go to 57
	    if (sname(ispec) .eq. 'hepv' .and. sname(jspec) .eq. 'fapv' .and. k .eq. 1) go to 57
	    if (sname(ispec) .eq. 'hlpv' .and. sname(jspec) .eq. 'fapv' .and. k .eq. 1) go to 57
c	    if (sname(ispec) .eq. 'wu  ' .and. sname(jspec) .eq. 'mag ' .and. k .eq. 2) go to 57
c	    if (r(icfe,ispec,k) .eq. r(icfe,jspec,k)) then
	    if (r(icfe,ispec,k) .ne. 0. .and. r(icfe,jspec,k) .ne. 0.) then
	     iastate(ispec,jspec,k) = .true.
	     iastate(jspec,ispec,k) = .true.
	     write(31,*) 'Found iron in different states on the same site of the same phase in two different species',ispec,jspec,k,
     &        sname(ispec)(1:4),sname(jspec)(1:4),iastate(ispec,jspec,k),iastate(jspec,ispec,k)
	    end if
57	   continue
56	  continue
55	 continue
54	continue
541	continue

C  Following assumes that the first site is completely ordered in the end-member and that it does not mix with any other 
C  components of any other species on that site.
	do 64 lph=1,nph
	 do 65 ispec=iphase(lph),iphase(lph)+mphase(lph)-1
	  if (phname(lph)(1:2) .eq. 'mw' .and. sname(ispec)(1:3) .eq. 'mag') then
	   write(31,*) 'Found completely ordered iron site.  set r_ijk=0',lph,icfe,ispec,ione,phname(lph)(1:4),sname(ispec)(1:4)
	   r(icfe,ispec,1) = 0.
	  end if
c	  if (phname(lph)(1:2) .eq. 'pv' .and. sname(ispec)(1:4) .eq. 'hepv') then
c	   write(31,*) 'Found completely ordered iron site.  set r_ijk=0',lph,icfe,ispec,ione,phname(lph)(1:4),sname(ispec)(1:4)
c	   r(icfe,ispec,1) = 0.
c	  end if
c	  if (phname(lph)(1:2) .eq. 'pv' .and. sname(ispec)(1:4) .eq. 'hlpv') then
c	   write(31,*) 'Found completely ordered iron site.  set r_ijk=0',lph,icfe,ispec,ione,phname(lph)(1:4),sname(ispec)(1:4)
c	   r(icfe,ispec,1) = 0.
c	  end if
65	 continue
64	continue
	    
C  Assumes that the species with the extra site is not the first species of the phase
        do 28 lph=1,nph
         nmin = 1000
         do 29 ispec=iphase(lph),iphase(lph)+mphase(lph)-1
          nmin = min(nmin,nsitsp(ispec))
          if (nsitsp(ispec) .ne. nmin .and. ispec .gt. iphase(lph)) then
           write(31,*) 'WARNING: not all species have the same number of sites in phase',lph,ispec,nmin,nsitsp(ispec)
           nsitsp(ispec) = 1
          else
           nsitsp(ispec) = 0
          end if
29       continue
         nsite(lph) = nmin
28      continue
	print*, 'number of sites in phase 3 (opx)',nsite(3)
         
        nco = nc
        call sform(s,b,n1,q1,q2,nspec,nco,nc,ncs,nnull,nnulls,absents)
        call sitered(r,wreg,f,nsite,iphase,mophase,iastate,nph,nspec,nco)
c        call regread(wreg,comp,nsite,nco,nph,dirname,ndname)
        call regread(wreg,vreg,sname,nspec,mphase,nph,dirname,ndname)
	print*, 'after regread',wreg(1,1,1,2),wreg(1,1,2,1),vreg(1,1,1,2),vreg(1,1,2,1)

        print*, nph,' Phases'
        write(31,*) nph,' Phases'
        write(31,'(27a5)') (phname(iph),iph=1,nph)
        write(31,*) nspec,' Species'
        write(31,'(27a5)') (sname(ispec),ispec=1,nspec)
        write(31,*) 'Species Arrangement'
        do 24 iph=1,nph
         write(31,'(i2,1x,27a5)') iph,phname(iph),
     &    (sname(i),i=iphase(iph),iphase(iph)+mphase(iph)-1)
24      continue
        write(31,*) 'Species composition'
        write(31,'(7x,9(a6,12x))') 'Total ','site 1','site 2','site 3','site 4','site 5'
        write(31,'(5x,240a3)') (comp(j),j=1,nco),(comp(j),j=1,nco),
     &                         (comp(j),j=1,nco),(comp(j),j=1,nco),
     &                         (comp(j),j=1,nco),(comp(j),j=1,nco)
        do 26 ispec=1,nspec
         write(31,'(a5,240i3)') sname(ispec),(int(s(j,ispec)),j=1,nco),((int(r(j,ispec,k)),j=1,nco),k=1,nsite(lphase(ispec)))
26      continue
        write(31,*) 'Species parameters'
        write(31,*) 'Debye model, Birch-Murnaghan, no electronic terms, no zero-point'
        write(31,'(3x,a7,a6,18a9)') 'species','fn','Z','wm','T0','F0','V0',
     &     'K0','K0p','K0pp','theta0','gam0','q0','G0','G0p','eta0','Vbw1','Vbw2','Vsp1','Vsp2'
        do 25 ispec=1,nspec
         write(31,'(i2,1x,a5,19f9.2)') ispec,sname(ispec),
     &                        (apar(ispec,ipar),ipar=1,10),
     &                        (apar(ispec,ipar),ipar=26,27),
     &                        (apar(ispec,ipar),ipar=35,37),(apar(ispec,ipar),ipar=51,54)
25      continue
	print*, 'before do 27',wreg(1,1,1,2),wreg(1,1,2,1)
        write(31,*) 'Non-zero regular solution parameters'
        do 27 iph=1,nph
          do 27 k=1,nsite(iph)
           do 27 ia=1,nspec-1
            do 27 ib=ia+1,nspec
             if (wreg(iph,k,ia,ib) .ne. 0.) 
     &        write(31,'(a5,4i5,1x,2a5,3f14.5)') phname(iph),iph,k,ia,ib,
     &         sname(ia),sname(ib),wreg(iph,k,ia,ib),wreg(iph,k,ib,ia),vreg(iph,k,ia,ib)
27      continue

        write(31,*) 'Initial number of total constraints, compositional constraints, species constraints, total number of species'
        write(31,*) 'nc,nco,ncs,nspec = ',nc,nco,ncs,nspec
	write(31,*) 'ssmall,tradd = ',ssmall,tradd

        return
        end
