        subroutine cformula(ispec,nc,nsite,lox,comp,s,r,rep)

C  Reads in the Chemical Formula of species ispec
C  Identifies the Components Present in the Formula
C  Assigns Stoichiometric and Site Coefficients: s,r

        include 'P1'
	include 'const.inc'
	include 'elem.inc'

	integer ispec,nc,nsite,i,i1,i2,ic,ielem,iloc,ip1,isite,jatom,jelem,ksite,mult,n,ncall,ncaug,nchar,j,k
	integer idchar,nmchar
	double precision fn,stoich,wm
        logical lox,open,square
        character*1 blank
        character*2 comp(natomp),yatom
        character*80 form
	character*132 temp
	character*80 subs(100)
	double precision r(ncompp,nspecp,nsitep),s(ncompp,nspecp),rep(ncompp,nspecp,nsitep)
c	double precision rat(ncompp,nspecp,nsitep),sat(ncompp,nspecp)
	integer nchara(100)
        data ncall/0/
        ncall = ncall + 1

C  Initialize
	if (ncall .eq. 1) then
	 do 99 i=1,ncompp
	  do 99 j=1,nspecp
	   sat(i,j) = 0.
	   do 99 k=1,nsitep
	    rat(i,j,k) = 0.
	    rep(i,j,k) = 0.
99	 continue
	end if

C  Point elements to components
        if (ncall .eq. 1) then
	 do 6 ic=1,nc
          do 6 jelem=1,nelem
          if (comp(ic)(1:2) .eq. elnam(jelem)(1:2)) ipelem(ic) = jelem
6        continue
	end if

c  Assumes that bulk composition is given as oxides
	if (.not. lox) return
        blank = ' '
        
C  Read Formula and Count Number of Characters
        read(1,'(a80)') temp
	call parse(temp,subs,nchara,n,80)
	form = subs(1)
	nchar = nchara(1)
	print*, nchar,form(1:nchar)

        ksite = 0
        jatom = 0
        ncaug = 0
        mult = 0
	open = .false.
	square = .false.

C  Loop Through Characters in the Formula
	fn = 0.
	wm = 0.
	iloc = 1
	ksite = 0
10	continue
	i = iloc
        ip1 = i + 1
	if (idchar(form(i:i)) .eq. 5 .or. idchar(form(i:i)) .eq. 8) then
	 open = .true.
	 if (idchar(form(i:i)) .eq. 8) square = .true.
	 ksite = ksite + 1
	end if
	if (idchar(form(i:i)) .eq. 6 .or. idchar(form(i:i)) .eq. 9) open = .false.
C  Atom?
        if (idchar(form(i:i)) .eq. 3) then
         if (idchar(form(ip1:ip1)) .eq. 2) then
          yatom = form(i:ip1)
	  iloc = iloc + 3
         else
          yatom = form(i:i)//blank
	  iloc = iloc + 2
         end if
	 if (.not. open) ksite = ksite + 1
         do 3 jelem=1,nelem
          if (yatom(1:2) .eq. elnam(jelem)(1:2)) then
	   ielem = jelem
           print*, 'Element found',ielem,elnam(ielem),wat(ielem),ksite
          end if
3        continue
C  Find Stoichiometric Coefficient
	i1 = iloc
	ip1 = iloc + 1
        if (idchar(form(ip1:ip1)) .eq. 1) then
	 i2 = iloc + 1
	else
	 i2 = iloc
	end if
	stoich = nmchar(form,i1,i2)
	print*, 'stoich = ',stoich
	fn = fn + stoich
	wm = wm + stoich*wat(ielem)
C  Allow for repeated occurence of an element on a site
	if (rat(ielem,ispec,ksite) .ne. 0.) rep(ielem,ispec,ksite) = stoich/(stoich+rat(ielem,ispec,ksite))
	rat(ielem,ispec,ksite) = rat(ielem,ispec,ksite) + stoich
c	rat(ielem,ispec,ksite) = stoich
	if (square) then
	 rat(ielem,ispec,ksite) = -stoich
	 if (.not. open) square = .false.
	end if
	sat(ielem,ispec) = sat(ielem,ispec) + stoich
	else
	 iloc = iloc + 1
	 if (iloc .le. nchar) go to 10
        end if
	if (iloc .le. nchar) go to 10
	print '(a7,2f12.5)', 'fn,wm= ',fn,wm
	print '(100i3)', (nint(sat(jelem,ispec)),jelem=1,36)
	do 5 isite=1,ksite
	 print '(100i3)', (nint(rat(jelem,ispec,isite)),jelem=1,36)
5	continue

C  Assign r and s
	nsite = ksite
	do 7 ic=1,nc
	 s(ic,ispec) = sat(ipelem(ic),ispec)
	 do 8 isite=1,nsite
	  r(ic,ispec,isite) = rat(ipelem(ic),ispec,isite)
8	 continue
7	continue
	 
	return
	end
