        subroutine formula(ispec,nc,nsite,lox,comp,s,r)

C  Reads in the Chemical Formula of species ispec
C  Identifies the Components Present in the Formula
C  Assigns Stoichiometric and Site Coefficients: s,r

        include 'P1'

	integer ispec,nc,nsite,i,im1,ip1,j,j1,j2,jatom,k,ksite,m,mult,ncaug,nchar
        logical begin,lox,close,assign,ladd
        character*1 number(10),c,blank,underscore,right,left,comma
        character*2 comp(natomp),xatom,yatom1,yatom2,ox
        character*80 form,temp,blank80
        double precision r(ncompp,nspecp,nsitep),s(ncompp,nspecp)
        double precision coeff(ncompp)
        integer jatoma(ncompp)
        data number/'0','1','2','3','4','5','6','7','8','9'/
c  Assumes that bulk composition is given as oxides
	if (.not. lox) return
        ox = 'O '
        blank = ' '
        blank80 = ' '
        underscore='_'
        left = '('
        right = ')'
        comma = ','
        begin = .false.
        ladd = .false.
        
C  Read Formula and Count Number of Characters
        nchar = 0
        temp = blank80
        form = blank80
        read(1,'(a80)') temp
        do 1 i=1,30
         c = temp(i:i)
         if (c .ne. blank) then
          begin=.true.
          nchar = nchar + 1
          form(nchar:nchar) = c
         end if
         if (begin .and. c .eq. blank) go to 10
1       continue
10      continue

        ksite = 0
        jatom = 0
        ncaug = 0
        mult = 0
        close = .true.
        assign = .false.
CMAC    do 12 j=1,ncomp
CMAC
        do 12 j=1,ncompp
CMAC
         jatoma(j) = 0
         coeff(j) = 0.
12      continue

C  Loop Through Characters in the Formula
        do 2 i=1,nchar
         im1 = i - 1
         ip1 = i + 1
         yatom1 = form(i:ip1)
         yatom2 = form(i:i)//blank
         if (form(i:i) .eq. left) close = .false.
         if (form(i:i) .eq. right) close = .true.

C  Identify Chemical Components Present in the Formula
C  Remember that a stoich. coeff. must be assigned to the atom 
C  (ladd=.true.)
         do 3 j=1,nc
          xatom = comp(j)
          if (yatom1 .eq. xatom) then
           mult = mult + 1
           jatoma(mult) = j
           ladd = .true.
           go to 2
          else if (yatom2 .eq. xatom) then
           mult = mult + 1
           jatoma(mult) = j
           ladd = .true.
           go to 2
          end if
3        continue

C  Unidentifiable Chemical Component in Formula
C  Augment Number of Components, nc
c        do 9 j=1,natomp
c         xatom = atom(j)
c         if (yatom1 .eq. xatom) then
c          if (lox .and. xatom .eq. ox) go to 2
c          ncaug = ncaug + 1
c          jatom = nc + ncaug
c          go to 2
c         else if (yatom2 .eq. xatom) then
c          if (lox .and. xatom .eq. ox) go to 2
c          ncaug = ncaug + 1
c          jatom = nc + ncaug
c          go to 2
c         end if
c9       continue

C  Identify Stoichiometric Coefficients
C  Remember that stoich. coeff(s). must be assigned to all atoms on this site
C  (assign=.true.)
         if (ladd) then
          do 4 j1=1,9
           if (form(i:i) .eq. number(j1) .and. 
     &         form(i-1:i-1) .eq. underscore) then
            coeff(mult) = j1 - 1
            do 5 j2=1,10
             if (form(i+1:i+1) .eq. number(j2)) then
              coeff(mult) = coeff(mult) + 9.*(j1-1) + (j2-1)
             end if
5           continue
            assign = .true.
            ladd = .false.
           end if
4         continue
         end if

C  Assign Stoichiometric Coefficients to s,r matrices
         if (jatoma(1) .ne. 0 .and. assign .and. close) then
          ksite = ksite + 1
          do 41 m=1,mult
           jatom = jatoma(m)
           r(jatom,ispec,ksite) = coeff(m)
           jatoma(m) = 0
           coeff(m) = 0
41        continue
          jatom = 0
          mult = 0
          assign = .false.
         end if
2       continue
        nsite = ksite
c       nc = nc + ncaug

        do 6 j=1,nc
         s(j,ispec) = 0.0
         do 6 k=1,nsite
          s(j,ispec) = s(j,ispec) + r(j,ispec,k)
6       continue

c       print*, 'Stoichiometric Coefficients for species',ispec
c       print*, (s(j,ispec),j=1,nc)

c       print*, 'Site Coefficients for species',ispec
c       do 71 j=1,nc
c71     print*, (r(j,ispec,k),k=1,nsite)

        return
        end
