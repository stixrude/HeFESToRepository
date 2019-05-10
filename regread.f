C  Read in regular solution parameters from phase files in parameter sub-directory phase
C  5 Sept. 2018.  Added ability to read pressure dependence of regular solution parameters from the same file
C  This addition is backwards compatible: if the volume block is not present in the phase file, the routine 
C  assumes that the volume of solution is ideal.
         subroutine regread(wreg,vreg,sname,nspec,mphase,nph,dirname,ndname)

        include 'P1'
	integer nspec,nph,ndname,i,ia,ib,ic,id,ii,j,k,l,lph,ncmax,nfirst,nlast,nspeci
        character*1 blank,label
        character*4 xatom(natomp)
        character*80 dirname,fname,newfname,subs(100)
        character*132 header
        character*80 sname(nspecp)
        integer mphase(nphasep),nchar(100)
	integer ip(nspecp)
        double precision wreg(nphasep,nsitep,nspecp,nspecp),vreg(nphasep,nsitep,nspecp,nspecp)
        double precision temp(nspecp,nspecp),vtemp(nspecp,nspecp)
        blank = ' '
        dirname = dirname(1:ndname)//'/'//'PHASE/'
        ndname = ndname + 7
c        ncmax = 6

        do 1 i=1,nphasep
         do 1 j=1,nsitep
          do 1 k=1,nspecp
           do 1 l=1,nspecp
            wreg(i,j,k,l) = 0.0
            vreg(i,j,k,l) = 0.0
            temp(k,l) = 0.0
	    vtemp(k,l) = 0.0
1       continue

        lph = 0
	nspeci = 1
        rewind (51)
2	continue
         do 11 k=1,nspec
          do 11 l=1,nspec
           temp(k,l) = 0.0
	   vtemp(k,l) = 0.0
11       continue
         if (lph .eq. nph) go to 50
C  Search for phase
         read(51,'(a)',end=50) fname
         if (fname(1:5) .ne. 'phase') go to 2
         lph = lph + 1
         nfirst = 7
         nlast = 7
         if (fname(nfirst:nlast) .eq. 'blank') print*, 'Blank phase name'
         do 3 j=8,50
          if (fname(j:j) .eq. blank) go to 4
          nlast = nlast + 1
3        continue
4        continue
         newfname = fname(nfirst:nlast)
         newfname = dirname(1:ndname)//newfname
         open(108+lph,file=newfname,status='old',err=20)
         go to 30
20       print*, 'Regular solution parameter file not found.  Assuming ideal mixing.',newfname
30       continue
C  Found new phase
	 if (mphase(lph) .le. 1) go to 2
	 read(108+lph,'(a132)') header
	 call parse(header,subs,nchar,ncmax,132)
	 backspace (108+lph)
         read(108+lph,*) (xatom(ii),ii=1,ncmax)
	 print*, header,ncmax,(xatom(ii),ii=1,ncmax)
         do 12 ia=1,ncmax
          read(108+lph,*) (temp(ia,ib),ib=1,ncmax)
	  print*, (temp(ia,ib),ib=1,ncmax)
12       continue
         read(108+lph,'(a)',err=1210,end=1210) label
         do 121 ia=1,ncmax
          read(108+lph,*) (vtemp(ia,ib),ib=1,ncmax)
	  print*, (vtemp(ia,ib),ib=1,ncmax)
121      continue
	 go to 1220
1210	 print*, 'Volume of solution is ideal'
1220	 continue
	 do 21 ic=1,ncmax
	  ip(ic) = nspecp
c	  nspeci = 1
c	  do 22 id=1,nspec
c	   if (xatom(ic) .eq. sname(id)(1:4)) ip(ic) = id
c	   print*, ic,id,ip(ic),xatom(ic),sname(id)
c22	  continue
	  do 22 id=nspeci,nspec
	   if (xatom(ic) .eq. sname(id)(1:4)) then
	    ip(ic) = id
	    nspeci = id+1
	    go to 21
	   end if
	   print*, 'Test',ic,id,ip(ic),xatom(ic),sname(id)
22	  continue
21	 continue
	 do 23 ia=1,ncmax
	  do 24 ib=1,ncmax
	   if (ia .eq. ib) go to 24
	   if (ip(ib) .eq. nspecp) go to 24
	   if (ip(ia) .eq. nspecp) go to 23
	   wreg(lph,1,ip(ia),ip(ib)) = temp(ia,ib)
	   vreg(lph,1,ip(ia),ip(ib)) = vtemp(ia,ib)
	   print '(a17,5i5,1x,4a5,6f12.5)', 'Found w parameter',lph,ia,ib,ip(ia),ip(ib),sname(ip(ia)),sname(ip(ib))
     &      ,xatom(ia),xatom(ib),wreg(lph,1,ip(ia),ip(ib)),vreg(lph,1,ip(ia),ip(ib))
24	  continue
23	 continue

C  For backwards compatibility 19 Sept. 2015: if W_i<j=W_j>i set O_ij=W_j>i=0.
	do 231 ia=1,ncmax-1
	 do 231 ib=ia+1,ncmax
	  if (ip(ib) .eq. nspecp) go to 231
	  if (ip(ia) .eq. nspecp) go to 231
	  if (wreg(lph,1,ip(ia),ip(ib)) .eq. wreg(lph,1,ip(ib),ip(ia)) .and. wreg(lph,1,ip(ia),ip(ib)) .ne. 0.) then
	   wreg(lph,1,ip(ib),ip(ia)) = 0.
	   vreg(lph,1,ip(ib),ip(ia)) = 0.
	   print '(a17,5i5,4a5,6f12.5)', 'Set o to zero',lph,ia,ib,ip(ia),ip(ib),sname(ip(ia)),sname(ip(ib))
     &      ,xatom(ia),xatom(ib),wreg(lph,1,ip(ia),ip(ib)),wreg(lph,1,ip(ib),ip(ia))
	   print '(a17,5i5,4a5,6f12.5)', 'Set ov to zero',lph,ia,ib,ip(ia),ip(ib),sname(ip(ia)),sname(ip(ib))
     &      ,xatom(ia),xatom(ib),vreg(lph,1,ip(ia),ip(ib)),vreg(lph,1,ip(ib),ip(ia))
	  else
	   print '(a17,5i5,1x,4a5,6f12.5)', 'W and O',lph,ia,ib,ip(ia),ip(ib),sname(ip(ia)),sname(ip(ib))
     &      ,xatom(ia),xatom(ib),wreg(lph,1,ip(ia),ip(ib)),wreg(lph,1,ip(ib),ip(ia)),
     &                           vreg(lph,1,ip(ia),ip(ib)),vreg(lph,1,ip(ib),ip(ia))
	  end if
231	continue

c         do 21 ia=1,nspec
c          do 22 ic=1,ncmax
c           if (sname(ia) .eq. xatom(ic)) go to 40
c22        continue
c40        continue
c          do 21 ib=1,nspec
c           do 23 id=1,ncmax
c            if (sname(ib) .eq. xatom(id)) go to 50
c23         continue
c50         continue
c           wreg(lph,1,ia,ib) = temp(ic,id)
c	   print '(a17,4i5,4a5,f12.5)', 'Found w parameter',ia,ib,ic,id,sname(ia),sname(ib),xatom(ic),xatom(id),wreg(lph,1,ia,ib)
c21       continue
         read(108+lph,'(a)',err=10,end=10) fname
10       continue
         close (108+lph)
	 go to 2
50	continue
        close (108+lph)

        return
        end
