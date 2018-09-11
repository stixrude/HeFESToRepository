        subroutine tform(n,iphase,mphase,itrmin,itrsmin,nspec,nph,absent,absents)

        include 'P1'
	include 'const.inc'
	include 'numpar.inc'

	integer itrmin,itrsmin,nspec,nph,i,iph,ispec
	double precision dnrm2
        logical absents(nspecp),absent(nphasep)
        integer iphase(nphasep),mphase(nphasep)
        double precision n(nspecp)
        double precision ntrmin,nsum(nphasep)

C  Identify trace species as those that are: 1) not absent and 2) have abs(n_i) < ssmall
C  Note that dnrm2 is a LAPACK routine that returns the euclidean norm of a vector

        itrsmin = 0
        ntrmin = +1.e15
        do 8 i=1,nspec
         if (absents(i)) go to 8
         if (abs(n(i)) .lt. ssmall) then
          ntrmin = min(ntrmin,abs(n(i)))
          if (ntrmin .eq. abs(n(i))) itrsmin = i
c          write(31,*) 'Identified trace species',i,n(i)
         end if
8       continue
        if (itrsmin .ne. 0) write(31,*) 'Smallest trace species',itrsmin,n(itrsmin)

C  Identify trace phases as those that are: 1) not absent and 2) have only trace or absent species

        itrmin = 0
        ntrmin = +1.e15
        do 10 iph=1,nph
         if (absent(iph)) go to 10
         do 11 ispec=iphase(iph),iphase(iph)+mphase(iph)-1
          if (.not. absents(ispec)) go to 10
11       continue
	 nsum(iph) = dnrm2(mphase(iph),n(iphase(iph)),ione)/float(mphase(iph))
         ntrmin = min(ntrmin,nsum(iph))
         if (ntrmin .eq. nsum(iph)) itrmin = iph
c         write(31,*) 'Identified trace phase',iph,(n(i),i=iphase(iph),iphase(iph)+mphase(iph)-1)
10      continue
c        if (itrmin .ne. 0) write(31,*) 'Smallest trace phase',itrmin,nsum(itrmin)

        return
        end
