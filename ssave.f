        subroutine ssave(fret,qual,nnew,fretsav,qualsav,absentsav,absentssav,nnewsav,ftol,succes,ires)

        include 'P1'
	include 'numpar.inc'
        include 'chem.inc'
        include 'absent.inc'
        include 'const.inc'

	integer i,ires,ispec
	double precision, parameter :: qreltest = 1.e-1
	double precision fret,qual,fretsav,qualsav,ftol,vsum,fnagg,Ti,Pi,apar
        double precision nnew(nspecp)
        double precision nnewsav(nspecp)
        logical absentsav(nspecp),absentssav(nspecp),succes,valid
        common /state/ apar(nspecp,nparp),Ti,Pi
	succes = .false.

c	call nform(nnew,n,n1,q2,nspec,nnull)
        fnagg = 0.
        do 12 ispec=1,nspec
         fnagg = fnagg + apar(ispec,1)*n(ispec)
12      continue

	if (.not. valid(vsum,nnew)) then
	 if (vsum .gt. ssmall) then
	  write(31,*) 'ssave returning without considering invalid solution',fret,fretsav,qual,qualsav,nnew(nnull+1),ftol,vsum
	  return
	 end if
	end if
        call qcalc(nnew,qual)

	write(31,'(a42,99e16.9)') 'ssave fret,fretsav,qual,qualsav,ftol,fnagg',fret,fretsav,qual,qualsav,ftol,fnagg
c     &   ftol*fnagg,fret-fretsav,abs(qual-qualsav)/qualsav
C  Only from gibmin to force trace removal.
	if (ftol .eq. 10.) then
	 write(31,*) 'ssave Accept to remove trace species'
	 go to 30
	end if
C  Energy is better and quality is better!
	if (fret .lt. fretsav .and. qual .lt. qualsav) then
	 write(31,*) 'ssave Energy is better and quality is better!'
	 go to 30
	end if
C  Energy is better and quality is worse but acceptable
        if (fret .lt. fretsav .and. qual .lt. qtest) then
	 write(31,*) 'ssave Energy is better and quality is worse but acceptable'
	 go to 30
	end if
cC  Energy is better and quality is not acceptable but only slightly worse
c        if (fret .lt. fretsav .and. qual .gt. qtest .and. abs(qual-qualsav)/qualsav .lt. qreltest) then
c	 write(31,*) 'ssave Energy is better and quality is not acceptable but only slightly worse'
c	 go to 30
c	end if
C  Energy is worse but acceptable and quality is better
	if (fret .lt. fretsav+fnagg*ftol .and. qual .lt. qualsav) then
	 write(31,*) 'ssave Energy is worse but acceptable and quality is better'
	 go to 30
	end if
cC  Energy is worse but acceptable and quality not acceptable but only slightly worse
c	if (fret .lt. fretsav+fnagg*ftol .and. abs(qual-qualsav)/qualsav .lt. qreltest) then
c	 write(31,*) 'ssave Energy is worse but acceptable and quality is only slightly worse'
c	 go to 30
c	end if
C  Energy is much better regardless of quality
	if (fret .lt. fretsav-fnagg*ftol) then
	 write(31,*) 'ssave Energy is much better regardless of quality'
	 go to 30
	end if
	write(31,*) 'ssave returning without saving'
	return
30	continue
	succes = .true.
	write(31,*) 'ssave saving new solution'
        fretsav = fret
        qualsav = qual
        do 32 i=1,nspec
         absentssav(i) = absents(i)
         nnewsav(i) = nnew(i)
32      continue
        do 36 i=1,nph
36      absentsav(i) = absent(i)

        return
        end
