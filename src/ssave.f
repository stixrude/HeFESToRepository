        subroutine ssave(fret,qual,nnew,fretsav,qualsav,absentsav,absentssav,nnewsav,ftol,succes,ires)

        include 'P1'
	include 'numpar.inc'
        include 'chem.inc'
        include 'absent.inc'
        include 'const.inc'

	integer i,ires
	double precision fret,qual,fretsav,qualsav,ftol,vsum
        double precision nnew(nspecp)
        double precision nnewsav(nspecp)
        logical absentsav(nspecp),absentssav(nspecp),succes,valid
	succes = .false.

	if (.not. valid(vsum,nnew)) then
	 if (vsum .gt. ssmall) then
	  write(31,*) 'ssave returning without considering invalid solution',fret,fretsav,qual,qualsav,nnew(nnull+1),ftol,vsum
	  return
	 end if
	end if
c	if (ires .lt. 0) then
c	 write(31,*) 'ssave returning without considering invalid solution',fret,fretsav,qual,qualsav,nnew(nnull+1),ftol,ires
c	 return
c	end if
        call qcalc(nnew,qual)

	write(31,*) 'ssave fret,fretsav,qual,qualsav',fret,fretsav,qual,qualsav
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
C  Energy is worse but acceptable and quality is better
	if (fret .lt. fretsav+ftol .and. qual .lt. qualsav) then
	 write(31,*) 'ssave Energy is worse but acceptable and quality is better'
	 go to 30
	end if
C  Energy is much better regardless of quality
	if (fret .lt. fretsav-ftol) then
	 write(31,*) 'ssave Energy is much better regardless of quality'
	 go to 30
	end if
	write(31,*) 'ssave returning without saving',fret,fretsav,qual,qualsav,nnew(nnull+1),ftol
	return
30	continue
	succes = .true.
	write(31,*) 'ssave saving new solution: fret,fretsav,ftol',fret,fretsav,ftol,nnew(nnull+1)
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
