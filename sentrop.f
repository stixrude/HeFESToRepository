	subroutine sentrop(Ttry,entrop,dsdtmol,wmagg)

	include 'P1'
	include 'chem.inc'
	include 'absent.inc'
	integer ispec
	double precision sspeco(nspecp),chempot,rsum,smixi,smag,volsum
	double precision apar,Ti,Pi,Ttry,entrop,gspec,dsdtmol
	double precision volve,entve,cpve,bkve,wmagg,wm
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /volent/ volve,entve,cpve,bkve
	Ti = Ttry

	entrop = 0.
	dsdtmol = 0.
	wmagg = 0.
	do 1 ispec=1,nspec
	 if (absents(ispec)) go to 1
	 wm = apar(ispec,3)
         gspeca(ispec) = gspec(ispec)
	 sspeco(ispec) = entve
         cspeca(ispec) = cpve
         call cp(ispec,n,chempot,rsum,volsum,smixi,smag)
         sspeca(ispec) = sspeco(ispec) + smixi + smag
	 entrop = entrop + n(ispec)*sspeca(ispec)
	 dsdtmol = dsdtmol + n(ispec)*cspeca(ispec)
	 wmagg = wmagg + n(ispec)*wm
c	 write(31,*) 'sentrop',Pi,Ti,ispec,n(ispec),sspeca(ispec),entrop
1	continue
	dsdtmol = dsdtmol/Ttry
	 
	return
	end
