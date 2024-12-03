        subroutine lagscomp(s,n,lags,lagc,gspeca,lphase,nspec,nco,absent,absents)

C  Compute the chemical affinities of absent species in present phases lags
C  lags_j = mu_j - s_ij lagc_i 
C  See Harvie, Greenberg, Weare, GCA, 51, 1045, eq. 28.  

        include 'P1'
	include 'numpar.inc'

	integer nspec,nco,i,ispec
	double precision apar,chempot,cpa,fn,Pi,rsum,smag,smixi,Ti,volsum
	logical absents(nspecp),absent(nphasep)
        logical spinod(nspecp),spinph(nphasep)
	integer lphase(nspecp)
	double precision gspeca(nspecp),cpcomp(nspecp)
	double precision s(ncompp,nspecp),n(nspecp),lags(nspecp)
	double precision nlags(nspecp)
	double precision lagc(ncompp)
        character*80 phname(nphasep),sname(nspecp)
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /names/ phname,sname
        common /spinc/ spinod,spinph

	call lagcomp(s,cpcomp,lagc,nspec,nco,absents)

        do 5 i=1,nspecp
         lags(i) = 0.0
	 nlags(i) = n(i)
5	continue

c	write(31,*) 'Check for species addition or subtraction'
c        write(31,'(a7,2x,3a12)') 'species','affinity'
	do 1 ispec=1,nspec
	 fn = apar(ispec,1)
	 if (spinod(ispec)) go to 1
	 if (absent(lphase(ispec))) go to 1
	 if (.not. absents(ispec)) go to 1
	 nlags(ispec) = tradd
	 call cp(ispec,nlags,chempot,rsum,volsum,smixi,smag)
	 nlags(ispec) = n(ispec)
	 lags(ispec) = gspeca(ispec)/1000. + chempot/1000.
	 cpa = lags(ispec)
	 lags(ispec) = lags(ispec) - cpcomp(ispec)
	 lags(ispec) = lags(ispec)/fn
c         write(31,'(i3,1x,a5,10f12.4)') ispec,sname(ispec),lags(ispec)
1	continue

        return
        end
