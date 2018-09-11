        subroutine setup(binit,dbulk,P1,dP,nP,T1,dT,nT,superad,tfreeze,nbulk,adiabat,chcalc,adcalc,hucalc,frozen,nvet,nvep)

        include 'P1'
        include 'chem.inc'
	include 'absent.inc'
        include 'lag.inc'

	integer np,nt,nbulk,i,ispec,iensem,nvet,nvep,iextra
	double precision P1,dp,T1,dt,superad,tfreeze,P2,Pi,T2,Ti
        logical adiabat,chcalc,adcalc,hucalc,frozen
        logical isochor
        double precision nnew(nspecp)
        double precision binit(ncompp),bfinal(ncompp),dbulk(ncompp),zero
        common /chor/ ispec,isochor
        adiabat = .false.
        isochor = .false.
	adcalc = .false.
	chcalc = .false.
	hucalc = .false.

C  Set universal constants including identity matrix
        call const

	superad = 0.0
	tfreeze = 0.0
	nP = 0
	nT = 0
        open(51,file='control',status='old')
        read(51,*) P1,P2,nP,T1,T2,nT,iensem,superad,tfreeze
        close (51)
        call readin(binit,bfinal,nbulk)
	print*, 'back from readin'
        nnull = nspec - nc

        if (nc .gt. nspec) then
         print*, '********** ERROR **********'
         print*, 'Number of constraints',nc
         print*, 'exceeds number of species',nspec
         print*, 'Try reducing number of components or '
         print*, 'increasing number of phases present in initial guess'
         print*, '********************'
         stop
        end if

        call sform(s,b,n1,q1,q2,nspec,nco,nc,ncs,nnull,nnulls,absents)
        if (nphpres .gt. nco) then
         print*, '********** ERROR **********'
         print*, 'Phase Rule Violated'
         print*, 'Number of phases',nphpres
         print*, 'exceeds number of components',nco
         print*, 'Try eliminating one or more phases'
         print*, '********************'
         stop
        end if
        if (nnulls .gt. nnull) then
         print*, '********** ERROR **********'
         print*, 'Number of linearly independent components',nco - nnulls + nnull
         print*, 'is less than the total number of components',nco
         print*, 'Try eliminating one or more components'
         print*, '********************'
         stop
        end if

        do 3 i=1,nnull
3       nnew(i) = 0.
	do 4 ispec=1,nspec
	 dn(ispec) = 0.
	 cpcomp(ispec) = 0.
4	continue

        dP = 0.0
        dT = 0.0
	iextra = 0
	nvet = 0
	nvep = 0
c        if (nT .ge. 0 .and. nP .eq. -1) then
c         isochor = .true.
c         nP = 0
c        end if
c        if (nT .eq. -1 .and. nP .ge. 0) then
	if (iensem .eq. -1.) then
         adiabat = .true.
         nT = 0
        end if
c        if (nT .eq. -2 .and. nP .ge. 0) then
        if (iensem .eq. -2) then
         adcalc = .true.
	 iextra = iextra + 1
	 nvet = iextra
c         nT = 0
        end if
c        if (nT .eq. -3 .and. nP .ge. 0) then
        if (iensem .eq. -3) then
         hucalc = .true.
         nT = 0
        end if
        if (iensem .eq. -4) then
         chcalc = .true.
	 iextra = iextra + 1
	 nvep = iextra
	 print*, 'isochore'
        end if
        if (iensem .eq. -5) then
	 adcalc = .true.
         chcalc = .true.
	 iextra = iextra + 2
	 nvet = 1
	 nvep = 2
        end if
	if (iensem .eq. -6) then
	 frozen = .true.
	end if
        if (nP .lt. 0) stop 'nP less than zero'
        if (nT .lt. 0) stop 'nT less than zero'
        if (nP .ne. 0) dP = (P2 - P1)/float(nP)
        if (nT .ne. 0) dT = (T2 - T1)/float(nT)
        Pi = P1
        Ti = T1

        do 22 i=1,nco
         dbulk(i) = 0.0 
         if (nbulk .eq. 0) go to 22
         dbulk(i) = (bfinal(i) - binit(i))/float(nbulk)
22      continue

        call dossetup
        call nform(nnew,n,n1,q2,nspec,nnull)

	zero = 0.
        if (adiabat) print*, 'Read in a temperature profile, the following increment will be added = ',superad
        if (chcalc)  print*, 'Density and temperature are the independent variables.  Initial guess at the pressure = ',zero
        if (adcalc)  print*, 'Pressure and entropy are the independent variables.  Initial guess at the temperature = ',superad
        if (adcalc .and. chcalc)  print*, 'Density and entropy are the independent variables.  
     &     Initial guess at the temperature = ',superad

        return
        end
