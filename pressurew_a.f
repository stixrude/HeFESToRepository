	double precision function pressurew(Vi)

	include 'P1'

	logical isochor
	integer ispec
	double precision apar,Ti,Pi,Vi
	double precision tcr,rhocr,pcr,rcnstw,xmcapw
	double precision RHO,T,PnkT,FNkT,UNkT,CV,CHIT,CHIR,PMbar,USPEC

        common /state/ apar(nspecp,nparp),Ti,Pi
        common /chor/ ispec,isochor

      DATA tcr / 647.096d+00 /, rhocr / 322.0d+00 /,
     $ pcr / 22.064d0 /, rcnstw / 0.46151805d+00 /,
     $ xmcapw / 0.018015268d0 /

	RHO = 1000.*xmcapw/Vi		! g/cm^3
	T = Ti

         call H2OFIT(RHO,T,PnkT,FNkT,UNkT,CV,CHIT,CHIR,PMbar,USPEC)

	pressurew = PMbar*100.		! GPa

	pressurew = pressurew - Pi

	return
	end
