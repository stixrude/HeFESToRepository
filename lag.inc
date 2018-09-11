        double precision lags(nspecp),nphase(nphasep),nlag(nspecp),cpcomp(nspecp)
	double precision lagc(ncompp),rep(ncompp,nspecp,nsitep)
        integer iphase(nphasep),mphase(nphasep),lphase(nspecp)
        integer iophase(nphasep),mophase(nphasep)
	integer lph
	logical iastate(nspecp,nspecp,nsitep)
        common /lagcom/ cpcomp,lags,lagc,nphase,nlag,rep,lph,iphase,mphase,lphase,iophase,mophase,iastate
