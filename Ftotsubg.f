        subroutine Ftotsubg(ispec,Vi,Ftot)

        include 'P1'
        include 'const.inc'
        include 'theory.inc'
	
	logical sfix,vofix,cvfix
	integer ispec,i,j,nbm,nobm,noth,mfit,maxij,lineart,jelem,noln
	double precision vi,volnl,cp,cv,gamma,alp,ftot,ph,ent,deltas,tcal,zeta,gsh,uth,uto,thet,q,etas
	double precision pzp,akt,aktel,aktig,aktxs,apar,cvel,cvig,cvxs,d2fdv2,d2tdt2,dgdt
	double precision dfdv,dtdt,e,eel,eig,exs,f,fac,fel
	double precision fig,fmth,fn,fxs,Pi,pig,sel,sig,sxs,theta,Ti,to,vo,dfac
        double precision K,Ks,Kxs,Kig,Kel,Kelp,P,Pxs,Pel
        double precision daktdvel,betael
	double precision aliq(nparp,nparp),aliqc,cliq0,cliq1,cliq2,tee,dteedt,acof,bcof
	double precision betaig,daktdvig,Kigp
	double precision volve,entve,cpve,bkve
        double precision acp,bcp,ccp,dcp,ecp,fcp,gcp,hcp,H,Href,S,Sref,So
        double precision, parameter :: Tglass = 1500.0
        common /volent/ volve,entve,cpve,bkve
        common /state/ apar(nspecp,nparp),Ti,Pi
	common /liqc/ aliqc(nparp,nparp),mfit,nobm,noth,sfix,vofix,cvfix,maxij,lineart,nbm,noln
	double precision, parameter :: fsmall=1.e-12
        double precision, parameter :: P1bar = 1.e-4

        Vo = apar(ispec,6)
        To = apar(ispec,4)
        So = apar(ispec,10)
        acp = apar(ispec,11)
        bcp = apar(ispec,12)
        ccp = apar(ispec,13)
        dcp = apar(ispec,14)
        ecp = apar(ispec,15)
        fcp = 0.
        gcp = 0.
        hcp = 0.

        Cp = acp + bcp*Ti + ccp/Ti**2 + dcp/sqrt(Ti) + ecp*Ti**2
        Cv = Cp - fn*Rgas*Ti
        E = Cv*Ti                                   ! J/mol

        S =    acp*log(Ti) + bcp*Ti - 0.5*ccp*Ti**(-2) - 2.*dcp*Ti**(-0.5)
     &    + 0.5*ecp*Ti**2 - 1./3.*fcp*Ti**(-3)
     &    - gcp*Ti**(-1) + 0.5*hcp*(log(Ti))**2
        Sref =  acp*log(To) + bcp*To - 0.5*ccp*To**(-2) - 2.*dcp*To**(-0.5)
     &    + 0.5*ecp*To**2 - 1./3.*fcp*To**(-3)
     &    - gcp*To**(-1) + 0.5*hcp*(log(To))**2
        ent = S - Sref + So

        H =    acp*Ti + 0.5*bcp*Ti**2 - ccp*Ti**(-1) + 2.*dcp*Ti**(+0.5)
     &    + 1./3.*ecp*Ti**3 - 0.5*fcp*Ti**(-2) + gcp*log(Ti)
     &    + hcp*(Ti*log(Ti) - Ti)
        Href = acp*To + 0.5*bcp*To**2 - ccp*To**(-1) + 2.*dcp*To**(+0.5)
     &    + 1./3.*ecp*To**3 - 0.5*fcp*To**(-2) + gcp*log(To)
     &    + hcp*(To*log(To) - To)
        H = H - Href

        Ftot = H - Ti*ent + To*So + Rgas*Ti*log(Pi/P1bar)
	
        return
        end
