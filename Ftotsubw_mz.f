        subroutine Ftotsubw(ispec,Vi,Ftot)

        include 'P1'

	integer ispec
      REAL(8) tcr, rhocr, pcr, rcnstw,xmcapw
        double precision Vi,volnl,Cp,Cv,gamma,K,Ks,alp,Ftot,ph,ent,deltas,
     &                   tcal,zeta,Gsh,uth,uto,thet,q,etas,dGdT,pzp,Sel,Eel,Pel,Cvel,Eig,Pig,P,E,w,SNk
	double precision apar,Ti,Pi,volve,entve,cpve,bkve,Fpv
        double precision RHO,T,PnkT,FNkT,UNkT,CHIT,CHIR,PMbar,USPEC

        common /state/ apar(nspecp,nparp),Ti,Pi
	common /volent/ volve,entve,cpve,bkve

        double precision, parameter :: Sconst = 86.808 + 46.315042475974771             !  Recovers JANAF value at 372.78 K in J/mol/K
        double precision, parameter :: Fconst = -236.839 - (-19.399324932220778)        !  Recovers JANAF value of DG_f at 300 K in J/mo

      DATA tcr / 647.096d+00 /, rhocr / 322.0d+00 /,
     $ pcr / 22.064d0 /, rcnstw / 0.46151805d+00 /,
     $ xmcapw / 0.018015268d0 /

        RHO = 1000.*xmcapw/Vi           ! g/cm^3
        T = Ti

         call H2OFIT(RHO,T,PnkT,FNkT,UNkT,CV,CHIT,CHIR,PMbar,USPEC)

	P = PMbar*100.				! GPa
	K = CHIR*P				! GPa
	alp = CHIT/CHIR/Ti			! K-1
	CV = 3.*rcnstw*xmcapw*1000.*CV		! J/mol/K
	gamma = 1000.*alp*K*Vi/CV		! -
	CP = CV*(1. + alp*gamma*Ti)		! J/mol/K
	KS = K*(1. + alp*gamma*Ti)		! GPa
	w = 1000.*sqrt(KS/rho)			! m/s
	SNk = UNkT - FNkT			! -
	ent = 3.*rcnstw*xmcapw*1000.*SNk + Sconst   	! J/mol/K
	Ftot = 3.*rcnstw*xmcapw*1000.*Ti*FNkT - Sconst*Ti + 1000.*Fconst	! J/mol/K
	E = 3.*rcnstw*xmcapw*1000.*Ti*UNkT	! J/mol/K

        Fpv = 1000.*Pi*Vi
 	Ftot = Ftot + Fpv

        volve = Vi
        entve = ent
        cpve = Cp
        bkve = K

c	print*, 'V (cm3/mol)',Vi
c	print*, 'P (GPa)',P
c	print*, 'KT (GPa)',K
c	print*, 'KS (GPa)',KS
c	print*, 'w (m/s)',w
c	print*, 'alp (K-1)',alp
c	print*, 'CV (J/mol/K) (kJ/kg/K)',CV,CV/(xmcapw*1000.)
c	print*, 'CP (J/mol/K) (kJ/kg/K)',CP,CP/(xmcapw*1000.)
c	print*, 'Entropy (J/molK) (kJ/kg/K)',ent,ent/(xmcapw*1000.)
c	print*, 'gamma',gamma
c	print*, 'Helmholtz Free Energy (J/molK) (kJ/kg/K)',Ftot,Ftot/(xmcapw*1000.)
c	print*, 'Internal Energy (J/molK) (kJ/kg/K)',E,E/(xmcapw*1000.)

	return
	end
