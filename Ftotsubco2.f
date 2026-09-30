	subroutine Ftotsubco2(ispec,Vi,Ftot)
       USE def_constants
       USE properties
	include 'P1'
	integer ispec,ncall
        double precision Vi,volnl,Cp,Cv,gamma,K,Ks,alp,Ftot,ph,ent,deltas,
     &                   tcal,zeta,Gsh,uth,uto,thet,q,etas,dGdT,pzp,Sel,Eel,Pel,Cvel,Eig,Pig,P,E
	double precision T,v,c,s,h_helmho,g,KT,enthalpy,deriv,alpha
	double precision apar,Ti,Pi,volve,entve,cpve,bkve,Fpv
        double precision, parameter :: Sconst = 214.025d0 - 0.13522408887259521d0            !  Recovers JANAF value at 300 K in J/mol/K
        double precision, parameter :: Fconst = -394.394d0 - (-64.17945442266272d0)          !  Recovers JANAF value of DG_f at 300 K in kJ/mol/K
        common /volent/ volve,entve,cpve,bkve
        common /state/ apar(nspecp,nparp),Ti,Pi
        data ncall/0/
        ncall = ncall + 1
        if (ncall .eq. 1) write(31,*) "INFORMATION: CO2 FTR OF SPAN AND WAGNER 1996"
	ispec = 1

        v = Vi/W/1.d6                           ! m^3/kg

	CALL inter_energy(Ti,v,e)
	CALL heat_cap_v(Ti,v,cv)
	CALL heat_cap_p(Ti,v,cp)
	CALL sound_speed(Ti,v,c)
	KS = c**2/v/1.d9
	KT = cv/cp*KS
	CALL entropy(Ti,v,s)
	CALL helmho(Ti,v,h_helmho)
	CALL gibbs(Ti,v,g)
	enthalpy = g + Ti*s
	call dpdT_v(Ti,v,deriv)
	alpha = deriv/KT/1.d9
        gamma = alpha/KT/(cv/v)*1.e9   

        Cp = cp*W                       ! J/mol/K
        Cv = cv*W                       ! J/mol/K
        K = KT                          ! GPa
        KS = KS                         ! GPa
        alp = alpha                     ! 1/K
        Ftot = 1000.*Fconst + h_helmho*W - Sconst*Ti               ! J/mol
        ent = s*W + Sconst              ! J/mol/K
        P = Pi                          ! GPa
        E = e*W                         ! J/mol
        Fpv = 1000.*Pi*Vi               ! J/mol
        Ftot = Ftot + Fpv               ! J/mol
        E = Ftot + ent*Ti               ! J/mol

        volve = Vi
        entve = ent
        cpve = Cp
        bkve = K

	return
	end
