	subroutine hillert(Ti,Tc,p,smax,glan,slan,cplan)

C  Expression for the magnetic terms according to Hillert's book, pg. 523.  See also dorogokupetsetal_17 Eqs. 20,21.
C  p is the fraction of the total magnetic enthalpy that is absorbed at T>T_c.

	double precision Tc,p,smax,glan,slan,cplan,t,Ti,Pi,apar
	double precision gfunc,dgdt,d2gdt2
	double precision afe,bfe,cfe,denom

	t = Ti/Tc
C  jacobsschmidfetzer_10 (Eqs. 20a,20b) provide approximate numerical values of coefficients by evaluating the exact expressions given in the references above and assuming p=0.40.
C  Recast so that the Gibbs free energy is zero at zero temperature.
        afe = 0.9053
	bfe = 0.91805
        cfe = 0.64173
	denom = 518./1125. + 11692./15975.*(1./p - 1.)
	afe = 79./(140.*p)/denom
	bfe = 474./497.*(1./p - 1.)/denom
	cfe = 1./denom

        if (t .le. 1.0) then
         gfunc =  -bfe*(t**3/6. + t**9/135. +    t**15/600.)
         dgdt =   -bfe*(t**2/2. + t**8/15. +     t**14/40.)
         d2gdt2 = -bfe*(t +    8.*t**7/15. + 14.*t**13/40.)
        else
         gfunc =    -cfe*(t**(-5)/10. + t**(-15)/315. +    t**(-25)/1500.) - 1. + afe/t
         dgdt =      cfe*(t**(-6)/2. +  t**(-16)/21 +      t**(-26)/60.) - afe/t**2
         d2gdt2 = -cfe*(3*t**(-7) + 16.*t**(-17)/21. + 26.*t**(-27)/60.) + 2.*afe/t**3
        end if
        glan = smax*Ti*gfunc
        slan = -smax*gfunc -smax*t*dgdt
        cplan = -smax*(2.*t*dgdt + t**2*d2gdt2)


	return
	end
