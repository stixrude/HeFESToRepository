        subroutine qr19(z,Ti,qs,qp)
        include 'P1'
        include 'const.inc'
	include 'adqref.inc'
C  Romanowicz, B., A global tomographic model of shear attenuation 
C  in the upper mantle, Journal of Geophysical Research, 100, B7, 
C  12375-12394, 1995.
C  Lower mantle value from PREM
	integer i,jlo,klo,mint,ncall
	double precision z,Ti,qs,qp,dq,dt,pad,qs00,tref
        double precision d(13),q(13)
        character*80 fname
        double precision, parameter :: qlarge=1.e15, age=100., Eact=424000., alpha=0.26, qslarge=9999.
	integer, parameter :: nd=13, nadmax=2000
        data d/2891.,648.,647.,511.,422.,370.,310.,250.,200., 160.
     &       ,100., 81.5,0./
        data q/32.1, 32.1,67.1,66.2,65.1,65.2,67.8,76.6,145.3,154.9
     &       ,153.5,17.1,17.1/
        data ncall/0/
        ncall = ncall + 1
        if (ncall .eq. 1) then
	 write(31,*) 'Base Q model: QR19'
         write(31,*) 'Age for Q computation (Ma) = ',age
         write(31,*) 'Reference temperature profile for Q computation = ',parsetname
         write(31,*) 'Potential temperature for Q computation = ',tad(1)
        end if

        mint = 2
        call bserch(d,nd,z,jlo)
        klo = min(max(jlo-(mint-1)/2,1),nd+1-mint)
        call neville(d(klo),q(klo),mint,z,qs00,dq)
        qs00 = 10000./qs00

        call bserch(dad,nad,z,jlo)
        klo = min(max(jlo-(mint-1)/2,1),nad+1-mint)
        call neville(dad(klo),tad(klo),mint,z,tref,dt)

C  Approximate account of lithospheric age for Pacific after:
C  Romanowicz, B., Attenuation tomography of the upper mantle: A review of
C  current status, Pure and Applied Geophysics, 153, 257-272, 1995 (Figure 3).
C  by
C  Stixrude & Lithgow-Bertelloni, Mineralogy and elasticity of the
C  oceanic upper mantle: Origin of the low velocity zone, Journal of Geophysical
C  Research, 2005. (Figure 9 caption)

        qs = 1000./(1000./qs00 + 3.*(1. - age/100.))

C  Approximate account of temperature dependence of qs added 31 July, 2007
        qs = qs*exp(alpha*Eact/Rgas*(1./Ti - 1./tref))
        qs = min(qs,qslarge)

C  Assume Poisson solid with no bulk attenuation
        qp = 9./4.*qs

C  Test
c       qs = qlarge
c       qp = qlarge

        return
        end
