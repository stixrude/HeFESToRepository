        double precision function volumel(ispec,x1)
        include 'P1'

	integer ispec,ires,jspec,nb,nseg
	double precision x1,apar,fret,Pi,pressurel,Ti,vlan,vlow,vo,vsplow,vspupp,vupp,x2,xx,vsp,plow,pupp
	double precision zeroin,p1,p2
        logical isochor
        double precision xb1(10),xb2(10)
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /chor/ jspec,isochor
        external pressurel
        double precision, parameter :: tol=1.e-12, vsmall=1.0, fracv = 1.e-2
	jspec = ispec
	vlan = apar(jspec,40)
        vlow = apar(jspec,51)
        vupp = apar(jspec,52)
        vsplow = apar(jspec,53)
        vspupp = apar(jspec,54)
C  For the liquid apar(j,51) apar(j,52) are not relevant as these are derived from instabilities in the expression for the vibrational frequency
C  Use the spinodal limits instead, derived from the reference isotherm (bulk modulus->0)
	vlow = vsplow
	vupp = vspupp

        Vo = apar(ispec,6)
C  Fix volume to Vo and return if isochoric conditions have been chosen
        if (isochor) then
         volumel = Vo
         Pi = pressurel(Vo)
         return
        end if
C  Find volume by cage using the last succesfully found volume as guess (x1)
        x2 = x1*(1.0 + fracv)
c	print*, 'In volumel entering zbrac and pressurel'
c	print*, 'In volumel entering cage and pressurel',vlow,vupp
c        call zbrac(pressurel,x1,x2,succes)
	call cage(pressurel,x1,x2,vlow,vupp,ires)
c	print*, 'back from zbrac in volumel',succes,x1,x2
	p1 = pressurel(x1)
	p2 = pressurel(x2)
	plow = pressurel(vlow)
	pupp = pressurel(vupp)
c	print*, 'back from cage in volumel',x1,x2,vlow,vupp,p1,p2
c        if (succes) then
        if (ires .eq. 1) then
c         volumel = zbrent(pressurel,x1,x2,tol)
         volumel = zeroin(x1,x2,pressurel,tol)
c	 print*, 'Found volume',volumel
        else
C  Zbrac failed.  The following logic assumes that if a solution Vsol exists,
C  it satisfies Vo < Vsol < Vsp where Vsp is the spinodal limit.
C  Find spinodal volume at this temperature by finding the volume at which the pressure is a minimum.
C  Assume that T>T_0 and that Vsp(T)<Vsp(T_0)
         print *, 'volume Failed to find V cage',ispec,Ti,Pi,x1,x2,p1,p2,vlow,vupp,plow,pupp
         x1 = Vo - vlan
         x2 = vspupp
c         x1 = Vo - vlan
c         x2 = min(vsp,vupp)
         xx = x2*(1. - tol)
c         call nlmin_VL(xx,x1,x2,fret)
         call nlmin_VL(xx,vlow,vupp,fret)
         vsp = xx
	  print*, 'volumel after nlmin_VL',x1,x2,xx,Vo,vlan,vsp,vlow,vupp,fret
C  Value of function pressurel is positive at Vsp: no solution.
         if (fret .gt. 0.) then
          volumel = -1
	  print*, 'volumel Failed after nlmin_VL',x1,x2,xx,Vo,vlan,vsp,vupp,fret
          return
         end if
C  Find volume by zbrak.  Search between V=Vo and V=Vsp
         nb = 10
	 nseg = 10
c         call cages(pressurel,x1,x2,nseg,xb1,xb2,nb)
         call cages(pressurel,vlow,vsp,nseg,xb1,xb2,nb)
c         call zbrak(pressurel,x1,x2,nseg,xb1,xb2,nb)
c         if (jspec .eq. 8) print*, 'volume zbrak',nb
         if (nb .ge. 1) then
c          volumel = zbrent(pressurel,xb1(1),xb2(1),tol)
          volumel = zeroin(xb1(1),xb2(1),pressurel,tol)
	  print*, 'volumel found',volumel
         else
C  Zbrak failed.
c         if (jspec .eq. 8) print '(a20,i5,5f12.5)', 'volume zbrak failed ',ispec,Ti,Pi,xb1(1),xb2(1),volume
          volumel = -2
	  print*, 'volumel Failed after zbrak',x1,x2,xx,Vo,vlan,vsp,vupp,fret
         end if
        end if

        return
        end
