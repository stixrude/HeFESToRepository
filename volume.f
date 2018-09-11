        function volume(ispec,x1)
        include 'P1'

C  Find the volume of species ispec, given a pressure Pi, and temperature Ti
	double precision x1,apar,fret,Pi,pressure,Ti,vlan,vlow,Vo,vsplow,vspupp,vupp,vsp
	double precision x2,xx,volume,zeroin,pxb1,pxb2,pv
	integer ispec,ires,jspec,nb,nseg
        logical isochor
        double precision xb1(10),xb2(10)
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /chor/ jspec,isochor
        external pressure
        double precision, parameter :: tol=1.e-15, vsmall=1.0

	jspec = ispec
	vlan = apar(jspec,40)
	vlow = apar(jspec,51)
	vupp = apar(jspec,52)
	vsplow = apar(jspec,53)
	vspupp = apar(jspec,54)

        Vo = apar(ispec,6)
c	write(31,*) 'vibrational and spinodal limits = ',ispec,vlow,vupp,vsplow,vspupp,x1
C  Fix volume to Vo and return if isochoric conditions have been chosen
        if (isochor) then
         volume = Vo
         Pi = pressure(Vo)
         return
        end if
C  Find volume by cage using the last succesfully found volume as guess (x1)
        x2 = x1*(1.0 + tol)
c	write(31,*) 'Calling cage the first time',vlow,vupp
	call cage(pressure,x1,x2,vlow,vupp,ires)
	if (ires .eq. 1) then
         volume = zeroin(x1,x2,pressure,tol)
        else
C  Zbrac failed.  The following logic assumes that if a solution Vsol exists,
C  it satisfies Vo < Vsol < min(Vsp,Vupp) where Vsp is the spinodal and Vupp is the vibrational limit.
C  Find spinodal volume at this temperature.
C  Assume that T>T_0 and that Vsp(T)<Vsp(T_0)
c         write(31,*) 'volume Failed to find V cage',ispec,Pi,Ti,x1,x2
         x1 = Vo - vlan
         x2 = min(vspupp,vupp)
	 xx = x2*(1. - tol)
	 call nlmin_V(xx,x1,x2,fret,ires)
	 vsp = xx
         write(31,*) 'volume spinodal',ispec,Pi,Ti,x1,x2,vsp,fret,ires
C  Value of function pressure is positive at Vsp: no solution.
         if (fret .gt. 0.) then
	  x1 = vsp
          volume = -1
          return
         end if
C  Find volume by zbrak.  Search between V=Vo and V=Vsp
         nb = 10
	 nseg = 10
	 call cages(pressure,x1,x2,nseg,xb1,xb2,nb)
         if (nb .ge. 1) then
          volume = zeroin(xb1(1),xb2(1),pressure,tol)
	  pxb1 = pressure(xb1(1))
	  pxb2 = pressure(xb2(1))
	  pv = pressure(volume)
c          write(31,*) 'volume zeroin',ispec,Pi,Ti,x1,x2,xb1(1),xb2(1),pxb1,pxb2,volume,pv
         else
C  Zbrak failed.
          volume = -2
         end if
        end if

        return
        end
