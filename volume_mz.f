        double precision function volumew(ispec,x1)
        include 'P1'

	integer ispec,ires,jspec,nb,nseg
	double precision x1,apar,fret,Pi,pressurew,Ti,vlan,vlow,vo,vsplow,vspupp,vupp,x2,xx,vsp,plow,pupp
	double precision zeroin,p1,p2
        logical isochor
        double precision xb1(10),xb2(10)
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /chor/ jspec,isochor
        external pressurew
        double precision, parameter :: tol=1.e-12, vsmall=1.0, fracv = 1.e-1
	jspec = ispec

C  Volume at the critical point
	vupp = 55.948757763975155
	vupp = 60.0
C  Volume at 5270 K, 275 GPa
	vlow = 4.0

c        Vo = apar(ispec,6)
C  Fix volume to Vo and return if isochoric conditions have been chosen
        if (isochor) then
         volumew = Vo
         Pi = pressurew(Vo)
         return
        end if
C  Find volume by cage using the last succesfully found volume as guess (x1)
        x2 = x1*(1.0 + fracv)
c	print*, 'In volumew entering zbrac and pressurew'
c	print*, 'In volumew entering cage and pressurew',vlow,vupp,Ti,Pi
c        call zbrac(pressurew,x1,x2,succes)
	call cage(pressurew,x1,x2,vlow,vupp,ires)
	p1 = pressurew(x1)
	p2 = pressurew(x2)
	plow = pressurew(vlow)
	pupp = pressurew(vupp)
c	print*, 'back from cage in volumew',x1,x2,p1,p2
c        if (succes) then
        if (ires .eq. 1) then
c         volumew = zbrent(pressurew,x1,x2,tol)
         volumew = zeroin(x1,x2,pressurew,tol)
c	 print*, 'Found volume',volumew,pressurew(volumew)
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
c         call nlmin_VW(xx,x1,x2,fret)
         call nlmin_VW(xx,vlow,vupp,fret)
         vsp = xx
	  print*, 'volumew after nlmin_VW',x1,x2,xx,Vo,vlan,vsp,vlow,vupp,fret
C  Value of function pressurew is positive at Vsp: no solution.
         if (fret .gt. 0.) then
          volumew = -1
	  print*, 'volumew Failed after nlmin_VW',x1,x2,xx,Vo,vlan,vsp,vupp,fret
          return
         end if
C  Find volume by zbrak.  Search between V=Vo and V=Vsp
         nb = 10
	 nseg = 10
c         call cages(pressurew,x1,x2,nseg,xb1,xb2,nb)
         call cages(pressurew,vlow,vsp,nseg,xb1,xb2,nb)
c         call zbrak(pressurew,x1,x2,nseg,xb1,xb2,nb)
c         if (jspec .eq. 8) print*, 'volume zbrak',nb
         if (nb .ge. 1) then
c          volumew = zbrent(pressurew,xb1(1),xb2(1),tol)
          volumew = zeroin(xb1(1),xb2(1),pressurew,tol)
	  print*, 'volumew found',volumew
         else
C  Zbrak failed.
c         if (jspec .eq. 8) print '(a20,i5,5f12.5)', 'volume zbrak failed ',ispec,Ti,Pi,xb1(1),xb2(1),volume
          volumew = -2
	  print*, 'volumew Failed after zbrak',x1,x2,xx,Vo,vlan,vsp,vupp,fret
         end if
        end if

        return
        end
