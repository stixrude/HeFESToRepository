        double precision function volumel(ispec,x1)
        include 'P1'

	integer ispec,ires,jspec,nb,nseg
	double precision x1,apar,fret,Pi,pressurel,Ti,vlan,vlow,vo,vsplow,vspupp,vupp,x2,xx,vsp,plow,pupp
	double precision zeroin,p1,p2,fn,videal,videalgas
	double precision htl,pv,pxb1,pxb2
	double precision vlowgraph,vuppgraph,dv,v,one,vsp2
	integer iv,nv
        logical isochor,igraph
        double precision xb1(10),xb2(10)
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /chor/ jspec,isochor
        external pressurel
        double precision, parameter :: tol=1.e-12, vsmall=1.0, fracv = 1.e-2
	one = 1.0
	igraph = .false.

C  Set bounds on the volume
	jspec = ispec
	vlan = apar(jspec,40)
        vlow = apar(jspec,51)
        vupp = apar(jspec,52)
        vsplow = apar(jspec,53)
        vspupp = apar(jspec,54)
C  For the liquid apar(j,51) apar(j,52) are not relevant as these are derived from instabilities in the expression for the vibrational frequency
C  Use the spinodal limits instead, derived from the reference isotherm (bulk modulus->0)

	fn = apar(jspec,1)
	htl = apar(jspec,31)
	videalgas = videal(fn,Pi,Ti)
C  --> Testing: Make graph of p vs. v
	if (igraph) then
	 vlowgraph = vlow
	 vuppgraph = 10.0*videalgas
	 open(121,file='vgraph.txt',status='unknown')
	 nv = 1000
	 dv = (log10(vuppgraph) - log10(vlowgraph))/float(nv)
	 do 1 iv=1,nv+1
	  v = vlowgraph*10**(dv*float(iv-1))
	  write(121,*) v,pressurel(v)+Pi,Pi,Ti,videalgas
1	 continue
	 close (121)
	end if
C  <-- 
	if (htl .eq. 4.0) then
	 vlow = vspupp
	 vupp = max(3.0*videalgas,vspupp)
	 if (x1 .gt. vupp .or. x1 .lt. vlow) then
	  x1 = videalgas
	  if (x1 .gt. vupp .or. x1 .lt. vlow) x1 = (vlow + vupp)/2.
	 end if
	end if
	if (htl .eq. 1.0) then
	 vlow = vsplow
	 vupp = 3.0*vspupp
	 if (x1 .gt. vupp .or. x1 .lt. vlow) x1 = Vo
	end if
c	vlow = vsplow
c	vupp = max(3.0*videalgas,vspupp)

        Vo = apar(ispec,6)
C  Fix volume to Vo and return if isochoric conditions have been chosen
        if (isochor) then
         volumel = Vo
         Pi = pressurel(Vo)
         return
        end if
C  Find volume by cage using the last succesfully found volume as guess (x1)
        x2 = x1*(1.0 + fracv)
c        write(31,'(a,i5,99e12.5)') 'In volumel entering cage and pressure',ispec,Pi,Ti,vlow,vupp,x1
c        print '(a,i5,99e12.5)', 'In volumel entering cage and pressure',ispec,Pi,Ti,vlow,vupp,x1
	if (htl .eq. 4.0 .and. Pi .le. 0.) then
	 volumel = -1
	 go to 10
	end if
	call cage(pressurel,x1,x2,vlow,vupp,ires)
	p1 = pressurel(x1)
	p2 = pressurel(x2)
	plow = pressurel(vlow)
	pupp = pressurel(vupp)
        if (ires .eq. 1) then
         volumel = zeroin(x1,x2,pressurel,tol)
c         write(31,'(a,2i5,99e12.5)') 'Found volumel by zbrac',ispec,ires,Pi,Ti,volumel,x1,x2,p1,p2,vlow,vupp,plow,pupp
c         print '(a,2i5,99e12.5)', 'Found volumel by zbrac',ispec,ires,Pi,Ti,volumel,x1,x2,p1,p2,vlow,vupp,plow,pupp
	 return
	end if
C  Zbrac failed.  Try zbrak, first seach for extrema.
c        write(31,'(a,2i5,99e12.5)') 'volumel failed to find V cage',ispec,ires,Pi,Ti,x1,x2,p1,p2,vlow,vupp,plow,pupp
c        print '(a,2i5,99e12.5)', 'volumel failed to find V cage',ispec,ires,Pi,Ti,x1,x2,p1,p2,vlow,vupp,plow,pupp
C  Search for the asbolute minimum in the pressure along the isotherm T=Ti
        x1 = vsplow
        x2 = 3.*vspupp
        xx = x1*(1. + tol)
        call nlmin_VL(xx,x1,x2,fret,one,ires)
        vsp = xx
c        write(31,'(a,2i5,99e12.5)') 'volumel Pmin',ispec,ires,Pi,Ti,vsp,fret+Pi,x1,x2,fret
c        print '(a,2i5,99e12.5)', 'volumel Pmin',ispec,ires,Pi,Ti,vsp,fret+Pi,x1,x2,fret,vsp-x1,vsp-x2
	if (ires .lt. 0 .and. ires .ne. -4) then
c	 print '(a,i5)', 'volumel failed to find minimum pressure: nlmin_VL failed',ires
c	 volumel = -2
c	 go to 10
	end if
	if (ires .gt. 0 .or. ires .eq. -4) then
	 if (abs(vsp-x2) .le. tol) then
c	  print '(a,99e12.5)', 'volumel failed to find minimum pressure: solution equals upper bound',vsp,fret+Pi,x2
	  if (htl .eq. 1.0) then
c	   volumel = -1
c	   go to 10
	  end if
	 end if
	 if (abs(vsp-x1) .le. tol) then
c	  print '(a,99e12.5)', 'volumel failed to find minimum pressure: solution equals lower bound',vsp,fret+Pi,x1
	 end if
         if (Pi .lt. 0. .and. fret .gt. 0.) then
c	  print '(a,99e12.5)', 'volumel: Pi<0 and Pi<Pmin. spinodal instability.',vsp,Pi,fret+Pi
c	  volumel = -1
c	  go to 10
         end if
	 if (Pi .gt. 0. .and. fret .gt. 0.) then
c	  print '(a,99e12.5)', 'volumel: Pi>0 and Pi<Pmin. liquid spinodally unstable.',vsp,Pi,fret+Pi
	  if (htl .eq. 1.0) then
c	   volumel = -1
c	   go to 10
	  end if
	 end if
	end if
C  Search for the local maximum along the isotherm T=Ti
        x1 = (1. + 0.01)*vsp
        x2 = 2.*videalgas
        xx = x1*(1. + tol)
        call nlmin_VL(xx,x1,x2,fret,-one,ires)
	fret = pressurel(xx)
	vsp2 = xx
c        write(31,*) 'volumel Pmax',ispec,Pi,Ti,vsp2,fret+Pi,x1,x2,fret,ires
c        print '(a,2i5,99e12.5)', 'volumel Pmax',ispec,ires,Pi,Ti,vsp2,fret+Pi,x1,x2,fret
	if (ires .lt. 0 .and. ires .ne. -4) then
c	 print '(a,i5)', 'volumel failed to find maximum pressure: nlmin_VL failed',ires
c	 volumel = -2
c	 go to 10
	end if
	if (ires .gt. 0 .or. ires .eq. -4) then
	 if (fret+Pi .lt. 0.) then
c	  print '(a,99e12.5)', 'volumel failed to find maximum pressure: Pmax<0',vsp2,fret+Pi
c	  volumel = -2
c	  go to 10
	 end if
	 if (abs(vsp2-x1) .le. tol) then
c	  print '(a,99e12.5)', 'volumel failed to find maximum pressure: solution equals lower bound',vsp2,fret+Pi,x1
c	  volumel = -2
c	  go to 10
	 end if
	 if (abs(vsp2-x2) .le. tol) then
c	  print '(a,99e12.5)', 'volumel failed to find maximum pressure: solution equals upper bound',vsp2,fret+Pi,x2
c	  if (htl .eq. 1.0) x1 = Vo
	  if (htl .eq. 1.0) continue
	  if (htl .eq. 4.0) then
c	   volumel = -2
c	   go to 10
	  end if
	 end if
	 if (fret .lt. 0.) then
c	  print '(a,99e12.5)', 'volumel: Pi>Pmax. vapor spinodally unstable.',vsp2,Pi,fret+Pi
	  if (htl .eq. 4.0) then
c	   volumel = -1
c	   go to 10
	  end if
	 end if
	end if
C  Find volume by zbrak.  Search between vlow and vsp for liquid.  Search between vsp2 and vupp for vapor.
C  Logic assumes that Pmin<Pi<Pmax
c        write(31,'(a,i5,99e12.5)') 'starting zbrak',ispec,Pi,Ti,volumel,x1,x2,xb1(1),xb2(1),pxb1,pxb2,pv
c        print '(a,i5,99e12.5)', 'starting zbrak',ispec,Pi,Ti,volumel,x1,x2,xb1(1),xb2(1),pxb1,pxb2,pv
	if (htl .eq. 1) then
	 x1 = vlow
	 x2 = vsp
	end if
	if (htl .eq. 4) then
	 x1 = vsp2
	 x2 = vupp
	end if
	volumel = 0
c        write(31,'(a,i5,99e12.5)') 'starting zbrak',ispec,Pi,Ti,x1,x2
c        print '(a,i5,99e12.5)', 'starting zbrak',ispec,Pi,Ti,x1,x2
        nb = 10
	nseg = 10
        call cages(pressurel,x1,x2,nseg,xb1,xb2,nb)
        if (nb .ge. 1) then
         volumel = zeroin(xb1(1),xb2(1),pressurel,tol)
         pxb1 = pressurel(xb1(1))
         pxb2 = pressurel(xb2(1))
         pv = pressurel(volumel)
c         write(31,'(a,i5,99e12.5)') 'volumel by zbrak',ispec,Pi,Ti,volumel,x1,x2,xb1(1),xb2(1),pxb1,pxb2,pv
c         print '(a,i5,99e12.5)', 'volumel by zbrak',ispec,Pi,Ti,volumel,x1,x2,xb1(1),xb2(1),pxb1,pxb2,pv
	 return
        end if
        if (nb .lt. 1) then
C  Zbrak failed.
c         write(31,'(a,i5,99e12.5)') 'WARNING: Volumel failed to find V zbrak',ispec,Pi,Ti,volumel,x1,x2,xb1(1),xb2(1),pxb1,pxb2,pv
c         print '(a,i5,99e12.5)', 'WARNING: volumel failed to find V zbrak',ispec,Pi,Ti,volumel,x1,x2,xb1(1),xb2(1),pxb1,pxb2,pv
         volumel = -2
	 go to 10
        end if

	return

10	continue
	if (volumel .eq. -1 .or. volumel .eq. -2) then
	 if (htl .eq. 1.0) x1 = Vo
	 if (htl .eq. 4.0) x1 = videalgas
	end if

        return
        end
